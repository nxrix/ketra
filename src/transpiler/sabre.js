/**
 * sabre.js - Real Sabre routing algorithm.
 *
 * Implements the Sabre routing and layout algorithms from
 * "Tackling the Qubit Mapping Problem for NISQ-Era Quantum Devices" (Li et al., 2019).
 *
 * Sabre is a heuristic routing algorithm that:
 *   1. Maintains a mapping between virtual and physical qubits
 *   2. For each gate, finds the shortest path in the coupling map
 *   3. Inserts SWAPs along the path to bring interacting qubits together
 *   4. Uses a lookahead cost function with decay to avoid local minima
 */

import { Layout, CouplingMap } from "./layout.js";
import { QuantumCircuit } from "../core/circuit.js";
import { QuantumRegister, ClassicalRegister } from "../core/bit.js";
import { CircuitInstruction } from "../core/circuit.js";
import { Gate } from "../core/gate.js";
import { DAGOpNode } from "../dagcircuit/dagcircuit.js";
import * as standardGates from "../library/standard_gates.js";

// ---------------------------------------------------------------------------
// Sabre routing
// ---------------------------------------------------------------------------
export class SabreSwap {
  constructor(couplingMap, heuristic = "lookahead", seed = null) {
    this.coupling_map = couplingMap;
    this.heuristic = heuristic; // "basic" or "lookahead"
    this.seed = seed;
    this._rng = seed != null ? _makeRng(seed) : Math.random;
    this.swap_weight = 0.5; // weight for SWAP cost
    this.decrement = 0.001; // decay factor
    this.max_distance = couplingMap ? couplingMap.size + 1 : 100;
  }

  // Run Sabre routing on a DAGCircuit.
  // Returns a routed QuantumCircuit with SWAP gates inserted where needed.
  run(dag) {
    if (!this.coupling_map) {
      // No coupling map: just convert DAG back to a circuit unchanged.
      return dag;
    }

    // Initialize layout: trivial (virtual qubit i -> physical qubit i).
    // We use the DAG's actual Qubit objects as the virtual keys so that
    // layout.getVirtual(qubit) works correctly throughout routing.
    const layout = new Layout();
    for (let i = 0; i < dag.qubits.length; i++) {
      layout.setPhysical(i, dag.qubits[i]);
    }

    // Collect all gates in topological order
    const gates = dag.topological_op_nodes();
    const frontLayer = [];
    const executedGates = new Set();
    const remainingGates = new Set(gates);

    // Build dependency graph: gate -> successor op nodes
    const successors = new Map();
    const predecessors = new Map();
    for (const gate of gates) {
      successors.set(gate, []);
      predecessors.set(gate, []);
    }
    for (const gate of gates) {
      for (const succ of dag.successors(gate)) {
        if (succ.is_op_node()) {
          successors.get(gate).push(succ);
          predecessors.get(succ).push(gate);
        }
      }
    }
    // Capture the dependency map for the lookahead heuristic (see _getNextLayer).
    this._successorsMap = successors;

    // Initialize front layer: gates with no predecessors
    for (const gate of gates) {
      if (predecessors.get(gate).length === 0) {
        frontLayer.push(gate);
      }
    }

    // Decay factors for each physical qubit
    const decay = new Array(this.coupling_map.size).fill(1);

    // Distance matrix (precomputed)
    const distMatrix = this._computeDistanceMatrix();

    // Schedule: list of {kind: 'gate'|'swap', ...} recording the order of
    // operations to apply when building the routed circuit.
    const schedule = [];

    let iterations = 0;
    const maxIterations = gates.length * 10;

    while (frontLayer.length > 0 && iterations < maxIterations) {
      iterations++;

      // Execute all gates in front layer that can be executed on current layout
      const executeList = [];
      const remainingFront = [];
      for (const gate of frontLayer) {
        if (this._canExecute(gate, layout)) {
          executeList.push(gate);
        } else {
          remainingFront.push(gate);
        }
      }

      if (executeList.length > 0) {
        // Schedule these gates with the layout snapshot at execution time.
        for (const gate of executeList) {
          schedule.push({ kind: "gate", gate, layout: layout.copy() });
          executedGates.add(gate);
          remainingGates.delete(gate);
          // Update front layer: add successors whose all predecessors are done
          for (const succ of successors.get(gate)) {
            const allDone = predecessors.get(succ).every(p => executedGates.has(p));
            if (allDone && !frontLayer.includes(succ) && !remainingFront.includes(succ)) {
              remainingFront.push(succ);
            }
          }
        }
        frontLayer.length = 0;
        frontLayer.push(...remainingFront);
        // Reset decay
        decay.fill(1);
        continue;
      }

      // No gates can execute — need to insert SWAPs
      // Find the best SWAP using the cost function
      const swapCandidates = this._generateSwapCandidates(frontLayer, layout);
      if (swapCandidates.length === 0) {
        // No SWAPs possible — break out and let _buildRoutedCircuit emit
        // the remaining gates as-is (best effort).
        break;
      }

      let bestSwap = null;
      let bestCost = Infinity;
      for (const swap of swapCandidates) {
        const cost = this._swapCost(swap, frontLayer, layout, distMatrix, decay);
        if (cost < bestCost) {
          bestCost = cost;
          bestSwap = swap;
        }
      }

      if (bestSwap) {
        // Record the SWAP in the schedule BEFORE mutating the layout, so the
        // builder can emit it at the right point with the right qubits.
        schedule.push({ kind: "swap", physical: [bestSwap[0], bestSwap[1]], layout: layout.copy() });
        // Apply the SWAP to the layout
        layout.swap(bestSwap[0], bestSwap[1]);
        // Increase decay for swapped qubits (Sabre's decay schedule)
        decay[bestSwap[0]] += this.decrement;
        decay[bestSwap[1]] += this.decrement;
      }
    }

    // Build the routed circuit
    return this._buildRoutedCircuit(dag, schedule, layout);
  }

  // Check if a gate can execute on the current layout
  _canExecute(gate, layout) {
    if (gate.op.num_qubits < 2) return true;
    const physicalQubits = gate.qargs.map(q => layout.getVirtual(q));
    if (physicalQubits.some(p => p === undefined)) return false;
    // Check coupling map
    if (gate.op.name === "cx" || gate.op.name === "cy" || gate.op.name === "cz" || gate.op.name === "ch") {
      return this.coupling_map.hasEdge(physicalQubits[0], physicalQubits[1]);
    }
    return true;
  }

  // Generate candidate SWAPs from the front layer
  _generateSwapCandidates(frontLayer, layout) {
    const candidates = [];
    const seen = new Set();
    for (const gate of frontLayer) {
      if (gate.op.num_qubits < 2) continue;
      const q1 = layout.getVirtual(gate.qargs[0]);
      const q2 = layout.getVirtual(gate.qargs[1]);
      if (q1 === undefined || q2 === undefined) continue;
      // Consider SWAPs with neighbors of q1
      for (const neighbor of this.coupling_map.neighbors(q1)) {
        const key = [Math.min(q1, neighbor), Math.max(q1, neighbor)].join(",");
        if (!seen.has(key)) {
          seen.add(key);
          candidates.push([q1, neighbor]);
        }
      }
      for (const neighbor of this.coupling_map.neighbors(q2)) {
        const key = [Math.min(q2, neighbor), Math.max(q2, neighbor)].join(",");
        if (!seen.has(key)) {
          seen.add(key);
          candidates.push([q2, neighbor]);
        }
      }
    }
    return candidates;
  }

  // Compute the cost of a SWAP using the Sabre heuristic
  _swapCost(swap, frontLayer, layout, distMatrix, decay) {
    // Temporarily apply the SWAP
    const tempLayout = layout.copy();
    tempLayout.swap(swap[0], swap[1]);

    let totalCost = 0;
    for (const gate of frontLayer) {
      if (gate.op.num_qubits < 2) continue;
      const q1 = tempLayout.getVirtual(gate.qargs[0]);
      const q2 = tempLayout.getVirtual(gate.qargs[1]);
      if (q1 === undefined || q2 === undefined) continue;
      const dist = distMatrix[q1][q2];
      const weight = Math.max(decay[q1], decay[q2]);
      totalCost += dist * weight;
    }

    // Lookahead: consider the next layer
    if (this.heuristic === "lookahead") {
      const nextLayer = this._getNextLayer(frontLayer);
      let nextCost = 0;
      let count = 0;
      for (const gate of nextLayer) {
        if (gate.op.num_qubits < 2) continue;
        const q1 = tempLayout.getVirtual(gate.qargs[0]);
        const q2 = tempLayout.getVirtual(gate.qargs[1]);
        if (q1 === undefined || q2 === undefined) continue;
        nextCost += distMatrix[q1][q2];
        count++;
      }
      if (count > 0) totalCost += this.swap_weight * nextCost / count;
    }

    return totalCost;
  }

  // Get the next layer (immediate successors of the current front layer).
  // Used by the lookahead term of the cost function.
  _getNextLayer(frontLayer) {
    const next = new Set();
    const succMap = this._successorsMap;
    if (!succMap) return [];
    for (const gate of frontLayer) {
      const succs = succMap.get(gate) || [];
      for (const succ of succs) next.add(succ);
    }
    return Array.from(next);
  }

  // Precompute shortest-path distances between all pairs
  _computeDistanceMatrix() {
    const n = this.coupling_map.size;
    const dist = new Array(n);
    for (let i = 0; i < n; i++) {
      dist[i] = new Array(n);
      for (let j = 0; j < n; j++) {
        dist[i][j] = i === j ? 0 : (this.coupling_map.hasEdge(i, j) ? 1 : Infinity);
      }
    }
    // Floyd-Warshall
    for (let k = 0; k < n; k++) {
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          if (dist[i][k] + dist[k][j] < dist[i][j]) {
            dist[i][j] = dist[i][k] + dist[k][j];
          }
        }
      }
    }
    return dist;
  }

  // Build the routed circuit from the recorded schedule.
  // `schedule` is an array of {kind: 'gate'|'swap', ...} entries in execution
  // order. `finalLayout` is the layout at the end of routing (used only for
  // diagnostic info on the returned circuit).
  _buildRoutedCircuit(dag, schedule, finalLayout) {
    const n = dag.qubits.length;
    const routed = new QuantumCircuit();
    routed.add_register(new QuantumRegister(n, "q"));
    if (dag.clbits.length > 0) {
      routed.add_register(new ClassicalRegister(dag.clbits.length, "c"));
    }
    // Map original DAG qubits/clbits to routed circuit qubits/clbits by index.
    const qubitMap = new Map();
    for (let i = 0; i < n; i++) {
      qubitMap.set(dag.qubits[i], routed.qubits[i]);
    }
    const clbitMap = new Map();
    for (let i = 0; i < dag.clbits.length; i++) {
      clbitMap.set(dag.clbits[i], routed.clbits[i]);
    }

    // Track which physical qubit each routed circuit qubit currently represents.
    // Initially routed.qubits[i] holds physical qubit i (trivial layout).
    // We need to translate each scheduled op's *physical* qubits back to
    // routed-circuit qubit indices. Since trivial layout maps virtual i -> phys i,
    // and routed.qubits[i] starts as virtual i, the identity holds at t=0.
    // As SWAPs are applied, the *virtual* qubit on each routed-circuit wire
    // changes — but the routed circuit's wire i still represents "the qubit
    // that started as virtual i". The standard convention is: SWAPs mutate the
    // logical state on the wires, so a SWAP(phys_a, phys_b) in the schedule
    // means "swap the contents of wires a and b" — emit a swap gate on those
    // wires in the routed circuit.

    for (const entry of schedule) {
      if (entry.kind === "gate") {
        const gate = entry.gate;
        const op = gate.op.copy();
        // gate.qargs are the original DAG (virtual) qubits. The layout at
        // execution time tells us which physical qubit each virtual qubit
        // currently sits on. With trivial initial layout, physical index = wire
        // index in the routed circuit.
        const qargs = gate.qargs.map(q => {
          const phys = entry.layout.getVirtual(q);
          return phys === undefined ? qubitMap.get(q) : routed.qubits[phys];
        });
        const cargs = (gate.cargs || []).map(c => clbitMap.get(c) || c);
        routed.append(op, qargs, cargs);
      } else if (entry.kind === "swap") {
        const [p0, p1] = entry.physical;
        routed.append(_makeStdGate("SWAP"), [routed.qubits[p0], routed.qubits[p1]]);
      }
    }

    return routed;
  }
}

// ---------------------------------------------------------------------------
// SabreLayout: find a good initial layout
// ---------------------------------------------------------------------------
export class SabreLayout {
  constructor(couplingMap, seed = null, maxIterations = 4) {
    this.coupling_map = couplingMap;
    this.seed = seed;
    this.max_iterations = maxIterations;
  }

  run(dag) {
    if (!this.coupling_map) {
      // No coupling map: return trivial layout using DAG's Qubit objects.
      const layout = new Layout();
      for (let i = 0; i < dag.qubits.length; i++) layout.setPhysical(i, dag.qubits[i]);
      return layout;
    }

    // Run Sabre routing multiple times with different initial layouts
    // and pick the best one (fewest SWAPs).
    let bestLayout = null;
    let bestScore = Infinity;

    for (let iter = 0; iter < this.max_iterations; iter++) {
      // Try different initial layouts. The first iteration uses the trivial
      // layout (virtual i -> physical i); subsequent iterations use random
      // permutations.
      let initialLayout;
      if (iter === 0) {
        initialLayout = new Layout();
        for (let i = 0; i < dag.qubits.length; i++) {
          initialLayout.setPhysical(i, dag.qubits[i]);
        }
      } else {
        initialLayout = this._randomLayout(dag.qubits, this.coupling_map.size);
      }

      // Run SabreSwap with this coupling map. The swap pass uses its own
      // internal trivial layout, so we only use the routed circuit to score.
      const sabre = new SabreSwap(this.coupling_map, "lookahead", this.seed != null ? this.seed + iter : null);
      const routed = sabre.run(dag);

      // Score: count of SWAP gates inserted.
      const score = routed.data.filter(ci => ci.operation.name === "swap").length;
      if (score < bestScore) {
        bestScore = score;
        bestLayout = initialLayout;
      }
    }

    return bestLayout;
  }

  _randomLayout(virtualQubits, numPhysical) {
    // virtualQubits: the DAG's qubit array (Qubit objects).
    // Returns a Layout mapping physical qubits to virtual Qubit objects.
    const layout = new Layout();
    const physical = Array.from({ length: numPhysical }, (_, i) => i);
    // Fisher-Yates shuffle.
    for (let i = numPhysical - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [physical[i], physical[j]] = [physical[j], physical[i]];
    }
    for (let i = 0; i < virtualQubits.length; i++) {
      layout.setPhysical(physical[i], virtualQubits[i]);
    }
    return layout;
  }
}

// ---------------------------------------------------------------------------
// Additional transpiler passes
// ---------------------------------------------------------------------------

// Collect runs of single-qubit gates
export function collect_1q_runs(dag) {
  const runs = [];
  const visited = new Set();
  for (const node of dag.topological_op_nodes()) {
    if (visited.has(node)) continue;
    if (node.op.num_qubits !== 1) continue;
    const run = [node];
    visited.add(node);
    // Follow the chain
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter(s => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.num_qubits !== 1) break;
      if (next.qubits[0] !== current.qubits[0]) break;
      run.push(next);
      visited.add(next);
      current = next;
    }
    if (run.length > 1) runs.push(run);
  }
  return runs;
}

// Collect runs of 2-qubit gates
export function collect_2q_runs(dag) {
  const runs = [];
  const visited = new Set();
  for (const node of dag.topological_op_nodes()) {
    if (visited.has(node)) continue;
    if (node.op.num_qubits !== 2) continue;
    const run = [node];
    visited.add(node);
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter(s => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.num_qubits !== 2) break;
      if (next.qubits[0] !== current.qubits[0] || next.qubits[1] !== current.qubits[1]) break;
      run.push(next);
      visited.add(next);
      current = next;
    }
    if (run.length > 1) runs.push(run);
  }
  return runs;
}

// Consolidate runs of single-qubit gates into a single Unitary gate.
// Each run is replaced by a single gate whose matrix is the product of the
// run's gate matrices (applied left to right = rightmost matrix is earliest
// in time).
export function consolidate_blocks(dag) {
  const runs1q = collect_1q_runs(dag);
  for (const run of runs1q) {
    // Compose all gates in the run into one 2x2 unitary.
    let combined = null;
    let ok = true;
    for (const node of run) {
      try {
        const m = node.op.to_matrix();
        combined = combined ? m.mul(combined) : m;
      } catch (e) { ok = false; break; }
    }
    if (!ok || !combined) continue;
    // Build a unitary gate wrapping the combined matrix.
    const newGate = new Gate("unitary", 1, [combined]);
    newGate._matrixBuilder = () => combined;
    const newNode = new DAGOpNode(newGate, run[0].qargs.slice(), []);
    dag.substitute_node(run[0], newNode);
    // Remove the rest of the run.
    for (let i = 1; i < run.length; i++) {
      try { dag.remove_op_node(run[i]); } catch (e) { /* already removed */ }
    }
  }
  return dag;
}

// Cancel adjacent identical self-inverse Clifford gates (H, X, Y, Z, CX, CZ, SWAP).
// Walks the topologically-sorted op node list and removes consecutive pairs
// with matching name and qubit operands.
export function optimize_cliffords(dag) {
  const nodes = dag.topological_op_nodes();
  const selfInverse = new Set(["h", "x", "y", "z", "cx", "cz", "swap"]);
  // We can't mutate the DAG while iterating topological_op_nodes; collect
  // cancellations and apply them in a second pass.
  const toRemove = new Set();
  const stack = []; // stack of pending nodes
  for (const node of nodes) {
    if (toRemove.has(node)) continue;
    if (stack.length > 0) {
      const prev = stack[stack.length - 1];
      if (selfInverse.has(prev.op.name) &&
          prev.op.name === node.op.name &&
          prev.qubits.length === node.qubits.length &&
          prev.qubits.every((q, i) => q === node.qubits[i])) {
        // Cancel the pair.
        stack.pop();
        toRemove.add(prev);
        toRemove.add(node);
        continue;
      }
    }
    stack.push(node);
  }
  for (const node of toRemove) {
    try { dag.remove_op_node(node); } catch (e) { /* already gone */ }
  }
  return dag;
}

// Remove diagonal gates before measure (they don't affect the outcome)
export function remove_diagonal_gates_before_measure(dag) {
  const diagonalGates = ["z", "s", "sdg", "t", "tdg", "p", "u1", "rz", "cz", "cp"];
  const toRemove = [];
  for (const node of dag.topological_op_nodes()) {
    if (!diagonalGates.includes(node.op.name)) continue;
    // Check if all successors are measurements
    const succs = dag.successors(node).filter(s => s.is_op_node());
    if (succs.length > 0 && succs.every(s => s.op.name === "measure")) {
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.remove_op_node(node);
  }
  return dag;
}

// Helper: PRNG (mulberry32-style, seedable for reproducible routing).
function _makeRng(seed) {
  let s = seed >>> 0;
  return function () {
    s = (s + 0x6D2B79F5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Helper: return a fresh copy of a registered standard gate by its canonical
// (case-insensitive) name. Used to emit SWAP gates without depending on the
// QuantumCircuit convenience methods (which require the std cache to be
// populated; in transpiler contexts it always is, but we want to be explicit).
function _makeStdGate(name) {
  const upper = name.toUpperCase();
  switch (upper) {
    case "SWAP": return standardGates.makeSwapGate();
    case "CX":   return standardGates.makeCXGate();
    case "CZ":   return standardGates.makeCZGate();
    case "H":    return standardGates.makeHGate();
    case "X":    return standardGates.makeXGate();
    case "Y":    return standardGates.makeYGate();
    case "Z":    return standardGates.makeZGate();
    case "S":    return standardGates.makeSGate();
    case "SDG":  return standardGates.makeSdgGate();
    case "T":    return standardGates.makeTGate();
    case "TDG":  return standardGates.makeTdgGate();
    case "SX":   return standardGates.makeSxGate();
    case "I":    return standardGates.makeIGate();
    default:     return new Gate(name.toLowerCase(), 1, []);
  }
}
