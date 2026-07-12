/**
 * passes.js - Transpiler passes for layout, routing, and optimization.
 *
 * Layout passes: TrivialLayout, DenseLayout, SabreLayout, ApplyLayout.
 * Routing passes: BasicSwap, LookaheadSwap, StochasticSwap, SabreSwap.
 *   Note: the real Sabre implementation lives in ./sabre.js. The SabreSwap
 *   exported here delegates to it so callers importing from this module get
 *   the real algorithm, not a stub.
 * Optimization passes: CommutativeCancellation, CXCancellation, Optimize1qGates,
 *   OptimizeSwapBeforeMeasure, RemoveBarriers, RemoveResetInZeroState,
 *   DAGFixedPointPass.
 */

import { Layout, CouplingMap, AnalysisPass, TransformationPass } from "./layout.js";
import { DAGCircuit, DAGOpNode } from "./../dagcircuit/dagcircuit.js";
import { QuantumRegister, ClassicalRegister } from "./../core/bit.js";
import { Gate, ControlledGate, Instruction } from "./../core/gate.js";
import * as standardGates from "./../library/standard_gates.js";
import * as generalizedGates from "./../library/generalized_gates.js";
import {
  SabreSwap as RealSabreSwap,
  SabreLayout as RealSabreLayout,
} from "./sabre.js";

// ---------------------------------------------------------------------------
// Layout passes
// ---------------------------------------------------------------------------
export class TrivialLayout extends TransformationPass {
  constructor(couplingMap = null) {
    super();
    this.coupling_map = couplingMap;
  }

  run(dag) {
    // Map virtual qubit i -> physical qubit i (identity layout).
    const layout = Layout.trivial(dag.qubits.length);
    dag._layout = layout;
    return dag;
  }
}

export class DenseLayout extends TransformationPass {
  constructor(couplingMap = null) {
    super();
    this.coupling_map = couplingMap;
  }

  run(dag) {
    if (!this.coupling_map) {
      return new TrivialLayout().run(dag);
    }
    // Dense layout: pick the `n` most-connected physical qubits from the
    // coupling map and assign them to the circuit's virtual qubits. This is
    // a greedy heuristic — qiskit's version also considers error rates.
    const n = dag.qubits.length;
    const degree = new Map();
    for (let p = 0; p < this.coupling_map.size; p++) {
      degree.set(p, this.coupling_map.neighbors(p).length);
    }
    const sorted = Array.from(degree.entries())
      .sort((a, b) => b[1] - a[1])
      .slice(0, n)
      .map(([p]) => p);
    const layout = new Layout();
    for (let i = 0; i < n; i++) {
      layout.setPhysical(sorted[i], i);
    }
    dag._layout = layout;
    return dag;
  }
}

export class SabreLayout extends TransformationPass {
  constructor(couplingMap = null, seed = null, maxIterations = 3) {
    super();
    this.coupling_map = couplingMap;
    this.seed = seed;
    this.max_iterations = maxIterations;
  }

  run(dag) {
    if (!this.coupling_map) {
      return new TrivialLayout().run(dag);
    }
    const real = new RealSabreLayout(this.coupling_map, this.seed, this.max_iterations);
    dag._layout = real.run(dag);
    return dag;
  }
}

export class ApplyLayout extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    if (!dag._layout) return dag;
    // Apply the layout by expanding to the full backend size.
    const layout = dag._layout;
    const maxPhys = Math.max(...Array.from(layout._p2v.keys())) + 1;
    if (maxPhys === dag.qubits.length) return dag;
    // Expand the DAG to add ancilla qubits for unused physical positions.
    const newQubitsNeeded = maxPhys - dag.qubits.length;
    if (newQubitsNeeded > 0) {
      const ancReg = new QuantumRegister(newQubitsNeeded, "ancilla");
      dag.add_qreg(ancReg);
    }
    return dag;
  }
}

// ---------------------------------------------------------------------------
// Routing passes
// ---------------------------------------------------------------------------
export class BasicSwap extends TransformationPass {
  constructor(couplingMap = null) {
    super();
    this.coupling_map = couplingMap;
  }

  run(dag) {
    if (!this.coupling_map) return dag;
    // Greedy SWAP insertion: for each 2-qubit gate whose qubits are not
    // adjacent in the coupling map, insert SWAP gates along the shortest
    // path so the gate becomes adjacent, then insert the reverse SWAPs
    // after the gate to restore the layout. This is the standard
    // "BasicSwap" algorithm.
    //
    // We rebuild the DAG from scratch with the SWAPs inserted (rather than
    // trying to splice into the existing DAG), which is simpler and
    // avoids relying on internal DAG mutation methods.
    const newDag = new DAGCircuit();
    for (const [name, reg] of dag.qregs.entries()) {
      newDag.add_qreg(new QuantumRegister(reg.bits.length, name));
    }
    for (const [name, reg] of dag.cregs.entries()) {
      newDag.add_creg(new ClassicalRegister(reg.bits.length, name));
    }
    newDag.name = dag.name;
    newDag.global_phase = dag.global_phase;
    newDag.metadata = dag.metadata;

    for (const node of dag.topological_op_nodes()) {
      const op = node.op;
      const qargs = node.qargs;
      const cargs = node.cargs;
      if (op.num_qubits !== 2) {
        newDag.apply_operation(op.copy(), qargs.map(q => _mapQubit(q, dag, newDag)), cargs.map(c => _mapClbit(c, dag, newDag)));
        continue;
      }
      const q0 = dag.qubits.indexOf(qargs[0]);
      const q1 = dag.qubits.indexOf(qargs[1]);
      if (this.coupling_map.hasEdge(q0, q1)) {
        newDag.apply_operation(op.copy(), qargs.map(q => _mapQubit(q, dag, newDag)), cargs.map(c => _mapClbit(c, dag, newDag)));
        continue;
      }
      const path = this.coupling_map.shortestPath(q0, q1);
      if (!path || path.length < 2) {
        throw new Error(`BasicSwap: no path between ${q0} and ${q1} in coupling map`);
      }
      // Insert forward SWAPs along the path. After these SWAPs:
      //   - the qubit originally at path[0] (the gate's first operand) is
      //     now at path[last],
      //   - the qubit originally at path[last] (the gate's second operand,
      //     IF path[last] === q1) is now at path[last-1].
      // So the gate should be applied to (path[last-1], path[last]) — but
      // note the order is reversed from the original (q0, q1) because the
      // qubits have been physically swapped. For symmetric gates (SWAP,
      // CZ) this doesn't matter; for asymmetric gates (CX, CY, CRX, …)
      // the control/target ordering is preserved relative to the qubits'
      // new positions.
      const forwardSwaps = [];
      for (let i = 0; i + 1 < path.length; i++) {
        const a = path[i];
        const b = path[i + 1];
        forwardSwaps.push([a, b]);
        const swapOp = standardGates.SwapGate.copy();
        const qa = newDag.qubits[a];
        const qb = newDag.qubits[b];
        newDag.apply_operation(swapOp, [qa, qb], []);
      }
      // Apply the original gate on the qubits at the end of the path
      // (which are now adjacent). The operand ordering matches the
      // original gate's (control, target) — the SWAP chain moves the
      // first operand to path[last] and the second to path[last-1].
      const lastIdx = path.length - 1;
      const newQargs = [
        newDag.qubits[path[lastIdx]],
        newDag.qubits[path[lastIdx - 1]],
      ];
      newDag.apply_operation(op.copy(), newQargs, cargs.map(c => _mapClbit(c, dag, newDag)));
      // Insert reverse SWAPs to restore the layout.
      for (let i = forwardSwaps.length - 1; i >= 0; i--) {
        const [a, b] = forwardSwaps[i];
        const swapOp = standardGates.SwapGate.copy();
        const qa = newDag.qubits[a];
        const qb = newDag.qubits[b];
        newDag.apply_operation(swapOp, [qa, qb], []);
      }
    }
    return newDag;
  }
}

// Map a qubit from the old DAG to the corresponding qubit in the new DAG
// (by index, since both DAGs have the same register structure).
function _mapQubit(q, oldDag, newDag) {
  const idx = oldDag.qubits.indexOf(q);
  return newDag.qubits[idx];
}
function _mapClbit(c, oldDag, newDag) {
  const idx = oldDag.clbits.indexOf(c);
  return newDag.clbits[idx];
}

export class LookaheadSwap extends TransformationPass {
  constructor(couplingMap = null, searchDepth = 5) {
    super();
    this.coupling_map = couplingMap;
    this.search_depth = searchDepth;
  }

  run(dag) {
    if (!this.coupling_map) return dag;
    // LookaheadSwap uses the real Sabre router with the "lookahead" heuristic,
    // which considers the next layer of gates when choosing SWAPs.
    const sabre = new RealSabreSwap(this.coupling_map, "lookahead", null);
    return sabre.run(dag);
  }
}

export class StochasticSwap extends TransformationPass {
  constructor(couplingMap = null, seed = null, trials = 20) {
    super();
    this.coupling_map = couplingMap;
    this.seed = seed;
    this.trials = trials;
  }

  run(dag) {
    if (!this.coupling_map) return dag;
    // StochasticSwap runs Sabre multiple times with different seeds and picks
    // the routing with the fewest SWAPs.
    let best = null;
    let bestSwaps = Infinity;
    for (let t = 0; t < this.trials; t++) {
      const sabre = new RealSabreSwap(this.coupling_map, "lookahead", this.seed != null ? this.seed + t : null);
      const routed = sabre.run(dag);
      const numSwaps = routed.data.filter(ci => ci.operation.name === "swap").length;
      if (numSwaps < bestSwaps) {
        bestSwaps = numSwaps;
        best = routed;
      }
    }
    return best || dag;
  }
}

export class SabreSwap extends TransformationPass {
  constructor(couplingMap = null, heuristic = "lookahead", seed = null) {
    super();
    this.coupling_map = couplingMap;
    this.heuristic = heuristic;
    this.seed = seed;
  }

  run(dag) {
    if (!this.coupling_map) return dag;
    // Delegate to the real Sabre router from sabre.js.
    const sabre = new RealSabreSwap(this.coupling_map, this.heuristic, this.seed);
    return sabre.run(dag);
  }
}

// ---------------------------------------------------------------------------
// Optimization passes
// ---------------------------------------------------------------------------
export class CXCancellation extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    // Cancel adjacent identical CNOTs: CNOT(a,b) CNOT(a,b) -> identity.
    const opNodes = dag.topological_op_nodes();
    const toRemove = new Set();
    for (const node of opNodes) {
      if (toRemove.has(node)) continue;
      if (node.op.name !== "cx") continue;
      for (const succ of dag.successors(node)) {
        if (succ.is_op_node() && succ.op.name === "cx" &&
            succ.qargs.length === 2 &&
            succ.qargs[0] === node.qargs[0] &&
            succ.qargs[1] === node.qargs[1]) {
          toRemove.add(node);
          toRemove.add(succ);
          break;
        }
      }
    }
    for (const node of toRemove) {
      try { dag.remove_op_node(node); } catch (e) { /* already removed */ }
    }
    return dag;
  }
}

export class CommutativeCancellation extends TransformationPass {
  constructor(basisGates = ["cx", "rz", "sx", "x"]) {
    super();
    this.basis_gates = basisGates;
  }

  run(dag) {
    // Cancel adjacent commuting gates on the same qubits when they are
    // mutual inverses (e.g. RZ(a) RZ(-a) -> identity). The current
    // implementation handles the most common cases: consecutive RZ gates
    // on a single qubit, and consecutive single-qubit gates with the same
    // name and opposite parameters. The full CommutativeCancellation pass
    // would also handle cross-qubit commuting pairs (e.g. CNOTs on
    // disjoint qubit sets), which requires a more expensive commutativity
    // analysis.
    const opNodes = dag.topological_op_nodes();
    const toRemove = new Set();
    const stack = [];
    for (const node of opNodes) {
      if (toRemove.has(node)) continue;
      if (stack.length > 0) {
        const prev = stack[stack.length - 1];
        // Cancel RZ(a) followed by RZ(-a) on the same qubit.
        if (prev.op.name === "rz" && node.op.name === "rz" &&
            prev.qargs[0] === node.qargs[0] &&
            prev.op.params.length > 0 && node.op.params.length > 0) {
          const a = typeof prev.op.params[0] === "number" ? prev.op.params[0] : 0;
          const b = typeof node.op.params[0] === "number" ? node.op.params[0] : 0;
          if (Math.abs(a + b) < 1e-12) {
            stack.pop();
            toRemove.add(prev);
            toRemove.add(node);
            continue;
          }
        }
      }
      stack.push(node);
    }
    for (const node of toRemove) {
      try { dag.remove_op_node(node); } catch (e) { /* already removed */ }
    }
    return dag;
  }
}

export class Optimize1qGates extends TransformationPass {
  constructor(basisGates = ["u1", "u2", "u3", "cx"]) {
    super();
    this.basis_gates = basisGates;
  }

  run(dag) {
    // Merge consecutive single-qubit gates on the same qubit into a single
    // U3 gate via Euler Z-Y-Z decomposition.
    const opNodes = dag.topological_op_nodes();
    for (const node of opNodes) {
      if (node.op.num_qubits !== 1) continue;
      const qarg = node.qargs[0];
      let current = node;
      const chain = [node];
      while (true) {
        const succs = dag.successors(current).filter(s => s.is_op_node());
        if (succs.length !== 1) break;
        const next = succs[0];
        if (next.op.num_qubits !== 1 || next.qargs[0] !== qarg) break;
        chain.push(next);
        current = next;
      }
      if (chain.length > 1) {
        // Compose all gates' matrices into one 2x2 unitary.
        let combined = null;
        let ok = true;
        for (const n of chain) {
          try {
            const mat = n.op.to_matrix();
            combined = combined ? mat.mul(combined) : mat;
          } catch (e) { ok = false; break; }
        }
        if (!ok || !combined) continue;
        // Decompose combined into U3 Euler angles (Z-Y-Z convention).
        const euler = _decomposeZYZ(combined);
        const newGate = generalizedGates.makeU3Gate(euler.theta, euler.phi, euler.lambda);
        const newNode = new DAGOpNode(newGate, [qarg], []);
        dag.substitute_node(chain[0], newNode);
        // Remove the rest of the chain.
        for (let i = 1; i < chain.length; i++) {
          try { dag.remove_op_node(chain[i]); } catch (e) {}
        }
      }
    }
    return dag;
  }
}

// Z-Y-Z Euler decomposition of a 2x2 unitary.
// U = e^{i alpha} RZ(phi) RY(theta) RZ(lambda).
// We drop the global phase alpha (irrelevant for quantum states).
function _decomposeZYZ(m) {
  const a = m.get(0, 0);
  const b = m.get(0, 1);
  const c = m.get(1, 0);
  const d = m.get(1, 1);
  // theta = 2 * atan2(|b|, |a|)
  const theta = 2 * Math.atan2(b.abs(), a.abs());
  // lambda = arg(a) - arg(c_conj) ... use the standard formula:
  //   lambda = arg(a) + arg(b)  (up to sign; we verify via round-trip in tests)
  //   phi    = arg(b) - arg(a) + lambda  ... etc.
  // The robust derivation (see Nielsen & Chuang exercise 4.11):
  //   if sin(theta/2) != 0 and cos(theta/2) != 0:
  //     lambda = arg(a) - arg(c) + pi   (mod 2pi)  when using RY(theta) = exp(-i theta Y/2)
  //     phi    = arg(b) + arg(c) - pi   (mod 2pi)
  // We use the symmetric form that matches makeU3Gate's matrix definition:
  //   U3(theta, phi, lambda) = [[cos(theta/2), -e^{i lambda} sin(theta/2)],
  //                              [e^{i phi} sin(theta/2), e^{i (phi+lambda)} cos(theta/2)]]
  // Matching: a = cos(theta/2), b = -e^{i lambda} sin(theta/2),
  //           c = e^{i phi} sin(theta/2), d = e^{i (phi+lambda)} cos(theta/2)
  const halfTheta = theta / 2;
  let phi, lambda;
  if (Math.abs(Math.sin(halfTheta)) < 1e-12) {
    // theta ~ 0: U is diagonal, phi is irrelevant.
    phi = 0;
    lambda = a.arg() - d.arg() >= 0 ? 0 : 0; // both phases equal up to global
    lambda = 2 * a.arg();
  } else if (Math.abs(Math.cos(halfTheta)) < 1e-12) {
    // theta ~ pi: U is off-diagonal.
    phi = c.arg();
    lambda = b.arg() + Math.PI;
  } else {
    phi = c.arg();
    lambda = b.arg() + Math.PI;
  }
  return { theta, phi, lambda };
}

export class OptimizeSwapBeforeMeasure extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    // Remove SWAPs that immediately precede measurements by relabeling which
    // qubit is measured. We find SWAP nodes whose only successors are
    // measurements, remove the SWAP, and swap the measurement target qubits.
    const opNodes = dag.topological_op_nodes();
    for (const node of opNodes) {
      if (node.op.name !== "swap") continue;
      const succs = dag.successors(node).filter(s => s.is_op_node());
      if (succs.length === 2 && succs.every(s => s.op.name === "measure")) {
        // Swap the classical targets of the two measurements.
        const m0 = succs[0], m1 = succs[1];
        const tmp = m0.cargs[0];
        m0.cargs[0] = m1.cargs[0];
        m1.cargs[0] = tmp;
        try { dag.remove_op_node(node); } catch (e) {}
      }
    }
    return dag;
  }
}

export class RemoveBarriers extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    const opNodes = dag.op_nodes();
    for (const node of opNodes) {
      if (node.op.name === "barrier") {
        dag.remove_op_node(node);
      }
    }
    return dag;
  }
}

export class RemoveResetInZeroState extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    // Remove reset gates on qubits that are still in their initial |0> state
    // (i.e. no prior gate has touched them).
    const opNodes = dag.topological_op_nodes();
    const zeroQubits = new Set();
    for (const q of dag.qubits) zeroQubits.add(q);
    for (const node of opNodes) {
      if (node.op.name === "reset") {
        if (zeroQubits.has(node.qargs[0])) {
          dag.remove_op_node(node);
        }
      } else {
        for (const q of node.qargs) zeroQubits.delete(q);
      }
    }
    return dag;
  }
}

export class DAGFixedPointPass extends AnalysisPass {
  constructor() { super(); }

  run(dag) {
    // No-op marker pass; used as a wrapper for iterative pass application.
    return dag;
  }
}
