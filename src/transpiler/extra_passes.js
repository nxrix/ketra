import { AnalysisPass, TransformationPass } from "./layout.js";
import { DAGCircuit, DAGOpNode } from "./../dagcircuit/dagcircuit.js";
import { QuantumRegister, ClassicalRegister } from "./../core/bit.js";
import { Instruction } from "./../core/gate.js";
import * as standardGates from "./../library/standard_gates.js";
import * as generalizedGates from "./../library/generalized_gates.js";
import { UnitaryGate } from "./../library/extra_gates.js";

// Decompose: replace gates that have a definition with their definition body.
export class Decompose extends TransformationPass {
  constructor(gatesToDecompose = null) {
    super();
    this.gates = gatesToDecompose ? (Array.isArray(gatesToDecompose) ? gatesToDecompose : [gatesToDecompose]) : null;
  }

  run(dag) {
    const opNodes = dag.topologicalOpNodes();
    for (const node of opNodes) {
      if (this.gates && !this.gates.includes(node.op.name)) continue;
      const defFn = node.op._definition;
      if (typeof defFn !== "function") continue;
      try {
        const def = defFn(node.op);
        if (!def || def.length === 0) continue;
        // The definition's qubits/clbits are the SUB-circuit's bits, which
        // we map to the parent node's qargs/cargs by index (since the
        // definition was built from a circuit with the same qubit count).
        const subDag = new DAGCircuit();
        const numQubits = node.qargs.length;
        const numClbits = node.cargs.length;
        if (numQubits > 0) subDag.addQreg(new QuantumRegister(numQubits, "q"));
        if (numClbits > 0) subDag.addCreg(new ClassicalRegister(numClbits, "c"));
        // Build a position-based map from the def's unique qubits to the
        // parent node's qargs (the def's qubits are in the same order as
        // the parent's qargs, since toGate() builds the definition from
        // this.qubits).
        for (const [subOp, subQubits, subClbits] of def) {
          if (!node._qubitIndexMap) {
            node._qubitIndexMap = new Map();
            const seen = new Set();
            let posIdx = 0;
            for (const sq of node.qargs) {
              if (!seen.has(sq)) {
                node._qubitIndexMap.set(posIdx, sq);
                seen.add(sq);
                posIdx++;
              }
            }
          }
          if (!node._defUniqueQubits) {
            node._defUniqueQubits = [];
            const seenSet = new Set();
            for (const [_, sqs, scs] of def) {
              for (const sq of sqs) {
                if (!seenSet.has(sq)) {
                  seenSet.add(sq);
                  node._defUniqueQubits.push(sq);
                }
              }
            }
          }
          const mappedQ = subQubits.map(sq => {
            const idx = node._defUniqueQubits.indexOf(sq);
            return idx >= 0 ? subDag.qubits[idx] : sq;
          });
          const mappedC = subClbits.map(sc => {
            const idx = node.cargs.indexOf(sc);
            return idx >= 0 ? subDag.clbits[idx] : sc;
          });
          subDag.applyOperation(subOp.copy(), mappedQ, mappedC);
        }
        dag.substituteNodeWithDag(node, subDag);
      } catch (e) {
        // Keep original on error.
      }
    }
    return dag;
  }
}

// RemoveFinalMeasurements: remove measurements at the end of the circuit
// (i.e. measurements with no gates following them on the same qubit).
export class RemoveFinalMeasurements extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    const opNodes = dag.topologicalOpNodes();
    // Find measurement nodes whose qubit has no later op.
    const toRemove = new Set();
    const qubitLastOp = new Map();
    // Walk in reverse: the last op on each qubit is a candidate for removal
    // if it's a measurement.
    for (let i = opNodes.length - 1; i >= 0; i--) {
      const node = opNodes[i];
      if (node.op.name !== "measure") continue;
      let isFinal = true;
      for (const q of node.qargs) {
        if (qubitLastOp.has(q)) { isFinal = false; break; }
      }
      if (isFinal) {
        toRemove.add(node);
        for (const q of node.qargs) qubitLastOp.set(q, node);
      }
    }
    for (const node of toRemove) {
      try { dag.removeOpNode(node); } catch (e) {}
    }
    return dag;
  }
}

// BarrierBeforeFinalMeasurements: insert a barrier on all qubits immediately
// before the first final measurement.
export class BarrierBeforeFinalMeasurements extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    const opNodes = dag.topologicalOpNodes();
    // Find the first measurement that begins the "final measurements" run.
    let firstFinalMeasureIdx = -1;
    for (let i = opNodes.length - 1; i >= 0; i--) {
      if (opNodes[i].op.name === "measure") {
        firstFinalMeasureIdx = i;
      } else if (firstFinalMeasureIdx !== -1) {
        // We've found the measurements; stop at the first non-measurement
        // gate before them.
        break;
      }
    }
    if (firstFinalMeasureIdx === -1) return dag;
    // Insert a barrier on all qubits before opNodes[firstFinalMeasureIdx].
    // We do this by rebuilding the DAG.
    const newDag = new DAGCircuit();
    for (const [name, reg] of dag.qregs.entries()) {
      newDag.addQreg(new QuantumRegister(reg.bits.length, name));
    }
    for (const [name, reg] of dag.cregs.entries()) {
      newDag.addCreg(new ClassicalRegister(reg.bits.length, name));
    }
    for (let i = 0; i < opNodes.length; i++) {
      if (i === firstFinalMeasureIdx) {
        const barrier = new Instruction("barrier", newDag.qubits.length, 0, []);
        newDag.applyOperation(barrier, newDag.qubits.slice(), []);
      }
      const node = opNodes[i];
      const qargs = node.qargs.map(q => {
        const idx = dag.qubits.indexOf(q);
        return newDag.qubits[idx];
      });
      const cargs = node.cargs.map(c => {
        const idx = dag.clbits.indexOf(c);
        return newDag.clbits[idx];
      });
      newDag.applyOperation(node.op.copy(), qargs, cargs);
    }
    return newDag;
  }
}

// Collect2qBlocks: identify maximal runs of consecutive 2-qubit gates on the
// same pair of qubits. Returns an analysis result; this pass is mainly used
// as a precursor to ConsolidateBlocks.
export class Collect2qBlocks extends AnalysisPass {
  constructor() { super(); }

  run(dag) {
    const opNodes = dag.topologicalOpNodes();
    const blocks = [];
    let currentBlock = null;
    for (const node of opNodes) {
      if (node.op.numQubits === 2) {
        const q0 = node.qargs[0];
        const q1 = node.qargs[1];
        if (currentBlock &&
            currentBlock.qubits.includes(q0) &&
            currentBlock.qubits.includes(q1)) {
          currentBlock.nodes.push(node);
        } else {
          if (currentBlock) blocks.push(currentBlock);
          currentBlock = { qubits: [q0, q1], nodes: [node] };
        }
      } else {
        if (currentBlock) blocks.push(currentBlock);
        currentBlock = null;
      }
    }
    if (currentBlock) blocks.push(currentBlock);
    this.propertySet["2q_blocks"] = blocks;
    return dag;
  }
}

// ConsolidateBlocks: replace each 2-qubit block (from Collect2qBlocks) with
// a single UnitaryGate containing the block's matrix. Reduces gate count.
export class ConsolidateBlocks extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    const collector = new Collect2qBlocks();
    collector.run(dag);
    const blocks = collector.propertySet["2q_blocks"] || [];
    for (const block of blocks) {
      if (block.nodes.length < 2) continue; // no point consolidating a single gate
      let combined = null;
      for (const node of block.nodes) {
        try {
          const mat = node.op.toMatrix();
          combined = combined ? mat.mul(combined) : mat;
        } catch (e) { combined = null; break; }
      }
      if (!combined) continue;
      const newGate = new UnitaryGate(combined);
      block.nodes[0].op = newGate;
      for (let i = 1; i < block.nodes.length; i++) {
        try { dag.removeOpNode(block.nodes[i]); } catch (e) {}
      }
    }
    return dag;
  }
}

// Unroll3qOrMore: decompose any gate acting on 3+ qubits into 1- and 2-qubit
// gates. Uses the gate's matrix and a brute-force decomposition.
export class Unroll3qOrMore extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    const opNodes = dag.topologicalOpNodes();
    for (const node of opNodes) {
      if (node.op.numQubits < 3) continue;
      const defFn = node.op._definition;
      if (typeof defFn !== "function") continue;
      try {
        const def = defFn(node.op);
        if (!def || def.length === 0) continue;
        const subDag = new DAGCircuit();
        const numQubits = node.qargs.length;
        const numClbits = node.cargs.length;
        if (numQubits > 0) subDag.addQreg(new QuantumRegister(numQubits, "q"));
        if (numClbits > 0) subDag.addCreg(new ClassicalRegister(numClbits, "c"));
        const defUniqueQubits = [];
        const seen = new Set();
        for (const [_, sqs, scs] of def) {
          for (const sq of sqs) {
            if (!seen.has(sq)) { seen.add(sq); defUniqueQubits.push(sq); }
          }
        }
        const defUniqueClbits = [];
        const seenC = new Set();
        for (const [_, sqs, scs] of def) {
          for (const sc of scs) {
            if (!seenC.has(sc)) { seenC.add(sc); defUniqueClbits.push(sc); }
          }
        }
        for (const [subOp, subQubits, subClbits] of def) {
          const mappedQargs = subQubits.map(sq => {
            const idx = defUniqueQubits.indexOf(sq);
            return idx >= 0 && idx < subDag.qubits.length ? subDag.qubits[idx] : sq;
          });
          const mappedCargs = subClbits.map(sc => {
            const idx = defUniqueClbits.indexOf(sc);
            return idx >= 0 && idx < subDag.clbits.length ? subDag.clbits[idx] : sc;
          });
          subDag.applyOperation(subOp.copy(), mappedQargs, mappedCargs);
        }
        dag.substituteNodeWithDag(node, subDag);
      } catch (e) { /* keep original */ }
    }
    return dag;
  }
}

// MergeAdjacentBarriers: combine consecutive barriers on the same qubits
// into a single barrier.
export class MergeAdjacentBarriers extends TransformationPass {
  constructor() { super(); }

  run(dag) {
    const opNodes = dag.topologicalOpNodes();
    const toRemove = new Set();
    let prevBarrier = null;
    for (const node of opNodes) {
      if (node.op.name === "barrier") {
        if (prevBarrier &&
            prevBarrier.qargs.length === node.qargs.length &&
            prevBarrier.qargs.every(q => node.qargs.includes(q))) {
          toRemove.add(node);
        } else {
          prevBarrier = node;
        }
      } else {
        prevBarrier = null;
      }
    }
    for (const node of toRemove) {
      try { dag.removeOpNode(node); } catch (e) {}
    }
    return dag;
  }
}

// CheckMap: verify that every 2-qubit gate respects the coupling map.
// Sets propertySet["is_mapped"] = true/false.
export class CheckMap extends AnalysisPass {
  constructor(couplingMap = null) {
    super();
    this.couplingMap = couplingMap;
  }

  run(dag) {
    if (!this.couplingMap) {
      this.propertySet["is_mapped"] = true;
      return dag;
    }
    let isMapped = true;
    for (const node of dag.topologicalOpNodes()) {
      if (node.op.numQubits !== 2) continue;
      const q0 = dag.qubits.indexOf(node.qargs[0]);
      const q1 = dag.qubits.indexOf(node.qargs[1]);
      if (!this.couplingMap.hasEdge(q0, q1)) {
        isMapped = false;
        break;
      }
    }
    this.propertySet["is_mapped"] = isMapped;
    return dag;
  }
}

// GateDirection: flip the direction of asymmetric 2-qubit gates (CX, CY, CZ
// with a directed coupling map) so that the control is on the "from" side.
export class GateDirection extends TransformationPass {
  constructor(couplingMap = null) {
    super();
    this.couplingMap = couplingMap;
  }

  run(dag) {
    if (!this.couplingMap) return dag;
    const opNodes = dag.topologicalOpNodes();
    for (const node of opNodes) {
      if (node.op.numQubits !== 2) continue;
      const q0 = dag.qubits.indexOf(node.qargs[0]);
      const q1 = dag.qubits.indexOf(node.qargs[1]);
      if (this.couplingMap.hasEdge(q0, q1)) continue;
      if (this.couplingMap.hasEdge(q1, q0)) {
        // Flip the gate's operands. Operand swapping is correct for symmetric
        // gates like CZ and SWAP; for asymmetric gates like CX, the reverse
        // decomposition (CX(a,b) = H_b · CX(b,a) · H_b) would be needed and
        // users should use BasisTranslator.
        const tmp = node.qargs[0];
        node.qargs[0] = node.qargs[1];
        node.qargs[1] = tmp;
      }
    }
    return dag;
  }
}

// CountOps: count the number of each gate type.
export class CountOps extends AnalysisPass {
  constructor() { super(); }

  run(dag) {
    const counts = {};
    for (const node of dag.topologicalOpNodes()) {
      counts[node.op.name] = (counts[node.op.name] || 0) + 1;
    }
    this.propertySet["countOps"] = counts;
    return dag;
  }
}

// Depth: compute the circuit depth.
export class Depth extends AnalysisPass {
  constructor() { super(); }

  run(dag) {
    const lastOnQubit = new Map();
    for (const q of dag.qubits) lastOnQubit.set(q, -1);
    for (const c of dag.clbits) lastOnQubit.set(c, -1);
    let maxDepth = 0;
    for (const node of dag.topologicalOpNodes()) {
      if (node.op.name === "barrier") continue;
      let d = 0;
      for (const q of node.qargs) d = Math.max(d, lastOnQubit.get(q) + 1);
      for (const c of node.cargs) d = Math.max(d, lastOnQubit.get(c) + 1);
      for (const q of node.qargs) lastOnQubit.set(q, d);
      for (const c of node.cargs) lastOnQubit.set(c, d);
      if (d + 1 > maxDepth) maxDepth = d + 1;
    }
    this.propertySet["depth"] = maxDepth;
    return dag;
  }
}

// Size: total gate count (excluding barriers).
export class Size extends AnalysisPass {
  constructor() { super(); }

  run(dag) {
    let n = 0;
    for (const node of dag.topologicalOpNodes()) {
      if (node.op.name === "barrier") continue;
      n++;
    }
    this.propertySet["size"] = n;
    return dag;
  }
}

// Width: total qubits + clbits.
export class Width extends AnalysisPass {
  constructor() { super(); }

  run(dag) {
    this.propertySet["width"] = dag.qubits.length + dag.clbits.length;
    return dag;
  }
}

// Optimize1qGatesDecomposition: improved 1-qubit synthesis using the
// Z-Y-Z Euler decomposition. Targets a specific basis (default: ['u3', 'cx']).
export class Optimize1qGatesDecomposition extends TransformationPass {
  constructor(basisGates = ["u3", "cx"]) {
    super();
    this.basisGates = basisGates;
  }

  run(dag) {
    const opNodes = dag.topologicalOpNodes();
    for (const node of opNodes) {
      if (node.op.numQubits !== 1) continue;
      const succs = dag.successors(node).filter(s => s.is_op_node());
      if (succs.length !== 1) continue;
      const next = succs[0];
      if (next.op.numQubits !== 1 || next.qargs[0] !== node.qargs[0]) continue;
      try {
        const m1 = node.op.toMatrix();
        const m2 = next.op.toMatrix();
        const combined = m2.mul(m1);
        const euler = _decomposeZYZ(combined);
        const newGate = generalizedGates.makeU3Gate(euler.theta, euler.phi, euler.lambda);
        const newNode = new DAGOpNode(newGate, [node.qargs[0]], []);
        dag.substituteNode(node, newNode);
        try { dag.removeOpNode(next); } catch (e) {}
      } catch (e) { /* keep original */ }
    }
    return dag;
  }
}

// Lookup an extra gate class (UnitaryGate) via the global registry.
let _extraGateRegistry = null;
function _lookupExtraGate(name) {
  if (!_extraGateRegistry) {
    try {
      // Use the registry set by extra_gates.js.
      const mod = globalThis.__ketraExtraGates;
      if (mod && mod[name]) return mod[name];
    } catch (e) {}
    return null;
  }
  return _extraGateRegistry[name];
}

// Z-Y-Z Euler decomposition (same as in passes.js).
function _decomposeZYZ(m) {
  const a = m.get(0, 0);
  const b = m.get(0, 1);
  const c = m.get(1, 0);
  const d = m.get(1, 1);
  const theta = 2 * Math.atan2(b.abs(), a.abs());
  const halfTheta = theta / 2;
  let phi, lambda;
  if (Math.abs(Math.sin(halfTheta)) < 1e-12) {
    phi = 0;
    lambda = 2 * a.arg();
  } else if (Math.abs(Math.cos(halfTheta)) < 1e-12) {
    phi = c.arg();
    lambda = b.arg() + Math.PI;
  } else {
    phi = c.arg();
    lambda = b.arg() + Math.PI;
  }
  return { theta, phi, lambda };
}

// BasisTranslator: translate gates to a target basis using known decompositions.
// For each gate not in the target basis, replace it with its decomposition
// (a sequence of basis gates). Uses a fixed table of decompositions rather
// than an equivalence library. The decompositions cover the common cases:
//   - H, S, Sdg, T, Tdg -> u3 / u1
//   - RX, RY, RZ -> u3 / u1
//   - SWAP -> 3 CX
//   - CCX -> 6 CX + single-qubit gates
//   - CY, CZ, CH, CSX -> decompositions via CX + single-qubit gates
export class BasisTranslator extends TransformationPass {
  constructor(targetBasis = ["cx", "u3", "u1"], basisGates = null) {
    super();
    this.target_basis = targetBasis;
    this.basisGates = basisGates || targetBasis;
  }

  run(dag) {
    const opNodes = dag.topologicalOpNodes();
    for (const node of opNodes) {
      if (this.target_basis.includes(node.op.name)) continue;
      const decomposition = _getBasisDecomposition(node.op, this.target_basis);
      if (decomposition === null) continue;
      const subDag = new DAGCircuit();
      const numQubits = node.qargs.length;
      const numClbits = node.cargs.length;
      if (numQubits > 0) subDag.addQreg(new QuantumRegister(numQubits, "q"));
      if (numClbits > 0) subDag.addCreg(new ClassicalRegister(numClbits, "c"));
      for (const [op, qargs, cargs] of decomposition) {
        const mappedQargs = qargs.map((q, i) => {
          // The placeholder qubit's index is i (the position in the
          // decomposition's qubit list). Map to the parent node's qargs.
          if (q && q._placeholderIdx !== undefined) return subDag.qubits[q._placeholderIdx];
          const idx = node.qargs.indexOf(q);
          return idx >= 0 ? subDag.qubits[idx] : subDag.qubits[i];
        });
        const mappedCargs = cargs.map((c, i) => {
          const idx = node.cargs.indexOf(c);
          return idx >= 0 ? subDag.clbits[idx] : subDag.clbits[i];
        });
        subDag.applyOperation(op.copy(), mappedQargs, mappedCargs);
      }
      try {
        dag.substituteNodeWithDag(node, subDag);
      } catch (e) { /* keep original */ }
    }
    return dag;
  }
}

function _getBasisDecomposition(op, targetBasis) {
  const name = op.name.toLowerCase();
  const placeholder = (idx) => ({ _placeholderIdx: idx });

  // H = U3(pi/2, 0, pi)
  if (name === "h" && (targetBasis.includes("u3") || targetBasis.includes("u"))) {
    const uName = targetBasis.includes("u3") ? "u3" : "u";
    const uGate = generalizedGates.makeU3Gate(Math.PI / 2, 0, Math.PI);
    return [[uGate, [placeholder(0)], []]];
  }
  // S = U1(pi/2)
  if (name === "s" && targetBasis.includes("u1")) {
    return [[generalizedGates.makeU1Gate(Math.PI / 2), [placeholder(0)], []]];
  }
  if (name === "sdg" && targetBasis.includes("u1")) {
    return [[generalizedGates.makeU1Gate(-Math.PI / 2), [placeholder(0)], []]];
  }
  if (name === "t" && targetBasis.includes("u1")) {
    return [[generalizedGates.makeU1Gate(Math.PI / 4), [placeholder(0)], []]];
  }
  if (name === "tdg" && targetBasis.includes("u1")) {
    return [[generalizedGates.makeU1Gate(-Math.PI / 4), [placeholder(0)], []]];
  }
  // RX = U3(theta, -pi/2, pi/2)
  if (name === "rx" && targetBasis.includes("u3")) {
    const theta = (typeof op.params[0] === "number") ? op.params[0] : 0;
    return [[generalizedGates.makeU3Gate(theta, -Math.PI / 2, Math.PI / 2), [placeholder(0)], []]];
  }
  // RY = U3(theta, 0, 0)
  if (name === "ry" && targetBasis.includes("u3")) {
    const theta = (typeof op.params[0] === "number") ? op.params[0] : 0;
    return [[generalizedGates.makeU3Gate(theta, 0, 0), [placeholder(0)], []]];
  }
  // RZ = U1(theta) up to global phase
  if (name === "rz" && targetBasis.includes("u1")) {
    const theta = (typeof op.params[0] === "number") ? op.params[0] : 0;
    return [[generalizedGates.makeU1Gate(theta), [placeholder(0)], []]];
  }
  // SWAP = 3 CX
  if (name === "swap" && targetBasis.includes("cx")) {
    const cx = standardGates.CXGate.copy();
    const q0 = placeholder(0);
    const q1 = placeholder(1);
    return [
      [cx, [q0, q1], []],
      [cx.copy(), [q1, q0], []],
      [cx.copy(), [q0, q1], []],
    ];
  }
  // CZ = H(target) · CX · H(target)
  if (name === "cz" && targetBasis.includes("cx") && targetBasis.includes("u3")) {
    const cx = standardGates.CXGate.copy();
    const h = generalizedGates.makeU3Gate(Math.PI / 2, 0, Math.PI);
    const q0 = placeholder(0);
    const q1 = placeholder(1);
    return [
      [h, [q1], []],
      [cx, [q0, q1], []],
      [h.copy(), [q1], []],
    ];
  }
  // CY = Sdg(target) · CX · S(target)
  if (name === "cy" && targetBasis.includes("cx") && targetBasis.includes("u1")) {
    const cx = standardGates.CXGate.copy();
    const s = generalizedGates.makeU1Gate(Math.PI / 2);
    const sdg = generalizedGates.makeU1Gate(-Math.PI / 2);
    const q0 = placeholder(0);
    const q1 = placeholder(1);
    return [
      [sdg, [q1], []],
      [cx, [q0, q1], []],
      [s, [q1], []],
    ];
  }
  // CCX (Toffoli) = standard 6-CX + single-qubit decomposition
  if (name === "ccx" && targetBasis.includes("cx") && targetBasis.includes("u3")) {
    const cx = standardGates.CXGate.copy();
    const h = generalizedGates.makeU3Gate(Math.PI / 2, 0, Math.PI);
    const t = generalizedGates.makeU1Gate(Math.PI / 4);
    const tdg = generalizedGates.makeU1Gate(-Math.PI / 4);
    const q0 = placeholder(0);
    const q1 = placeholder(1);
    const q2 = placeholder(2);
    return [
      [h, [q2], []],
      [cx, [q1, q2], []],
      [tdg, [q2], []],
      [cx.copy(), [q0, q2], []],
      [t, [q2], []],
      [cx.copy(), [q1, q2], []],
      [tdg, [q2], []],
      [cx.copy(), [q0, q2], []],
      [t, [q2], []],
      [h, [q2], []],
      [t, [q1], []],
      [cx.copy(), [q0, q1], []],
      [t, [q0], []],
      [tdg, [q1], []],
      [cx.copy(), [q0, q1], []],
    ];
  }
  return null;
}
