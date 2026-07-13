import { DAGOpNode } from "../dagcircuit/dagcircuit.js";
import { Gate } from "../core/gate.js";
import * as standardGates from "../library/standard_gates.js";

// Collect runs of single-qubit gates on the same qubit
export function collect1qRuns(dag) {
  const runs = [];
  const visited = new Set();
  for (const node of dag.topologicalOpNodes()) {
    if (visited.has(node)) continue;
    if (node.op.numQubits !== 1) continue;
    if (!node.qargs || node.qargs.length === 0) continue;
    if (["measure", "reset", "barrier"].includes(node.op.name)) continue;
    const run = [node];
    visited.add(node);
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter(s => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.numQubits !== 1) break;
      if (!next.qargs || next.qargs.length === 0) break;
      if (next.qargs[0] !== current.qargs[0]) break;
      if (["measure", "reset", "barrier"].includes(next.op.name)) break;
      run.push(next);
      visited.add(next);
      current = next;
    }
    if (run.length > 1) runs.push(run);
  }
  return runs;
}

// Collect runs of 2-qubit gates on the same pair of qubits
export function collect2qRuns(dag) {
  const runs = [];
  const visited = new Set();
  for (const node of dag.topologicalOpNodes()) {
    if (visited.has(node)) continue;
    if (node.op.numQubits !== 2) continue;
    if (["measure", "reset", "barrier"].includes(node.op.name)) continue;
    const run = [node];
    visited.add(node);
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter(s => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.numQubits !== 2) break;
      if (next.qargs[0] !== current.qargs[0] || next.qargs[1] !== current.qargs[1]) break;
      run.push(next);
      visited.add(next);
      current = next;
    }
    if (run.length > 1) runs.push(run);
  }
  return runs;
}

// Consolidate blocks of single-qubit gates into a single unitary
export function consolidateBlocks(dag) {
  const runs = collect1qRuns(dag);
  for (const run of runs) {
    // Compose all gates in the run
    let combined = null;
    for (const node of run) {
      try {
        const m = node.op.toMatrix();
        combined = combined ? m.mul(combined) : m;
      } catch (e) { break; }
    }
    if (combined) {
      // Replace the run with a single unitary gate
      const newGate = new Gate("unitary", 1, [combined]);
      newGate._matrixBuilder = () => combined;
      const newNode = new DAGOpNode(newGate, run[0].qargs.slice(), []);
      dag.substituteNode(run[0], newNode);
      // Remove the rest of the run
      for (let i = 1; i < run.length; i++) {
        try { dag.removeOpNode(run[i]); } catch (e) {}
      }
    }
  }
  return dag;
}

// Optimize Clifford gates: cancel adjacent identical self-inverse Cliffords
// (H, X, Y, Z, CX, CZ, SWAP, SX, SXDG) and simplify S*S into SDG.
export function optimizeCliffords(dag) {
  const selfInverse = new Set(["h", "x", "y", "z", "cx", "cz", "swap", "sx", "sxdg"]);
  const nodes = dag.topologicalOpNodes();
  const toRemove = new Set();
  const sSimplifications = []; // {replace: node, with: DAGOpNode}
  const stack = [];
  for (const node of nodes) {
    if (toRemove.has(node)) continue;
    if (stack.length > 0) {
      const prev = stack[stack.length - 1];
      if (prev.op.name === node.op.name &&
          prev.qargs.length === node.qargs.length &&
          prev.qargs.every((q, i) => q === node.qargs[i])) {
        if (selfInverse.has(prev.op.name)) {
          stack.pop();
          toRemove.add(prev);
          toRemove.add(node);
          continue;
        }
        // S * S = SDG
        if (prev.op.name === "s") {
          stack.pop();
          const newGate = new Gate("sdg", 1, []);
          newGate._matrixBuilder = () => standardGates.makeSdgGate().toMatrix();
          sSimplifications.push({ replace: prev, with: new DAGOpNode(newGate, node.qargs.slice(), []) });
          toRemove.add(node);
          continue;
        }
      }
    }
    stack.push(node);
  }
  for (const { replace, with: newNode } of sSimplifications) {
    dag.substituteNode(replace, newNode);
  }
  for (const node of toRemove) {
    try { dag.removeOpNode(node); } catch (e) { /* already gone */ }
  }
  return dag;
}

// Remove diagonal gates that immediately precede measurements
export function removeDiagonalGatesBeforeMeasure(dag) {
  const diagonalGates = ["z", "s", "sdg", "t", "tdg", "p", "u1", "rz", "cz", "cp"];
  const toRemove = [];
  for (const node of dag.topologicalOpNodes()) {
    if (!diagonalGates.includes(node.op.name)) continue;
    const succs = dag.successors(node).filter(s => s.is_op_node());
    if (succs.length > 0 && succs.every(s => s.op.name === "measure")) {
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.removeOpNode(node);
  }
  return dag;
}

// Elide permutations: SWAP gates immediately before measurement can be
// removed by relabeling which qubit each measurement reads.
export function elidePermutations(dag) {
  const toRemove = [];
  for (const node of dag.topologicalOpNodes()) {
    if (node.op.name !== "swap") continue;
    const succs = dag.successors(node).filter(s => s.is_op_node());
    if (succs.length === 2 && succs.every(s => s.op.name === "measure")) {
      // Swap the classical targets of the two measurements so the SWAP
      // becomes a no-op on the logical state.
      const m0 = succs[0], m1 = succs[1];
      const tmp = m0.cargs[0];
      m0.cargs[0] = m1.cargs[0];
      m1.cargs[0] = tmp;
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.removeOpNode(node);
  }
  return dag;
}

// Remove redundant gates: identity gates and zero-angle rotations.
export function removeRedundantGates(dag) {
  const toRemove = [];
  for (const node of dag.topologicalOpNodes()) {
    if (node.op.name === "id") toRemove.push(node);
    if (["rz", "rx", "ry", "p", "u1"].includes(node.op.name) &&
        node.op.params.length > 0 &&
        typeof node.op.params[0] === "number" &&
        Math.abs(node.op.params[0]) < 1e-15) {
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.removeOpNode(node);
  }
  return dag;
}

// Commutative cancellation: cancel adjacent identical self-inverse gates.
export function commutativeCancellation(dag) {
  const nodes = dag.topologicalOpNodes();
  const selfInverse = new Set(["h", "x", "y", "z", "cx", "cz", "swap"]);
  const toRemove = new Set();
  const stack = [];
  for (const node of nodes) {
    if (toRemove.has(node)) continue;
    if (stack.length > 0) {
      const prev = stack[stack.length - 1];
      if (selfInverse.has(prev.op.name) &&
          prev.op.name === node.op.name &&
          prev.qargs.length === node.qargs.length &&
          prev.qargs.every((q, i) => q === node.qargs[i])) {
        stack.pop();
        toRemove.add(prev);
        toRemove.add(node);
        continue;
      }
    }
    stack.push(node);
  }
  for (const node of toRemove) {
    try { dag.removeOpNode(node); } catch (e) { /* already gone */ }
  }
  return dag;
}

// Template optimization: replace common gate patterns with simpler equivalents.
// H Z H -> X, H X H -> Z, H Sdg S H = H (no-op S Sdg) H -> identity, etc.
export function templateOptimization(dag) {
  const nodes = dag.topologicalOpNodes();
  const replacements = []; // {node: oldNode, gate: newGate}
  const toRemove = [];
  for (let i = 0; i < nodes.length - 2; i++) {
    const n1 = nodes[i], n2 = nodes[i + 1], n3 = nodes[i + 2];
    if (n1.op.name === "h" && n2.op.name === "z" && n3.op.name === "h" &&
        n1.qargs[0] === n2.qargs[0] && n2.qargs[0] === n3.qargs[0]) {
      // H Z H = X
      const newGate = new Gate("x", 1, []);
      newGate._matrixBuilder = () => standardGates.makeXGate().toMatrix();
      replacements.push({ node: n1, gate: newGate });
      toRemove.push(n2, n3);
    } else if (n1.op.name === "h" && n2.op.name === "x" && n3.op.name === "h" &&
               n1.qargs[0] === n2.qargs[0] && n2.qargs[0] === n3.qargs[0]) {
      // H X H = Z
      const newGate = new Gate("z", 1, []);
      newGate._matrixBuilder = () => standardGates.makeZGate().toMatrix();
      replacements.push({ node: n1, gate: newGate });
      toRemove.push(n2, n3);
    }
  }
  for (const { node, gate } of replacements) {
    dag.substituteNode(node, new DAGOpNode(gate, node.qargs.slice(), []));
  }
  for (const node of toRemove) {
    try { dag.removeOpNode(node); } catch (e) {}
  }
  return dag;
}
