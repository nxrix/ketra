/**
 * dag_passes.js - Additional DAG analysis and transformation passes.
 *
 * collect_1q_runs, collect_2q_runs, ConsolidateBlocks,
 * OptimizeCliffords, RemoveDiagonalGatesBeforeMeasure,
 * ElidePermutations, RemoveRedundantGates.
 */

import { ComplexMatrix } from "../math/linalg.js";
import { DAGOpNode } from "../dagcircuit/dagcircuit.js";
import { Gate } from "../core/gate.js";
import * as standardGates from "../library/standard_gates.js";

// Collect runs of single-qubit gates on the same qubit
export function collect_1q_runs(dag) {
  const runs = [];
  const visited = new Set();
  for (const node of dag.topological_op_nodes()) {
    if (visited.has(node)) continue;
    if (node.op.num_qubits !== 1) continue;
    if (!node.qubits || node.qubits.length === 0) continue;
    if (["measure", "reset", "barrier"].includes(node.op.name)) continue;
    const run = [node];
    visited.add(node);
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter(s => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.num_qubits !== 1) break;
      if (!next.qubits || next.qubits.length === 0) break;
      if (next.qubits[0] !== current.qubits[0]) break;
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
export function collect_2q_runs(dag) {
  const runs = [];
  const visited = new Set();
  for (const node of dag.topological_op_nodes()) {
    if (visited.has(node)) continue;
    if (node.op.num_qubits !== 2) continue;
    if (["measure", "reset", "barrier"].includes(node.op.name)) continue;
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

// Consolidate blocks of single-qubit gates into a single unitary
export function consolidate_blocks(dag) {
  const runs = collect_1q_runs(dag);
  for (const run of runs) {
    // Compose all gates in the run
    let combined = null;
    for (const node of run) {
      try {
        const m = node.op.to_matrix();
        combined = combined ? m.mul(combined) : m;
      } catch (e) { break; }
    }
    if (combined) {
      // Replace the run with a single unitary gate
      const newGate = new Gate("unitary", 1, [combined]);
      newGate._matrixBuilder = () => combined;
      const newNode = new DAGOpNode(newGate, run[0].qargs.slice(), []);
      dag.substitute_node(run[0], newNode);
      // Remove the rest of the run
      for (let i = 1; i < run.length; i++) {
        try { dag.remove_op_node(run[i]); } catch (e) {}
      }
    }
  }
  return dag;
}

// Optimize Clifford gates: cancel adjacent identical self-inverse Cliffords
// (H, X, Y, Z, CX, CZ, SWAP, SX, SXDG) and simplify S*S into SDG.
export function optimize_cliffords(dag) {
  const selfInverse = new Set(["h", "x", "y", "z", "cx", "cz", "swap", "sx", "sxdg"]);
  const nodes = dag.topological_op_nodes();
  const toRemove = new Set();
  const sSimplifications = []; // {replace: node, with: DAGOpNode}
  const stack = [];
  for (const node of nodes) {
    if (toRemove.has(node)) continue;
    if (stack.length > 0) {
      const prev = stack[stack.length - 1];
      if (prev.op.name === node.op.name &&
          prev.qubits.length === node.qubits.length &&
          prev.qubits.every((q, i) => q === node.qubits[i])) {
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
          newGate._matrixBuilder = () => standardGates.makeSdgGate().to_matrix();
          sSimplifications.push({ replace: prev, with: new DAGOpNode(newGate, node.qargs.slice(), []) });
          toRemove.add(node);
          continue;
        }
      }
    }
    stack.push(node);
  }
  for (const { replace, with: newNode } of sSimplifications) {
    dag.substitute_node(replace, newNode);
  }
  for (const node of toRemove) {
    try { dag.remove_op_node(node); } catch (e) { /* already gone */ }
  }
  return dag;
}

// Remove diagonal gates that immediately precede measurements
export function remove_diagonal_gates_before_measure(dag) {
  const diagonalGates = ["z", "s", "sdg", "t", "tdg", "p", "u1", "rz", "cz", "cp"];
  const toRemove = [];
  for (const node of dag.topological_op_nodes()) {
    if (!diagonalGates.includes(node.op.name)) continue;
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

// Elide permutations: SWAP gates immediately before measurement can be
// removed by relabeling which qubit each measurement reads.
export function elide_permutations(dag) {
  const toRemove = [];
  for (const node of dag.topological_op_nodes()) {
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
    dag.remove_op_node(node);
  }
  return dag;
}

// Remove redundant gates: identity gates and zero-angle rotations.
export function remove_redundant_gates(dag) {
  const toRemove = [];
  for (const node of dag.topological_op_nodes()) {
    if (node.op.name === "id") toRemove.push(node);
    if (["rz", "rx", "ry", "p", "u1"].includes(node.op.name) &&
        node.op.params.length > 0 &&
        typeof node.op.params[0] === "number" &&
        Math.abs(node.op.params[0]) < 1e-15) {
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.remove_op_node(node);
  }
  return dag;
}

// Commutative cancellation: cancel adjacent identical self-inverse gates.
export function commutative_cancellation(dag) {
  const nodes = dag.topological_op_nodes();
  const selfInverse = new Set(["h", "x", "y", "z", "cx", "cz", "swap"]);
  const toRemove = new Set();
  const stack = [];
  for (const node of nodes) {
    if (toRemove.has(node)) continue;
    if (stack.length > 0) {
      const prev = stack[stack.length - 1];
      if (selfInverse.has(prev.op.name) &&
          prev.op.name === node.op.name &&
          prev.qubits.length === node.qubits.length &&
          prev.qubits.every((q, i) => q === node.qubits[i])) {
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

// Template optimization: replace common gate patterns with simpler equivalents.
// H Z H -> X, H X H -> Z, H Sdg S H = H (no-op S Sdg) H -> identity, etc.
export function template_optimization(dag) {
  const nodes = dag.topological_op_nodes();
  const replacements = []; // {node: oldNode, gate: newGate}
  const toRemove = [];
  for (let i = 0; i < nodes.length - 2; i++) {
    const n1 = nodes[i], n2 = nodes[i + 1], n3 = nodes[i + 2];
    if (n1.op.name === "h" && n2.op.name === "z" && n3.op.name === "h" &&
        n1.qubits[0] === n2.qubits[0] && n2.qubits[0] === n3.qubits[0]) {
      // H Z H = X
      const newGate = new Gate("x", 1, []);
      newGate._matrixBuilder = () => standardGates.makeXGate().to_matrix();
      replacements.push({ node: n1, gate: newGate });
      toRemove.push(n2, n3);
    } else if (n1.op.name === "h" && n2.op.name === "x" && n3.op.name === "h" &&
               n1.qubits[0] === n2.qubits[0] && n2.qubits[0] === n3.qubits[0]) {
      // H X H = Z
      const newGate = new Gate("z", 1, []);
      newGate._matrixBuilder = () => standardGates.makeZGate().to_matrix();
      replacements.push({ node: n1, gate: newGate });
      toRemove.push(n2, n3);
    }
  }
  for (const { node, gate } of replacements) {
    dag.substitute_node(node, new DAGOpNode(gate, node.qargs.slice(), []));
  }
  for (const node of toRemove) {
    try { dag.remove_op_node(node); } catch (e) {}
  }
  return dag;
}
