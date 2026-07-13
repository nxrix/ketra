import { QuantumCircuit } from "./../core/circuit.js";
import { Parameter } from "./../core/parameter.js";
import * as generalizedGates from "./generalized_gates.js";

// QuantumCircuit.prototype.cp/crx/cry/crz are defined in core/circuit.js.

export function bellState() {
  const qc = new QuantumCircuit(2);
  qc.h(0);
  qc.cx(0, 1);
  return qc;
}

export function ghzState(numQubits) {
  if (numQubits < 2) throw new Error("GHZ requires at least 2 qubits");
  const qc = new QuantumCircuit(numQubits);
  qc.h(0);
  for (let i = 0; i < numQubits - 1; i++) {
    qc.cx(i, i + 1);
  }
  return qc;
}

export function qft(numQubits) {
  const qc = new QuantumCircuit(numQubits);
  for (let i = numQubits - 1; i >= 0; i--) {
    qc.h(i);
    for (let j = i - 1; j >= 0; j--) {
      const lam = Math.PI / Math.pow(2, i - j);
      qc.cp(lam, j, i);
    }
  }
  for (let i = 0; i < Math.floor(numQubits / 2); i++) {
    qc.swap(i, numQubits - 1 - i);
  }
  return qc;
}

export function qftInverse(numQubits) {
  return qft(numQubits).inverse();
}

export function quantumVolume(numQubits, depth = null, seed = null) {
  const d = depth || numQubits;
  const qc = new QuantumCircuit(numQubits);
  const rng = seed != null ? _makeRng(seed) : Math.random;
  for (let layer = 0; layer < d; layer++) {
    const perm = _randomPermutation(numQubits, rng);
    for (let i = 0; i + 1 < numQubits; i += 2) {
      const q1 = perm[i];
      const q2 = perm[i + 1];
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q1);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q2);
      qc.cx(q1, q2);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q1);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q2);
      qc.cx(q1, q2);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q1);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q2);
      qc.cx(q1, q2);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q1);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q2);
    }
    qc.barrier();
  }
  return qc;
}

export function realAmplitudes(numQubits, reps = 3, entanglement = "full") {
  const qc = new QuantumCircuit(numQubits);
  for (let q = 0; q < numQubits; q++) {
    qc.ry(new Parameter(`theta[0][${q}]`), q);
  }
  qc.barrier();
  for (let r = 1; r <= reps; r++) {
    _applyEntanglement(qc, numQubits, entanglement, "cx");
    qc.barrier();
    for (let q = 0; q < numQubits; q++) {
      qc.ry(new Parameter(`theta[${r}][${q}]`), q);
    }
    qc.barrier();
  }
  return qc;
}

export function efficientSU2(numQubits, reps = 3, entanglement = "full") {
  const qc = new QuantumCircuit(numQubits);
  for (let q = 0; q < numQubits; q++) {
    qc.ry(new Parameter(`theta[0][${q}]`), q);
    qc.rz(new Parameter(`phi[0][${q}]`), q);
  }
  qc.barrier();
  for (let r = 1; r <= reps; r++) {
    _applyEntanglement(qc, numQubits, entanglement, "cx");
    qc.barrier();
    for (let q = 0; q < numQubits; q++) {
      qc.ry(new Parameter(`theta[${r}][${q}]`), q);
      qc.rz(new Parameter(`phi[${r}][${q}]`), q);
    }
    qc.barrier();
  }
  return qc;
}

export function pauliTwoDesign(numQubits, reps = 3, seed = 12345) {
  const qc = new QuantumCircuit(numQubits);
  const rng = _makeRng(seed);
  for (let q = 0; q < numQubits; q++) {
    qc.ry(new Parameter(`p[0][${q}]`), q);
  }
  const paulis = ["x", "y", "z", "xx", "yy", "zz"];
  for (let r = 0; r < reps; r++) {
    const layer = paulis[Math.floor(rng() * paulis.length)];
    if (layer === "x") qc.x(r % numQubits);
    else if (layer === "y") qc.y(r % numQubits);
    else if (layer === "z") qc.z(r % numQubits);
    else if (layer === "xx") qc.cx(r % numQubits, (r + 1) % numQubits);
    else if (layer === "yy") qc.cy(r % numQubits, (r + 1) % numQubits);
    else if (layer === "zz") qc.cz(r % numQubits, (r + 1) % numQubits);
    for (let q = 0; q < numQubits; q++) {
      qc.rz(new Parameter(`p[${r + 1}][${q}]`), q);
    }
    qc.barrier();
  }
  return qc;
}

export function twoLocal(numQubits, rotationBlocks = ["ry"], entanglementBlocks = ["cx"], reps = 3, entanglement = "full") {
  // TwoLocal builds an ansatz with `reps + 1` rotation layers separated by
  // `reps` entanglement layers. Each (layer, qubit, block) triple must get
  // its own uniquely-named Parameter so that VQE/QAOA optimizers can treat
  // them as independent variables. The naming matches qiskit's convention:
  //   `theta[layer][qubit]` for the per-block angle (with a block suffix
  //   when there are multiple rotation blocks).
  const qc = new QuantumCircuit(numQubits);
  const rotName = (block) => block.toLowerCase();
  const paramFor = (layer, q, block) => {
    if (rotationBlocks.length === 1) {
      return new Parameter(`theta[${layer}][${q}]`);
    }
    return new Parameter(`${rotName(block)}[${layer}][${q}]`);
  };
  const applyRotation = (q, layer, block) => {
    const p = paramFor(layer, q, block);
    if (block === "rx") qc.rx(p, q);
    else if (block === "ry") qc.ry(p, q);
    else if (block === "rz") qc.rz(p, q);
    else throw new Error(`twoLocal: unknown rotation block ${block}`);
  };
  // Initial rotation layer (layer 0).
  for (let q = 0; q < numQubits; q++) {
    for (const block of rotationBlocks) applyRotation(q, 0, block);
  }
  // Repeated entanglement + rotation layers.
  for (let r = 1; r <= reps; r++) {
    for (const block of entanglementBlocks) {
      _applyEntanglement(qc, numQubits, entanglement, block.toLowerCase());
    }
    for (let q = 0; q < numQubits; q++) {
      for (const block of rotationBlocks) applyRotation(q, r, block);
    }
  }
  return qc;
}

// GraphState: prepare |+> on all qubits, then entangle via graph edges
export function graphState(numQubits, edges) {
  const qc = new QuantumCircuit(numQubits);
  for (let q = 0; q < numQubits; q++) qc.h(q);
  for (const [i, j] of edges) qc.cz(i, j);
  return qc;
}

// PauliEvolutionGate: e^{-i t P} where P is a Pauli string.
// `time` may be a number or a ParameterExpression (for VQE/QAOA ansaetze).
export function pauliEvolution(pauliLabel, time) {
  const n = pauliLabel.length;
  const qc = new QuantumCircuit(n);
  // Decomposition: rotate each non-I qubit to the Z basis, CNOT-chain them
  // onto the last qubit, apply RZ(2t) there, then undo the chain and basis.
  // Reference: Nielsen & Chuang, section 4.2 (e^{-i θ Z⊗Z⊗...⊗Z} circuit).
  const chars = pauliLabel.toUpperCase().split("");
  // First, change basis: X -> H, Y -> Sdg H, Z -> nothing, I -> nothing
  for (let i = 0; i < n; i++) {
    const ch = chars[n - 1 - i]; // qubit 0 is rightmost char
    if (ch === "X") qc.h(i);
    else if (ch === "Y") { qc.sdg(i); qc.h(i); }
  }
  // CNOT chain: connect all non-I qubits
  const nonI = [];
  for (let i = 0; i < n; i++) {
    if (chars[n - 1 - i] !== "I") nonI.push(i);
  }
  for (let k = 0; k + 1 < nonI.length; k++) {
    qc.cx(nonI[k], nonI[k + 1]);
  }
  // RZ(2 * t) on the last qubit of the chain (only if there is one)
  if (nonI.length > 0) {
    // Use Parameter-aware multiplication: if t is a ParameterExpression,
    // t.mul(2) keeps it bindable; if it's a number, 2*t also works.
    const angle = (time && typeof time.mul === "function") ? time.mul(2) : (2 * time);
    qc.rz(angle, nonI[nonI.length - 1]);
  }
  // Undo CNOT chain
  for (let k = nonI.length - 2; k >= 0; k--) {
    qc.cx(nonI[k], nonI[k + 1]);
  }
  // Undo basis change
  for (let i = 0; i < n; i++) {
    const ch = chars[n - 1 - i];
    if (ch === "X") qc.h(i);
    else if (ch === "Y") { qc.h(i); qc.s(i); }
  }
  return qc;
}

// HiddenLinearFunction
export function hiddenLinearFunction(numQubits, matrix) {
  // matrix is a 2D array of 0/1
  const qc = new QuantumCircuit(numQubits);
  for (let q = 0; q < numQubits; q++) qc.h(q);
  for (let i = 0; i < numQubits; i++) {
    for (let j = i + 1; j < numQubits; j++) {
      if (matrix[i][j]) qc.cz(i, j);
    }
  }
  for (let q = 0; q < numQubits; q++) qc.h(q);
  return qc;
}

// IQFT (alias for qftInverse)
export const iqft = qftInverse;

// Helpers
function _applyEntanglement(qc, n, entanglement, gate) {
  // Apply a 2-qubit entangling gate between qubits i and j.
  const apply2q = (i, j) => {
    switch (gate) {
      case "cx": qc.cx(i, j); break;
      case "cz": qc.cz(i, j); break;
      case "cy": qc.cy(i, j); break;
      case "ch": qc.ch(i, j); break;
      case "swap": qc.swap(i, j); break;
      case "iswap": qc.iswap(i, j); break;
      case "dcx": qc.dcx(i, j); break;
      case "csx": qc.csx(i, j); break;
      default: throw new Error(`Unknown entanglement block: ${gate}`);
    }
  };
  if (entanglement === "full") {
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) apply2q(i, j);
    }
  } else if (entanglement === "linear") {
    for (let i = 0; i + 1 < n; i++) apply2q(i, i + 1);
  } else if (entanglement === "circular") {
    for (let i = 0; i + 1 < n; i++) apply2q(i, i + 1);
    if (n > 1) apply2q(n - 1, 0);
  } else if (entanglement === "sca") {
    // SCA (sorted circular anticommuting): reverse direction per layer
    for (let i = 0; i + 1 < n; i++) apply2q(i, i + 1);
    if (n > 1) apply2q(n - 1, 0);
  } else if (entanglement === "pairwise") {
    // Adjacent disjoint pairs: (0,1), (2,3), ...
    for (let i = 0; i + 1 < n; i += 2) apply2q(i, i + 1);
  } else {
    throw new Error(`Unknown entanglement type: ${entanglement}`);
  }
}

function _randomPermutation(n, rng) {
  const arr = Array.from({ length: n }, (_, i) => i);
  for (let i = n - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    const tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
  }
  return arr;
}

function _makeRng(seed) {
  let s = seed >>> 0;
  return function () {
    s |= 0; s = (s + 0x6D2B79F5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
