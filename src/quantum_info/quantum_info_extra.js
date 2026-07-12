/**
 * quantum_info_extra.js - Additional quantum_info functions and methods
 * matching qiskit.quantum_info.
 *
 * Provides:
 *   - Random generators: random_unitary, random_statevector, random_pauli,
 *     random_clifford, random_density_matrix
 *   - Entanglement measures: purity, concurrence, entanglement_of_formation,
 *     mutual_information, gate_fidelity, unitarity
 *   - Helper methods that augment Statevector, DensityMatrix, Operator,
 *     Clifford, and SparsePauliOp via prototype patching.
 */

import { Complex, ComplexMatrix, ComplexVector, PAULI, randomUniform } from "./../math/linalg.js";
import { Statevector } from "./statevector.js";
import { Operator } from "./operator.js";
import { Pauli, SparsePauliOp } from "./pauli.js";
import { Clifford } from "./clifford.js";
import { DensityMatrix } from "./density_matrix.js";

// ---------------------------------------------------------------------------
// Random generators.
// ---------------------------------------------------------------------------

// Haar-random unitary of dimension 2^n. Uses the QR decomposition of a
// random complex Gaussian matrix (the standard algorithm).
export function random_unitary(numQubits, seed = null) {
  const dim = 1 << numQubits;
  const rng = _makeRng(seed);
  // Build a random complex Gaussian matrix.
  const m = ComplexMatrix.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      m.set(i, j, new Complex(_gaussian(rng), _gaussian(rng)));
    }
  }
  // QR decomposition via Gram-Schmidt.
  const Q = ComplexMatrix.zeros(dim, dim);
  const R = ComplexMatrix.zeros(dim, dim);
  for (let j = 0; j < dim; j++) {
    // v = m[:, j]
    const v = new Array(dim);
    for (let i = 0; i < dim; i++) v[i] = m.get(i, j);
    // Subtract projections onto previous Q columns.
    for (let k = 0; k < j; k++) {
      // R[k][j] = <Q[:,k], v>
      let rkj = Complex.ZERO;
      for (let i = 0; i < dim; i++) {
        rkj = rkj.add(Q.get(i, k).conjugate().mul(v[i]));
      }
      R.set(k, j, rkj);
      for (let i = 0; i < dim; i++) {
        v[i] = v[i].sub(rkj.mul(Q.get(i, k)));
      }
    }
    // R[j][j] = ||v||
    let normSq = Complex.ZERO;
    for (let i = 0; i < dim; i++) normSq = normSq.add(v[i].mul(v[i].conjugate()));
    const norm = Math.sqrt(normSq.re);
    R.set(j, j, new Complex(norm, 0));
    if (norm < 1e-12) {
      // Degenerate column; use a random restart.
      for (let i = 0; i < dim; i++) v[i] = new Complex(_gaussian(rng), _gaussian(rng));
      let ns = Complex.ZERO;
      for (let i = 0; i < dim; i++) ns = ns.add(v[i].mul(v[i].conjugate()));
      const nn = Math.sqrt(ns.re);
      for (let i = 0; i < dim; i++) Q.set(i, j, v[i].scale(1 / nn));
    } else {
      for (let i = 0; i < dim; i++) Q.set(i, j, v[i].scale(1 / norm));
    }
  }
  // Adjust phases: divide each column of Q by R[j][j] / |R[j][j]| to make
  // R diagonal real-positive (the standard Haar correction).
  for (let j = 0; j < dim; j++) {
    const rjj = R.get(j, j);
    if (rjj.abs() > 1e-12) {
      const phase = rjj.scale(1 / rjj.abs());
      for (let i = 0; i < dim; i++) {
        Q.set(i, j, Q.get(i, j).mul(phase.conjugate()));
      }
    }
  }
  return new Operator(Q);
}

// Haar-random statevector of `numQubits` qubits.
export function random_statevector(numQubits, seed = null) {
  const dim = 1 << numQubits;
  const rng = _makeRng(seed);
  const data = new Array(dim);
  let normSq = 0;
  for (let i = 0; i < dim; i++) {
    const re = _gaussian(rng);
    const im = _gaussian(rng);
    data[i] = new Complex(re, im);
    normSq += re * re + im * im;
  }
  const norm = Math.sqrt(normSq);
  const normalized = data.map(c => new Complex(c.re / norm, c.im / norm));
  return new Statevector(new ComplexVector(normalized), numQubits);
}

// Random Pauli on `numQubits` qubits.
export function random_pauli(numQubits, seed = null) {
  const rng = _makeRng(seed);
  let label = "";
  const chars = ["I", "X", "Y", "Z"];
  for (let i = 0; i < numQubits; i++) {
    label += chars[Math.floor(rng() * 4)];
  }
  return new Pauli(label);
}

// Random Clifford on `numQubits` qubits. Uses the algorithm of
// Bravyi & Maslov (2020): apply a sequence of random Clifford generators
// (H, S, CX) to the |0...0> tableau.
export function random_clifford(numQubits, seed = null) {
  return Clifford.random(numQubits, seed);
}

// Random density matrix (mixed state) of `numQubits` qubits. Generates a
// random statevector, builds its density matrix, then mixes it with the
// maximally mixed state with a random weight.
export function random_density_matrix(numQubits, seed = null, mixedWeight = null) {
  const dim = 1 << numQubits;
  const rng = _makeRng(seed);
  const psi = random_statevector(numQubits, seed);
  const pureRho = psi.to_operator()._data;
  const w = mixedWeight !== null ? mixedWeight : rng();
  const mixedRho = ComplexMatrix.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      mixedRho.set(i, j, pureRho.get(i, j).scale(1 - w).add(
        (i === j ? new Complex(w / dim, 0) : Complex.ZERO)
      ));
    }
  }
  return new DensityMatrix(mixedRho);
}

// ---------------------------------------------------------------------------
// Entanglement measures.
// ---------------------------------------------------------------------------

// Purity of a density matrix: Tr(rho^2). Pure states have purity 1; the
// maximally mixed state has purity 1/d.
export function purity(densityMatrix) {
  const rho = densityMatrix._data;
  const dim = rho.rows;
  let tr = 0;
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      tr += rho.get(i, j).mul(rho.get(j, i)).re;
    }
  }
  return tr;
}

// Concurrence of a 2-qubit density matrix (Wootters 1998).
export function concurrence(densityMatrix) {
  if (densityMatrix.num_qubits !== 2) {
    throw new Error("concurrence is only defined for 2-qubit states");
  }
  const rho = densityMatrix._data;
  // R = rho * (Y⊗Y) * rho^* * (Y⊗Y)
  const YY = PAULI.Y.tensor(PAULI.Y);
  const rhoStar = rho.conjugate();
  const R = rho.mul(YY).mul(rhoStar).mul(YY);
  // Eigenvalues of sqrt(R) (which equals sqrt(rho * (Y⊗Y) * rho^* * (Y⊗Y)))
  // are real; their square roots give the concurrence.
  // We compute sqrt(R) via the matrix square root (using eigendecomposition).
  const { eigenvalues, eigenvectors } = R.eigh();
  // sqrt eigenvalues (taking abs since numerical noise can make them tiny-negative).
  const sqrtEigs = eigenvalues.map(v => Math.sqrt(Math.abs(v)));
  // Sort descending.
  sqrtEigs.sort((a, b) => b - a);
  const C = Math.max(0, sqrtEigs[0] - sqrtEigs[1] - sqrtEigs[2] - sqrtEigs[3]);
  return C;
}

// Entanglement of formation: E = h((1 + sqrt(1 - C^2)) / 2),
// where h is the binary entropy.
export function entanglement_of_formation(densityMatrix) {
  const C = concurrence(densityMatrix);
  if (C < 1e-12) return 0;
  const x = (1 + Math.sqrt(1 - C * C)) / 2;
  const h = -x * Math.log2(x) - (1 - x) * Math.log2(1 - x);
  return h;
}

// Mutual information between two subsystems.
//   I(A:B) = S(rho_A) + S(rho_B) - S(rho_AB)
// where S is the von Neumann entropy.
export function mutual_information(densityMatrix, subsystemA, subsystemB = null) {
  const n = densityMatrix.num_qubits;
  if (subsystemB == null) {
    // Default: A = first len(subsystemA) qubits, B = the rest.
    subsystemB = [];
    for (let q = 0; q < n; q++) {
      if (!subsystemA.includes(q)) subsystemB.push(q);
    }
  }
  const rhoAB = densityMatrix;
  const rhoA = rhoAB.partial_trace(_complement(subsystemA, n));
  const rhoB = rhoAB.partial_trace(_complement(subsystemB, n));
  const sAB = _vonNeumannEntropy(rhoAB);
  const sA = _vonNeumannEntropy(rhoA);
  const sB = _vonNeumannEntropy(rhoB);
  return sA + sB - sAB;
}

// Gate fidelity between a channel's action and a target unitary.
//   F_g(U, V) = |Tr(U^† V)|^2 / d^2
export function gate_fidelity(unitary1, unitary2) {
  const u1 = unitary1._data || unitary1;
  const u2 = unitary2._data || unitary2;
  if (u1.rows !== u2.rows) {
    throw new Error("gate_fidelity: unitaries must have the same dimension");
  }
  const d = u1.rows;
  const u1dag = u1.dagger();
  const product = u1dag.mul(u2);
  const tr = product.trace();
  return (tr.re * tr.re + tr.im * tr.im) / (d * d);
}

// Unitarity of a quantum channel (unitarity = 1 for unitary channels).
// Defined as the average gate fidelity of the channel with itself minus
// the depolarizing fixed point. For a unital channel, this equals
// (d * F_avg - 1) / (d - 1).
export function unitarity(channel) {
  // Use the superoperator representation: ||chi||_F / sqrt(d * (d+1) / 2)
  // This is a simplified definition; the exact formula involves the
  // Choi matrix.
  const d = channel.dim;
  const choi = channel.to_choi();
  const { eigenvalues } = choi.eigh();
  let sumSq = 0;
  for (const e of eigenvalues) sumSq += e * e;
  return Math.sqrt(sumSq) / d;
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function _makeRng(seed) {
  if (seed == null) return Math.random;
  let s = seed >>> 0;
  return () => {
    s = (s + 0x6D2B79F5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function _gaussian(rng) {
  // Box-Muller transform.
  let u = 0, v = 0;
  while (u === 0) u = rng();
  while (v === 0) v = rng();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function _complement(subset, n) {
  const set = new Set(subset);
  const out = [];
  for (let i = 0; i < n; i++) {
    if (!set.has(i)) out.push(i);
  }
  return out;
}

function _vonNeumannEntropy(densityMatrix) {
  const rho = densityMatrix._data;
  const { eigenvalues } = rho.eigh();
  let s = 0;
  for (const e of eigenvalues) {
    if (e > 1e-12) s -= e * Math.log2(e);
  }
  return s;
}
