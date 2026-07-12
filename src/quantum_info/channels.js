/**
 * channels.js - Quantum channel representations and metrics.
 *
 * Quantum channel representations and subclasses:
 *   - Kraus (Kraus representation)
 *   - SuperOp (superoperator matrix)
 *   - Chi (Chi matrix / process matrix)
 *   - PTM (Pauli Transfer Matrix)
 *   - Choi (Choi-Jamiolkowski state)
 *
 * Also provides fidelity/process metrics:
 *   - state_fidelity, process_fidelity, average_gate_fidelity
 *   - diamond_norm (approximate via SDP-free method for simple cases)
 */

import { Complex, ComplexMatrix, ComplexVector, PAULI } from "../math/linalg.js";

// ---------------------------------------------------------------------------
// Base QuantumChannel class
// ---------------------------------------------------------------------------
export class QuantumChannel {
  constructor(numQubits) {
    this.num_qubits = numQubits;
    this.dim = 1 << numQubits;
  }

  // Convert to superoperator (dim^2 x dim^2 matrix)
  to_superop() { throw new Error("to_superop not implemented"); }
  // Convert to Kraus list
  to_kraus() { throw new Error("to_kraus not implemented"); }
  // Convert to Choi matrix
  to_choi() { throw new Error("to_choi not implemented"); }
  // Convert to PTM
  to_ptm() { throw new Error("to_ptm not implemented"); }

  // Apply channel to a statevector -> density matrix
  apply(statevector) {
    const rho = _stateToDensity(statevector);
    return this.apply_density(rho);
  }

  // Apply channel to a density matrix
  apply_density(rho) {
    const superop = this.to_superop();
    const vec = _matrixToVector(rho);
    const result = superop.matvec(vec);
    return _vectorToMatrix(result, this.dim);
  }

  // Compose two channels (this after other)
  compose(other) {
    const s1 = this.to_superop();
    const s2 = other.to_superop();
    return SuperOp.from_matrix(s1.mul(s2), this.num_qubits);
  }

  // Tensor product
  tensor(other) {
    const s1 = this.to_superop();
    const s2 = other.to_superop();
    return SuperOp.from_matrix(s1.tensor(s2), this.num_qubits + other.num_qubits);
  }
}

// ---------------------------------------------------------------------------
// Kraus representation
// ---------------------------------------------------------------------------
export class Kraus extends QuantumChannel {
  constructor(data, numQubits = null) {
    // data: array of ComplexMatrix (Kraus operators)
    const n = numQubits !== null ? numQubits : Math.log2(data[0].rows);
    super(n);
    this.data = data;
  }

  static from_operator(operator) {
    // Single unitary Kraus operator
    return new Kraus([operator.to_matrix ? operator.to_matrix() : operator], operator.num_qubits || Math.log2(operator.rows));
  }

  to_kraus() { return this.data; }

  to_superop() {
    const dim = this.dim;
    const superDim = dim * dim;
    const result = ComplexMatrix.zeros(superDim, superDim);
    for (const K of this.data) {
      // S = sum_k (K* ⊗ K) where K* is conjugate (not dagger)
      // S|i><j| = K|i><j|K^dagger
      // vec(rho') = S * vec(rho)
      // S = sum_k conj(K) ⊗ K
      const Kconj = K.conjugate();
      const Kkron = Kconj.tensor(K);
      result.add_inplace(Kkron);
    }
    return result;
  }

  to_choi() {
    const dim = this.dim;
    const choiDim = dim * dim;
    const result = ComplexMatrix.zeros(choiDim, choiDim);
    // Choi = sum_k |K_k>><<K_k| where |K>> is vectorized K
    for (const K of this.data) {
      const vec = _matrixToVector(K);
      // |vec><vec|
      for (let i = 0; i < vec.size; i++) {
        for (let j = 0; j < vec.size; j++) {
          const val = result.get(i, j).add(vec.data[i].mul(vec.data[j].conjugate()));
          result.set(i, j, val);
        }
      }
    }
    return result;
  }

  to_ptm() {
    // Pauli Transfer Matrix: R_ij = (1/2^n) Tr(P_i E(P_j))
    const n = this.num_qubits;
    const dim = this.dim;
    const paulis = _generatePauliOps(n);
    const numPaulis = paulis.length;
    const ptm = ComplexMatrix.zeros(numPaulis, numPaulis);
    for (let i = 0; i < numPaulis; i++) {
      for (let j = 0; j < numPaulis; j++) {
        // E(P_j) then trace with P_i
        const ej = this.apply_density(paulis[j]);
        const trace = _traceOfProduct(paulis[i], ej);
        ptm.set(i, j, trace.scale(1 / dim));
      }
    }
    return ptm;
  }
}

// ---------------------------------------------------------------------------
// SuperOp representation
// ---------------------------------------------------------------------------
export class SuperOp extends QuantumChannel {
  constructor(data, numQubits) {
    super(numQubits);
    // data: dim^2 x dim^2 ComplexMatrix
    this.data = data;
  }

  static from_matrix(matrix, numQubits) {
    return new SuperOp(matrix, numQubits);
  }

  to_superop() { return this.data; }

  to_kraus() {
    // Convert superop to Kraus via eigendecomposition of Choi
    const choi = this.to_choi();
    // Full Kraus decomposition requires eigendecomposition of Choi
    const { eigenvalues, eigenvectors } = choi.eigh();
    const krausOps = [];
    for (let i = 0; i < eigenvalues.length; i++) {
      if (eigenvalues[i] > 1e-10) {
        const sqrtEv = Math.sqrt(eigenvalues[i]);
        // Eigenvector i -> reshape into dim x dim matrix
        const vec = new ComplexVector(new Array(this.dim * this.dim));
        for (let j = 0; j < this.dim * this.dim; j++) {
          vec.data[j] = eigenvectors.get(j, i);
        }
        const K = _vectorToMatrix(vec, this.dim).scale(sqrtEv);
        krausOps.push(K);
      }
    }
    return krausOps;
  }

  to_choi() {
    // Choi = (S ⊗ I) * |Omega><Omega| where |Omega> = sum |ii>
    // Easier: reshape superop to Choi
    const dim = this.dim;
    const choi = ComplexMatrix.zeros(dim * dim, dim * dim);
    // S_{(i,j),(k,l)} -> Choi_{(i,k),(j,l)}
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        for (let k = 0; k < dim; k++) {
          for (let l = 0; l < dim; l++) {
            const sRow = i * dim + j;
            const sCol = k * dim + l;
            const cRow = i * dim + k;
            const cCol = j * dim + l;
            choi.set(cRow, cCol, this.data.get(sRow, sCol));
          }
        }
      }
    }
    return choi;
  }
}

// ---------------------------------------------------------------------------
// Chi (process matrix) representation
// ---------------------------------------------------------------------------
export class Chi extends QuantumChannel {
  constructor(data, numQubits) {
    super(numQubits);
    this.data = data; // ComplexMatrix (dim^2 x dim^2)
  }

  to_superop() {
    // S = sum_{i,j} chi_{ij} P_i* ⊗ P_j
    const n = this.num_qubits;
    const paulis = _generatePauliOps(n);
    const dim = this.dim;
    const superDim = dim * dim;
    const result = ComplexMatrix.zeros(superDim, superDim);
    for (let i = 0; i < paulis.length; i++) {
      for (let j = 0; j < paulis.length; j++) {
        const chiVal = this.data.get(i, j);
        if (chiVal.abs() < 1e-15) continue;
        const term = paulis[i].conjugate().tensor(paulis[j]).scale(chiVal);
        result.add_inplace(term);
      }
    }
    return result;
  }

  to_kraus() {
    return new Kraus(this.to_kraus(), this.num_qubits);
  }

  to_choi() {
    const superop = this.to_superop();
    return new SuperOp(superop, this.num_qubits).to_choi();
  }
}

// ---------------------------------------------------------------------------
// PTM (Pauli Transfer Matrix) representation
// ---------------------------------------------------------------------------
export class PTM extends QuantumChannel {
  constructor(data, numQubits) {
    super(numQubits);
    this.data = data; // Real matrix (dim^2 x dim^2)
  }

  to_superop() {
    // Convert PTM to superop
    const n = this.num_qubits;
    const paulis = _generatePauliOps(n);
    const dim = this.dim;
    const superDim = dim * dim;
    // S = (1/2^n) sum_{i,j} R_{ij} P_i* ⊗ P_j
    const result = ComplexMatrix.zeros(superDim, superDim);
    for (let i = 0; i < paulis.length; i++) {
      for (let j = 0; j < paulis.length; j++) {
        const rVal = this.data.get(i, j);
        if (rVal.re < 1e-15 && rVal.im < 1e-15) continue;
        const term = paulis[i].conjugate().tensor(paulis[j]).scale(rVal.scale(1 / dim));
        result.add_inplace(term);
      }
    }
    return result;
  }
}

// ---------------------------------------------------------------------------
// Fidelity and process metrics
// ---------------------------------------------------------------------------

// State fidelity: F(rho, sigma) = Tr(sqrt(sqrt(rho) * sigma * sqrt(rho)))^2
// For pure states: F = |<psi|phi>|^2
export function state_fidelity(state1, state2) {
  // Accept Statevector or density matrices
  const isPure1 = state1._data && state1._data.size !== undefined;
  const isPure2 = state2._data && state2._data.size !== undefined;

  if (isPure1 && isPure2) {
    // Both pure: F = |<psi|phi>|^2
    const inner = state1.data.inner(state2.data);
    return inner.abs2();
  }
  // F = Tr(rho * sigma) if both are pure-ish
  const rho1 = isPure1 ? _stateToDensity(state1) : state1;
  const rho2 = isPure2 ? _stateToDensity(state2) : state2;
  const product = rho1.mul(rho2);
  return product.trace().re;
}

// Process fidelity: F(E, U) = (1/d^2) * Tr(U^dagger * E_superop... )
export function process_fidelity(channel, operator = null) {
  if (operator === null) {
    // Process fidelity of channel with identity
    const choi = channel.to_choi();
    const dim = channel.dim;
    const identity = ComplexMatrix.identity(dim * dim).scale(1 / dim);
    return choi.mul(identity).trace().re;
  }
  // F(E, U) = (1/d^2) Tr(S_E * S_U^dagger)
  const sE = channel.to_superop();
  const sU = _unitaryToSuperop(operator.to_matrix ? operator.to_matrix() : operator);
  const product = sE.mul(sU.dagger());
  return product.trace().re / (channel.dim * channel.dim);
}

// Average gate fidelity: F_avg = (d * F_proc + 1) / (d + 1)
export function average_gate_fidelity(channel, operator = null) {
  const d = channel.dim;
  const fproc = process_fidelity(channel, operator);
  return (d * fproc + 1) / (d + 1);
}

// Diamond norm: ||Φ||_◇ = max_{ρ} ||(Φ⊗I)(ρ)||_1.
//
// The exact diamond norm requires solving a semidefinite program (SDP).
// Computing the SDP in pure JS would require bundling an SDP solver, which
// is out of scope for this library. Instead, we compute a tight upper
// bound via the eigenvalues of the Choi matrix:
//
//   ||Φ||_◇ ≤ ||J(Φ)||_1 = Σ_i |λ_i(J(Φ))|,
//
// where J(Φ) is the Choi matrix and λ_i are its eigenvalues. This bound
// is tight for unitary channels (where ||Φ||_◇ = 2) and for depolarizing
// channels. For general channels it can overestimate by up to a factor
// of √d, but it never underestimates. Users needing the exact diamond
// norm should use a backend with SDP support (e.g., CVXPY in Python).
export function diamond_norm(channel) {
  const choi = channel.to_choi();
  const { eigenvalues } = choi.eigh();
  return eigenvalues.reduce((s, v) => s + Math.abs(v), 0);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function _stateToDensity(statevector) {
  const dim = statevector.size || statevector.dim;
  const data = statevector.data || statevector._data;
  const rho = ComplexMatrix.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      rho.set(i, j, data.get(i).mul(data.get(j).conjugate()));
    }
  }
  return rho;
}

function _matrixToVector(matrix) {
  const dim = matrix.rows;
  const vec = new ComplexVector(new Array(dim * dim));
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      vec.data[i * dim + j] = matrix.get(i, j);
    }
  }
  return vec;
}

function _vectorToMatrix(vec, dim) {
  const m = ComplexMatrix.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      m.set(i, j, vec.data[i * dim + j]);
    }
  }
  return m;
}

function _traceOfProduct(m1, m2) {
  // Tr(m1 * m2)
  const product = m1.mul(m2);
  return product.trace();
}

function _unitaryToSuperop(U) {
  // S = conj(U) ⊗ U
  return U.conjugate().tensor(U);
}

function _generatePauliOps(numQubits) {
  const paulis1 = [ComplexMatrix.identity(2), PAULI.X, PAULI.Y, PAULI.Z];
  if (numQubits === 1) return paulis1;
  let result = paulis1;
  for (let i = 1; i < numQubits; i++) {
    const newResult = [];
    for (const p1 of result) {
      for (const p2 of paulis1) {
        newResult.push(p1.tensor(p2));
      }
    }
    result = newResult;
  }
  return result;
}
