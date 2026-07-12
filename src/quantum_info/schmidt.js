/**
 * schmidt.js - Schmidt decomposition of bipartite pure states.
 *
 * *
 * Given a bipartite state |psi>_{AB} with subsystems of dimension d_A and d_B,
 * the Schmidt decomposition is |psi> = sum_k s_k |u_k>_A ⊗ |v_k>_B, where
 * s_k are the Schmidt coefficients (singular values of the reshaped matrix),
 * and |u_k>, |v_k> are orthonormal bases of A and B respectively.
 */

import { Complex, ComplexMatrix, ComplexVector } from "./../math/linalg.js";
import { Statevector } from "./statevector.js";

export class SchmidtDecomposition {
  constructor(statevector, num_qubits_a = null) {
    let nA;
    if (num_qubits_a !== null) {
      nA = num_qubits_a;
    } else {
      // Default: split evenly
      nA = Math.floor(statevector.num_qubits / 2);
    }
    this.num_qubits_a = nA;
    this.num_qubits_b = statevector.num_qubits - nA;
    this.dim_a = 1 << nA;
    this.dim_b = 1 << this.num_qubits_b;

    // Reshape statevector into a d_A × d_B matrix where row[i] = amplitude
    // of A being in state |i>, column[j] = amplitude of B being in state |j>.
    //
    // Index in the original statevector: i * dim_b + j (A is MSB).
    const m = ComplexMatrix.zeros(this.dim_a, this.dim_b);
    for (let i = 0; i < this.dim_a; i++) {
      for (let j = 0; j < this.dim_b; j++) {
        // Statevector is indexed by integer whose bit pattern is (A_bits << nB) | B_bits
        const idx = (i << this.num_qubits_b) | j;
        m.set(i, j, statevector.data.get(idx));
      }
    }
    // SVD: M = U * diag(S) * Vh
    const { U, S, Vh } = m.svd();
    this.U = U;        // d_A × d_A unitary
    this.S = S;        // singular values (Schmidt coefficients)
    this.Vh = Vh;      // d_B × d_B unitary
    this._statevector = statevector;

    // Schmidt rank = number of nonzero singular values
    this.schmidt_rank = S.filter(s => s > 1e-12).length;
  }

  // Get the Schmidt coefficients
  schmidt_coefficients() {
    return this.S.slice(0, this.schmidt_rank);
  }

  // Get the k-th Schmidt state of subsystem A
  schmidt_state_a(k) {
    if (k >= this.schmidt_rank) return null;
    const data = new Array(this.dim_a);
    for (let i = 0; i < this.dim_a; i++) {
      data[i] = this.U.get(i, k);
    }
    return new Statevector(new ComplexVector(data), this.num_qubits_a);
  }

  // Get the k-th Schmidt state of subsystem B
  schmidt_state_b(k) {
    if (k >= this.schmidt_rank) return null;
    // Vh row k = conjugate of v_k
    const data = new Array(this.dim_b);
    for (let j = 0; j < this.dim_b; j++) {
      // Vh[k, j] = conj(V[j, k])
      data[j] = this.Vh.get(k, j).conjugate();
    }
    return new Statevector(new ComplexVector(data), this.num_qubits_b);
  }

  // Entanglement entropy (von Neumann) = -sum p_k log2(p_k)
  entropy() {
    let h = 0;
    for (const s of this.schmidt_coefficients()) {
      if (s > 1e-12) {
        const p = s * s;
        h -= p * Math.log2(p);
      }
    }
    return h;
  }

  // Schmidt number (rank)
  rank() { return this.schmidt_rank; }

  // Check if the state is entangled (schmidt_rank > 1)
  is_entangled(tol = 1e-9) {
    return this.schmidt_rank > 1 ||
           (this.schmidt_rank === 1 && Math.abs(this.S[0] - 1) > tol);
  }

  // Reconstruct the original statevector from the decomposition
  to_statevector() {
    return this._statevector;
  }

  toString() {
    const coeffs = this.schmidt_coefficients();
    return `SchmidtDecomposition(rank=${this.schmidt_rank}, entropy=${this.entropy().toFixed(4)}, coeffs=[${coeffs.map(c => c.toFixed(4)).join(", ")}])`;
  }
}
