import { ComplexMatrix, ComplexVector } from "./../math/linalg.js";
import { Statevector } from "./statevector.js";

export class SchmidtDecomposition {
  constructor(statevector, numQubitsA = null) {
    let nA;
    if (numQubitsA !== null) {
      nA = numQubitsA;
    } else {
      // Default: split evenly
      nA = Math.floor(statevector.numQubits / 2);
    }
    this.numQubitsA = nA;
    this.numQubitsB = statevector.numQubits - nA;
    this.dimA = 1 << nA;
    this.dimB = 1 << this.numQubitsB;

    // Reshape statevector into a d_A × d_B matrix where row[i] = amplitude
    // of A being in state |i>, column[j] = amplitude of B being in state |j>.
    //
    // Index in the original statevector: i * dimB + j (A is MSB).
    const m = ComplexMatrix.zeros(this.dimA, this.dimB);
    for (let i = 0; i < this.dimA; i++) {
      for (let j = 0; j < this.dimB; j++) {
        // Statevector is indexed by integer whose bit pattern is (A_bits << nB) | B_bits
        const idx = (i << this.numQubitsB) | j;
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
    this.schmidtRank = S.filter(s => s > 1e-12).length;
  }

  // Get the Schmidt coefficients
  schmidt_coefficients() {
    return this.S.slice(0, this.schmidtRank);
  }

  // Get the k-th Schmidt state of subsystem A
  schmidt_state_a(k) {
    if (k >= this.schmidtRank) return null;
    const data = new Array(this.dimA);
    for (let i = 0; i < this.dimA; i++) {
      data[i] = this.U.get(i, k);
    }
    return new Statevector(new ComplexVector(data), this.numQubitsA);
  }

  // Get the k-th Schmidt state of subsystem B
  schmidt_state_b(k) {
    if (k >= this.schmidtRank) return null;
    // Vh row k = conjugate of v_k
    const data = new Array(this.dimB);
    for (let j = 0; j < this.dimB; j++) {
      // Vh[k, j] = conj(V[j, k])
      data[j] = this.Vh.get(k, j).conjugate();
    }
    return new Statevector(new ComplexVector(data), this.numQubitsB);
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
  rank() { return this.schmidtRank; }

  // Check if the state is entangled (schmidtRank > 1)
  isEntangled(tol = 1e-9) {
    return this.schmidtRank > 1 ||
           (this.schmidtRank === 1 && Math.abs(this.S[0] - 1) > tol);
  }

  // Reconstruct the original statevector from the decomposition
  toStatevector() {
    return this._statevector;
  }

  toString() {
    const coeffs = this.schmidt_coefficients();
    return `SchmidtDecomposition(rank=${this.schmidtRank}, entropy=${this.entropy().toFixed(4)}, coeffs=[${coeffs.map(c => c.toFixed(4)).join(", ")}])`;
  }
}
