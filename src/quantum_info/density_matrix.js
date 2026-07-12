/**
 * density_matrix.js - DensityMatrix class for mixed quantum states.
 *
 * Supports construction from
 * Statevector, Operator, circuit, or label. Provides evolve, probabilities,
 * purity, trace, partial_trace, expectation_value, etc.
 */

import { Complex, ComplexMatrix, ComplexVector, sampleDistribution } from "./../math/linalg.js";
import { Statevector, _registerDensityMatrixClass } from "./statevector.js";
import { Operator } from "./operator.js";

export class DensityMatrix {
  constructor(data, numQubits = null) {
    if (data instanceof ComplexMatrix) {
      this._data = data;
    } else if (Array.isArray(data) && Array.isArray(data[0])) {
      this._data = ComplexMatrix.fromRows(data);
    } else {
      throw new TypeError("DensityMatrix requires a ComplexMatrix or 2D array");
    }
    if (this._data.rows !== this._data.cols) {
      throw new Error("DensityMatrix must be square");
    }
    const nq = Math.log2(this._data.rows);
    if (!Number.isInteger(nq)) {
      throw new Error(`DensityMatrix dimension must be 2^n, got ${this._data.rows}`);
    }
    this._numQubits = numQubits !== null ? numQubits : nq;
  }

  static fromStatevector(statevector) {
    const dim = statevector.dim;
    const data = ComplexMatrix.zeros(dim, dim);
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        const a = statevector.data.get(i);
        const b = statevector.data.get(j);
        data.set(i, j, a.mul(b.conjugate()));
      }
    }
    return new DensityMatrix(data, statevector.num_qubits);
  }

  static fromOperator(operator) {
    // Treat an Operator as the density matrix directly
    return new DensityMatrix(operator.data, operator.num_qubits);
  }

  static fromCircuit(circuit, initState = null) {
    // Evolve |0>...|0> through the circuit, then convert to density matrix
    const sv = initState ? initState : Statevector.zero(circuit.num_qubits);
    let dm = DensityMatrix.fromStatevector(sv);
    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (op.num_qubits === 0) continue;
      if (typeof op.to_matrix !== "function") continue;
      // rho' = U rho U^dagger
      const gateMatrix = op.to_matrix();
      const qubitIndices = ci.qubits.map(q => circuit._qubit_index.get(q));
      const fullOp = _embedGate(gateMatrix, circuit.num_qubits, qubitIndices);
      const newRho = fullOp.mul(dm._data).mul(fullOp.dagger());
      dm = new DensityMatrix(newRho, circuit.num_qubits);
    }
    return dm;
  }

  static fromLabel(label) {
    return DensityMatrix.fromStatevector(Statevector.fromLabel(label));
  }

  static zero(numQubits) {
    return new DensityMatrix(ComplexMatrix.zeros(1 << numQubits), numQubits);
  }

  static identity(numQubits) {
    // Maximally mixed state: I / 2^n
    return new DensityMatrix(ComplexMatrix.identity(1 << numQubits).scale(1 / (1 << numQubits)), numQubits);
  }

  get data() { return this._data; }
  get dim() { return this._data.rows; }
  get num_qubits() { return this._numQubits; }

  copy() {
    return new DensityMatrix(new ComplexMatrix(this._data.rows, this._data.cols, this._data.data.slice()), this._numQubits);
  }

  // Trace of density matrix (should be 1 for valid quantum states)
  trace() { return this._data.trace(); }

  // Purity: Tr(rho^2). Pure states have purity 1, mixed < 1.
  purity() {
    const sq = this._data.mul(this._data);
    const tr = sq.trace();
    return tr.re; // should be real for Hermitian rho
  }

  is_valid(tol = 1e-9) {
    // Hermitian, trace 1, positive semi-definite
    if (!this._data.isHermitian(tol)) return false;
    const tr = this.trace();
    if (Math.abs(tr.re - 1) > tol || Math.abs(tr.im) > tol) return false;
    // Check positive semi-definite via eigenvalues
    try {
      const { eigenvalues } = this._data.eigh();
      return eigenvalues.every(ev => ev >= -tol);
    } catch (e) {
      return false;
    }
  }

  is_pure(tol = 1e-9) {
    return Math.abs(this.purity() - 1) < tol;
  }

  // Evolve under a unitary: rho' = U rho U^dagger
  evolve(other) {
    if (other && other._data && other._data instanceof ComplexMatrix) {
      if (other.num_qubits !== this._numQubits) {
        throw new Error("evolve: operator qubit count mismatch");
      }
      return new DensityMatrix(other._data.mul(this._data).mul(other._data.dagger()), this._numQubits);
    }
    if (other && typeof other.to_matrix === "function") {
      if (other.num_qubits !== this._numQubits) {
        throw new Error("evolve: gate qubit count mismatch");
      }
      const m = other.to_matrix();
      return new DensityMatrix(m.mul(this._data).mul(m.dagger()), this._numQubits);
    }
    throw new TypeError("evolve: unsupported operand type");
  }

  // Tensor product with another DensityMatrix
  tensor(other) {
    return new DensityMatrix(this._data.tensor(other._data), this._numQubits + other._numQubits);
  }

  expand(other) {
    return new DensityMatrix(other._data.tensor(this._data), this._numQubits + other._numQubits);
  }

  // Add another density matrix (for mixture)
  add(other) {
    if (this._numQubits !== other._numQubits) {
      throw new Error("add: qubit count mismatch");
    }
    return new DensityMatrix(this._data.add(other._data), this._numQubits);
  }

  scale(c) {
    return new DensityMatrix(this._data.scale(c), this._numQubits);
  }

  // Partial trace over a qubit
  partial_trace(qubits) {
    const qList = Array.isArray(qubits) ? qubits : [qubits];
    let result = this;
    // Sort descending so we trace out higher-index qubits first (keeps indices stable)
    const sorted = qList.slice().sort((a, b) => b - a);
    let currentQubits = this._numQubits;
    let currentData = this._data;
    for (const q of sorted) {
      currentData = currentData.partialTrace(currentQubits, q);
      currentQubits--;
    }
    return new DensityMatrix(currentData, currentQubits);
  }

  // Probabilities of measuring each computational basis state
  probabilities(qargs = null) {
    const diag = new Array(this._data.rows);
    for (let i = 0; i < this._data.rows; i++) {
      diag[i] = this._data.get(i, i).re;
    }
    if (qargs == null) return diag;
    const indices = Array.isArray(qargs) ? qargs : [qargs];
    return _marginalProbabilitiesDM(diag, this._numQubits, indices);
  }

  // Sample measurement outcomes
  sample(counts = 1024, qargs = null, rng) {
    const probs = this.probabilities(qargs);
    const results = {};
    const nq = qargs == null ? this._numQubits : (Array.isArray(qargs) ? qargs.length : 1);
    for (let i = 0; i < counts; i++) {
      const idx = sampleDistribution(probs, rng);
      const bits = new Array(nq).fill("0");
      for (let q = 0; q < nq; q++) bits[q] = ((idx >> q) & 1).toString();
      const key = bits.reverse().join("");
      results[key] = (results[key] || 0) + 1;
    }
    return results;
  }

  // Expectation value <O> = Tr(O rho)
  expectation_value(operator) {
    const m = (operator instanceof ComplexMatrix) ? operator : operator._data;
    const prod = m.mul(this._data);
    const tr = prod.trace();
    return tr;
  }

  // Convert to a Statevector if pure (returns null if not pure)
  to_statevector() {
    if (!this.is_pure()) return null;
    // Find the eigenvalue closest to 1 and use the corresponding eigenvector
    const { eigenvalues, eigenvectors } = this._data.eigh();
    let maxIdx = 0;
    let maxVal = -Infinity;
    for (let i = 0; i < eigenvalues.length; i++) {
      if (eigenvalues[i] > maxVal) {
        maxVal = eigenvalues[i];
        maxIdx = i;
      }
    }
    const dim = this._data.rows;
    const data = new Array(dim);
    for (let i = 0; i < dim; i++) {
      data[i] = eigenvectors.get(i, maxIdx);
    }
    return new Statevector(new ComplexVector(data), this._numQubits);
  }

  equals(other, tol) {
    return other instanceof DensityMatrix && this._data.equals(other._data, tol);
  }

  to_matrix() { return this._data; }

  to_dict() {
    const rows = this._data.toRows();
    return rows.map(r => r.map(c => ({ re: c.re, im: c.im })));
  }

  toString() {
    return `DensityMatrix(num_qubits=${this._numQubits})\n${this._data.toString()}`;
  }
}

function _embedGate(gateMatrix, numQubits, qubitIndices) {
  const dim = 1 << numQubits;
  const out = ComplexMatrix.zeros(dim, dim);
  const k = qubitIndices.length;
  for (let row = 0; row < dim; row++) {
    for (let col = 0; col < dim; col++) {
      let subRow = 0, subCol = 0, restMatches = true;
      for (let qi = 0; qi < k; qi++) {
        const rb = (row >> qubitIndices[qi]) & 1;
        const cb = (col >> qubitIndices[qi]) & 1;
        subRow |= (rb << qi);
        subCol |= (cb << qi);
      }
      for (let q = 0; q < numQubits; q++) {
        if (qubitIndices.indexOf(q) === -1) {
          const rb = (row >> q) & 1;
          const cb = (col >> q) & 1;
          if (rb !== cb) { restMatches = false; break; }
        }
      }
      if (!restMatches) continue;
      out.set(row, col, gateMatrix.get(subRow, subCol));
    }
  }
  return out;
}

function _marginalProbabilitiesDM(diag, numQubits, indices) {
  const outDim = 1 << indices.length;
  const out = new Array(outDim).fill(0);
  for (let i = 0; i < diag.length; i++) {
    let subIdx = 0;
    for (let qi = 0; qi < indices.length; qi++) {
      const bit = (i >> indices[qi]) & 1;
      subIdx |= (bit << qi);
    }
    out[subIdx] += diag[i];
  }
  return out;
}

// Register the DensityMatrix class so Statevector.partial_trace() can find it.
_registerDensityMatrixClass(DensityMatrix);
