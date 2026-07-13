/**
 * operator.js - Quantum Operator class.
 *
 * */

import { Complex, ComplexMatrix, kronMatrices, PAULI } from "./../math/linalg.js";
import { _registerOperatorClass } from "./../core/circuit.js";

export class Operator {
  constructor(data) {
    if (data instanceof ComplexMatrix) {
      this._data = data;
    } else if (Array.isArray(data) && Array.isArray(data[0])) {
      this._data = ComplexMatrix.fromRows(data);
    } else {
      throw new TypeError("Operator requires a ComplexMatrix or 2D array");
    }
    if (this._data.rows !== this._data.cols) {
      throw new Error("Operator must be square");
    }
    this._numQubits = Math.log2(this._data.rows);
    if (!Number.isInteger(this._numQubits)) {
      throw new Error("Operator dimension must be 2^n");
    }
  }

  static fromMatrix(matrix) { return new Operator(matrix); }

  static fromLabel(label) {
    const c = (re, im = 0) => new Complex(re, im);
    const gates = {
      "I": ComplexMatrix.identity(2),
      "X": PAULI.X, "Y": PAULI.Y, "Z": PAULI.Z,
      "0": ComplexMatrix.fromRows([[c(1), c(0)], [c(0), c(0)]]),
      "1": ComplexMatrix.fromRows([[c(0), c(0)], [c(0), c(1)]]),
      "+": ComplexMatrix.fromRows([[c(0.5), c(0.5)], [c(0.5), c(0.5)]]),
      "-": ComplexMatrix.fromRows([[c(0.5), c(-0.5)], [c(-0.5), c(0.5)]]),
      "r": ComplexMatrix.fromRows([[c(0.5, 0.5), c(0)], [c(0.5, -0.5), c(0)]]),
      "l": ComplexMatrix.fromRows([[c(0.5, -0.5), c(0)], [c(0.5, 0.5), c(0)]]),
    };
    const mats = [];
    for (let i = label.length - 1; i >= 0; i--) {
      const ch = label[i].toUpperCase();
      if (!gates[ch]) throw new Error(`Unknown label character: ${ch}`);
      mats.push(gates[ch]);
    }
    return new Operator(kronMatrices(mats));
  }

  static fromCircuit(circuit) {
    const op = Operator.identity(circuit.numQubits);
    for (const ci of circuit.data) {
      if (ci.operation.name === "barrier" || ci.operation.name === "measure" || ci.operation.name === "reset") continue;
      if (ci.operation.numQubits === 0) continue;
      const gateOp = Operator.fromGate(ci.operation);
      const fullOp = gateOp.embedIntoCircuit(circuit.numQubits, ci.qubits.map(q => circuit._qubit_index.get(q)));
      op._data = fullOp._data.mul(op._data);
    }
    return op;
  }

  static fromGate(gate) {
    if (typeof gate.toMatrix === "function") {
      return new Operator(gate.toMatrix());
    }
    throw new TypeError("Cannot construct Operator from gate without matrix");
  }

  static identity(numQubits) {
    return new Operator(ComplexMatrix.identity(1 << numQubits));
  }

  static zero(numQubits) {
    return new Operator(ComplexMatrix.zeros(1 << numQubits));
  }

  static fromPhase(phase, numQubits = 0) {
    // ScalarOp equivalent: just multiply identity by phase
    return new Operator(ComplexMatrix.identity(1 << numQubits).scale(phase));
  }

  get data() { return this._data; }
  get dim() { return this._data.rows; }
  get numQubits() { return this._numQubits; }

  compose(other) {
    if (this._numQubits !== other._numQubits) {
      throw new Error("compose: qubit count mismatch");
    }
    return new Operator(this._data.mul(other._data));
  }

  tensor(other) { return new Operator(this._data.tensor(other._data)); }
  expand(other) { return new Operator(other._data.tensor(this._data)); }
  adjoint() { return new Operator(this._data.dagger()); }
  conjugate() { return new Operator(this._data.conjugate()); }
  transpose() { return new Operator(this._data.transpose()); }
  trace() { return this._data.trace(); }
  det() { return this._data.det(); }

  isUnitary(tol) { return this._data.isUnitary(tol); }
  isHermitian(tol) { return this._data.isHermitian(tol); }

  pow(n) {
    if (n === 0) return Operator.identity(this._numQubits);
    let result = this;
    for (let i = 1; i < n; i++) result = result.compose(this);
    return result;
  }

  exp() { return new Operator(this._data.expm()); }

  equals(other, tol) {
    return other instanceof Operator && this._data.equals(other._data, tol);
  }

  equalsUpToPhase(other, tol) {
    const eps = (typeof tol === "number") ? tol : 1e-9;
    if (!(other instanceof Operator) || this._numQubits !== other._numQubits) return false;
    let phase = null;
    for (let i = 0; i < this._data.rows; i++) {
      for (let j = 0; j < this._data.cols; j++) {
        const a = this._data.get(i, j);
        const b = other._data.get(i, j);
        if (a.abs() > eps) {
          if (b.abs() < eps) return false;
          const r = b.div(a);
          if (phase === null) phase = r;
          else if (!phase.equals(r, eps)) return false;
        } else if (b.abs() > eps) {
          return false;
        }
      }
    }
    return true;
  }

  embedIntoCircuit(numQubits, qubitIndices) {
    if (qubitIndices.length !== this._numQubits) {
      throw new Error("embedIntoCircuit: qubit count mismatch");
    }
    const dim = 1 << numQubits;
    const out = ComplexMatrix.zeros(dim, dim);
    const indices = qubitIndices;
    for (let row = 0; row < dim; row++) {
      for (let col = 0; col < dim; col++) {
        let subRow = 0;
        let subCol = 0;
        let restMatches = true;
        for (let qi = 0; qi < indices.length; qi++) {
          const bit = (row >> indices[qi]) & 1;
          subRow |= (bit << qi);
        }
        for (let qi = 0; qi < indices.length; qi++) {
          const bit = (col >> indices[qi]) & 1;
          subCol |= (bit << qi);
        }
        for (let q = 0; q < numQubits; q++) {
          if (indices.indexOf(q) === -1) {
            const rb = (row >> q) & 1;
            const cb = (col >> q) & 1;
            if (rb !== cb) { restMatches = false; break; }
          }
        }
        if (!restMatches) continue;
        const elem = this._data.get(subRow, subCol);
        out.set(row, col, elem);
      }
    }
    return new Operator(out);
  }

  // Apply to a statevector
  applyToVector(vector) { return this._data.matvec(vector); }

  // Partial trace: trace out the qubits in `qubitsToTraceOut`, returning an
  // Operator on the remaining qubits.
  partialTrace(qubitsToTraceOut) {
    const traceOut = (Array.isArray(qubitsToTraceOut) ? qubitsToTraceOut : [qubitsToTraceOut]).slice();
    if (traceOut.length === 0) return this;
    // Trace out one qubit at a time, using the existing ComplexMatrix.partialTrace.
    // We need to renumber the remaining qubits after each trace.
    const n = this._numQubits;
    let current = this._data;
    let currentN = n;
    // Sort descending so we trace out higher-indexed qubits first (this
    // keeps the bit positions of lower-indexed qubits stable).
    traceOut.sort((a, b) => b - a);
    for (const q of traceOut) {
      current = current.partialTrace(currentN, q);
      currentN -= 1;
    }
    return new Operator(current);
  }

  // Sum of two operators: returns this + other.
  sum(other) {
    if (this._numQubits !== other._numQubits) {
      throw new Error("sum: qubit count mismatch");
    }
    const result = ComplexMatrix.zeros(this._data.rows, this._data.cols);
    for (let i = 0; i < this._data.rows; i++) {
      for (let j = 0; j < this._data.cols; j++) {
        result.set(i, j, this._data.get(i, j).add(other._data.get(i, j)));
      }
    }
    return new Operator(result);
  }

  // Subtract: returns this - other.
  subtract(other) {
    if (this._numQubits !== other._numQubits) {
      throw new Error("subtract: qubit count mismatch");
    }
    const result = ComplexMatrix.zeros(this._data.rows, this._data.cols);
    for (let i = 0; i < this._data.rows; i++) {
      for (let j = 0; j < this._data.cols; j++) {
        result.set(i, j, this._data.get(i, j).sub(other._data.get(i, j)));
      }
    }
    return new Operator(result);
  }

  // Convert to a measurement channel: project the input state onto the operator
  toMatrix() { return this._data; }

  // Eigenvalues (for Hermitian operators)
  eigvals() {
    if (!this.isHermitian()) {
      throw new Error("eigvals only supports Hermitian operators (use svd for general)");
    }
    return this._data.eigh().eigenvalues;
  }

  // SVD decomposition
  svd() { return this._data.svd(); }

  toDict() {
    const rows = this._data.toRows();
    return {
      rows: this._data.rows,
      cols: this._data.cols,
      data: rows.map(r => r.map(c => ({ re: c.re, im: c.im }))),
    };
  }

  toString() {
    return `Operator(numQubits=${this._numQubits})\n${this._data.toString()}`;
  }
}

// Register the Operator class so QuantumCircuit.toGate()/toInstruction() can find it.
_registerOperatorClass(Operator);
