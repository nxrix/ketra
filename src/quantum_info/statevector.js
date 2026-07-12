/**
 * statevector.js - Quantum Statevector class.
 *
 * */

import { Complex, ComplexVector, ComplexMatrix, kronVectors, sampleDistribution, PAULI } from "./../math/linalg.js";
import { Operator } from "./operator.js";
import { _registerStatevectorClass } from "./../core/circuit.js";

// DensityMatrix is imported lazily to avoid a circular import (density_matrix.js
// imports Statevector from this file). We register the class below.
let _DensityMatrixClass = null;
export function _registerDensityMatrixClass(cls) { _DensityMatrixClass = cls; }

const ZERO = new ComplexVector([new Complex(1, 0), new Complex(0, 0)]);
const ONE  = new ComplexVector([new Complex(0, 0), new Complex(1, 0)]);

const LABEL_VECTORS = {
  "0": ZERO, "1": ONE,
  "+": new ComplexVector([new Complex(1 / Math.SQRT2, 0), new Complex(1 / Math.SQRT2, 0)]),
  "-": new ComplexVector([new Complex(1 / Math.SQRT2, 0), new Complex(-1 / Math.SQRT2, 0)]),
  "r": new ComplexVector([new Complex(1 / Math.SQRT2, 0), new Complex(0, 1 / Math.SQRT2)]),
  "l": new ComplexVector([new Complex(1 / Math.SQRT2, 0), new Complex(0, -1 / Math.SQRT2)]),
};

export class Statevector {
  constructor(data, numQubits = null) {
    if (data instanceof ComplexVector) {
      this._data = data;
    } else if (Array.isArray(data)) {
      this._data = new ComplexVector(data);
    } else {
      throw new TypeError("Statevector requires a ComplexVector or array");
    }
    const dim = this._data.size;
    const nq = Math.log2(dim);
    if (!Number.isInteger(nq)) {
      throw new Error(`Statevector dimension must be 2^n, got ${dim}`);
    }
    this._numQubits = numQubits !== null ? numQubits : nq;
  }

  static fromLabel(label) {
    const parts = [];
    for (let i = 0; i < label.length; i++) {
      const ch = label[i].toLowerCase();
      if (!LABEL_VECTORS[ch]) throw new Error(`Unknown state label: ${ch}`);
      parts.push(LABEL_VECTORS[ch]);
    }
    if (parts.length === 0) return new Statevector(new ComplexVector([Complex.ONE]), 0);
    return new Statevector(kronVectors(parts), parts.length);
  }

  static zero(numQubits) {
    const data = ComplexVector.zeros(1 << numQubits);
    data.data[0] = Complex.ONE;
    return new Statevector(data, numQubits);
  }

  static one(numQubits) {
    const data = ComplexVector.zeros(1 << numQubits);
    data.data[(1 << numQubits) - 1] = Complex.ONE;
    return new Statevector(data, numQubits);
  }

  static fromInstruction(instruction) {
    const nq = instruction.num_qubits;
    const sv = Statevector.zero(nq);
    if (typeof instruction.to_matrix === "function") {
      const m = instruction.to_matrix();
      return new Statevector(m.matvec(sv._data), nq);
    }
    throw new Error("Instruction has no matrix");
  }

  static fromCircuit(circuit, initState = null) {
    let sv = initState ? initState.copy() : Statevector.zero(circuit.num_qubits);
    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier") continue;
      if (op.name === "measure" || op.name === "reset") continue;
      if (op.num_qubits === 0) continue;
      if (typeof op.to_matrix !== "function") {
        throw new Error(`Cannot evolve Statevector under non-unitary ${op.name}`);
      }
      const fullOp = _embedGate(op, circuit.num_qubits, ci.qubits.map(q => circuit._qubit_index.get(q)));
      sv = new Statevector(fullOp.matvec(sv._data), circuit.num_qubits);
    }
    return sv;
  }

  static from_int(index, numQubits) {
    const data = ComplexVector.zeros(1 << numQubits);
    data.data[index] = Complex.ONE;
    return new Statevector(data, numQubits);
  }

  get data() { return this._data; }
  get dim() { return this._data.size; }
  get num_qubits() { return this._numQubits; }

  copy() {
    return new Statevector(new ComplexVector(this._data.data.slice()), this._numQubits);
  }

  conjugate() { return new Statevector(this._data.conjugate(), this._numQubits); }

  inner(other) { return this._data.inner(other._data); }
  dot(other) { return this._data.inner(other._data); }

  norm() { return this._data.norm(); }

  is_unitary() { return false; }

  probabilities(qargs = null) {
    const probs = this._data.probabilities();
    if (qargs == null) return probs;
    const indices = Array.isArray(qargs) ? qargs : [qargs];
    return _marginalProbabilities(probs, this._numQubits, indices);
  }

  evolve(other) {
    if (other instanceof Statevector) {
      return new Statevector(this._data.tensor(other._data), this._numQubits + other._numQubits);
    }
    if (other && other._data && other._data instanceof ComplexMatrix) {
      if (other.num_qubits !== this._numQubits) {
        throw new Error("evolve: operator qubit count mismatch");
      }
      return new Statevector(other._data.matvec(this._data), this._numQubits);
    }
    if (other && typeof other.to_matrix === "function") {
      if (other.num_qubits !== this._numQubits) {
        throw new Error("evolve: gate qubit count mismatch");
      }
      return new Statevector(other.to_matrix().matvec(this._data), this._numQubits);
    }
    if (other && other.data && other.data instanceof Array) {
      if (typeof Statevector.fromCircuit !== "function") {
        throw new Error("Cannot evolve by circuit without circuit module loaded");
      }
      return Statevector.fromCircuit(other, this);
    }
    throw new TypeError("evolve: unsupported operand type");
  }

  tensor(other) {
    return new Statevector(this._data.tensor(other._data), this._numQubits + other._numQubits);
  }

  expand(other) {
    return new Statevector(other._data.tensor(this._data), this._numQubits + other._numQubits);
  }

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

  sample_memory(counts = 1024, qargs = null, rng) {
    const probs = this.probabilities(qargs);
    const nq = qargs == null ? this._numQubits : (Array.isArray(qargs) ? qargs.length : 1);
    const results = new Array(counts);
    for (let i = 0; i < counts; i++) {
      const idx = sampleDistribution(probs, rng);
      const bits = new Array(nq).fill("0");
      for (let q = 0; q < nq; q++) bits[q] = ((idx >> q) & 1).toString();
      results[i] = bits.reverse().join("");
    }
    return results;
  }

  expectation_value(operator) {
    const m = (operator instanceof ComplexMatrix) ? operator : operator._data;
    const op = m.matvec(this._data);
    return this._data.inner(op);
  }

  // Compute the density matrix rho = |psi><psi|
  to_operator() {
    const dim = this._data.size;
    const m = ComplexMatrix.zeros(dim, dim);
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        m.set(i, j, this._data.get(i).mul(this._data.get(j).conjugate()));
      }
    }
    return new Operator(m);
  }

  equals(other, tol) {
    if (!(other instanceof Statevector)) return false;
    if (other._numQubits !== this._numQubits) return false;
    return this._data.equals(other._data, tol);
  }

  // Check equality up to global phase
  equiv(other, tol = 1e-9) {
    if (!(other instanceof Statevector) || other._numQubits !== this._numQubits) return false;
    let phase = null;
    for (let i = 0; i < this._data.size; i++) {
      const a = this._data.get(i);
      const b = other._data.get(i);
      if (a.abs() > tol) {
        if (b.abs() < tol) return false;
        const r = b.div(a);
        if (phase === null) phase = r;
        else if (!phase.equals(r, tol)) return false;
      } else if (b.abs() > tol) {
        return false;
      }
    }
    return true;
  }

  to_dict() {
    const out = {};
    for (let i = 0; i < this._data.size; i++) {
      const bits = i.toString(2).padStart(this._numQubits, "0");
      out[bits] = this._data.get(i);
    }
    return out;
  }

  // Measure a single qubit (or all qubits if no argument) and collapse
  // the state. Returns a { bit, statevector } pair.
  measure(qubit = null, rng) {
    if (qubit == null) {
      // Measure all qubits, returning a bitstring and a collapsed (basis) state.
      const probs = this.probabilities();
      const idx = sampleDistribution(probs, rng);
      const bits = new Array(this._numQubits).fill("0");
      for (let q = 0; q < this._numQubits; q++) bits[q] = ((idx >> q) & 1).toString();
      const bitstring = bits.reverse().join("");
      // Collapsed state: the basis state |idx>.
      const collapsed = new Array(this._data.size).fill(Complex.ZERO);
      collapsed[idx] = Complex.ONE;
      return { bits: bitstring, statevector: new Statevector(new ComplexVector(collapsed), this._numQubits) };
    }
    // Measure a single qubit: compute marginal probability of |1> on `qubit`,
    // sample the outcome, then collapse.
    const probs = this.probabilities([qubit]);
    const prob0 = probs[0];
    const r = rng ? rng() : Math.random();
    // Outcome is 0 if r < prob0, else 1. (Matches qiskit's convention
    // where the cumulative-distribution sample picks the lowest outcome
    // when r falls in the first probability bin.)
    const outcome = r < prob0 ? 0 : 1;
    // Collapse: zero out the amplitudes inconsistent with `outcome`, then renormalize.
    const dim = this._data.size;
    const newData = new Array(dim);
    let normSq = 0;
    for (let i = 0; i < dim; i++) {
      if (((i >> qubit) & 1) === outcome) {
        newData[i] = this._data.get(i);
        normSq += newData[i].re * newData[i].re + newData[i].im * newData[i].im;
      } else {
        newData[i] = Complex.ZERO;
      }
    }
    const norm = Math.sqrt(normSq);
    for (let i = 0; i < dim; i++) {
      newData[i] = new Complex(newData[i].re / norm, newData[i].im / norm);
    }
    return { bit: outcome, statevector: new Statevector(new ComplexVector(newData), this._numQubits) };
  }

  // Reset a qubit to |0>: measure it, and if outcome is |1>, apply X.
  reset(qubit) {
    const { bit, statevector } = this.measure(qubit);
    if (bit === 0) return statevector;
    // Apply X to flip the qubit back to |0>.
    const X = PAULI.X;
    const dim = statevector._data.size;
    const newData = new Array(dim);
    for (let i = 0; i < dim; i++) {
      const j = i ^ (1 << qubit);
      newData[i] = statevector._data.get(j);
    }
    return new Statevector(new ComplexVector(newData), this._numQubits);
  }

  // Partial trace: trace out the qubits NOT in `keep`, returning a
  // DensityMatrix on the kept qubits.
  partial_trace(qargsToKeep) {
    const keep = Array.isArray(qargsToKeep) ? qargsToKeep : [qargsToKeep];
    // Build the density matrix of the full state, then partial-trace it.
    const rho = this.to_operator(); // Operator on all qubits
    const traceOut = _complement(keep, this._numQubits);
    const reducedOp = rho.partial_trace(traceOut);
    // Convert to a DensityMatrix (lazy to avoid circular import).
    if (!_DensityMatrixClass) {
      throw new Error("DensityMatrix class not registered. Import quantum_info/density_matrix.js first.");
    }
    return new _DensityMatrixClass(reducedOp._data);
  }

  // Expand dimensions: tensor with |0...0> on `numQubits` additional qubits
  // at the front (high-order bits).
  expand_dims(numQubits) {
    const newDim = this._data.size << numQubits;
    const newData = new Array(newDim).fill(Complex.ZERO);
    for (let i = 0; i < this._data.size; i++) {
      newData[i << numQubits] = this._data.get(i);
    }
    return new Statevector(new ComplexVector(newData), this._numQubits + numQubits);
  }

  toString() {
    const terms = [];
    for (let i = 0; i < this._data.size; i++) {
      const amp = this._data.get(i);
      if (amp.abs() < 1e-10) continue;
      const bits = i.toString(2).padStart(this._numQubits, "0");
      terms.push(`${amp.toString()}|${bits}>`);
    }
    return terms.join(" + ");
  }
}

function _embedGate(gate, numQubits, qubitIndices) {
  if (typeof gate.to_matrix !== "function") {
    throw new Error("Cannot embed gate without matrix");
  }
  const gateMatrix = gate.to_matrix();
  if (gate.num_qubits !== qubitIndices.length) {
    throw new Error("embedGate: qubit count mismatch");
  }
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

function _marginalProbabilities(probs, numQubits, indices) {
  const sortedIndices = indices.slice().sort((a, b) => a - b);
  const outDim = 1 << indices.length;
  const out = new Array(outDim).fill(0);
  for (let i = 0; i < probs.length; i++) {
    let subIdx = 0;
    for (let qi = 0; qi < indices.length; qi++) {
      const bit = (i >> indices[qi]) & 1;
      subIdx |= (bit << qi);
    }
    out[subIdx] += probs[i];
  }
  return out;
}

// Register the Statevector class so QuantumCircuit.initialize() can find it.
_registerStatevectorClass(Statevector);

// Return the complement of a subset of qubit indices (the qubits NOT in `subset`).
function _complement(subset, n) {
  const set = new Set(subset);
  const out = [];
  for (let i = 0; i < n; i++) {
    if (!set.has(i)) out.push(i);
  }
  return out;
}
