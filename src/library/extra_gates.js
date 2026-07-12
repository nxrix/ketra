/**
 * extra_gates.js - Additional gate types matching qiskit.extensions and
 * qiskit.circuit.library.
 *
 * Implements:
 *   - UnitaryGate      (arbitrary 2^n x 2^n unitary)
 *   - DiagonalGate     (diagonal unitary, given as a list of phases)
 *   - PermutationGate  (qubit permutation)
 *   - HamiltonianGate  (e^{-iHt} for a Hermitian H)
 *   - Initialize       (state preparation from amplitudes or labels)
 *   - MCRYGate / MCRXGate / MCRZGate / MCPhaseGate (multi-controlled rotations)
 *   - MCMTGate         (multi-control multi-target)
 */

import { Complex, ComplexMatrix } from "./../math/linalg.js";
import { Gate, ControlledGate, Instruction } from "./../core/gate.js";
import { _registerParamBuilder, _registerStd, _registerExtraGateClass } from "./../core/circuit.js";
import { Statevector } from "./../quantum_info/statevector.js";
import { SparsePauliOp } from "./../quantum_info/pauli.js";
import { Operator } from "./../quantum_info/operator.js";

const r = (re, im = 0) => new Complex(re, im);

// ---------------------------------------------------------------------------
// UnitaryGate: an arbitrary 2^n × 2^n unitary matrix.
// ---------------------------------------------------------------------------
export class UnitaryGate extends Gate {
  constructor(matrix, label = null) {
    // Accept a ComplexMatrix, an Operator, or a 2D array of Complex / numbers.
    let mat;
    if (matrix instanceof ComplexMatrix) {
      mat = matrix;
    } else if (matrix && matrix._data instanceof ComplexMatrix) {
      // Operator
      mat = matrix._data;
    } else if (Array.isArray(matrix)) {
      mat = ComplexMatrix.fromRows(matrix.map(row =>
        row.map(v => v instanceof Complex ? v : r(v))
      ));
    } else {
      throw new TypeError("UnitaryGate requires a ComplexMatrix, Operator, or 2D array");
    }
    const n = Math.log2(mat.rows);
    if (!Number.isInteger(n)) {
      throw new Error(`UnitaryGate: matrix dimension ${mat.rows} is not 2^n`);
    }
    super("unitary", n, []);
    this._matrix = mat;
    this.label = label;
    this._matrixBuilder = () => mat;
  }

  copy() {
    const g = new UnitaryGate(this._matrix, this.label);
    return g;
  }

  // Allow a unitary gate to be parameter-bound (it has no parameters, so
  // this is a no-op, but bind_parameters calls it).
  bind() { return this; }

  control(numCtrl = 1) {
    return new ControlledGate(this, numCtrl);
  }
}

// ---------------------------------------------------------------------------
// DiagonalGate: a diagonal unitary U = diag(e^{i*diag[0]}, ..., e^{i*diag[d-1]}).
// ---------------------------------------------------------------------------
export class DiagonalGate extends Gate {
  constructor(diag) {
    const n = Math.log2(diag.length);
    if (!Number.isInteger(n)) {
      throw new Error(`DiagonalGate: diag length ${diag.length} is not 2^n`);
    }
    super("diagonal", n, diag);
    this._diag = diag.slice();
    this._matrixBuilder = (params) => {
      const d = (params && params.length ? params : this._diag);
      const dim = 1 << n;
      const m = ComplexMatrix.zeros(dim, dim);
      for (let i = 0; i < dim; i++) {
        m.set(i, i, r(Math.cos(d[i]), Math.sin(d[i])));
      }
      return m;
    };
  }

  copy() {
    return new DiagonalGate(this._diag);
  }
}

// ---------------------------------------------------------------------------
// PermutationGate: permutes qubits according to a pattern.
//   pattern = [2, 0, 1] means qubit 0 → position 2, qubit 1 → position 0,
//   qubit 2 → position 1. (qiskit's convention: pattern[i] = which qubit
//   ends up at position i.)
// ---------------------------------------------------------------------------
export class PermutationGate extends Gate {
  constructor(pattern) {
    const n = pattern.length;
    super("permutation", n, []);
    this._pattern = pattern.slice();
    this._matrixBuilder = () => {
      const dim = 1 << n;
      const m = ComplexMatrix.zeros(dim, dim);
      for (let i = 0; i < dim; i++) {
        // qubit j of input goes to position pattern[j] in output.
        let j = 0;
        for (let q = 0; q < n; q++) {
          if ((i >> q) & 1) {
            j |= (1 << pattern.indexOf(q));
          }
        }
        m.set(j, i, r(1));
      }
      return m;
    };
  }

  copy() { return new PermutationGate(this._pattern); }
}

// ---------------------------------------------------------------------------
// HamiltonianGate: U = e^{-i * time * H}, where H is a Hermitian operator.
// ---------------------------------------------------------------------------
export class HamiltonianGate extends Gate {
  constructor(operator, time, label = null) {
    let mat;
    if (operator instanceof ComplexMatrix) {
      mat = operator;
    } else if (operator && operator._data instanceof ComplexMatrix) {
      mat = operator._data;
    } else if (operator instanceof SparsePauliOp) {
      mat = operator.to_matrix();
    } else if (Array.isArray(operator)) {
      mat = ComplexMatrix.fromRows(operator.map(row =>
        row.map(v => v instanceof Complex ? v : r(v))
      ));
    } else {
      throw new TypeError("HamiltonianGate requires a matrix, Operator, or SparsePauliOp");
    }
    const n = Math.log2(mat.rows);
    if (!Number.isInteger(n)) {
      throw new Error(`HamiltonianGate: dimension ${mat.rows} is not 2^n`);
    }
    super("hamiltonian", n, [time]);
    this._hamiltonian = mat;
    this._time = time;
    this.label = label;
    this._matrixBuilder = (params) => {
      const t = (params && params.length ? params[0] : this._time);
      // e^{-i t H} via the matrix exponential: -i t H = scale(-i t) · H.
      const scaledH = mat.scale(new Complex(0, -t));
      return scaledH.expm();
    };
  }

  copy() {
    return new HamiltonianGate(this._hamiltonian, this._time, this.label);
  }
}

// ---------------------------------------------------------------------------
// Initialize: state preparation. U|0...0> = |psi>, where |psi> is given as
// amplitudes or as a state label.
// ---------------------------------------------------------------------------
export class Initialize extends Instruction {
  constructor(amplitudes, numQubits = null) {
    let amps;
    let n;
    if (typeof amplitudes === "string") {
      const sv = Statevector.fromLabel(amplitudes);
      amps = sv._data.data;
      n = amplitudes.length;
    } else if (Array.isArray(amplitudes)) {
      amps = amplitudes.map(v => v instanceof Complex ? v : r(v));
      n = Math.log2(amps.length);
      if (!Number.isInteger(n)) {
        throw new Error(`Initialize: amplitudes length ${amps.length} is not 2^n`);
      }
    } else {
      throw new TypeError("Initialize requires an array of amplitudes or a label string");
    }
    if (numQubits !== null) n = numQubits;
    super("initialize", n, 0, amps);
    this._amplitudes = amps;
    this._matrixBuilder = () => {
      const dim = amps.length;
      const m = ComplexMatrix.identity(dim);
      for (let i = 0; i < dim; i++) m.set(i, 0, amps[i]);
      for (let j = 1; j < dim; j++) m.set(0, j, r(0));
      return m;
    };
  }

  copy() {
    return new Initialize(this._amplitudes);
  }
}

// ---------------------------------------------------------------------------
// Multi-controlled rotation gates.
// ---------------------------------------------------------------------------
export class MCPhaseGate extends Gate {
  constructor(lambda, numCtrlQubits) {
    super("mcphase", numCtrlQubits + 1, [lambda]);
    this._lambda = lambda;
    this._numCtrl = numCtrlQubits;
    this._matrixBuilder = (params) => {
      const l = (params && params.length ? params[0] : this._lambda);
      const numCtrl = this._numCtrl;
      const dim = 1 << (numCtrl + 1);
      const m = ComplexMatrix.identity(dim);
      // Set the |1...1> diagonal element to e^{i*lambda}.
      m.set(dim - 1, dim - 1, r(Math.cos(l), Math.sin(l)));
      return m;
    };
  }
  copy() { return new MCPhaseGate(this._lambda, this._numCtrl); }
}

export class MCRXGate extends Gate {
  constructor(theta, numCtrlQubits) {
    super("mcrx", numCtrlQubits + 1, [theta]);
    this._theta = theta;
    this._numCtrl = numCtrlQubits;
    this._matrixBuilder = (params) => {
      const t = (params && params.length ? params[0] : this._theta);
      const numCtrl = this._numCtrl;
      const dim = 1 << (numCtrl + 1);
      const m = ComplexMatrix.identity(dim);
      // Build RX(t) on the last 2x2 block (the |1...1> subspace).
      const ct = Math.cos(t / 2);
      const st = Math.sin(t / 2);
      const idx00 = dim - 2;
      const idx11 = dim - 1;
      m.set(idx00, idx00, r(ct));
      m.set(idx00, idx11, r(0, -st));
      m.set(idx11, idx00, r(0, -st));
      m.set(idx11, idx11, r(ct));
      return m;
    };
  }
  copy() { return new MCRXGate(this._theta, this._numCtrl); }
}

export class MCRYGate extends Gate {
  constructor(theta, numCtrlQubits) {
    super("mcry", numCtrlQubits + 1, [theta]);
    this._theta = theta;
    this._numCtrl = numCtrlQubits;
    this._matrixBuilder = (params) => {
      const t = (params && params.length ? params[0] : this._theta);
      const numCtrl = this._numCtrl;
      const dim = 1 << (numCtrl + 1);
      const m = ComplexMatrix.identity(dim);
      const ct = Math.cos(t / 2);
      const st = Math.sin(t / 2);
      const idx00 = dim - 2;
      const idx11 = dim - 1;
      m.set(idx00, idx00, r(ct));
      m.set(idx00, idx11, r(-st));
      m.set(idx11, idx00, r(st));
      m.set(idx11, idx11, r(ct));
      return m;
    };
  }
  copy() { return new MCRYGate(this._theta, this._numCtrl); }
}

export class MCRZGate extends Gate {
  constructor(theta, numCtrlQubits) {
    super("mcrz", numCtrlQubits + 1, [theta]);
    this._theta = theta;
    this._numCtrl = numCtrlQubits;
    this._matrixBuilder = (params) => {
      const t = (params && params.length ? params[0] : this._theta);
      const numCtrl = this._numCtrl;
      const dim = 1 << (numCtrl + 1);
      const m = ComplexMatrix.identity(dim);
      const ct = Math.cos(t / 2);
      const st = Math.sin(t / 2);
      const idx00 = dim - 2;
      const idx11 = dim - 1;
      m.set(idx00, idx00, r(ct, -st));
      m.set(idx11, idx11, r(ct, st));
      return m;
    };
  }
  copy() { return new MCRZGate(this._theta, this._numCtrl); }
}

// ---------------------------------------------------------------------------
// MCMTGate: Multi-Control Multi-Target gate.
//   Applies `base_gate` to each target qubit when all control qubits are |1>.
// ---------------------------------------------------------------------------
export class MCMTGate extends Gate {
  constructor(baseGate, numCtrlQubits, numTargetQubits) {
    super("mcmt", numCtrlQubits + numTargetQubits, baseGate.params.slice());
    this._baseGate = baseGate;
    this._numCtrl = numCtrlQubits;
    this._numTarget = numTargetQubits;
    this._matrixBuilder = () => {
      const baseMat = baseGate.to_matrix();
      const baseDim = baseMat.rows; // 2^(base num_qubits), typically 2
      if (numTargetQubits === 1 && baseDim === 2) {
        // Single-target MCMT is just a controlled version of the base gate.
        const cg = new ControlledGate(baseGate, numCtrlQubits);
        return cg.to_matrix();
      }
      // For multi-target, build the full unitary by stacking controlled
      // versions of the base gate on each target.
      // The full Hilbert space is 2^(numCtrl + numTarget).
      const totalDim = 1 << (numCtrlQubits + numTargetQubits);
      const m = ComplexMatrix.identity(totalDim);
      // For each target qubit position, when all controls are |1>, apply
      // the base gate to that target.
      const ctrlOn = (1 << numCtrlQubits) - 1;
      // The block where all controls are |1> has the target qubits in
      // positions [numCtrl, numCtrl+1, ..., numCtrl+numTarget-1].
      // We apply I⊗...⊗I⊗base to each target individually.
      // For simplicity (and matching qiskit's behavior for the common
      // case of single-qubit base gates), we apply the base gate
      // simultaneously to all targets via tensor product.
      const targetBlockDim = 1 << numTargetQubits;
      // Build the target-block unitary: base^{⊗ numTarget}.
      let targetU = baseMat;
      for (let i = 1; i < numTargetQubits; i++) {
        targetU = targetU.tensor(baseMat);
      }
      // Place targetU in the |1...1> control block.
      for (let bi = 0; bi < targetBlockDim; bi++) {
        for (let bj = 0; bj < targetBlockDim; bj++) {
          const row = ctrlOn + (bi << numCtrlQubits);
          const col = ctrlOn + (bj << numCtrlQubits);
          m.set(row, col, targetU.get(bi, bj));
        }
      }
      return m;
    };
  }
  copy() {
    return new MCMTGate(this._baseGate.copy(), this._numCtrl, this._numTarget);
  }
}

// Register a generic builder so QuantumCircuit.unitary(...) can look it up.
_registerParamBuilder("UNITARY", (matrix) => new UnitaryGate(matrix));
_registerParamBuilder("DIAGONAL", (diag) => new DiagonalGate(diag));
_registerParamBuilder("PERMUTATION", (pattern) => new PermutationGate(pattern));
_registerParamBuilder("HAMILTONIAN", (op, t) => new HamiltonianGate(op, t));

// Register the extra gate classes so QuantumCircuit methods can find them
// (avoids circular imports — circuit.js can't statically import this file
// because we import _registerParamBuilder from it above).
_registerExtraGateClass("UnitaryGate", UnitaryGate);
_registerExtraGateClass("DiagonalGate", DiagonalGate);
_registerExtraGateClass("PermutationGate", PermutationGate);
_registerExtraGateClass("HamiltonianGate", HamiltonianGate);
_registerExtraGateClass("Initialize", Initialize);
_registerExtraGateClass("MCPhaseGate", MCPhaseGate);
_registerExtraGateClass("MCRXGate", MCRXGate);
_registerExtraGateClass("MCRYGate", MCRYGate);
_registerExtraGateClass("MCRZGate", MCRZGate);
_registerExtraGateClass("MCMTGate", MCMTGate);
