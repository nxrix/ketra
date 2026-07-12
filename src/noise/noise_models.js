/**
 * noise_models.js - Quantum noise models and error channels.
 *
 * Quantum noise models and error channels.
 * Provides:
 *   - QuantumError: a list of (Kraus operators, probability) tuples
 *   - NoiseModel: maps errors to gates/qubits
 *   - Standard channels: depolarizing, bit_flip, phase_flip, amplitude_damping,
 *     phase_damping, pauli_x/y/z, kraus, readout_error
 *   - Helpers: combine_errors, mixed_unitary_error
 */

import { Complex, ComplexMatrix, ComplexVector, PAULI } from "../math/linalg.js";

// ---------------------------------------------------------------------------
// QuantumError
// ---------------------------------------------------------------------------
export class QuantumError {
  constructor(terms = []) {
    // terms: array of { operators: ComplexMatrix[], probability: number }
    this.terms = terms;
    // Verify probabilities sum to 1 (or close)
    const total = terms.reduce((s, t) => s + t.probability, 0);
    if (terms.length > 0 && Math.abs(total - 1) > 1e-9) {
      // Renormalize
      for (const t of this.terms) t.probability /= total;
    }
  }

  get size() { return this.terms.length; }
  get num_qubits() {
    if (this.terms.length === 0) return 0;
    const rows = this.terms[0].operators[0].rows;
    const n = Math.log2(rows);
    if (!Number.isInteger(n)) {
      throw new Error(`QuantumError: Kraus operator has non-power-of-2 dimension ${rows}`);
    }
    return n;
  }

  // Apply this error to a density matrix (Kraus representation)
  apply(densityMatrix) {
    const dim = densityMatrix.rows;
    const result = ComplexMatrix.zeros(dim, dim);
    for (const term of this.terms) {
      for (const K of term.operators) {
        // rho' += p * K * rho * K^dagger
        const Krho = K.mul(densityMatrix);
        const KrhoKd = Krho.mul(K.dagger());
        result.add_inplace(KrhoKd.scale(term.probability));
      }
    }
    return result;
  }

  // Compose two errors (sequential application)
  compose(other) {
    const newTerms = [];
    for (const t1 of this.terms) {
      for (const t2 of other.terms) {
        const newOps = [];
        for (const K1 of t1.operators) {
          for (const K2 of t2.operators) {
            newOps.push(K2.mul(K1));
          }
        }
        newTerms.push({ operators: newOps, probability: t1.probability * t2.probability });
      }
    }
    return new QuantumError(newTerms);
  }

  // Tensor product of two errors
  tensor(other) {
    const newTerms = [];
    for (const t1 of this.terms) {
      for (const t2 of other.terms) {
        const newOps = [];
        for (const K1 of t1.operators) {
          for (const K2 of t2.operators) {
            newOps.push(K1.tensor(K2));
          }
        }
        newTerms.push({ operators: newOps, probability: t1.probability * t2.probability });
      }
    }
    return new QuantumError(newTerms);
  }

  copy() {
    return new QuantumError(this.terms.map(t => ({
      operators: t.operators.map(K => new ComplexMatrix(K.rows, K.cols, K.data.slice())),
      probability: t.probability,
    })));
  }
}

// ---------------------------------------------------------------------------
// Standard error channels
// ---------------------------------------------------------------------------

// Depolarizing error: with prob p, replace state with I/2^n
export function depolarizing_error(prob, numQubits) {
  // For 1 qubit: E(rho) = (1-p)*rho + p*I/2
  // Kraus: K0 = sqrt(1-p)*I, K1 = sqrt(p/3)*X, K2 = sqrt(p/3)*Y, K3 = sqrt(p/3)*Z
  // For n qubits: standard depolarizing
  const dim = 1 << numQubits;
  const identity = ComplexMatrix.identity(dim);
  if (numQubits === 1) {
    const p = prob;
    const K0 = identity.scale(Math.sqrt(1 - p));
    const K1 = PAULI.X.scale(Math.sqrt(p / 3));
    const K2 = PAULI.Y.scale(Math.sqrt(p / 3));
    const K3 = PAULI.Z.scale(Math.sqrt(p / 3));
    return new QuantumError([
      { operators: [K0], probability: 1 - p },
      { operators: [K1], probability: p / 3 },
      { operators: [K2], probability: p / 3 },
      { operators: [K3], probability: p / 3 },
    ]);
  }
  // Multi-qubit: use Pauli group
  const paulis = _generatePauliGroup(numQubits);
  const numNonId = paulis.length - 1;
  const terms = [{ operators: [identity], probability: 1 - prob }];
  for (let i = 1; i < paulis.length; i++) {
    terms.push({ operators: [paulis[i]], probability: prob / numNonId });
  }
  return new QuantumError(terms);
}

// Bit flip: with prob p, apply X
export function bit_flip_error(prob) {
  const K0 = ComplexMatrix.identity(2).scale(Math.sqrt(1 - prob));
  const K1 = PAULI.X.scale(Math.sqrt(prob));
  return new QuantumError([
    { operators: [K0], probability: 1 - prob },
    { operators: [K1], probability: prob },
  ]);
}

// Phase flip: with prob p, apply Z
export function phase_flip_error(prob) {
  const K0 = ComplexMatrix.identity(2).scale(Math.sqrt(1 - prob));
  const K1 = PAULI.Z.scale(Math.sqrt(prob));
  return new QuantumError([
    { operators: [K0], probability: 1 - prob },
    { operators: [K1], probability: prob },
  ]);
}

// Pauli X error: with prob p, apply X (alias for bit_flip)
export function pauli_x_error(prob) { return bit_flip_error(prob); }
export function pauli_y_error(prob) {
  const K0 = ComplexMatrix.identity(2).scale(Math.sqrt(1 - prob));
  const K1 = PAULI.Y.scale(Math.sqrt(prob));
  return new QuantumError([
    { operators: [K0], probability: 1 - prob },
    { operators: [K1], probability: prob },
  ]);
}
export function pauli_z_error(prob) { return phase_flip_error(prob); }

// Amplitude damping: |1> -> |0> with prob gamma
export function amplitude_damping_error(gamma) {
  const K0 = ComplexMatrix.fromRows([
    [new Complex(1, 0), new Complex(0, 0)],
    [new Complex(0, 0), new Complex(Math.sqrt(1 - gamma), 0)],
  ]);
  const K1 = ComplexMatrix.fromRows([
    [new Complex(0, 0), new Complex(Math.sqrt(gamma), 0)],
    [new Complex(0, 0), new Complex(0, 0)],
  ]);
  return new QuantumError([
    { operators: [K0, K1], probability: 1 },
  ]);
}

// Phase damping: loses phase info without losing energy
export function phase_damping_error(gamma) {
  const K0 = ComplexMatrix.fromRows([
    [new Complex(1, 0), new Complex(0, 0)],
    [new Complex(0, 0), new Complex(Math.sqrt(1 - gamma), 0)],
  ]);
  const K1 = ComplexMatrix.fromRows([
    [new Complex(0, 0), new Complex(0, 0)],
    [new Complex(0, 0), new Complex(Math.sqrt(gamma), 0)],
  ]);
  return new QuantumError([
    { operators: [K0, K1], probability: 1 },
  ]);
}

// Reset error: with prob p, reset to |0>
export function reset_error(prob, numQubits = 1) {
  const dim = 1 << numQubits;
  const proj0 = ComplexMatrix.zeros(dim, dim);
  proj0.set(0, 0, new Complex(1, 0));
  const identity = ComplexMatrix.identity(dim);
  const K0 = identity.scale(Math.sqrt(1 - prob));
  const K1 = proj0.scale(Math.sqrt(prob));
  return new QuantumError([
    { operators: [K0], probability: 1 - prob },
    { operators: [K1], probability: prob },
  ]);
}

// General Kraus error
export function kraus_error(krausOps) {
  // krausOps: array of ComplexMatrix
  // Single term with all Kraus operators, probability 1
  return new QuantumError([
    { operators: krausOps, probability: 1 },
  ]);
}

// Mixed unitary error: list of (unitary, probability) pairs
export function mixed_unitary_error(errors) {
  // errors: [{ unitary: ComplexMatrix, probability: number }]
  return new QuantumError(errors.map(e => ({
    operators: [e.unitary],
    probability: e.probability,
  })));
}

// Readout error: prob of measuring wrong bit value
export function ReadoutError(probabilities) {
  // probabilities: [[p00, p01], [p10, p11]] where pij = P(measured=i | actual=j)
  // For 1 qubit: [[p(0|0), p(0|1)], [p(1|0), p(1|1)]]
  this.probabilities = probabilities;
  this.num_qubits = Math.log2(probabilities.length);

  this.apply = function(classicalBits) {
    // classicalBits: array of 0/1
    const measured = [];
    for (let i = 0; i < classicalBits.length; i++) {
      const actual = classicalBits[i];
      const r = Math.random();
      const probActual = this.probabilities[actual]; // [p(0|actual), p(1|actual)]
      measured.push(r < probActual[0] ? 0 : 1);
    }
    return measured;
  };

  this.copy = function() {
    return new ReadoutError(this.probabilities.map(row => row.slice()));
  };
}

// ---------------------------------------------------------------------------
// NoiseModel
// ---------------------------------------------------------------------------
export class NoiseModel {
  constructor() {
    // Map from gate name -> [QuantumError]
    this._local_quantum_errors = new Map();
    // Map from gate name -> {qubit: QuantumError} (per-qubit errors)
    this._local_readout_errors = new Map();
    // Map from gate name -> QuantumError (applied to all qubits)
    this._basis_gate_errors = new Map();
    this._readout_errors = [];
  }

  get basis_gates() {
    const gates = new Set();
    for (const g of this._local_quantum_errors.keys()) gates.add(g);
    for (const g of this._basis_gate_errors.keys()) gates.add(g);
    return Array.from(gates);
  }

  // Add a quantum error for a specific gate on specific qubits
  add_quantum_error(error, gates, qubits = null) {
    const gateList = Array.isArray(gates) ? gates : [gates];
    for (const gate of gateList) {
      const key = qubits ? `${gate}:${qubits.join(",")}` : gate;
      if (qubits) {
        if (!this._local_quantum_errors.has(gate)) {
          this._local_quantum_errors.set(gate, new Map());
        }
        this._local_quantum_errors.get(gate).set(qubits.join(","), error);
      } else {
        this._basis_gate_errors.set(gate, error);
      }
    }
    return this;
  }

  // Add a readout error for specific qubits
  add_readout_error(error, qubits = null) {
    if (qubits) {
      this._readout_errors.push({ error, qubits });
    } else {
      this._readout_errors.push({ error, qubits: "all" });
    }
    return this;
  }

  // Get the quantum error for a gate on specific qubits (or null)
  get_quantum_error(gate, qubits) {
    // Check local first
    if (this._local_quantum_errors.has(gate)) {
      const local = this._local_quantum_errors.get(gate);
      const key = qubits.join(",");
      if (local.has(key)) return local.get(key);
    }
    // Then basis gate
    if (this._basis_gate_errors.has(gate)) {
      return this._basis_gate_errors.get(gate);
    }
    return null;
  }

  // Check if a gate has an error
  has_quantum_error(gate, qubits) {
    return this.get_quantum_error(gate, qubits) !== null;
  }

  // Get all readout errors
  get_readout_errors() {
    return this._readout_errors;
  }

  copy() {
    const nm = new NoiseModel();
    nm._local_quantum_errors = new Map(this._local_quantum_errors);
    nm._basis_gate_errors = new Map(this._basis_gate_errors);
    nm._readout_errors = this._readout_errors.slice();
    return nm;
  }

  is_empty() {
    return this._local_quantum_errors.size === 0 &&
           this._basis_gate_errors.size === 0 &&
           this._readout_errors.length === 0;
  }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function _generatePauliGroup(numQubits) {
  // Generate all n-qubit Pauli operators (tensor products of I, X, Y, Z)
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

// Combine multiple errors into one (applied in sequence)
export function combine_errors(...errors) {
  if (errors.length === 0) return new QuantumError([]);
  let result = errors[0];
  for (let i = 1; i < errors.length; i++) {
    result = result.compose(errors[i]);
  }
  return result;
}
