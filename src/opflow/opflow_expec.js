/**
 * opflow_expec.js - Expectation value evaluators for opflow.
 *
 * PauliExpectation (uses Pauli measurements + statevector),
 * MatrixExpectation (exact matrix multiplication),
 * CircuitSampler (samples a state using a circuit).
 */

import { Complex, ComplexMatrix } from "../math/linalg.js";
import { Statevector } from "../quantum_info/statevector.js";
import { SparsePauliOp, Pauli } from "../quantum_info/pauli.js";
import { simulate } from "../simulator/statevector_simulator.js";

// ---------------------------------------------------------------------------
// MatrixExpectation: exact computation via matrix multiplication
// ---------------------------------------------------------------------------
export class MatrixExpectation {
  constructor() {}

  // Compute <ψ|O|ψ> for a CircuitStateFn and operator
  compute(stateFn, operator) {
    const sv = stateFn.primitive || stateFn;
    const statevector = sv instanceof Statevector ? sv : Statevector.fromCircuit(sv);
    const op = operator.primitive || operator;
    const matrix = op instanceof SparsePauliOp ? op.to_matrix() : op.to_matrix();
    const opResult = matrix.matvec(statevector.data);
    const result = statevector.data.inner(opResult);
    return result;
  }

  // Compute expectation for a list of (state, operator) pairs
  compute_list(stateFns, operators) {
    return stateFns.map((sf, i) => this.compute(sf, operators[i]));
  }
}

// ---------------------------------------------------------------------------
// PauliExpectation: decompose into Pauli measurements
// ---------------------------------------------------------------------------
// For each Pauli term in the operator, compute <ψ|P|ψ> separately.
// This is how real quantum computers compute expectation values.
export class PauliExpectation {
  constructor(grouping = false) {
    this.grouping = grouping; // if true, group commuting Paulis
  }

  compute(stateFn, operator) {
    const sv = stateFn.primitive || stateFn;
    const statevector = sv instanceof Statevector ? sv : Statevector.fromCircuit(sv);
    const op = operator.primitive || operator;
    const spo = op instanceof SparsePauliOp ? op : SparsePauliOp.from_list(op);

    // Sum <ψ|P_i|ψ> * coeff_i for each Pauli term
    let result = Complex.ZERO;
    for (let i = 0; i < spo.paulis.size; i++) {
      const pauli = spo.paulis.get(i);
      const coeff = spo.coeffs[i];
      const pauliMatrix = pauli.to_matrix();
      const opResult = pauliMatrix.matvec(statevector.data);
      const expP = statevector.data.inner(opResult);
      result = result.add(expP.mul(coeff));
    }
    return result;
  }

  compute_list(stateFns, operators) {
    return stateFns.map((sf, i) => this.compute(sf, operators[i]));
  }
}

// ---------------------------------------------------------------------------
// CircuitSampler: sample a state using a circuit
// ---------------------------------------------------------------------------
export class CircuitSampler {
  constructor(backend = null, shots = 1024) {
    this.backend = backend;
    this.shots = shots;
  }

  // Sample a state from a circuit (returns a Statevector or samples)
  sample(circuit, parameters = null) {
    let qc = circuit;
    if (parameters) {
      qc = circuit.bind_parameters(parameters);
    }
    if (this.backend) {
      const result = this.backend.run(qc, this.shots);
      return result;
    }
    return Statevector.fromCircuit(qc);
  }

  // Convert a CircuitStateFn to a DictStateFn (sampled state)
  convert(stateFn, parameters = null) {
    const circuit = stateFn.circuit || stateFn.primitive;
    return this.sample(circuit, parameters);
  }
}

// ---------------------------------------------------------------------------
// ExpectationFactory: pick the right expectation evaluator
// ---------------------------------------------------------------------------
export function get_expectation(operator, backend = null) {
  if (backend) {
    return new PauliExpectation();
  }
  return new MatrixExpectation();
}
