/**
 * opflow.js - Operator flow algebra (deprecated qiskit.opflow).
 *
 * Operator flow algebra: OperatorBase, PauliOp, PauliSumOp, MatrixOp,
 * CircuitOp, StateFn, CircuitStateFn, ListOp, EvolvedOp.
 *
 * This is qiskit's older high-level algebra system. While deprecated in
 * favor of qiskit.quantum_info + primitives, it's still widely used in
 * existing code.
 */

import { Complex, ComplexMatrix, ComplexVector } from "./../math/linalg.js";
import { Pauli, PauliList, SparsePauliOp } from "./../quantum_info/pauli.js";
import { Operator } from "./../quantum_info/operator.js";
import { Statevector } from "./../quantum_info/statevector.js";

// ---------------------------------------------------------------------------
// OperatorBase
// ---------------------------------------------------------------------------
export class OperatorBase {
  constructor() {
    this.num_qubits = 0;
    this.coeff = Complex.ONE;
  }

  add(other) { return new ListOp([this, other], "add"); }
  sub(other) { return new ListOp([this, other.scale(-1)], "add"); }
  compose(other) {
    // this ∘ other (matrix multiplication)
    return new ListOp([this, other], "compose");
  }
  tensor(other) { return new ListOp([this, other], "tensor"); }
  expand(other) { return new ListOp([other, this], "tensor"); }
  pow(n) { return new ListOp([this], "pow", { exponent: n }); }
  scale(c) {
    if (typeof c === "number") c = new Complex(c, 0);
    return new ListOp([this], "scale", { coeff: c });
  }
  adjoint() { return new ListOp([this], "adjoint"); }

  to_matrix() { throw new Error("OperatorBase.to_matrix not implemented"); }
  to_matrix_op() { return new MatrixOp(this.to_matrix()); }
  to_spmatrix_op() {
    return new PauliSumOp(SparsePauliOp.from_operator(new Operator(this.to_matrix())));
  }

  // Apply to a state
  apply_to(state) {
    const m = this.to_matrix();
    return m.matvec(state.data || state);
  }

  is_hermitian() { return this.to_matrix().isHermitian(); }
}

// ---------------------------------------------------------------------------
// PauliOp: A single Pauli operator
// ---------------------------------------------------------------------------
export class PauliOp extends OperatorBase {
  constructor(pauli, coeff = Complex.ONE) {
    super();
    this.primitive = pauli instanceof Pauli ? pauli : new Pauli(pauli);
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.num_qubits = this.primitive.num_qubits;
  }

  to_matrix() {
    return this.primitive.to_matrix().scale(this.coeff);
  }

  adjoint() {
    // Pauli^dagger = Pauli (up to sign for Y, but Y is self-adjoint up to phase)
    return new PauliOp(this.primitive, this.coeff.conjugate());
  }

  is_hermitian() {
    return Math.abs(this.coeff.im) < 1e-9;
  }

  toString() { return `${this.coeff.toString()} * ${this.primitive.label}`; }
}

// ---------------------------------------------------------------------------
// PauliSumOp: Sum of Pauli operators (wraps SparsePauliOp)
// ---------------------------------------------------------------------------
export class PauliSumOp extends OperatorBase {
  constructor(spmatrixop, coeff = Complex.ONE) {
    super();
    this.primitive = spmatrixop instanceof SparsePauliOp ? spmatrixop : SparsePauliOp.from_list(spmatrixop);
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.num_qubits = this.primitive.num_qubits;
  }

  static from_list(list) {
    return new PauliSumOp(SparsePauliOp.from_list(list));
  }

  to_matrix() {
    return this.primitive.to_matrix().scale(this.coeff);
  }

  adjoint() {
    return new PauliSumOp(this.primitive.adjoint(), this.coeff.conjugate());
  }

  is_hermitian() {
    return this.to_matrix().isHermitian();
  }

  // Reduce: combine like terms
  reduce() {
    return new PauliSumOp(this.primitive.simplify(), this.coeff);
  }

  // Convert to a list of [Pauli, coeff] pairs
  to_list() { return this.primitive.to_list(); }

  toString() { return this.primitive.toString(); }
}

// ---------------------------------------------------------------------------
// MatrixOp: An operator backed by a dense matrix
// ---------------------------------------------------------------------------
export class MatrixOp extends OperatorBase {
  constructor(matrix, coeff = Complex.ONE) {
    super();
    if (matrix instanceof ComplexMatrix) {
      this.primitive = matrix;
    } else if (matrix instanceof Operator) {
      this.primitive = matrix.data;
    } else {
      this.primitive = matrix;
    }
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    const n = Math.log2(this.primitive.rows);
    this.num_qubits = n;
  }

  to_matrix() {
    return this.primitive.scale(this.coeff);
  }

  adjoint() {
    return new MatrixOp(this.primitive.dagger(), this.coeff.conjugate());
  }

  is_hermitian() {
    return this.to_matrix().isHermitian();
  }

  toString() { return `MatrixOp(num_qubits=${this.num_qubits})`; }
}

// ---------------------------------------------------------------------------
// CircuitOp: An operator backed by a QuantumCircuit
// ---------------------------------------------------------------------------
export class CircuitOp extends OperatorBase {
  constructor(circuit, coeff = Complex.ONE) {
    super();
    this.primitive = circuit;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.num_qubits = circuit.num_qubits;
  }

  to_matrix() {
    const op = Operator.fromCircuit(this.primitive);
    return op.data.scale(this.coeff);
  }

  to_circuit() { return this.primitive; }

  adjoint() {
    return new CircuitOp(this.primitive.inverse(), this.coeff.conjugate());
  }

  toString() { return `CircuitOp(${this.primitive.name})`; }
}

// ---------------------------------------------------------------------------
// StateFn: A state functional (linear functional <psi| or vector |psi>)
// ---------------------------------------------------------------------------
export class StateFn extends OperatorBase {
  constructor(primitive, is_measurement = false, coeff = Complex.ONE) {
    super();
    if (primitive instanceof Statevector) {
      this.primitive = primitive;
    } else if (primitive instanceof ComplexVector) {
      this.primitive = new Statevector(primitive);
    } else if (typeof primitive === "string") {
      this.primitive = Statevector.fromLabel(primitive);
    } else {
      this.primitive = primitive;
    }
    this.is_measurement = is_measurement;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.num_qubits = this.primitive.num_qubits;
  }

  static from_label(label, is_measurement = false) {
    return new StateFn(Statevector.fromLabel(label), is_measurement);
  }

  to_matrix() {
    // If measurement: <psi| (row vector), else |psi> (column vector)
    if (this.is_measurement) {
      // Row vector = conjugate transpose of statevector
      const dim = this.primitive.dim;
      const m = ComplexMatrix.zeros(1, dim);
      for (let i = 0; i < dim; i++) {
        m.set(0, i, this.primitive.data.get(i).conjugate().mul(this.coeff));
      }
      return m;
    }
    // Column vector
    const dim = this.primitive.dim;
    const m = ComplexMatrix.zeros(dim, 1);
    for (let i = 0; i < dim; i++) {
      m.set(i, 0, this.primitive.data.get(i).mul(this.coeff));
    }
    return m;
  }

  adjoint() {
    return new StateFn(this.primitive, !this.is_measurement, this.coeff.conjugate());
  }

  // Apply an operator to the state
  compose(other) {
    if (other instanceof OperatorBase) {
      // |ψ> = O|ψ>
      const newState = this.primitive.evolve(other.to_matrix());
      return new StateFn(newState, this.is_measurement, this.coeff);
    }
    return super.compose(other);
  }

  toString() {
    const prefix = this.is_measurement ? "~" : "";
    return `${prefix}StateFn(${this.primitive.toString()})`;
  }
}

// ---------------------------------------------------------------------------
// CircuitStateFn: A state prepared by a circuit
// ---------------------------------------------------------------------------
export class CircuitStateFn extends StateFn {
  constructor(circuit, is_measurement = false, coeff = Complex.ONE) {
    const sv = Statevector.fromCircuit(circuit);
    super(sv, is_measurement, coeff);
    this.circuit = circuit;
  }

  to_circuit() { return this.circuit; }

  toString() { return `CircuitStateFn(${this.circuit.name})`; }
}

// ---------------------------------------------------------------------------
// ListOp: A list of OperatorBase with a combining operation
// ---------------------------------------------------------------------------
export class ListOp extends OperatorBase {
  constructor(oplist, comboFn = "add", auxFields = {}) {
    super();
    this.oplist = oplist;
    this.combo_fn = comboFn;
    this.aux_fields = auxFields;
    if (oplist.length > 0) this.num_qubits = oplist[0].num_qubits;
    this.coeff = auxFields.coeff || Complex.ONE;
  }

  to_matrix() {
    if (this.oplist.length === 0) return ComplexMatrix.identity(1);
    const matrices = this.oplist.map(op => op.to_matrix());
    let result;
    switch (this.combo_fn) {
      case "add":
        result = matrices[0];
        for (let i = 1; i < matrices.length; i++) result = result.add(matrices[i]);
        break;
      case "compose":
        result = matrices[0];
        for (let i = 1; i < matrices.length; i++) result = matrices[i].mul(result);
        break;
      case "tensor":
        result = matrices[0];
        for (let i = 1; i < matrices.length; i++) result = result.tensor(matrices[i]);
        break;
      case "adjoint":
        result = matrices[0].dagger();
        break;
      case "scale":
        result = matrices[0].scale(this.aux_fields.coeff || Complex.ONE);
        break;
      case "pow":
        result = matrices[0];
        const n = this.aux_fields.exponent || 1;
        for (let i = 1; i < n; i++) result = result.mul(matrices[0]);
        break;
      default:
        throw new Error(`ListOp: unknown combo function ${this.combo_fn}`);
    }
    return result.scale(this.coeff);
  }

  adjoint() {
    return new ListOp(
      this.oplist.map(op => op.adjoint()),
      this.combo_fn,
      Object.assign({}, this.aux_fields, { coeff: this.coeff.conjugate() })
    );
  }

  // Reduce: recursively simplify children
  reduce() {
    const reduced = this.oplist.map(op => op.reduce ? op.reduce() : op);
    return new ListOp(reduced, this.combo_fn, this.aux_fields);
  }

  toString() {
    const opStr = this.oplist.map(op => op.toString()).join(` ${this.combo_fn} `);
    return `(${opStr})`;
  }
}

// ---------------------------------------------------------------------------
// EvolvedOp: e^{-i t O} for an operator O (time evolution)
// ---------------------------------------------------------------------------
export class EvolvedOp extends OperatorBase {
  constructor(primitive, coeff = Complex.ONE, time = 1.0) {
    super();
    this.primitive = primitive;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.time = time;
    this.num_qubits = primitive.num_qubits;
  }

  to_matrix() {
    // e^{-i t O} = expm(-i t O)
    const O = this.primitive.to_matrix();
    const expMatrix = O.scale(new Complex(0, -this.time));
    return expMatrix.expm().scale(this.coeff);
  }

  adjoint() {
    return new EvolvedOp(this.primitive, this.coeff.conjugate(), -this.time);
  }

  toString() { return `EvolvedOp(e^(-i ${this.time} ${this.primitive}))`; }
}
