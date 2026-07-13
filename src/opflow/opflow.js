import { Complex, ComplexMatrix, ComplexVector } from "./../math/linalg.js";
import { Pauli, SparsePauliOp } from "./../quantum_info/pauli.js";
import { Operator } from "./../quantum_info/operator.js";
import { Statevector } from "./../quantum_info/statevector.js";

// OperatorBase
export class OperatorBase {
  constructor() {
    this.numQubits = 0;
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

  toMatrix() { throw new Error("OperatorBase.toMatrix not implemented"); }
  to_matrix_op() { return new MatrixOp(this.toMatrix()); }
  to_spmatrix_op() {
    return new PauliSumOp(SparsePauliOp.fromOperator(new Operator(this.toMatrix())));
  }

  // Apply to a state
  applyTo(state) {
    const m = this.toMatrix();
    return m.matvec(state.data || state);
  }

  isHermitian() { return this.toMatrix().isHermitian(); }
}

// PauliOp: A single Pauli operator
export class PauliOp extends OperatorBase {
  constructor(pauli, coeff = Complex.ONE) {
    super();
    this.primitive = pauli instanceof Pauli ? pauli : new Pauli(pauli);
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.numQubits = this.primitive.numQubits;
  }

  toMatrix() {
    return this.primitive.toMatrix().scale(this.coeff);
  }

  adjoint() {
    // Pauli^dagger = Pauli (up to sign for Y, but Y is self-adjoint up to phase)
    return new PauliOp(this.primitive, this.coeff.conjugate());
  }

  isHermitian() {
    return Math.abs(this.coeff.im) < 1e-9;
  }

  toString() { return `${this.coeff.toString()} * ${this.primitive.label}`; }
}

// PauliSumOp: Sum of Pauli operators (wraps SparsePauliOp)
export class PauliSumOp extends OperatorBase {
  constructor(spmatrixop, coeff = Complex.ONE) {
    super();
    this.primitive = spmatrixop instanceof SparsePauliOp ? spmatrixop : SparsePauliOp.fromList(spmatrixop);
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.numQubits = this.primitive.numQubits;
  }

  static fromList(list) {
    return new PauliSumOp(SparsePauliOp.fromList(list));
  }

  toMatrix() {
    return this.primitive.toMatrix().scale(this.coeff);
  }

  adjoint() {
    return new PauliSumOp(this.primitive.adjoint(), this.coeff.conjugate());
  }

  isHermitian() {
    return this.toMatrix().isHermitian();
  }

  // Reduce: combine like terms
  reduce() {
    return new PauliSumOp(this.primitive.simplify(), this.coeff);
  }

  // Convert to a list of [Pauli, coeff] pairs
  toList() { return this.primitive.toList(); }

  toString() { return this.primitive.toString(); }
}

// MatrixOp: An operator backed by a dense matrix
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
    this.numQubits = n;
  }

  toMatrix() {
    return this.primitive.scale(this.coeff);
  }

  adjoint() {
    return new MatrixOp(this.primitive.dagger(), this.coeff.conjugate());
  }

  isHermitian() {
    return this.toMatrix().isHermitian();
  }

  toString() { return `MatrixOp(numQubits=${this.numQubits})`; }
}

// CircuitOp: An operator backed by a QuantumCircuit
export class CircuitOp extends OperatorBase {
  constructor(circuit, coeff = Complex.ONE) {
    super();
    this.primitive = circuit;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.numQubits = circuit.numQubits;
  }

  toMatrix() {
    const op = Operator.fromCircuit(this.primitive);
    return op.data.scale(this.coeff);
  }

  toCircuit() { return this.primitive; }

  adjoint() {
    return new CircuitOp(this.primitive.inverse(), this.coeff.conjugate());
  }

  toString() { return `CircuitOp(${this.primitive.name})`; }
}

// StateFn: A state functional (linear functional <psi| or vector |psi>)
export class StateFn extends OperatorBase {
  constructor(primitive, isMeasurement = false, coeff = Complex.ONE) {
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
    this.isMeasurement = isMeasurement;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.numQubits = this.primitive.numQubits;
  }

  static fromLabel(label, isMeasurement = false) {
    return new StateFn(Statevector.fromLabel(label), isMeasurement);
  }

  toMatrix() {
    // If measurement: <psi| (row vector), else |psi> (column vector)
    if (this.isMeasurement) {
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
    return new StateFn(this.primitive, !this.isMeasurement, this.coeff.conjugate());
  }

  // Apply an operator to the state
  compose(other) {
    if (other instanceof OperatorBase) {
      // |ψ> = O|ψ>
      const newState = this.primitive.evolve(other.toMatrix());
      return new StateFn(newState, this.isMeasurement, this.coeff);
    }
    return super.compose(other);
  }

  toString() {
    const prefix = this.isMeasurement ? "~" : "";
    return `${prefix}StateFn(${this.primitive.toString()})`;
  }
}

// CircuitStateFn: A state prepared by a circuit
export class CircuitStateFn extends StateFn {
  constructor(circuit, isMeasurement = false, coeff = Complex.ONE) {
    const sv = Statevector.fromCircuit(circuit);
    super(sv, isMeasurement, coeff);
    this.circuit = circuit;
  }

  toCircuit() { return this.circuit; }

  toString() { return `CircuitStateFn(${this.circuit.name})`; }
}

// ListOp: A list of OperatorBase with a combining operation
export class ListOp extends OperatorBase {
  constructor(oplist, comboFn = "add", auxFields = {}) {
    super();
    this.oplist = oplist;
    this.comboFn = comboFn;
    this.auxFields = auxFields;
    if (oplist.length > 0) this.numQubits = oplist[0].numQubits;
    this.coeff = auxFields.coeff || Complex.ONE;
  }

  toMatrix() {
    if (this.oplist.length === 0) return ComplexMatrix.identity(1);
    const matrices = this.oplist.map(op => op.toMatrix());
    let result;
    switch (this.comboFn) {
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
        result = matrices[0].scale(this.auxFields.coeff || Complex.ONE);
        break;
      case "pow":
        result = matrices[0];
        const n = this.auxFields.exponent || 1;
        for (let i = 1; i < n; i++) result = result.mul(matrices[0]);
        break;
      default:
        throw new Error(`ListOp: unknown combo function ${this.comboFn}`);
    }
    return result.scale(this.coeff);
  }

  adjoint() {
    return new ListOp(
      this.oplist.map(op => op.adjoint()),
      this.comboFn,
      Object.assign({}, this.auxFields, { coeff: this.coeff.conjugate() })
    );
  }

  // Reduce: recursively simplify children
  reduce() {
    const reduced = this.oplist.map(op => op.reduce ? op.reduce() : op);
    return new ListOp(reduced, this.comboFn, this.auxFields);
  }

  toString() {
    const opStr = this.oplist.map(op => op.toString()).join(` ${this.comboFn} `);
    return `(${opStr})`;
  }
}

// EvolvedOp: e^{-i t O} for an operator O (time evolution)
export class EvolvedOp extends OperatorBase {
  constructor(primitive, coeff = Complex.ONE, time = 1.0) {
    super();
    this.primitive = primitive;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.time = time;
    this.numQubits = primitive.numQubits;
  }

  toMatrix() {
    // e^{-i t O} = expm(-i t O)
    const O = this.primitive.toMatrix();
    const expMatrix = O.scale(new Complex(0, -this.time));
    return expMatrix.expm().scale(this.coeff);
  }

  adjoint() {
    return new EvolvedOp(this.primitive, this.coeff.conjugate(), -this.time);
  }

  toString() { return `EvolvedOp(e^(-i ${this.time} ${this.primitive}))`; }
}
