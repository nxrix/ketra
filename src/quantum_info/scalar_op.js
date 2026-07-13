import { Complex, ComplexMatrix } from "./../math/linalg.js";
import { Operator } from "./operator.js";

export class ScalarOp {
  constructor(numQubits, coeff = null) {
    this._numQubits = numQubits;
    if (coeff instanceof Complex) {
      this._coeff = coeff;
    } else if (typeof coeff === "number") {
      this._coeff = new Complex(coeff, 0);
    } else if (coeff == null) {
      this._coeff = Complex.ONE;
    } else {
      throw new TypeError("ScalarOp coeff must be a number or Complex");
    }
  }

  static fromOperator(operator) {
    // Detect if operator is c*I
    const dim = operator.dim;
    const c0 = operator.data.get(0, 0);
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        if (i === j) {
          if (!operator.data.get(i, i).equals(c0, 1e-9)) {
            throw new Error("Operator is not a scalar multiple of identity");
          }
        } else if (operator.data.get(i, j).abs() > 1e-9) {
          throw new Error("Operator is not a scalar multiple of identity");
        }
      }
    }
    return new ScalarOp(operator.numQubits, c0);
  }

  get numQubits() { return this._numQubits; }
  get coeff() { return this._coeff; }
  get dim() { return 1 << this._numQubits; }

  toMatrix() {
    return ComplexMatrix.identity(this.dim).scale(this._coeff);
  }

  toOperator() {
    return new Operator(this.toMatrix());
  }

  compose(other) {
    if (other instanceof ScalarOp) {
      return new ScalarOp(this._numQubits, this._coeff.mul(other._coeff));
    }
    return other.scale(this._coeff);
  }

  tensor(other) {
    if (other instanceof ScalarOp) {
      return new ScalarOp(this._numQubits + other._numQubits, this._coeff.mul(other._coeff));
    }
    // ScalarOp ⊗ O = coeff * O (with extended identity)
    return other.scale(this._coeff);
  }

  expand(other) {
    return this.tensor(other);
  }

  adjoint() { return new ScalarOp(this._numQubits, this._coeff.conjugate()); }
  conjugate() { return new ScalarOp(this._numQubits, this._coeff.conjugate()); }
  transpose() { return this; }

  trace() { return this._coeff.scale(this.dim); }
  det() { return this._coeff.pow(this.dim); }

  isUnitary(tol = 1e-9) {
    return Math.abs(this._coeff.abs() - 1) < tol;
  }

  equals(other, tol) {
    if (!(other instanceof ScalarOp)) return false;
    if (this._numQubits !== other._numQubits) return false;
    return this._coeff.equals(other._coeff, tol);
  }

  applyToVector(vector) {
    return vector.scale(this._coeff);
  }

  expectationValue(statevector) {
    const norm2 = statevector.norm() ** 2;
    return this._coeff.scale(norm2);
  }

  toString() {
    return `ScalarOp(${this._coeff.toString()}, numQubits=${this._numQubits})`;
  }
}
