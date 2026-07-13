import { Operator } from "./../quantum_info/operator.js";
import { Statevector } from "./../quantum_info/statevector.js";
import { SparsePauliOp } from "./../quantum_info/pauli.js";

export class NumPyMinimumEigensolverResult {
  constructor(kwargs = {}) {
    this.eigenvalue = kwargs.eigenvalue || 0;
    this.eigenstate = kwargs.eigenstate || null;
    this.auxOperatorsEvaluated = kwargs.auxOperatorsEvaluated || [];
  }
}

export class NumPyMinimumEigensolver {
  constructor(options = {}) {
    this.auxOperators = options.auxOperators || null;
  }

  // Compute the minimum eigenvalue of `operator` (a SparsePauliOp, Operator,
  // or any object with a toMatrix() method).
  computeMinimumEigenvalue(operator, auxOperators = null) {
    let mat;
    if (operator instanceof SparsePauliOp) {
      mat = operator.toMatrix();
    } else if (operator instanceof Operator) {
      mat = operator._data;
    } else if (typeof operator.toMatrix === "function") {
      mat = operator.toMatrix();
    } else if (operator instanceof ComplexMatrix) {
      mat = operator;
    } else {
      throw new TypeError("NumPyMinimumEigensolver: operator must have toMatrix()");
    }
    // Diagonalize.
    const { eigenvalues, eigenvectors } = mat.eigh();
    let minIdx = 0;
    for (let i = 1; i < eigenvalues.length; i++) {
      if (eigenvalues[i] < eigenvalues[minIdx]) minIdx = i;
    }
    const eigenvalue = eigenvalues[minIdx];
    // Extract the eigenvector (column minIdx of eigenvectors).
    const dim = mat.rows;
    const nq = Math.log2(dim);
    const eigenvecData = new Array(dim);
    for (let i = 0; i < dim; i++) {
      eigenvecData[i] = eigenvectors.get(i, minIdx);
    }
    const eigenstate = new Statevector(
      // Lazy import to avoid circular deps.
      (function() {
        const ComplexVector = globalThis.__ketraComplexVector;
        if (ComplexVector) return new ComplexVector(eigenvecData);
        throw new Error("ComplexVector not registered");
      })(),
      nq,
    );
    const aux = auxOperators || this.auxOperators;
    const auxResults = [];
    if (aux) {
      for (const op of aux) {
        if (op instanceof SparsePauliOp) {
          auxResults.push([op, op.expectationValue(eigenstate)]);
        } else {
          auxResults.push([op, null]);
        }
      }
    }
    return new NumPyMinimumEigensolverResult({
      eigenvalue,
      eigenstate,
      auxOperatorsEvaluated: auxResults,
    });
  }
}

// Maximum eigensolver (alias for NumPyMinimumEigensolver with sign flip).
export class NumPyMaximumEigensolver extends NumPyMinimumEigensolver {
  computeMaximumEigenvalue(operator, auxOperators = null) {
    // Compute minimum of -operator, then negate.
    const negOp = operator instanceof SparsePauliOp
      ? new SparsePauliOp(operator.paulis, operator.coeffs.map(c => c.scale(-1)))
      : operator;
    const result = super.computeMinimumEigenvalue(negOp, auxOperators);
    return new NumPyMinimumEigensolverResult({
      eigenvalue: -result.eigenvalue,
      eigenstate: result.eigenstate,
      auxOperatorsEvaluated: result.auxOperatorsEvaluated,
    });
  }
}

// Lazy import of ComplexVector (avoid circular import at module load time).
import { ComplexVector } from "./../math/linalg.js";
globalThis.__ketraComplexVector = ComplexVector;
