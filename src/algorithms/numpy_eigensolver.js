/**
 * numpy_eigensolver.js - Exact classical eigensolver (qiskit.algorithms.NumPyMinimumEigensolver).
 *
 * Computes the exact minimum eigenvalue and eigenvector of a Hermitian
 * operator by diagonalizing its matrix representation. Used as a reference
 * for VQE benchmarking and for small-scale exact calculations.
 */

import { Operator } from "./../quantum_info/operator.js";
import { Statevector } from "./../quantum_info/statevector.js";
import { SparsePauliOp } from "./../quantum_info/pauli.js";

export class NumPyMinimumEigensolverResult {
  constructor(kwargs = {}) {
    this.eigenvalue = kwargs.eigenvalue || 0;
    this.eigenstate = kwargs.eigenstate || null;
    this.aux_operators_evaluated = kwargs.aux_operators_evaluated || [];
  }
}

export class NumPyMinimumEigensolver {
  constructor(options = {}) {
    this.aux_operators = options.aux_operators || null;
  }

  // Compute the minimum eigenvalue of `operator` (a SparsePauliOp, Operator,
  // or any object with a to_matrix() method).
  compute_minimum_eigenvalue(operator, auxOperators = null) {
    // Convert to a matrix.
    let mat;
    if (operator instanceof SparsePauliOp) {
      mat = operator.to_matrix();
    } else if (operator instanceof Operator) {
      mat = operator._data;
    } else if (typeof operator.to_matrix === "function") {
      mat = operator.to_matrix();
    } else if (operator instanceof ComplexMatrix) {
      mat = operator;
    } else {
      throw new TypeError("NumPyMinimumEigensolver: operator must have to_matrix()");
    }
    // Diagonalize.
    const { eigenvalues, eigenvectors } = mat.eigh();
    // Find the minimum eigenvalue.
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
        // ComplexVector constructor.
        const ComplexVector = globalThis.__ketraComplexVector;
        if (ComplexVector) return new ComplexVector(eigenvecData);
        // Fall back: require from linalg.
        throw new Error("ComplexVector not registered");
      })(),
      nq,
    );
    // Evaluate aux operators if provided.
    const aux = auxOperators || this.aux_operators;
    const auxResults = [];
    if (aux) {
      for (const op of aux) {
        if (op instanceof SparsePauliOp) {
          auxResults.push([op, op.expectation_value(eigenstate)]);
        } else {
          auxResults.push([op, null]);
        }
      }
    }
    return new NumPyMinimumEigensolverResult({
      eigenvalue,
      eigenstate,
      aux_operators_evaluated: auxResults,
    });
  }
}

// Maximum eigensolver (alias for NumPyMinimumEigensolver with sign flip).
export class NumPyMaximumEigensolver extends NumPyMinimumEigensolver {
  compute_maximum_eigenvalue(operator, auxOperators = null) {
    // Compute minimum of -operator, then negate.
    const negOp = operator instanceof SparsePauliOp
      ? new SparsePauliOp(operator.paulis, operator.coeffs.map(c => c.scale(-1)))
      : operator;
    const result = super.compute_minimum_eigenvalue(negOp, auxOperators);
    return new NumPyMinimumEigensolverResult({
      eigenvalue: -result.eigenvalue,
      eigenstate: result.eigenstate,
      aux_operators_evaluated: result.aux_operators_evaluated,
    });
  }
}

// Lazy import of ComplexVector (avoid circular import at module load time).
import { ComplexVector } from "./../math/linalg.js";
// Register so the closure above can find it.
globalThis.__ketraComplexVector = ComplexVector;
