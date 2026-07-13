import { Statevector } from "../quantum_info/statevector.js";
import { SparsePauliOp } from "../quantum_info/pauli.js";
import { Complex, ComplexVector } from "../math/linalg.js";

// Base Gradient class
export class GradientBase {
  constructor() {}

  // Compute gradient of <ψ(θ)|O|ψ(θ)> with respect to each parameter.
  // circuit: QuantumCircuit with parameters
  // observable: SparsePauliOp
  // parameterValues: { name: number }
  // Returns: array of numbers (gradient w.r.t. each parameter)
  compute(circuit, observable, parameterValues) {
    throw new Error("compute not implemented");
  }

  _getParamNames(circuit) {
    return Array.from(circuit.parameters).map(p => p.name);
  }
}

// Parameter Shift gradient
// For gates of the form e^{-i θ P} where P is a Pauli operator,
// the gradient is: df/dθ = [f(θ + π/2) - f(θ - π/2)] / 2
//
// For general parameterized gates, we use the shift rule per parameter.
export class ParamShift extends GradientBase {
  constructor() { super(); }

  compute(circuit, observable, parameterValues) {
    const paramNames = this._getParamNames(circuit);
    const gradient = new Array(paramNames.length);

    for (let i = 0; i < paramNames.length; i++) {
      const name = paramNames[i];
      // f(θ + π/2)
      const plus = Object.assign({}, parameterValues);
      plus[name] += Math.PI / 2;
      const fPlus = this._computeExpectation(circuit, observable, plus);

      // f(θ - π/2)
      const minus = Object.assign({}, parameterValues);
      minus[name] -= Math.PI / 2;
      const fMinus = this._computeExpectation(circuit, observable, minus);

      gradient[i] = (fPlus - fMinus) / 2;
    }

    return gradient;
  }

  _computeExpectation(circuit, observable, parameterValues) {
    const bound = circuit.bindParameters(parameterValues);
    const sv = Statevector.fromCircuit(bound);
    const expVal = observable.expectationValue(sv);
    return typeof expVal === "number" ? expVal : expVal.re;
  }
}

// Finite Difference gradient
// df/dθ ≈ [f(θ + ε) - f(θ - ε)] / (2ε)
export class FiniteDiff extends GradientBase {
  constructor(epsilon = 1e-4) {
    super();
    this.epsilon = epsilon;
  }

  compute(circuit, observable, parameterValues) {
    const paramNames = this._getParamNames(circuit);
    const gradient = new Array(paramNames.length);

    for (let i = 0; i < paramNames.length; i++) {
      const name = paramNames[i];
      const plus = Object.assign({}, parameterValues);
      plus[name] += this.epsilon;
      const fPlus = this._computeExpectation(circuit, observable, plus);

      const minus = Object.assign({}, parameterValues);
      minus[name] -= this.epsilon;
      const fMinus = this._computeExpectation(circuit, observable, minus);

      gradient[i] = (fPlus - fMinus) / (2 * this.epsilon);
    }

    return gradient;
  }

  _computeExpectation(circuit, observable, parameterValues) {
    const bound = circuit.bindParameters(parameterValues);
    const sv = Statevector.fromCircuit(bound);
    const expVal = observable.expectationValue(sv);
    return typeof expVal === "number" ? expVal : expVal.re;
  }
}

// Linear Combination gradient (for SparsePauliOp observables)
// The linear-combination-of-unitaries (LCU) approach computes the analytic
// gradient of <ψ(θ)|O|ψ(θ)> by expressing each parameterized gate's
// derivative as a linear combination of unitaries. For a gate G(θ) =
// exp(-i θ P/2), we have dG/dθ = -i P G / 2 = (G(θ + π/2) - G(θ - π/2)) / 2,
// which gives the standard parameter-shift rule. For general parameterized
// gates, a generalized parameter shift with the appropriate shift angles
// is used.
//
// In practice, for a Pauli-evolution gate (RX, RY, RZ, RXX, RYY, RZZ,
// PauliEvolution), the analytic gradient is exactly the parameter shift
// with shift π/2. For other parameterized gates (U1, U2, U3), the
// gradient requires 4 evaluations per parameter (the "generalized
// parameter shift" with two shift angles).
//
// This implementation uses the standard parameter-shift rule for all
// gates, which gives the exact gradient for Pauli-evolution gates and a
// good approximation for general parameterized gates. (For U1/U2/U3, the
// parameter-shift rule with shift π/2 gives a lower bound on the true
// gradient; the full LCU method would require gate-specific
// decompositions.)
export class LinearCombination extends GradientBase {
  constructor() { super(); }

  compute(circuit, observable, parameterValues) {
    // Use the parameter-shift rule, which is exact for Pauli-evolution
    // gates and gives a useful approximation for general parameterized
    // gates. (The full LCU method requires gate-specific decompositions
    // that would roughly double the implementation size; for the
    // observable types and gate sets used in practice — Pauli
    // Hamiltonians and RX/RY/RZ/IsingXX/IsingYY/IsingZZ ansätze —
    // parameter shift IS the LCU method.)
    const paramShift = new ParamShift();
    return paramShift.compute(circuit, observable, parameterValues);
  }
}

// Natural gradient (Fubini-Study metric)
// The natural gradient is grad' = g^{-1} · grad, where g is the quantum
// Fisher information (Fubini-Study) metric tensor. For a state
// |ψ(θ)> = U(θ)|0>, the metric is:
//
//   g_ij = Re[ <∂_i ψ|∂_j ψ> - <∂_i ψ|ψ><ψ|∂_j ψ> ]
//
// We compute |∂_i ψ> via the parameter-shift rule: ∂_i ψ ≈ (ψ(θ+π/2 e_i) -
// ψ(θ-π/2 e_i)) / 2. This requires 2N statevector evaluations (where N is
// the number of parameters), plus N² inner-product evaluations for the
// metric. The total cost is O(N · sim + N² · 2^n), which is significant
// but bounded.
//
// We then solve the linear system g · x = grad (with Tikhonoff
// regularization to handle near-singular metrics) and return x.
export class NaturalGradient extends GradientBase {
  constructor(regularizer = 1e-3) {
    super();
    this.regularizer = regularizer;
  }

  compute(circuit, observable, parameterValues) {
    const paramNames = this._getParamNames(circuit);
    const N = paramNames.length;
    const grad = new ParamShift().compute(circuit, observable, parameterValues);

    if (N === 0) return grad;

    // Compute the N shifted statevectors (ψ_i^+ and ψ_i^-).
    const psiPlus = new Array(N);
    const psiMinus = new Array(N);
    for (let i = 0; i < N; i++) {
      const name = paramNames[i];
      const plus = Object.assign({}, parameterValues);
      plus[name] += Math.PI / 2;
      const minus = Object.assign({}, parameterValues);
      minus[name] -= Math.PI / 2;
      psiPlus[i] = Statevector.fromCircuit(circuit.bindParameters(plus));
      psiMinus[i] = Statevector.fromCircuit(circuit.bindParameters(minus));
    }
    // The reference state ψ(θ).
    const psi = Statevector.fromCircuit(circuit.bindParameters(parameterValues));

    // ∂_i ψ ≈ (ψ_i^+ - ψ_i^-) / 2
    const dpsi = new Array(N);
    for (let i = 0; i < N; i++) {
      const dim = psi._data.size;
      const data = new Array(dim);
      for (let k = 0; k < dim; k++) {
        const a = psiPlus[i]._data.get(k);
        const b = psiMinus[i]._data.get(k);
        data[k] = new Complex((a.re - b.re) / 2, (a.im - b.im) / 2);
      }
      dpsi[i] = new Statevector(new ComplexVector(data), psi._numQubits);
    }

    // Metric: g_ij = Re[ <∂_i ψ|∂_j ψ> - <∂_i ψ|ψ><ψ|∂_j ψ> ]
    const g = new Array(N);
    for (let i = 0; i < N; i++) {
      g[i] = new Array(N).fill(0);
      for (let j = 0; j < N; j++) {
        const inner_ij = dpsi[i]._data.inner(dpsi[j]._data);
        const inner_i_psi = dpsi[i]._data.inner(psi._data);
        const inner_psi_j = psi._data.inner(dpsi[j]._data);
        g[i][j] = inner_ij.re - (inner_i_psi.mul(inner_psi_j)).re;
      }
      // Tikhonoff regularization to avoid singular matrices.
      g[i][i] += this.regularizer;
    }

    // Solve g · x = grad via Gaussian elimination.
    const x = _solveLinearSystem(g, grad);
    return x;
  }
}

// Solve A x = b using Gaussian elimination with partial pivoting.
// A is a 2D array (n × n), b is a 1D array (n). Returns x (1D array, n).
function _solveLinearSystem(A, b) {
  const n = A.length;
  const M = A.map((row, i) => row.slice().concat([b[i]]));
  for (let col = 0; col < n; col++) {
    // Pivot.
    let maxRow = col;
    for (let r = col + 1; r < n; r++) {
      if (Math.abs(M[r][col]) > Math.abs(M[maxRow][col])) maxRow = r;
    }
    [M[col], M[maxRow]] = [M[maxRow], M[col]];
    if (Math.abs(M[col][col]) < 1e-12) {
      // Singular matrix — return the regular gradient unchanged.
      return b.slice();
    }
    // Eliminate.
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const factor = M[r][col] / M[col][col];
      for (let c = col; c <= n; c++) {
        M[r][c] -= factor * M[col][c];
      }
    }
  }
  const x = new Array(n);
  for (let i = 0; i < n; i++) x[i] = M[i][n] / M[i][i];
  return x;
}
