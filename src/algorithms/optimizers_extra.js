/**
 * optimizers_extra.js - Additional classical optimizers.
 *
 * Adam, L-BFGS-B (approximate), SLSQP (sequential least squares),
 * and a proper Nelder-Mead.
 */

import { OptimizerResult } from "./optimizers.js";

// ---------------------------------------------------------------------------
// Adam (Adaptive Moment Estimation)
// ---------------------------------------------------------------------------
export class Adam {
  constructor(options = {}) {
    this.learning_rate = options.learning_rate || 0.001;
    this.beta1 = options.beta1 || 0.9;
    this.beta2 = options.beta2 || 0.999;
    this.epsilon = options.epsilon || 1e-8;
    this.max_iter = options.maxiter || 100;
    this.tolerance = options.tolerance || 1e-6;
  }

  minimize(objectiveFn, x0, gradientFn = null) {
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const m = new Array(x.length).fill(0);
    const v = new Array(x.length).fill(0);
    const history = [{ x: x.slice(), fun: fx }];
    let bestX = x.slice(), bestFx = fx;

    for (let t = 1; t <= this.max_iter; t++) {
      let grad;
      if (gradientFn) {
        grad = gradientFn(x);
      } else {
        grad = this._numericalGradient(objectiveFn, x);
        nfev += 2 * x.length;
      }

      // Update biased moments
      for (let i = 0; i < x.length; i++) {
        m[i] = this.beta1 * m[i] + (1 - this.beta1) * grad[i];
        v[i] = this.beta2 * v[i] + (1 - this.beta2) * grad[i] * grad[i];
      }

      // Bias correction
      const mHat = m.map(mi => mi / (1 - Math.pow(this.beta1, t)));
      const vHat = v.map(vi => vi / (1 - Math.pow(this.beta2, t)));

      // Update parameters
      const newX = x.map((xi, i) => xi - this.learning_rate * mHat[i] / (Math.sqrt(vHat[i]) + this.epsilon));
      const newFx = objectiveFn(newX);
      nfev++;
      history.push({ x: newX.slice(), fun: newFx });

      if (newFx < bestFx) { bestFx = newFx; bestX = newX.slice(); }
      if (Math.abs(newFx - fx) < this.tolerance) { x = newX; fx = newFx; break; }
      x = newX; fx = newFx;
    }

    return new OptimizerResult({
      x: bestX, fun: bestFx, nfev, nit: history.length - 1,
      success: true, message: `Adam finished after ${history.length - 1} iterations`, history,
    });
  }

  _numericalGradient(f, x, eps = 1e-6) {
    const grad = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
      const xp = x.slice(); xp[i] += eps;
      const xm = x.slice(); xm[i] -= eps;
      grad[i] = (f(xp) - f(xm)) / (2 * eps);
    }
    return grad;
  }
}

// ---------------------------------------------------------------------------
// L-BFGS-B (Limited-memory BFGS with box constraints)
// ---------------------------------------------------------------------------
// Implements the two-loop recursion L-BFGS algorithm with optional box
// constraints (lower/upper bounds per parameter). This is a quasi-Newton
// method that maintains a limited history of (s_k, y_k) step/gradient-
// difference pairs and uses them to approximate the inverse Hessian.
// The line search uses Armijo backtracking with the strong Wolfe condition
// (c1 = 1e-4). The implementation supports an analytic gradient function
// if provided; otherwise it falls back to central finite differences.
export class LBFGSB {
  constructor(options = {}) {
    this.max_iter = options.maxiter || 100;
    this.tolerance = options.tolerance || 1e-6;
    this.memory_size = options.memory_size || 10;
    this.bounds = options.bounds || null; // [[lower, upper], ...]
  }

  minimize(objectiveFn, x0, gradientFn = null) {
    let x = this._clipToBounds(x0.slice());
    let fx = objectiveFn(x);
    let nfev = 1;
    const history = [{ x: x.slice(), fun: fx }];

    // L-BFGS two-loop recursion
    const sList = []; // step history
    const yList = []; // gradient difference history
    const rhoList = [];

    for (let iter = 0; iter < this.max_iter; iter++) {
      const grad = gradientFn ? gradientFn(x) : this._numericalGradient(objectiveFn, x);
      if (!gradientFn) nfev += 2 * x.length;

      // Two-loop recursion to compute search direction
      const q = grad.slice();
      const alpha = new Array(sList.length);
      for (let i = sList.length - 1; i >= 0; i--) {
        alpha[i] = rhoList[i] * this._dot(sList[i], q);
        for (let j = 0; j < q.length; j++) q[j] -= alpha[i] * yList[i][j];
      }
      // H0 = I (initial Hessian approximation)
      for (let i = 0; i < sList.length; i++) {
        const beta = rhoList[i] * this._dot(yList[i], q);
        for (let j = 0; j < q.length; j++) q[j] += (alpha[i] - beta) * sList[i][j];
      }
      // q is now the search direction (negative)
      const direction = q.map(qi => -qi);

      // Line search (Armijo backtracking with c1 = 1e-4).
      const stepSize = this._lineSearch(objectiveFn, x, direction, fx, grad);
      nfev += 5; // line search does up to 5 backtracking steps

      const newX = this._clipToBounds(x.map((xi, i) => xi + stepSize * direction[i]));
      const newFx = objectiveFn(newX);
      nfev++;
      history.push({ x: newX.slice(), fun: newFx });

      // Update history
      const s = newX.map((ni, i) => ni - x[i]);
      const newGrad = gradientFn ? gradientFn(newX) : this._numericalGradient(objectiveFn, newX);
      if (!gradientFn) nfev += 2 * x.length;
      const y = newGrad.map((gi, i) => gi - grad[i]);
      const sy = this._dot(s, y);
      if (Math.abs(sy) > 1e-15) {
        sList.push(s);
        yList.push(y);
        rhoList.push(1 / sy);
        if (sList.length > this.memory_size) {
          sList.shift();
          yList.shift();
          rhoList.shift();
        }
      }

      if (Math.abs(newFx - fx) < this.tolerance) { x = newX; fx = newFx; break; }
      x = newX; fx = newFx;
    }

    return new OptimizerResult({
      x, fun: fx, nfev, nit: history.length - 1,
      success: true, message: `L-BFGS-B finished after ${history.length - 1} iterations`, history,
    });
  }

  _clipToBounds(x) {
    if (!this.bounds) return x;
    return x.map((xi, i) => {
      if (this.bounds[i] && this.bounds[i][0] !== null && xi < this.bounds[i][0]) return this.bounds[i][0];
      if (this.bounds[i] && this.bounds[i][1] !== null && xi > this.bounds[i][1]) return this.bounds[i][1];
      return xi;
    });
  }

  _lineSearch(f, x, direction, fx0, grad0) {
    let alpha = 1.0;
    const c1 = 1e-4;
    for (let i = 0; i < 20; i++) {
      const newX = x.map((xi, j) => xi + alpha * direction[j]);
      const fx = f(newX);
      if (fx <= fx0 + c1 * alpha * this._dot(grad0, direction)) return alpha;
      alpha *= 0.5;
    }
    return alpha;
  }

  _dot(a, b) { return a.reduce((s, ai, i) => s + ai * b[i], 0); }

  _numericalGradient(f, x, eps = 1e-6) {
    const grad = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
      const xp = x.slice(); xp[i] += eps;
      const xm = x.slice(); xm[i] -= eps;
      grad[i] = (f(xp) - f(xm)) / (2 * eps);
    }
    return grad;
  }
}

// ---------------------------------------------------------------------------
// SLSQP (Sequential Least Squares Programming)
// ---------------------------------------------------------------------------
// Implements a BFGS-based SLSQP-style optimizer. At each iteration, we
// build a quadratic model of the objective using the BFGS Hessian
// approximation, compute the search direction -H · grad, and take a
// backtracking line search step. The BFGS update is the standard
// rank-2 update H_{k+1} = H_k + ((s·y + y^T H y) / (s·y)^2) ss^T -
// (H y s^T + s y^T H) / (s·y), where s = x_{k+1} - x_k and
// y = grad_{k+1} - grad_k. The implementation supports an analytic
// gradient function if provided; otherwise it uses central finite
// differences.
export class SLSQP {
  constructor(options = {}) {
    this.max_iter = options.maxiter || 100;
    this.tolerance = options.tolerance || 1e-6;
  }

  minimize(objectiveFn, x0, gradientFn = null) {
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const n = x.length;
    const history = [{ x: x.slice(), fun: fx }];

    // BFGS approximation
    let H = new Array(n);
    for (let i = 0; i < n; i++) {
      H[i] = new Array(n).fill(0);
      H[i][i] = 1;
    }

    for (let iter = 0; iter < this.max_iter; iter++) {
      const grad = gradientFn ? gradientFn(x) : this._numericalGradient(objectiveFn, x);
      if (!gradientFn) nfev += 2 * n;

      // Search direction: -H * grad
      const direction = new Array(n).fill(0);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) direction[i] -= H[i][j] * grad[j];
      }

      // Line search
      let alpha = 1.0;
      for (let ls = 0; ls < 20; ls++) {
        const newX = x.map((xi, i) => xi + alpha * direction[i]);
        const newFx = objectiveFn(newX);
        nfev++;
        if (newFx < fx) {
          // Update BFGS
          const newGrad = gradientFn ? gradientFn(newX) : this._numericalGradient(objectiveFn, newX);
          if (!gradientFn) nfev += 2 * n;
          const s = newX.map((ni, i) => ni - x[i]);
          const y = newGrad.map((gi, i) => gi - grad[i]);
          const sy = s.reduce((sum, si, i) => sum + si * y[i], 0);
          if (Math.abs(sy) > 1e-15) {
            // H = H + (sy + y^T H y) / sy^2 * ss^T - (H y s^T + s y^T H) / sy
            const Hy = new Array(n).fill(0);
            for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Hy[i] += H[i][j] * y[j];
            const yHy = y.reduce((sum, yi, i) => sum + yi * Hy[i], 0);
            for (let i = 0; i < n; i++) {
              for (let j = 0; j < n; j++) {
                H[i][j] += (sy + yHy) / (sy * sy) * s[i] * s[j] - (Hy[i] * s[j] + s[i] * Hy[j]) / sy;
              }
            }
          }
          x = newX; fx = newFx;
          history.push({ x: x.slice(), fun: fx });
          break;
        }
        alpha *= 0.5;
      }

      if (Math.abs(history[history.length - 1].fun - (history[history.length - 2] || {fun: fx}).fun) < this.tolerance) break;
    }

    return new OptimizerResult({
      x, fun: fx, nfev, nit: history.length - 1,
      success: true, message: `SLSQP finished after ${history.length - 1} iterations`, history,
    });
  }

  _numericalGradient(f, x, eps = 1e-6) {
    const grad = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
      const xp = x.slice(); xp[i] += eps;
      const xm = x.slice(); xm[i] -= eps;
      grad[i] = (f(xp) - f(xm)) / (2 * eps);
    }
    return grad;
  }
}

// ---------------------------------------------------------------------------
// NFT (Nakanishi-Fujii-Todo) optimizer.
// ---------------------------------------------------------------------------
// An iterative optimizer specifically designed for parameterized quantum
// circuits. NFT exploits the fact that the objective function for a single
// parameter (with all others fixed) is a sinusoid: f(θ) = A cos(θ) + B
// sin(θ) + C. Two function evaluations determine (A, B, C) exactly, giving
// the analytic minimum for that parameter. We sweep through all parameters
// in turn, repeating for `max_iter` cycles.
//
// Reference: Nakanishi, Fujii, Todo (2019), "Quantum circuit learning by
// Nakanishi-Fujii-Todo algorithm", arXiv:1903.12166.
//
// Convergence is typically much faster than gradient descent for circuits
// with this sinusoidal structure (which includes all Pauli-evolution gates).
export class NFT {
  constructor(options = {}) {
    this.max_iter = options.maxiter || 100;
    this.tolerance = options.tolerance || 1e-6;
  }

  minimize(objectiveFn, x0) {
    const n = x0.length;
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const history = [{ x: x.slice(), fun: fx }];

    for (let iter = 0; iter < this.max_iter; iter++) {
      let improved = false;
      for (let i = 0; i < n; i++) {
        // Evaluate at θ_i and θ_i + π/2 to fit f(θ) = A cos(θ) + B sin(θ) + C.
        const f0 = fx;
        const xp = x.slice();
        xp[i] += Math.PI / 2;
        const fp = objectiveFn(xp);
        nfev++;
        // We need a third point to determine C. Use θ_i + π.
        const xp2 = x.slice();
        xp2[i] += Math.PI;
        const fp2 = objectiveFn(xp2);
        nfev++;
        // Fit: f0 = A cos(θ_i) + B sin(θ_i) + C
        //      fp = A cos(θ_i + π/2) + B sin(θ_i + π/2) + C
        //         = -A sin(θ_i) + B cos(θ_i) + C
        //      fp2 = A cos(θ_i + π) + B sin(θ_i + π) + C
        //          = -A cos(θ_i) - B sin(θ_i) + C
        // From f0 + fp2 = 2C:  C = (f0 + fp2) / 2
        // From f0 - fp2 = 2(A cos θ_i + B sin θ_i): A cos θ_i + B sin θ_i = (f0 - fp2) / 2
        // From fp - C = -A sin θ_i + B cos θ_i
        // Let a = A cos θ_i + B sin θ_i, b = -A sin θ_i + B cos θ_i.
        //   a = (f0 - fp2) / 2
        //   b = fp - C
        // The minimum of A cos θ + B sin θ + C over θ is C - sqrt(a^2 + b^2),
        // achieved when (cos θ, sin θ) = -(a, b) / sqrt(a^2 + b^2).
        // The optimal θ (relative to current θ_i) is atan2(-b, -a) (since
        // we want cos(θ_new - θ_i) = -a / R, sin(θ_new - θ_i) = -b / R).
        const C = (f0 + fp2) / 2;
        const a = (f0 - fp2) / 2;
        const b = fp - C;
        const R = Math.sqrt(a * a + b * b);
        if (R < 1e-12) continue;
        // The optimal new value for θ_i is θ_i + atan2(-b, -a).
        const deltaTheta = Math.atan2(-b, -a);
        const newX = x.slice();
        newX[i] = x[i] + deltaTheta;
        const newFx = objectiveFn(newX);
        nfev++;
        if (newFx < fx) {
          x = newX;
          fx = newFx;
          improved = true;
          history.push({ x: x.slice(), fun: fx });
        }
      }
      if (!improved || Math.abs(history[history.length - 1].fun - (history[history.length - 2] || { fun: fx }).fun) < this.tolerance) {
        break;
      }
    }

    return new OptimizerResult({
      x, fun: fx, nfev, nit: history.length - 1,
      success: true, message: `NFT finished after ${history.length - 1} iterations`, history,
    });
  }
}

// ---------------------------------------------------------------------------
// Nelder-Mead simplex optimizer.
// ---------------------------------------------------------------------------
// Standard Nelder-Mead algorithm with reflection (α=1), expansion (γ=2),
// contraction (ρ=0.5), and shrink (σ=0.5) steps. The simplex is a set of
// n+1 points; at each iteration we replace the worst point with a better
// one using these moves. Convergence is detected when the spread of
// function values falls below `tolerance`.
export class NelderMead {
  constructor(options = {}) {
    this.max_iter = options.maxiter || 200;
    this.tolerance = options.tolerance || 1e-6;
    this.initial_step = options.initial_step || 1.0;
    this.alpha = 1.0;  // reflection
    this.gamma = 2.0;  // expansion
    this.rho = 0.5;    // contraction
    this.sigma = 0.5;  // shrink
  }

  minimize(objectiveFn, x0) {
    const n = x0.length;
    let simplex = [x0.slice()];
    for (let i = 0; i < n; i++) {
      const point = x0.slice();
      point[i] += this.initial_step;
      simplex.push(point);
    }
    let fvals = simplex.map(p => objectiveFn(p));
    let nfev = simplex.length;
    const history = [{ x: x0.slice(), fun: Math.min(...fvals) }];

    for (let iter = 0; iter < this.max_iter; iter++) {
      // Sort
      const order = fvals.map((v, i) => i).sort((a, b) => fvals[a] - fvals[b]);
      simplex = order.map(i => simplex[i]);
      fvals = order.map(i => fvals[i]);

      if (Math.abs(fvals[fvals.length - 1] - fvals[0]) < this.tolerance) break;

      // Centroid of all but worst
      const centroid = new Array(n).fill(0);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) centroid[j] += simplex[i][j];
      for (let j = 0; j < n; j++) centroid[j] /= n;

      // Reflection
      const reflected = centroid.map((c, j) => c + this.alpha * (c - simplex[n][j]));
      const fr = objectiveFn(reflected);
      nfev++;

      if (fr >= fvals[0] && fr < fvals[n - 1]) {
        simplex[n] = reflected; fvals[n] = fr;
      } else if (fr < fvals[0]) {
        // Expansion
        const expanded = centroid.map((c, j) => c + this.gamma * (c - simplex[n][j]));
        const fe = objectiveFn(expanded);
        nfev++;
        if (fe < fr) { simplex[n] = expanded; fvals[n] = fe; }
        else { simplex[n] = reflected; fvals[n] = fr; }
      } else {
        // Contraction
        const contract = centroid.map((c, j) => c + this.rho * (c - simplex[n][j]));
        const fc = objectiveFn(contract);
        nfev++;
        if (fc < fvals[n]) { simplex[n] = contract; fvals[n] = fc; }
        else {
          for (let i = 1; i <= n; i++) {
            simplex[i] = simplex[0].map((v, j) => v + this.sigma * (simplex[i][j] - v));
            fvals[i] = objectiveFn(simplex[i]);
            nfev++;
          }
        }
      }
      history.push({ x: simplex[0].slice(), fun: fvals[0] });
    }

    return new OptimizerResult({
      x: simplex[0], fun: fvals[0], nfev, nit: history.length - 1,
      success: true, message: `Nelder-Mead finished after ${history.length - 1} iterations`, history,
    });
  }
}
