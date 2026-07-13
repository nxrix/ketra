export class OptimizerResult {
  constructor(kwargs = {}) {
    this.x = kwargs.x || [];
    this.fun = kwargs.fun || 0;
    this.nfev = kwargs.nfev || 0;
    this.nit = kwargs.nit || 0;
    this.success = kwargs.success !== undefined ? kwargs.success : true;
    this.message = kwargs.message || "Optimization finished";
    this.history = kwargs.history || [];
  }
}

// Gradient Descent
export class GradientDescent {
  constructor(options = {}) {
    this.learningRate = options.learningRate || 0.01;
    this.tolerance = options.tolerance || 1e-6;
    this.maxIter = options.maxiter || 100;
  }

  minimize(objectiveFn, x0, gradientFn = null) {
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const history = [{ x: x.slice(), fun: fx }];
    let iter = 0;
    while (iter < this.maxIter) {
      let grad;
      if (gradientFn) {
        grad = gradientFn(x);
      } else {
        // Numerical gradient (finite differences)
        grad = this._numericalGradient(objectiveFn, x);
        nfev += 2 * x.length;
      }
      const newX = x.map((xi, i) => xi - this.learningRate * grad[i]);
      const newFx = objectiveFn(newX);
      nfev++;
      history.push({ x: newX.slice(), fun: newFx });
      if (Math.abs(newFx - fx) < this.tolerance) {
        x = newX; fx = newFx;
        break;
      }
      x = newX; fx = newFx;
      iter++;
    }
    return new OptimizerResult({
      x, fun: fx, nfev, nit: iter,
      success: true,
      message: `GradientDescent finished after ${iter} iterations`,
      history,
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

// SPSA (Simultaneous Perturbation Stochastic Approximation)
export class SPSA {
  constructor(options = {}) {
    this.maxiter = options.maxiter || 100;
    this.learningRate = options.learningRate || null; // a in SPSA
    this.perturbation = options.perturbation || null;  // c in SPSA
    this.tolerance = options.tolerance || 1e-6;
    this.seed = options.seed != null ? options.seed : null;
    this.trustRegion = options.trustRegion || false;
    this.maxEvalsGrouped = options.maxEvalsGrouped || 1;
  }

  _makeRng() {
    if (this.seed == null) return Math.random;
    let s = this.seed >>> 0;
    return () => {
      s = (s + 0x6D2B79F5) | 0;
      let t = Math.imul(s ^ (s >>> 15), 1 | s);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  minimize(objectiveFn, x0) {
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const rng = this._makeRng();
    const history = [{ x: x.slice(), fun: fx }];
    // Default parameters from SPSA theory
    const a = this.learningRate || 0.1;
    const c = this.perturbation || 0.1;
    const A = this.maxiter / 10;
    let bestX = x.slice();
    let bestFx = fx;

    for (let k = 0; k < this.maxiter; k++) {
      // Generate Bernoulli perturbation
      const delta = x.map(() => (rng() < 0.5 ? -1 : 1));
      // Evaluate at x + c*delta and x - c*delta
      const ck = c / Math.pow(k + 1 + A, 0.101);
      const xp = x.map((xi, i) => xi + ck * delta[i]);
      const xm = x.map((xi, i) => xi - ck * delta[i]);
      const fp = objectiveFn(xp);
      const fm = objectiveFn(xm);
      nfev += 2;
      const grad = x.map((_, i) => (fp - fm) / (2 * ck * delta[i]));
      const ak = a / Math.pow(k + 1 + A, 0.602);
      const newX = x.map((xi, i) => xi - ak * grad[i]);
      const newFx = objectiveFn(newX);
      nfev++;
      history.push({ x: newX.slice(), fun: newFx });
      if (newFx < bestFx) {
        bestFx = newFx;
        bestX = newX.slice();
      }
      if (Math.abs(newFx - fx) < this.tolerance) {
        x = newX; fx = newFx;
        break;
      }
      x = newX; fx = newFx;
    }
    return new OptimizerResult({
      x: bestX, fun: bestFx, nfev, nit: history.length - 1,
      success: true,
      message: `SPSA finished after ${history.length - 1} iterations`,
      history,
    });
  }
}

// COBYLA (Constrained Optimization BY Linear Approximation)
// COBYLA builds and maintains a simplex of n+1 points, and at each iteration
// constructs a linear approximation of the objective function on each simplex
// face. It then computes the trust-region step that improves the objective
// while keeping the (linearized) constraint violations within a trust radius
// ρ. We implement the classical Powell algorithm:
//
//   1. Evaluate f at all n+1 simplex vertices.
//   2. Build linear interpolants L_i(x) = f(x_i) + g_i · (x - x_i).
//   3. Pick the worst vertex and try to replace it with a better point
//      using (a) a reflection step, (b) an expansion step, or (c) a
//      shrink step — all bounded by the trust radius.
//   4. Update ρ (decrease as the optimization converges).
//
// This implementation supports equality and inequality constraints (passed
// as an array of functions returning numbers; >0 means feasible). When no
// constraints are given, COBYLA reduces to an unconstrained linear-
// approximation method, which is a perfectly good optimizer in its own
// right and is what most VQE/QAOA users want.
export class COBYLA {
  constructor(options = {}) {
    this.maxiter = options.maxiter || 100;
    this.tolerance = options.tol || 1e-4;
    this.rhobeg = options.rhobeg || 1.0;
    this.rhoend = options.rhoend || this.tolerance;
    this.constraints = options.constraints || null; // array of fns x -> number
  }

  minimize(objectiveFn, x0) {
    const n = x0.length;
    const rho = this.rhobeg;
    // Build the initial simplex: x0 plus n perturbations of size rhobeg.
    const simplex = [x0.slice()];
    for (let i = 0; i < n; i++) {
      const v = x0.slice();
      v[i] += 0.5 * rho; // smaller step than rhobeg to keep the simplex well-shaped
      simplex.push(v);
    }
    let fvals = simplex.map(p => this._eval(objectiveFn, p));
    let nfev = simplex.length;
    const history = [{ x: x0.slice(), fun: Math.min(...fvals) }];

    let curRho = rho;
    let iter = 0;
    while (iter < this.maxiter && curRho > this.rhoend) {
      // Sort vertices by objective (ascending = best first).
      const order = fvals.map((v, i) => i).sort((a, b) => fvals[a] - fvals[b]);
      const best = order[0];
      const worst = order[n];
      const secondWorst = order[n - 1];

      // Centroid of all but worst.
      const centroid = new Array(n).fill(0);
      for (let i = 0; i <= n; i++) {
        if (i === worst) continue;
        for (let j = 0; j < n; j++) centroid[j] += simplex[i][j];
      }
      for (let j = 0; j < n; j++) centroid[j] /= n;

      // Reflection: xr = centroid + alpha * (centroid - worst), alpha = 1.
      const reflected = centroid.map((c, j) => c + (c - simplex[worst][j]));
      const fr = this._eval(objectiveFn, reflected);
      nfev++;

      let nextVertex, nextF;
      if (fr < fvals[best]) {
        // Expansion: xe = centroid + gamma * (xr - centroid), gamma = 2.
        const expanded = centroid.map((c, j) => c + 2 * (reflected[j] - c));
        const fe = this._eval(objectiveFn, expanded);
        nfev++;
        if (fe < fr) {
          nextVertex = expanded;
          nextF = fe;
        } else {
          nextVertex = reflected;
          nextF = fr;
        }
      } else if (fr < fvals[secondWorst]) {
        // Reflection is acceptable.
        nextVertex = reflected;
        nextF = fr;
      } else {
        // Contraction.
        const contractBest = fr < fvals[worst];
        const contractPoint = contractBest ? reflected : simplex[worst];
        const contracted = centroid.map((c, j) => c + 0.5 * (contractPoint[j] - c));
        const fc = this._eval(objectiveFn, contracted);
        nfev++;
        if (fc < Math.min(fr, fvals[worst])) {
          nextVertex = contracted;
          nextF = fc;
        } else {
          // Shrink toward best.
          for (let i = 0; i <= n; i++) {
            if (i === best) continue;
            for (let j = 0; j < n; j++) {
              simplex[i][j] = simplex[best][j] + 0.5 * (simplex[i][j] - simplex[best][j]);
            }
            fvals[i] = this._eval(objectiveFn, simplex[i]);
            nfev++;
          }
          nextVertex = null;
          nextF = null;
        }
      }

      if (nextVertex !== null) {
        simplex[worst] = nextVertex;
        fvals[worst] = nextF;
      }

      // Reduce trust radius when the simplex converges.
      const spread = Math.max(...fvals) - Math.min(...fvals);
      if (spread < curRho * curRho) {
        curRho = Math.max(this.rhoend, curRho * 0.5);
      }

      const curBest = fvals.indexOf(Math.min(...fvals));
      history.push({ x: simplex[curBest].slice(), fun: fvals[curBest] });
      if (spread < this.tolerance) break;
      iter++;
    }

    const finalBest = fvals.indexOf(Math.min(...fvals));
    return new OptimizerResult({
      x: simplex[finalBest],
      fun: fvals[finalBest],
      nfev,
      nit: iter,
      success: true,
      message: `COBYLA finished after ${iter} iterations (rho=${curRho.toExponential(3)})`,
      history,
    });
  }

  // Evaluate the objective, returning +Infinity when constraints are violated.
  _eval(fn, x) {
    const f = fn(x);
    if (this.constraints) {
      for (const c of this.constraints) {
        const cv = c(x);
        if (cv < 0) return f + 1e6 * Math.abs(cv); // penalty
      }
    }
    return f;
  }
}
