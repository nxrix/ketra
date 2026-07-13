import { SparsePauliOp } from "./../quantum_info/pauli.js";
import { Statevector } from "./../quantum_info/statevector.js";
import { Estimator } from "./../primitives/primitives.js";
import { SPSA } from "./optimizers.js";

export class VQEResult {
  constructor(kwargs = {}) {
    this.optimalParameters = kwargs.optimalParameters || {};
    this.optimalValue = kwargs.optimalValue || 0;
    this.optimalCircuit = kwargs.optimalCircuit || null;
    this.optimalState = kwargs.optimalState || null;
    this.costFunctionEvals = kwargs.costFunctionEvals || 0;
    this.optimizerResult = kwargs.optimizerResult || null;
    this.eigenvalue = this.optimalValue;
    this.eigenstate = this.optimalState;
  }
}

export class VQE {
  constructor(options = {}) {
    this.estimator = options.estimator || new Estimator();
    this.optimizer = options.optimizer || new SPSA({ maxiter: 100 });
    this.ansatz = options.ansatz || null;
    this.initialPoint = options.initialPoint || null;
    this.callback = options.callback || null;
    this.quantumInstance = options.quantumInstance || null;
  }

  // Compute the ground state of the given Hamiltonian.
  // operator: SparsePauliOp representing the Hamiltonian
  // ansatz: QuantumCircuit with parameters (overrides constructor's ansatz)
  // initialPoint: array of numbers, the initial parameter values
  computeMinimumEigenvalue(operator, options = {}) {
    const ansatz = options.ansatz || this.ansatz;
    if (!ansatz) throw new Error("VQE: ansatz is required");
    if (!operator) throw new Error("VQE: operator (Hamiltonian) is required");

    const paramSet = Array.from(ansatz.parameters);
    const paramNames = paramSet.map(p => p.name);

    let x0;
    if (options.initialPoint) {
      x0 = options.initialPoint.slice();
    } else if (this.initialPoint) {
      x0 = this.initialPoint.slice();
    } else {
      // Random initial point in [-pi, pi]
      x0 = paramNames.map(() => (Math.random() * 2 - 1) * Math.PI);
    }

    let evalCount = 0;
    const objectiveFn = (x) => {
      const params = {};
      for (let i = 0; i < paramNames.length; i++) {
        params[paramNames[i]] = x[i];
      }
      const boundAnsatz = ansatz.bindParameters(params);
      const result = this.estimator.run(boundAnsatz, operator);
      const value = result.values[0];
      evalCount++;
      if (this.callback) {
        this.callback(evalCount, x, value);
      }
      // Return real part (expectation value should be real for Hermitian ops)
      return typeof value === "number" ? value : value.re;
    };

    const optResult = this.optimizer.minimize(objectiveFn, x0);

    const optimalParams = {};
    for (let i = 0; i < paramNames.length; i++) {
      optimalParams[paramNames[i]] = optResult.x[i];
    }
    const optimalCircuit = ansatz.bindParameters(optimalParams);
    const optimalState = Statevector.fromCircuit(optimalCircuit);

    return new VQEResult({
      optimalParameters: optimalParams,
      optimalValue: optResult.fun,
      optimalCircuit: optimalCircuit,
      optimalState: optimalState,
      costFunctionEvals: evalCount,
      optimizerResult: optResult,
    });
  }
}
