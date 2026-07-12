/**
 * vqe.js - Variational Quantum Eigensolver.
 *
 * Uses a parameterized ansatz circuit and a
 * classical optimizer to find the ground state energy of a Hamiltonian.
 *
 * The Hamiltonian is provided as a SparsePauliOp. The ansatz is any
 * QuantumCircuit with parameters. The optimizer is one of the optimizers
 * from algorithms/optimizers.js.
 */

import { Complex } from "./../math/linalg.js";
import { SparsePauliOp } from "./../quantum_info/pauli.js";
import { Statevector } from "./../quantum_info/statevector.js";
import { Estimator } from "./../primitives/primitives.js";
import { SPSA, GradientDescent, OptimizerResult } from "./optimizers.js";

export class VQEResult {
  constructor(kwargs = {}) {
    this.optimal_parameters = kwargs.optimal_parameters || {};
    this.optimal_value = kwargs.optimal_value || 0;
    this.optimal_circuit = kwargs.optimal_circuit || null;
    this.optimal_state = kwargs.optimal_state || null;
    this.cost_function_evals = kwargs.cost_function_evals || 0;
    this.optimizer_result = kwargs.optimizer_result || null;
    this.eigenvalue = this.optimal_value;
    this.eigenstate = this.optimal_state;
  }
}

export class VQE {
  constructor(options = {}) {
    this.estimator = options.estimator || new Estimator();
    this.optimizer = options.optimizer || new SPSA({ maxiter: 100 });
    this.ansatz = options.ansatz || null;
    this.initial_point = options.initial_point || null;
    this.callback = options.callback || null;
    this.quantum_instance = options.quantum_instance || null;
  }

  // Compute the ground state of the given Hamiltonian.
  // operator: SparsePauliOp representing the Hamiltonian
  // ansatz: QuantumCircuit with parameters (overrides constructor's ansatz)
  // initial_point: array of numbers, the initial parameter values
  compute_minimum_eigenvalue(operator, options = {}) {
    const ansatz = options.ansatz || this.ansatz;
    if (!ansatz) throw new Error("VQE: ansatz is required");
    if (!operator) throw new Error("VQE: operator (Hamiltonian) is required");

    // Collect parameters from ansatz
    const paramSet = Array.from(ansatz.parameters);
    const paramNames = paramSet.map(p => p.name);

    // Initial point
    let x0;
    if (options.initial_point) {
      x0 = options.initial_point.slice();
    } else if (this.initial_point) {
      x0 = this.initial_point.slice();
    } else {
      // Random initial point in [-pi, pi]
      x0 = paramNames.map(() => (Math.random() * 2 - 1) * Math.PI);
    }

    // Objective function
    let evalCount = 0;
    const objectiveFn = (x) => {
      const params = {};
      for (let i = 0; i < paramNames.length; i++) {
        params[paramNames[i]] = x[i];
      }
      const boundAnsatz = ansatz.bind_parameters(params);
      const result = this.estimator.run(boundAnsatz, operator);
      const value = result.values[0];
      evalCount++;
      if (this.callback) {
        this.callback(evalCount, x, value);
      }
      // Return real part (expectation value should be real for Hermitian ops)
      return typeof value === "number" ? value : value.re;
    };

    // Optimize
    const optResult = this.optimizer.minimize(objectiveFn, x0);

    // Build the optimal state
    const optimalParams = {};
    for (let i = 0; i < paramNames.length; i++) {
      optimalParams[paramNames[i]] = optResult.x[i];
    }
    const optimalCircuit = ansatz.bind_parameters(optimalParams);
    const optimalState = Statevector.fromCircuit(optimalCircuit);

    return new VQEResult({
      optimal_parameters: optimalParams,
      optimal_value: optResult.fun,
      optimal_circuit: optimalCircuit,
      optimal_state: optimalState,
      cost_function_evals: evalCount,
      optimizer_result: optResult,
    });
  }
}
