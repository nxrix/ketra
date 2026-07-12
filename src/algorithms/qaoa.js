/**
 * qaoa.js - Quantum Approximate Optimization Algorithm.
 *
 * Solves combinatorial optimization problems by mapping them to a cost
 * Hamiltonian and applying alternating layers of cost and mixer unitaries.
 *
 * Reference: Farhi, Goldstone, Gutmann (2014).
 */

import { Complex } from "./../math/linalg.js";
import { QuantumCircuit } from "./../core/circuit.js";
import { Parameter } from "./../core/parameter.js";
import { SparsePauliOp, Pauli } from "./../quantum_info/pauli.js";
import { Statevector } from "./../quantum_info/statevector.js";
import { Estimator } from "./../primitives/primitives.js";
import { SPSA, OptimizerResult } from "./optimizers.js";
import { pauliEvolution } from "./../library/circuits.js";

export class QAOAResult {
  constructor(kwargs = {}) {
    this.optimal_parameters = kwargs.optimal_parameters || {};
    this.optimal_value = kwargs.optimal_value || 0;
    this.optimal_circuit = kwargs.optimal_circuit || null;
    this.optimal_state = kwargs.optimal_state || null;
    this.cost_function_evals = kwargs.cost_function_evals || 0;
    this.eigenvalue = this.optimal_value;
    this.eigenstate = this.optimal_state;
  }
}

export class QAOA {
  constructor(options = {}) {
    this.estimator = options.estimator || new Estimator();
    this.optimizer = options.optimizer || new SPSA({ maxiter: 100 });
    this.reps = options.reps || 1;  // p parameter (number of QAOA layers)
    this.initial_state = options.initial_state || null;
    this.mixer = options.mixer || null;
    this.callback = options.callback || null;
  }

  // Compute the minimum eigenvalue of `operator` (cost Hamiltonian).
  // operator: SparsePauliOp representing the cost Hamiltonian H_C
  compute_minimum_eigenvalue(operator, options = {}) {
    if (!operator) throw new Error("QAOA: cost operator is required");
    const numQubits = operator.num_qubits;

    // Default mixer: H_M = sum_i X_i
    let mixerOp;
    if (this.mixer) {
      mixerOp = this.mixer;
    } else {
      const mixerTerms = [];
      for (let i = 0; i < numQubits; i++) {
        const label = "I".repeat(i) + "X" + "I".repeat(numQubits - i - 1);
        mixerTerms.push([label, 1.0]);
      }
      mixerOp = SparsePauliOp.from_list(mixerTerms);
    }

    // Build QAOA ansatz with p layers.
    // For each cost term exp(-i gamma * coeff * P) we apply the corresponding
    // Pauli-evolution circuit (basis change + CNOT chain + RZ + undo).
    // Mixer layer exp(-i beta * X_q) is implemented as RX(2*beta) per qubit.
    const ansatz = new QuantumCircuit(numQubits);
    // Initial state: default is |+>^n (apply H to all)
    if (this.initial_state) {
      ansatz.compose(this.initial_state, null, null, false, true);
    } else {
      for (let q = 0; q < numQubits; q++) ansatz.h(q);
    }
    ansatz.barrier();
    // Parameters: gamma_1..gamma_p (cost) and beta_1..beta_p (mixer)
    const gammaParams = [];
    const betaParams = [];
    for (let p = 0; p < this.reps; p++) {
      gammaParams.push(new Parameter(`gamma_${p}`));
      betaParams.push(new Parameter(`beta_${p}`));
    }
    // Apply cost and mixer layers
    for (let p = 0; p < this.reps; p++) {
      // Cost: exp(-i gamma_p * H_C). Apply each Pauli term via the standard
      // Pauli-evolution decomposition (basis change + CNOT chain + RZ).
      for (const [label, coeff] of operator.to_list()) {
        if (label === "I".repeat(numQubits)) continue;
        const cr = (coeff instanceof Complex) ? coeff.re : Number(coeff);
        const ci = (coeff instanceof Complex) ? coeff.im : 0;
        if (Math.abs(cr) < 1e-15 && Math.abs(ci) < 1e-15) continue;
        // exp(-i gamma * coeff * P) -> RZ(2 * gamma * cr) on the last qubit
        // of the CNOT chain (with appropriate basis rotations). We build this
        // via pauliEvolution and compose it into the ansatz.
        // Use ParameterExpression for the time parameter so binding works.
        const timeExpr = gammaParams[p].mul(cr);
        const evolved = pauliEvolution(label, timeExpr);
        ansatz.compose(evolved, null, null, false, true);
      }
      // Mixer: exp(-i beta_p * H_M) = product over qubits of RX(2*beta) on each.
      for (let q = 0; q < numQubits; q++) {
        ansatz.rx(betaParams[p].mul(2), q);
      }
      ansatz.barrier();
    }

    // Initial point: gamma=0.5, beta=0.5
    const paramNames = gammaParams.concat(betaParams).map(p => p.name);
    const x0 = new Array(paramNames.length).fill(0.5);

    // Objective
    let evalCount = 0;
    const objectiveFn = (x) => {
      const params = {};
      for (let i = 0; i < paramNames.length; i++) params[paramNames[i]] = x[i];
      const boundAnsatz = ansatz.bind_parameters(params);
      const result = this.estimator.run(boundAnsatz, operator);
      const value = result.values[0];
      evalCount++;
      if (this.callback) this.callback(evalCount, x, value);
      return typeof value === "number" ? value : value.re;
    };

    const optResult = this.optimizer.minimize(objectiveFn, x0);

    const optimalParams = {};
    for (let i = 0; i < paramNames.length; i++) {
      optimalParams[paramNames[i]] = optResult.x[i];
    }
    const optimalCircuit = ansatz.bind_parameters(optimalParams);
    const optimalState = Statevector.fromCircuit(optimalCircuit);

    return new QAOAResult({
      optimal_parameters: optimalParams,
      optimal_value: optResult.fun,
      optimal_circuit: optimalCircuit,
      optimal_state: optimalState,
      cost_function_evals: evalCount,
    });
  }
}
