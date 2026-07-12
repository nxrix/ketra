/**
 * grover.js - Grover's search algorithm.
 *
 * Searches for marked items in an unstructured database using amplitude
 * amplification. Reference: Grover (1996).
 */

import { QuantumCircuit } from "../core/circuit.js";
import { Statevector } from "../quantum_info/statevector.js";
import { simulate } from "../simulator/statevector_simulator.js";

export class GroverResult {
  constructor(kwargs = {}) {
    this.top_measurement = kwargs.top_measurement || null;
    this.iterations = kwargs.iterations || 0;
    this.measurement = kwargs.measurement || null;
    this.circuit = kwargs.circuit || null;
  }
}

export class Grover {
  constructor(options = {}) {
    this.iterations = options.iterations || null; // auto-determine if null
    this.growth_rate = options.growth_rate || 2.0;
    this.sample_from_iterations = options.sample_from_iterations || false;
    this.sampler = options.sampler || null;
  }

  // Run Grover's algorithm.
  // oracle: QuantumCircuit that marks the solution(s) by flipping phase
  // numSolutions: number of marked items (for determining iterations)
  amplify(oracle, numSolutions = 1) {
    const numQubits = oracle.num_qubits;
    const N = 1 << numQubits;
    const M = numSolutions;

    // Determine optimal number of iterations
    let numIters;
    if (this.iterations !== null) {
      numIters = this.iterations;
    } else {
      // Optimal: pi/4 * sqrt(N/M)
      numIters = Math.floor(Math.PI / 4 * Math.sqrt(N / M));
      if (numIters < 1) numIters = 1;
    }

    // Build Grover circuit
    const circuit = new QuantumCircuit(numQubits, numQubits);

    // Step 1: Initialize to uniform superposition
    for (let q = 0; q < numQubits; q++) circuit.h(q);
    circuit.barrier();

    // Step 2: Repeat Grover iteration
    for (let i = 0; i < numIters; i++) {
      // Apply oracle (phase flip on marked items). compose() returns a new
      // circuit by default; use inplace=true to mutate `circuit`.
      circuit.compose(oracle, null, null, false, true);
      // Apply diffusion operator (Grover's diffusion)
      this._addDiffusion(circuit, numQubits);
      circuit.barrier();
    }

    // Step 3: Measure
    for (let q = 0; q < numQubits; q++) circuit.measure(q, q);

    // Simulate
    const result = simulate(circuit, 1024);
    const counts = result.get_counts().to_dict();

    // Find the most frequent measurement
    let bestKey = null, bestCount = -1;
    for (const [key, count] of Object.entries(counts)) {
      if (count > bestCount) { bestCount = count; bestKey = key; }
    }

    return new GroverResult({
      top_measurement: bestKey,
      iterations: numIters,
      measurement: counts,
      circuit,
    });
  }

  // Diffusion operator: 2|s><s| - I, where |s> = H^n |0...0>.
  // Decomposition: H^n · (2|0><0| - I) · H^n = H^n · X^n · (2|1...1><1...1| - I) · X^n · H^n
  // The middle operator is a multi-controlled Z (MCZ) up to a global phase.
  // For n=1: MCZ = Z. For n=2: MCZ = CZ. For n>=3: MCZ = H(target) · MCX(controls, target) · H(target).
  _addDiffusion(circuit, numQubits) {
    // H^n
    for (let q = 0; q < numQubits; q++) circuit.h(q);
    // X^n
    for (let q = 0; q < numQubits; q++) circuit.x(q);

    if (numQubits === 1) {
      // (2|1><1| - I) on a single qubit = Z
      circuit.z(0);
    } else if (numQubits === 2) {
      // MCZ with one control = CZ
      circuit.cz(0, 1);
    } else {
      // MCZ with n-1 controls: convert to MCX via H on the target.
      circuit.h(numQubits - 1);
      const controls = [];
      for (let q = 0; q < numQubits - 1; q++) controls.push(q);
      circuit.mcx(controls, numQubits - 1);
      circuit.h(numQubits - 1);
    }

    // X^n
    for (let q = 0; q < numQubits; q++) circuit.x(q);
    // H^n
    for (let q = 0; q < numQubits; q++) circuit.h(q);
  }
}
