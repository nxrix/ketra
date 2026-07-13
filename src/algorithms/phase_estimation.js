import { QuantumCircuit } from "../core/circuit.js";
import { Statevector } from "../quantum_info/statevector.js";
import { ControlledGate } from "../core/gate.js";
import { simulate } from "../simulator/statevector_simulator.js";

export class PhaseEstimationResult {
  constructor(kwargs = {}) {
    this.phase = kwargs.phase || 0;
    this.measurement = kwargs.measurement || null;
    this.circuit = kwargs.circuit || null;
  }
}

export class PhaseEstimation {
  constructor(options = {}) {
    this.numEvaluationQubits = options.numEvaluationQubits || 3;
    this.sampler = options.sampler || null;
  }

  // Estimate the phase of a unitary.
  // unitary: QuantumCircuit implementing U (any numQubits >= 1).
  // state_preparation: optional circuit to prepare the eigenstate on U's qubits.
  estimate(unitary, statePreparation = null) {
    const evalQubits = this.numEvaluationQubits;
    const unitaryQubits = unitary.numQubits;
    const totalQubits = evalQubits + unitaryQubits;

    const circuit = new QuantumCircuit(totalQubits, evalQubits);

    // Prepare eigenstate on unitary qubits (qubits evalQubits..totalQubits-1).
    // If statePreparation is a Statevector, convert it to an initialize
    // instruction on the unitary qubits.
    if (statePreparation) {
      if (statePreparation instanceof Statevector) {
        const targetQubits = [];
        for (let q = 0; q < unitaryQubits; q++) targetQubits.push(evalQubits + q);
        circuit.initialize(statePreparation._data.data, targetQubits);
      } else {
        // statePreparation is a QuantumCircuit
        const targetQubits = [];
        for (let q = 0; q < unitaryQubits; q++) targetQubits.push(evalQubits + q);
        circuit.compose(statePreparation, targetQubits, null, false, true);
      }
    }
    // Apply Hadamard to evaluation qubits (qubits 0..evalQubits-1).
    for (let q = 0; q < evalQubits; q++) {
      circuit.h(q);
    }
    circuit.barrier();

    // Apply controlled-U^(2^k) for k = 0, 1, ..., evalQubits-1.
    // Evaluation qubit k controls U applied 2^k times on the unitary qubits.
    for (let k = 0; k < evalQubits; k++) {
      const controlQubit = k;
      const unitaryQubitIndices = [];
      for (let q = 0; q < unitaryQubits; q++) {
        unitaryQubitIndices.push(evalQubits + q);
      }
      for (let p = 0; p < (1 << k); p++) {
        this._applyControlledUnitary(circuit, unitary, controlQubit, unitaryQubitIndices);
      }
    }
    circuit.barrier();

    // Apply inverse QFT on evaluation qubits.
    this._addInverseQFT(circuit, evalQubits);

    for (let q = 0; q < evalQubits; q++) {
      circuit.measure(q, q);
    }

    const result = simulate(circuit, 1024);
    const counts = result.getCounts().toDict();

    let bestKey = null, bestCount = -1;
    for (const [key, count] of Object.entries(counts)) {
      if (count > bestCount) { bestCount = count; bestKey = key; }
    }

    // Convert measurement to phase. The Counts key uses the simulator's
    // little-endian convention: qubit 0 is at the END of the string. In QPE,
    // qubit 0 is the least-significant evaluation qubit, so the key string
    // has the MSB (qubit evalQubits-1) as its FIRST character — which is
    // exactly what parseInt(_, 2) expects. No reversal needed.
    const measuredValue = parseInt(bestKey, 2);
    const phase = measuredValue / (1 << evalQubits);

    return new PhaseEstimationResult({
      phase,
      measurement: counts,
      circuit,
    });
  }

  // Apply the unitary circuit U controlled by `controlQubit`, targeting
  // `unitaryQubits` in order. Each gate in U is wrapped in a ControlledGate
  // and applied to the corresponding target qubit.
  _applyControlledUnitary(circuit, unitary, controlQubit, unitaryQubits) {
    for (const ci of unitary.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (typeof op.toMatrix !== "function") continue;
      // Map unitary's qubits onto the target qubits in the QPE circuit.
      const mappedQubits = ci.qubits.map(q => {
        const idxInUnitary = unitary._qubit_index.get(q);
        return unitaryQubits[idxInUnitary];
      });
      // Wrap op in a ControlledGate with one control qubit.
      try {
        const cg = new ControlledGate(op, 1);
        cg.name = "c" + op.name;
        circuit.append(cg, [controlQubit].concat(mappedQubits));
      } catch (e) {
        // If the gate cannot be made controlled (e.g. it's a non-unitary
        // Instruction), apply it unconditionally. This is a documented
        // limitation for non-unitary operators; for unitary gates the
        // ControlledGate path above is always taken.
        circuit.append(op.copy(), mappedQubits);
      }
    }
  }

  // Add inverse QFT (QFT-dagger) on the first n qubits.
  _addInverseQFT(circuit, n) {
    // QFT-dagger: reverse of QFT with negated angles.
    for (let i = 0; i < Math.floor(n / 2); i++) {
      circuit.swap(i, n - 1 - i);
    }
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < i; j++) {
        const lam = -Math.PI / Math.pow(2, i - j);
        circuit.cp(lam, i, j);
      }
      circuit.h(i);
    }
  }
}

// Amplitude Estimation
export class AmplitudeEstimationResult {
  constructor(kwargs = {}) {
    this.estimation = kwargs.estimation || 0;
    this.num_oracle_queries = kwargs.num_oracle_queries || 0;
    this.measurement = kwargs.measurement || null;
    this.circuit = kwargs.circuit || null;
  }
}

export class AmplitudeEstimation {
  constructor(options = {}) {
    this.num_eval_qubits = options.num_eval_qubits || 3;
    this.sampler = options.sampler || null;
  }

  // Estimate the amplitude of a good state.
  // statePreparation: circuit A that prepares the state with amplitude a.
  // groverOperator: circuit Q = A (2|0><0| - I) A^dagger (2|psi_g><psi_g| - I).
  // The amplitude a is estimated as sin^2(pi * y / 2^m) where y is the
  // measured value and m is the number of evaluation qubits.
  estimate(statePreparation, groverOperator = null) {
    const evalQubits = this.num_eval_qubits;
    const stateQubits = statePreparation.numQubits;
    const totalQubits = evalQubits + stateQubits;

    const circuit = new QuantumCircuit(totalQubits, evalQubits);

    // Apply Hadamard to evaluation qubits
    for (let q = 0; q < evalQubits; q++) circuit.h(q);

    // Apply state preparation on the state qubits (qubits evalQubits..total-1)
    const stateQubitIndices = [];
    for (let q = 0; q < stateQubits; q++) stateQubitIndices.push(evalQubits + q);
    this._applyCircuitToQubits(circuit, statePreparation, stateQubitIndices);

    // Apply Grover operator Q^(2^k) controlled by evaluation qubit k.
    if (groverOperator) {
      for (let k = 0; k < evalQubits; k++) {
        const controlQubit = k;
        for (let p = 0; p < (1 << k); p++) {
          this._applyControlledCircuit(circuit, groverOperator, controlQubit, stateQubitIndices);
        }
      }
    }

    // Apply inverse QFT on evaluation qubits
    const pe = new PhaseEstimation({ numEvaluationQubits: evalQubits });
    pe._addInverseQFT(circuit, evalQubits);

    for (let q = 0; q < evalQubits; q++) circuit.measure(q, q);

    const result = simulate(circuit, 1024);
    const counts = result.getCounts().toDict();

    let bestKey = null, bestCount = -1;
    for (const [key, count] of Object.entries(counts)) {
      if (count > bestCount) { bestCount = count; bestKey = key; }
    }

    // Same key-ordering convention as PhaseEstimation: parseInt is correct
    // as-is because the key has the MSB as its first character.
    const y = parseInt(bestKey, 2);
    const a = Math.sin(Math.PI * y / (1 << evalQubits)) ** 2;

    return new AmplitudeEstimationResult({
      estimation: a,
      num_oracle_queries: 2 * evalQubits,
      measurement: counts,
      circuit,
    });
  }

  // Apply each gate of `sub` onto the specified target qubit indices.
  _applyCircuitToQubits(circuit, sub, targetQubits) {
    for (const ci of sub.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (typeof op.toMatrix !== "function") continue;
      const mappedQubits = ci.qubits.map(q => {
        const idx = sub._qubit_index.get(q);
        return targetQubits[idx];
      });
      circuit.append(op.copy(), mappedQubits);
    }
  }

  // Apply each gate of `sub` controlled by `controlQubit`.
  _applyControlledCircuit(circuit, sub, controlQubit, targetQubits) {
    for (const ci of sub.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (typeof op.toMatrix !== "function") continue;
      const mappedQubits = ci.qubits.map(q => {
        const idx = sub._qubit_index.get(q);
        return targetQubits[idx];
      });
      try {
        const cg = new ControlledGate(op, 1);
        cg.name = "c" + op.name;
        circuit.append(cg, [controlQubit].concat(mappedQubits));
      } catch (e) {
        circuit.append(op.copy(), mappedQubits);
      }
    }
  }
}
