/**
 * statevector_simulator.js - A pure statevector quantum circuit simulator.
 *
 * Accepts a
 * QuantumCircuit, executes it on an internal statevector backend, returns a
 * Result object containing the final statevector and (optionally) sampled
 * measurement counts.
 */

import { Complex, ComplexVector, ComplexMatrix, sampleDistribution, randomUniform } from "./../math/linalg.js";
import { Result, Counts } from "./../result/result.js";

function _applyGateInPlace(state, gateMatrix, k, qubitIndices) {
  const n = state.size;
  const dimK = 1 << k;
  const processed = new Array(n).fill(false);
  const newState = new Array(n);
  for (let i = 0; i < n; i++) newState[i] = state.data[i];

  for (let i = 0; i < n; i++) {
    if (processed[i]) continue;
    const subs = new Array(dimK);
    const idxs = new Array(dimK);
    for (let s = 0; s < dimK; s++) {
      let idx = i;
      for (let qi = 0; qi < k; qi++) {
        const mask = 1 << qubitIndices[qi];
        idx &= ~mask;
      }
      for (let qi = 0; qi < k; qi++) {
        if ((s >> qi) & 1) idx |= (1 << qubitIndices[qi]);
      }
      idxs[s] = idx;
      subs[s] = state.data[idx];
      processed[idx] = true;
    }
    const newSubs = new Array(dimK);
    for (let rr = 0; rr < dimK; rr++) {
      let acc = Complex.ZERO;
      for (let cc = 0; cc < dimK; cc++) {
        acc = acc.add(gateMatrix.get(rr, cc).mul(subs[cc]));
      }
      newSubs[rr] = acc;
    }
    for (let s = 0; s < dimK; s++) newState[idxs[s]] = newSubs[s];
  }
  state.data = newState;
}

export class StatevectorSimulator {
  constructor(options = {}) {
    this.method = options.method || "statevector";
    this.precision = options.precision || "double";
    this.max_memory_mb = options.max_memory_mb || 8192;
    this.initial_statevector = options.initial_statevector || null;
  }

  run(circuit, shots = 1024, options = {}) {
    const experiments = Array.isArray(circuit) ? circuit : [circuit];
    const results = experiments.map(c => this._runOne(c, shots, options));
    return new Result({
      backend_name: "statevector_simulator",
      backend_version: "1.0.0",
      qobj_id: options.qobj_id || "qobj",
      job_id: options.job_id || `job_${Date.now()}`,
      success: true,
      results,
    });
  }

  _runOne(circuit, shots, options) {
    const n = circuit.num_qubits;
    const state = this.initial_statevector
      ? new ComplexVector(this.initial_statevector.data.slice())
      : ComplexVector.zeros(1 << n);
    if (!this.initial_statevector) state.data[0] = Complex.ONE;

    const finalMeasures = [];

    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier") continue;
      if (op.name === "delay") continue;
      if (op.name === "measure") {
        const qIdx = circuit._qubit_index.get(ci.qubits[0]);
        const cIdx = circuit._clbit_index.get(ci.clbits[0]);
        finalMeasures.push({ qIdx, cIdx });
        continue;
      }
      if (op.name === "reset") {
        const qIdx = circuit._qubit_index.get(ci.qubits[0]);
        const prob0 = _probOfZero(state, qIdx);
        if (randomUniform() < prob0) {
          _projectToZero(state, qIdx, Math.sqrt(prob0));
        } else {
          _projectToOne(state, qIdx, Math.sqrt(1 - prob0));
        }
        continue;
      }
      if (op.name === "if_else" || op.name === "while_loop") continue;
      if (typeof op.to_matrix !== "function") {
        throw new Error(`Cannot simulate ${op.name} (no matrix)`);
      }
      const gateMatrix = op.to_matrix();
      const qubitIndices = ci.qubits.map(q => circuit._qubit_index.get(q));
      _applyGateInPlace(state, gateMatrix, op.num_qubits, qubitIndices);
    }

    let counts = null;
    let memory = null;
    if (shots > 0) {
      const probs = state.probabilities();
      counts = {};
      memory = new Array(shots);
      const numClbits = circuit.num_clbits;
      const numQubits = circuit.num_qubits;
      for (let s = 0; s < shots; s++) {
        const idx = sampleDistribution(probs, options.rng);
        const bits = new Array(numClbits).fill("0");
        if (finalMeasures.length === 0) {
          for (let q = 0; q < numQubits && q < numClbits; q++) {
            bits[q] = ((idx >> q) & 1).toString();
          }
        } else {
          for (const m of finalMeasures) {
            const b = (idx >> m.qIdx) & 1;
            bits[m.cIdx] = b.toString();
          }
        }
        const key = bits.reverse().join("");
        counts[key] = (counts[key] || 0) + 1;
        memory[s] = key;
      }
    }

    return {
      success: true,
      shots,
      data: {
        statevector: new ComplexVector(state.data.slice()),
        counts: counts ? new Counts(counts) : null,
        memory,
      },
      status: "DONE",
      time_taken: 0,
    };
  }
}

function _probOfZero(state, qIdx) {
  let p0 = 0;
  for (let i = 0; i < state.size; i++) {
    if (((i >> qIdx) & 1) === 0) p0 += state.data[i].abs2();
  }
  return p0;
}

function _projectToZero(state, qIdx, norm) {
  for (let i = 0; i < state.size; i++) {
    if (((i >> qIdx) & 1) === 1) state.data[i] = Complex.ZERO;
    else state.data[i] = state.data[i].scale(1 / norm);
  }
}

function _projectToOne(state, qIdx, norm) {
  for (let i = 0; i < state.size; i++) {
    if (((i >> qIdx) & 1) === 0) state.data[i] = Complex.ZERO;
    else state.data[i] = state.data[i].scale(1 / norm);
  }
}

export function simulate(circuit, shots = 1024, options = {}) {
  const sim = new StatevectorSimulator(options);
  return sim.run(circuit, shots, options);
}
