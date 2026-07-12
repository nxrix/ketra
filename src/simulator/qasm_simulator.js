/**
 * qasm_simulator.js - Noisy quantum circuit simulator.
 *
 * Noisy quantum circuit simulator. Supports:
 *   - Statevector simulation with shot-based sampling
 *   - Noise models (depolarizing, bit flip, amplitude damping, etc.)
 *   - Readout errors
 *   - Multiple shots
 *
 * For noise: uses density matrix evolution (Kraus operators).
 */

import { Complex, ComplexMatrix, ComplexVector, sampleDistribution, randomUniform } from "../math/linalg.js";
import { Result, Counts } from "../result/result.js";
import { DensityMatrix } from "../quantum_info/density_matrix.js";

export class QasmSimulator {
  constructor(options = {}) {
    this.method = options.method || "automatic";
    this.noise_model = options.noise_model || null;
    this.basis_gates = options.basis_gates || null;
    this.max_memory_mb = options.max_memory_mb || 8192;
    this.seed = options.seed || null;
  }

  // Run one or more circuits
  run(circuits, shots = 1024, options = {}) {
    const circuitList = Array.isArray(circuits) ? circuits : [circuits];
    const noiseModel = options.noise_model || this.noise_model;
    const results = circuitList.map(c => this._runOne(c, shots, noiseModel, options));
    return new Result({
      backend_name: "qasm_simulator",
      backend_version: "1.0.0",
      qobj_id: options.qobj_id || "qobj",
      job_id: options.job_id || `job_${Date.now()}`,
      success: true,
      results,
    });
  }

  _runOne(circuit, shots, noiseModel, options) {
    const n = circuit.num_qubits;
    const hasNoise = noiseModel && !noiseModel.is_empty();

    if (hasNoise) {
      return this._runNoisy(circuit, shots, noiseModel, options);
    }
    return this._runIdeal(circuit, shots, options);
  }

  // Ideal simulation: statevector + sampling
  _runIdeal(circuit, shots, options) {
    const n = circuit.num_qubits;
    const state = ComplexVector.zeros(1 << n);
    state.data[0] = Complex.ONE;
    const finalMeasures = [];

    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "delay") continue;
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
      if (typeof op.to_matrix !== "function") continue;
      const gateMatrix = op.to_matrix();
      const qubitIndices = ci.qubits.map(q => circuit._qubit_index.get(q));
      _applyGateInPlace(state, gateMatrix, op.num_qubits, qubitIndices);
    }

    // Sample
    const probs = state.probabilities();
    const counts = {};
    const memory = new Array(shots);
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
          bits[m.cIdx] = ((idx >> m.qIdx) & 1).toString();
        }
      }
      const key = bits.reverse().join("");
      counts[key] = (counts[key] || 0) + 1;
      memory[s] = key;
    }

    return {
      success: true,
      shots,
      data: { counts: new Counts(counts), memory },
      status: "DONE",
    };
  }

  // Noisy simulation: density matrix evolution with Kraus operators
  _runNoisy(circuit, shots, noiseModel, options) {
    const n = circuit.num_qubits;
    // For each shot, we evolve a pure state and apply noise stochastically
    // (trajectory method). This is more efficient than density matrix for
    // large numbers of shots.
    const counts = {};
    const memory = new Array(shots);
    const finalMeasures = [];

    // Pre-collect measurements (qIdx / cIdx pairing for readout errors).
    for (const ci of circuit.data) {
      if (ci.operation.name === "measure") {
        finalMeasures.push({
          op: ci.operation,
          qIdx: circuit._qubit_index.get(ci.qubits[0]),
          cIdx: circuit._clbit_index.get(ci.clbits[0]),
          qubits: ci.qubits.map(q => circuit._qubit_index.get(q)),
          clbits: ci.clbits.map(c => circuit._clbit_index.get(c)),
        });
      }
    }

    for (let shot = 0; shot < shots; shot++) {
      const state = ComplexVector.zeros(1 << n);
      state.data[0] = Complex.ONE;
      let measuredBits = new Array(circuit.num_clbits).fill(0);

      for (const ci of circuit.data) {
        const op = ci.operation;
        if (op.name === "barrier" || op.name === "delay") continue;
        if (op.name === "measure") {
          const qIdx = circuit._qubit_index.get(ci.qubits[0]);
          const cIdx = circuit._clbit_index.get(ci.clbits[0]);
          // Measure qubit qIdx
          const prob1 = _probOfOne(state, qIdx);
          if (randomUniform() < prob1) {
            _projectToOne(state, qIdx, Math.sqrt(prob1));
            measuredBits[cIdx] = 1;
          } else {
            _projectToZero(state, qIdx, Math.sqrt(1 - prob1));
            measuredBits[cIdx] = 0;
          }
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
        if (typeof op.to_matrix !== "function") continue;

        // Apply the gate
        const gateMatrix = op.to_matrix();
        const qubitIndices = ci.qubits.map(q => circuit._qubit_index.get(q));
        _applyGateInPlace(state, gateMatrix, op.num_qubits, qubitIndices);

        // Apply noise if this gate has an error
        const error = noiseModel.get_quantum_error(op.name, qubitIndices);
        if (error) {
          _applyQuantumError(state, error, qubitIndices, n);
        }
      }

      // Apply readout errors. We track which classical bit came from which
      // qubit via finalMeasures so per-qubit readout errors apply only to the
      // corresponding classical bit. Readout errors registered with
      // qubits === "all" apply to every classical bit.
      for (const re of noiseModel.get_readout_errors()) {
        const targetQubits = re.qubits;
        for (const m of finalMeasures) {
          let appliesThisQubit = false;
          if (targetQubits === "all") appliesThisQubit = true;
          else if (Array.isArray(targetQubits) && targetQubits.includes(m.qIdx)) {
            appliesThisQubit = true;
          }
          if (!appliesThisQubit) continue;
          const r = randomUniform();
          // probabilities[actual][measured] — flip with prob[actual][1-actual]
          const errProb = re.error.probabilities[measuredBits[m.cIdx]][1 - measuredBits[m.cIdx]];
          if (r < errProb) measuredBits[m.cIdx] = 1 - measuredBits[m.cIdx];
        }
      }

      const key = measuredBits.reverse().join("");
      counts[key] = (counts[key] || 0) + 1;
      memory[shot] = key;
    }

    return {
      success: true,
      shots,
      data: { counts: new Counts(counts), memory },
      status: "DONE",
    };
  }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

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

function _probOfZero(state, qIdx) {
  let p0 = 0;
  for (let i = 0; i < state.size; i++) {
    if (((i >> qIdx) & 1) === 0) p0 += state.data[i].abs2();
  }
  return p0;
}

function _probOfOne(state, qIdx) {
  let p1 = 0;
  for (let i = 0; i < state.size; i++) {
    if (((i >> qIdx) & 1) === 1) p1 += state.data[i].abs2();
  }
  return p1;
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

// Apply a quantum error stochastically (trajectory method).
// `error` is a QuantumError whose Kraus operators act on the same qubits as
// the triggering gate. We sample one Kraus op per shot, weighted by the Born
// probability <psi|K^dagger K|psi>, and renormalize the post-Kraus state.
function _applyQuantumError(state, error, qubitIndices, totalQubits) {
  // Sample which error term occurs (some QuantumErrors group Kraus operators
  // into classical mixture terms; others put them all in one term).
  const r = randomUniform();
  let cum = 0;
  let chosenTerm = error.terms[0];
  for (const term of error.terms) {
    cum += term.probability;
    if (r < cum) { chosenTerm = term; break; }
  }
  if (chosenTerm.operators.length === 0) return;

  // Determine how many qubits each Kraus operator acts on.
  // For a k-qubit error the Kraus matrices are 2^k x 2^k.
  const k0 = chosenTerm.operators[0];
  const krausQubits = Math.log2(k0.rows);
  if (!Number.isInteger(krausQubits) || krausQubits < 1) {
    // Unknown shape — skip.
    return;
  }

  // Apply Kraus ops to the relevant qubits. The error's qubit count may be
  // smaller than the gate's qubit count (e.g. 1-qubit depolarizing on a
  // 2-qubit CX). In that case we iterate over each k-qubit subset of the
  // gate's qubitIndices — but for the standard noise model where the error
  // is registered for the full gate width, we just apply once on the first
  // k qubits.
  let krausOp;
  if (chosenTerm.operators.length === 1) {
    krausOp = chosenTerm.operators[0];
  } else {
    // Multiple Kraus ops: pick one stochastically based on Born probability
    // P(K) = <psi|K^dagger K|psi> restricted to the qubits K acts on.
    const probs = chosenTerm.operators.map(K => {
      // We compute <psi|K^dagger K|psi> over the FULL state. K is embedded
      // implicitly as K ⊗ I on the qubits it does not act on. The expectation
      // value factors through the marginal state on K's qubits, but for
      // simplicity we use the trace over the full state.
      const KdK = K.dagger().mul(K);
      // Apply KdK to the full state via _applyGateInPlace on the same qubits
      // (this gives us the marginal contribution).
      const tempState = { size: state.size, data: state.data.slice() };
      _applyGateInPlace(tempState, KdK, krausQubits, qubitIndices.slice(0, krausQubits));
      // <psi|KdK|psi> = sum |psi_i|^2 after applying KdK (which is Hermitian PSD).
      let p = 0;
      for (let i = 0; i < state.size; i++) {
        // tempState[i] = (KdK psi)_i; <psi|KdK|psi> = sum conj(psi_i) * (KdK psi)_i
        p += state.data[i].conjugate().mul(tempState.data[i]).re;
      }
      return Math.max(0, p);
    });
    const totalP = probs.reduce((a, b) => a + b, 0);
    if (totalP < 1e-15) return;
    const r2 = randomUniform() * totalP;
    let cum2 = 0;
    let idx = 0;
    for (let i = 0; i < probs.length; i++) {
      cum2 += probs[i];
      if (r2 < cum2) { idx = i; break; }
    }
    krausOp = chosenTerm.operators[idx];
  }

  // Apply Kraus op to the first `krausQubits` of qubitIndices.
  _applyGateInPlace(state, krausOp, krausQubits, qubitIndices.slice(0, krausQubits));

  // Renormalize the post-measurement state.
  const norm = state.norm();
  if (norm > 1e-15) {
    const inv = 1 / norm;
    for (let i = 0; i < state.size; i++) {
      state.data[i] = state.data[i].scale(inv);
    }
  }
}

// Convenience function
export function simulate_noisy(circuit, shots = 1024, options = {}) {
  const sim = new QasmSimulator(options);
  return sim.run(circuit, shots, options);
}
