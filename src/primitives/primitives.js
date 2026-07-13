import { Statevector } from "./../quantum_info/statevector.js";
import { SparsePauliOp, Pauli } from "./../quantum_info/pauli.js";
import { simulate as _simulate } from "./../simulator/statevector_simulator.js";
import { Result, Counts } from "./../result/result.js";

// Base classes
export class BaseEstimator {
  constructor(options = {}) {
    this.options = options;
    this._circuits = new Map();
  }

  run(circuits, observables, parameterValues = null, options = {}) {
    throw new Error("BaseEstimator.run not implemented");
  }
}

export class BaseSampler {
  constructor(options = {}) {
    this.options = options;
    this._circuits = new Map();
  }
}

// Estimator result
export class EstimatorResult {
  constructor(values, metadata = []) {
    this.values = values;  // array of complex expectation values
    this.metadata = metadata;
  }
}

export class SamplerResult {
  constructor(quasiDists, metadata = []) {
    this.quasiDists = quasiDists;
    this.metadata = metadata;
  }
}

// Estimator
export class Estimator extends BaseEstimator {
  constructor(options = {}) {
    super(options);
    this.shots = options.shots || 1024;
    this.backend = options.backend || "statevector";
  }

  // Compute expectation values for (circuit, observable) pairs.
  // circuits: QuantumCircuit or array of QuantumCircuit
  // observables: SparsePauliOp or Pauli or array of these (same length as circuits)
  // parameter_values: array of {paramName: number} dicts, or null
  run(circuits, observables, parameterValues = null, options = {}) {
    const circuitList = Array.isArray(circuits) ? circuits : [circuits];
    const obsList = Array.isArray(observables) ? observables : [observables];
    if (circuitList.length !== obsList.length) {
      throw new Error("Estimator: number of circuits must match number of observables");
    }

    const values = [];
    const metadata = [];

    for (let i = 0; i < circuitList.length; i++) {
      let circuit = circuitList[i];
      let observable = obsList[i];

      // Bind parameters if provided
      if (parameterValues && parameterValues[i]) {
        circuit = circuit.bindParameters(parameterValues[i]);
      }

      // Compute the statevector
      const sv = Statevector.fromCircuit(circuit);

      // Convert observable to SparsePauliOp
      let spo;
      if (observable instanceof SparsePauliOp) {
        spo = observable;
      } else if (observable instanceof Pauli) {
        spo = SparsePauliOp.fromList([[observable.label, 1.0]]);
      } else if (typeof observable === "string") {
        spo = SparsePauliOp.fromList([[observable, 1.0]]);
      } else {
        throw new TypeError("Observable must be Pauli, SparsePauliOp, or string");
      }

      // Compute <ψ|O|ψ>
      const expVal = spo.expectationValue(sv);
      values.push(expVal);
      metadata.push({ shots: 0 });
    }

    return new EstimatorResult(values, metadata);
  }
}

// Sampler
export class Sampler extends BaseSampler {
  constructor(options = {}) {
    super(options);
    this.shots = options.shots || 1024;
    this.backend = options.backend || "statevector";
  }

  // Sample measurement outcomes from circuits.
  run(circuits, parameterValues = null, options = {}) {
    const circuitList = Array.isArray(circuits) ? circuits : [circuits];
    const shots = options.shots || this.shots;

    const quasiDists = [];
    const metadata = [];

    for (let i = 0; i < circuitList.length; i++) {
      let circuit = circuitList[i];

      // Bind parameters if provided
      if (parameterValues && parameterValues[i]) {
        circuit = circuit.bindParameters(parameterValues[i]);
      }

      // If the circuit has measurements, run the simulator
      const hasMeasure = circuit.data.some(ci => ci.operation.name === "measure");
      if (hasMeasure) {
        const result = _simulate(circuit, shots);
        const counts = result.getCounts();
        // Convert to quasi-probability distribution
        const quasi = {};
        const total = counts.shots;
        for (const [key, val] of counts.items()) {
          quasi[key] = val / total;
        }
        quasiDists.push(quasi);
        metadata.push({ shots });
      } else {
        // Sample directly from the statevector
        const sv = Statevector.fromCircuit(circuit);
        const samples = sv.sample(shots);
        const quasi = {};
        for (const [key, val] of Object.entries(samples)) {
          quasi[key] = val / shots;
        }
        quasiDists.push(quasi);
        metadata.push({ shots });
      }
    }

    return new SamplerResult(quasiDists, metadata);
  }
}

// V2 primitives (qiskit 1.0+ API)
// The V2 API uses a "primitive job" pattern: run() returns a job object that
// can be awaited for results. Since JS is single-threaded here, we return
// the result synchronously wrapped in a minimal job-like interface.

export class PrimitiveResultV2 {
  constructor(results) {
    this.results = results;
  }
}

export class EstimatorResultV2 {
  constructor(evs, stds, metadata = {}) {
    this.evs = evs;       // expectation values (array of numbers)
    this.stds = stds;     // standard deviations (array of numbers)
    this.metadata = metadata;
  }
}

export class SamplerResultV2 {
  constructor(pubResults, metadata = {}) {
    this.pubResults = pubResults;
    this.metadata = metadata;
  }
}

export class PrimitivePubResult {
  constructor(data, metadata = {}) {
    this.data = data;
    this.metadata = metadata;
  }
}

// EstimatorV2: takes a list of "PUBs" (Primitive Unified Bundles), each of
// the form [circuit, observables, parameter_values]. Returns expectation
// values (and standard deviations) for each PUB.
export class EstimatorV2 extends BaseEstimator {
  constructor(options = {}) {
    super(options);
    this.defaultPrecision = options.defaultPrecision || 0.0;
    this.defaultShots = options.defaultShots || 1024;
  }

  run(pubs, options = {}) {
    if (!Array.isArray(pubs)) pubs = [pubs];
    const precision = options.defaultPrecision || this.defaultPrecision;
    const shots = options.defaultShots || this.defaultShots;
    const results = [];
    for (const pub of pubs) {
      // pub is [circuit, observables, parameter_values]
      const [circuit, observables, parameterValues] = pub;
      let boundCircuit = circuit;
      if (parameterValues) {
        boundCircuit = circuit.bindParameters(parameterValues);
      }
      const sv = Statevector.fromCircuit(boundCircuit);
      // observables can be a single SparsePauliOp or a list.
      const obsList = Array.isArray(observables) ? observables : [observables];
      const evs = [];
      const stds = [];
      for (const obs of obsList) {
        let spo;
        if (obs instanceof SparsePauliOp) spo = obs;
        else if (obs instanceof Pauli) spo = SparsePauliOp.fromList([[obs.label, 1.0]]);
        else if (typeof obs === "string") spo = SparsePauliOp.fromList([[obs, 1.0]]);
        else throw new TypeError("EstimatorV2: observable must be Pauli, SparsePauliOp, or string");
        const expVal = spo.expectationValue(sv);
        const ev = typeof expVal === "number" ? expVal : expVal.re;
        evs.push(ev);
        // For exact statevector simulation, the standard deviation is 0.
        stds.push(0);
      }
      results.push(new EstimatorResultV2(evs, stds, { shots: 0 }));
    }
    return new PrimitiveResultV2(results);
  }
}

// SamplerV2: takes a list of PUBs [circuit, parameter_values, shots].
// Returns measurement bit arrays and quasi-probability distributions.
export class SamplerV2 extends BaseSampler {
  constructor(options = {}) {
    super(options);
    this.defaultShots = options.defaultShots || 1024;
  }

  run(pubs, options = {}) {
    if (!Array.isArray(pubs)) pubs = [pubs];
    const defaultShots = options.defaultShots || this.defaultShots;
    const results = [];
    for (const pub of pubs) {
      // pub is [circuit, parameter_values, shots]
      const [circuit, parameterValues, shots] = pub;
      let boundCircuit = circuit;
      if (parameterValues) {
        boundCircuit = circuit.bindParameters(parameterValues);
      }
      const numShots = shots || defaultShots;
      // Run the simulator with measurements.
      const hasMeasure = boundCircuit.data.some(ci => ci.operation.name === "measure");
      let quasi;
      let measBits;
      if (hasMeasure) {
        const result = _simulate(boundCircuit, numShots);
        const counts = result.getCounts();
        quasi = {};
        const total = counts.shots;
        for (const [key, val] of counts.items()) {
          quasi[key] = val / total;
        }
        measBits = Object.keys(quasi);
      } else {
        const sv = Statevector.fromCircuit(boundCircuit);
        const samples = sv.sample(numShots);
        quasi = {};
        for (const [key, val] of Object.entries(samples)) {
          quasi[key] = val / numShots;
        }
        measBits = Object.keys(quasi);
      }
      results.push(new PrimitivePubResult(
        { counts: quasi, shots: numShots },
        { shots: numShots },
      ));
    }
    return new PrimitiveResultV2(results);
  }
}
