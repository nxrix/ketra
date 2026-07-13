/**
 * result.js - Result and Counts classes.
 *
 * */

export class Counts {
  constructor(data, hexKeys = false) {
    this._data = {};
    this.hexKeys = hexKeys;
    if (data) {
      for (const k in data) {
        this._data[k] = data[k];
      }
    }
  }

  get shots() {
    let total = 0;
    for (const k in this._data) total += this._data[k];
    return total;
  }

  get(key) { return this._data[key] || 0; }
  set(key, value) { this._data[key] = value; }
  keys() { return Object.keys(this._data); }
  values() { return Object.values(this._data); }
  items() { return Object.entries(this._data); }

  most_frequent() {
    let bestKey = null, bestCount = -1;
    for (const k in this._data) {
      if (this._data[k] > bestCount) { bestCount = this._data[k]; bestKey = k; }
    }
    return bestKey;
  }

  hex_outcomes() {
    const out = {};
    for (const k in this._data) {
      out["0x" + parseInt(k.split("").reverse().join(""), 2).toString(16)] = this._data[k];
    }
    return out;
  }

  toDict() { return Object.assign({}, this._data); }

  [Symbol.iterator]() {
    return Object.entries(this._data)[Symbol.iterator]();
  }

  marginal_counts(indices) {
    const idxs = Array.isArray(indices) ? indices : [indices];
    const out = {};
    for (const k in this._data) {
      let c = 0;
      for (let i = 0; i < k.length; i++) {
        const bit = parseInt(k[k.length - 1 - i], 10);
        c |= bit << i;
      }
      let new_c = 0;
      for (let i = 0; i < idxs.length; i++) {
        const bit = (c >> idxs[i]) & 1;
        new_c |= bit << i;
      }
      let newKey = "";
      for (let i = idxs.length - 1; i >= 0; i--) {
        newKey += ((new_c >> i) & 1).toString();
      }
      out[newKey] = (out[newKey] || 0) + this._data[k];
    }
    return new Counts(out, this.hexKeys);
  }

  clear() { this._data = {}; }
}

export class Result {
  constructor(kwargs) {
    this.backendName = kwargs.backendName || "unknown";
    this.backendVersion = kwargs.backendVersion || "0.0.0";
    this.qobjId = kwargs.qobjId || "qobj";
    this.jobId = kwargs.jobId || `job_${Date.now()}`;
    this.success = kwargs.success !== undefined ? kwargs.success : true;
    this.results = kwargs.results || [];
    this.date = kwargs.date || new Date().toISOString();
    this.status = kwargs.status || "COMPLETED";
  }

  getCounts(experiment_id = 0) {
    const exp = this.results[experiment_id];
    if (!exp) throw new Error(`No experiment at index ${experiment_id}`);
    return exp.data.counts;
  }

  get_statevector(experiment_id = 0) {
    const exp = this.results[experiment_id];
    if (!exp) throw new Error(`No experiment at index ${experiment_id}`);
    return exp.data.statevector;
  }

  get_memory(experiment_id = 0) {
    const exp = this.results[experiment_id];
    if (!exp) throw new Error(`No experiment at index ${experiment_id}`);
    return exp.data.memory;
  }

  get_experiment(experiment_id = 0) {
    return this.results[experiment_id];
  }

  toDict() {
    return {
      backendName: this.backendName,
      backendVersion: this.backendVersion,
      qobjId: this.qobjId,
      jobId: this.jobId,
      success: this.success,
      results: this.results,
      date: this.date,
      status: this.status,
    };
  }
}
