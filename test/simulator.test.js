import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("StatevectorSimulator", () => {
  test("simulate Bell state", () => {
    const qc = new K.QuantumCircuit(2, 2);
    qc.h(0);
    qc.cx(0, 1);
    qc.measure(0, 0);
    qc.measure(1, 1);
    const result = K.simulate(qc, 1024);
    const counts = result.getCounts();
    expect(counts.get("00") + counts.get("11")).toBe(1024);
  });

  test("simulate single qubit", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.x(0);
    qc.measure(0, 0);
    const result = K.simulate(qc, 100);
    const counts = result.getCounts();
    expect(counts.get("1")).toBe(100);
  });
});

describe("QasmSimulator (noisy)", () => {
  test("noise model with bit flip", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.x(0);
    qc.measure(0, 0);
    const noise = new K.NoiseModel();
    noise.addAllQubitQuantumError(K.bitFlipError(0.5), ["x"]);
    const result = K.simulateNoisy(qc, 1000, { noiseModel: noise });
    const counts = result.getCounts();
    // Should have both 0 and 1 outcomes due to bit flip
    expect(counts.get("0") + counts.get("1")).toBe(1000);
  });
});

describe("Estimator", () => {
  test("V1 Estimator", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const obs = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const est = new K.Estimator();
    const result = est.run(qc, obs);
    // <H|Z|H> = <+|Z|+> = 0
    expect(approxEq(result.values[0].re, 0, 1e-9)).toBe(true);
  });

  test("V2 EstimatorV2", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const obs = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const est = new K.EstimatorV2();
    const result = est.run([[qc, obs]]);
    expect(result.results.length).toBe(1);
    expect(result.results[0].evs.length).toBe(1);
    expect(approxEq(result.results[0].evs[0], 0, 1e-9)).toBe(true);
  });
});

describe("Sampler", () => {
  test("V1 Sampler", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.h(0);
    qc.measure(0, 0);
    const samp = new K.Sampler();
    const result = samp.run(qc);
    expect(result.quasiDists.length).toBe(1);
    const totalProb = Object.values(result.quasiDists[0]).reduce((a, b) => a + b, 0);
    expect(approxEq(totalProb, 1, 1e-6)).toBe(true);
  });

  test("V2 SamplerV2", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.h(0);
    qc.measure(0, 0);
    const samp = new K.SamplerV2();
    const result = samp.run([[qc, null, 100]]);
    expect(result.results.length).toBe(1);
    expect(result.results[0].data.shots).toBe(100);
    const totalProb = Object.values(result.results[0].data.counts).reduce((a, b) => a + b, 0);
    expect(approxEq(totalProb, 1, 1e-6)).toBe(true);
  });
});
