import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("DAG passes", () => {
  test("collect1qRuns", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    qc.t(0);
    qc.t(0);
    const dag = K.circuitToDag(qc);
    const runs = K.collect1qRuns(dag);
    expect(runs.length).toBe(1);
    expect(runs[0].length).toBe(3);
  });

  test("collect2qRuns", () => {
    const qc = new K.QuantumCircuit(2);
    qc.cx(0, 1);
    qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    const runs = K.collect2qRuns(dag);
    expect(runs.length).toBeGreaterThanOrEqual(1);
  });

  test("consolidateBlocks", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    qc.t(0);
    qc.t(0);
    const dag = K.circuitToDag(qc);
    const result = K.consolidateBlocks(dag);
    expect(result).toBeDefined();
  });

  test("optimizeCliffords", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    qc.h(0);
    const dag = K.circuitToDag(qc);
    const result = K.optimizeCliffords(dag);
    expect(result).toBeDefined();
  });

  test("removeDiagonalGatesBeforeMeasure", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.z(0);
    qc.measure(0, 0);
    const dag = K.circuitToDag(qc);
    const result = K.removeDiagonalGatesBeforeMeasure(dag);
    expect(result).toBeDefined();
  });

  test("elidePermutations", () => {
    const qc = new K.QuantumCircuit(2);
    qc.swap(0, 1);
    qc.swap(0, 1);
    const dag = K.circuitToDag(qc);
    const result = K.elidePermutations(dag);
    expect(result).toBeDefined();
  });

  test("removeRedundantGates", () => {
    const qc = new K.QuantumCircuit(1);
    qc.id(0);
    qc.id(0);
    const dag = K.circuitToDag(qc);
    const result = K.removeRedundantGates(dag);
    expect(result).toBeDefined();
  });

  test("commutativeCancellation", () => {
    const qc = new K.QuantumCircuit(1);
    qc.rz(0.5, 0);
    qc.rz(-0.5, 0);
    const dag = K.circuitToDag(qc);
    const result = K.commutativeCancellation(dag);
    expect(result).toBeDefined();
  });

  test("templateOptimization", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    const dag = K.circuitToDag(qc);
    const result = K.templateOptimization(dag);
    expect(result).toBeDefined();
  });
});

describe("Visualization", () => {
  test("drawCircuit text output", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.cx(0, 1);
    const text = K.drawCircuit(qc, "text");
    expect(typeof text).toBe("string");
    expect(text.length).toBeGreaterThan(0);
  });

  test("plot_histogram returns data", () => {
    const counts = { "00": 512, "11": 512 };
    const result = K.plot_histogram(counts);
    expect(result).toBeDefined();
  });

  test("plot_bloch_vector returns data", () => {
    const result = K.plot_bloch_vector([0, 0, 1]);
    expect(result).toBeDefined();
  });

  test("render_histogram_ascii", () => {
    const counts = { "00": 512, "11": 512 };
    const text = K.render_histogram_ascii(counts);
    expect(typeof text).toBe("string");
  });

  test("render_bloch_ascii", () => {
    const text = K.render_bloch_ascii([0, 0, 1]);
    expect(typeof text).toBe("string");
  });
});
