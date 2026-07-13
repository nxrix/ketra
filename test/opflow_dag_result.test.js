import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("Opflow classes", () => {
  test("PauliOp", () => {
    const p = new K.PauliOp("Z");
    expect(p.numQubits).toBe(1);
    const m = p.toMatrix();
    expect(m.rows).toBe(2);
  });

  test("PauliSumOp", () => {
    const p = K.PauliSumOp.fromList([["Z", 1.0], ["X", 0.5]]);
    expect(p.numQubits).toBe(1);
    const m = p.toMatrix();
    expect(m.rows).toBe(2);
  });

  test("MatrixOp", () => {
    const m = new K.MatrixOp(K.PAULI.X);
    expect(m.numQubits).toBe(1);
  });

  test("CircuitOp", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const c = new K.CircuitOp(qc);
    expect(c.numQubits).toBe(1);
  });

  test("StateFn", () => {
    const sf = new K.StateFn("0");
    expect(sf.numQubits).toBe(1);
    const m = sf.toMatrix();
    expect(m.rows).toBe(2);
  });

  test("CircuitStateFn", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const csf = new K.CircuitStateFn(qc);
    expect(csf.numQubits).toBe(1);
  });

  test("ListOp", () => {
    const p1 = new K.PauliOp("X");
    const p2 = new K.PauliOp("Z");
    const sum = p1.add(p2);
    expect(sum).toBeDefined();
    const m = sum.toMatrix();
    expect(m.rows).toBe(2);
  });

  test("EvolvedOp", () => {
    const p = new K.PauliOp("Z");
    const e = new K.EvolvedOp(p, Complex.ONE, 1.0);
    const m = e.toMatrix();
    expect(m.rows).toBe(2);
  });
});

describe("Opflow expectations", () => {
  test("MatrixExpectation", () => {
    const me = new K.MatrixExpectation();
    expect(me).toBeDefined();
  });

  test("PauliExpectation", () => {
    const pe = new K.PauliExpectation();
    expect(pe).toBeDefined();
  });
});

import { Complex } from "../src/math/linalg.js";

describe("DAG circuit", () => {
  test("DAGCircuit construction", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    expect(dag.qubits.length).toBe(2);
    const opNodes = dag.topologicalOpNodes();
    expect(opNodes.length).toBe(2);
  });

  test("dagToCircuit round-trip", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    const back = K.dagToCircuit(dag);
    expect(back.numQubits).toBe(2);
    expect(back.data.length).toBe(2);
  });

  test("DAG passes", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    expect(K.collect1qRuns(dag)).toBeDefined();
    expect(K.collect2qRuns(dag)).toBeDefined();
  });
});

describe("Result and Counts", () => {
  test("Counts", () => {
    const c = new K.Counts({ "00": 512, "11": 512 });
    expect(c.get("00")).toBe(512);
    expect(c.get("01")).toBe(0);
  });

  test("Result", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.x(0);
    qc.measure(0, 0);
    const result = K.simulate(qc, 100);
    const counts = result.getCounts();
    expect(counts.get("1")).toBe(100);
  });
});
