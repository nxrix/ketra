import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("Opflow expectations", () => {
  test("MatrixExpectation", () => {
    const me = new K.MatrixExpectation();
    expect(me).toBeDefined();
    const op = K.PauliSumOp.fromList([["Z", 1.0]]);
    const state = new K.StateFn("0");
    // <0|Z|0> = 1
    const result = me.compute(state, op);
    expect(result).toBeDefined();
  });

  test("PauliExpectation", () => {
    const pe = new K.PauliExpectation();
    expect(pe).toBeDefined();
    const op = K.PauliSumOp.fromList([["Z", 1.0]]);
    const state = new K.StateFn("0");
    const result = pe.compute(state, op);
    expect(result).toBeDefined();
  });

  test("CircuitSampler", () => {
    const cs = new K.CircuitSampler();
    expect(cs).toBeDefined();
  });

  test("get_expectation", () => {
    const fn = K.get_expectation;
    expect(typeof fn).toBe("function");
  });
});

describe("Opflow ListOp operations", () => {
  test("ListOp add", () => {
    const p1 = new K.PauliOp("X");
    const p2 = new K.PauliOp("Z");
    const sum = p1.add(p2);
    const m = sum.toMatrix();
    expect(m.rows).toBe(2);
  });

  test("ListOp compose", () => {
    const p1 = new K.PauliOp("X");
    const p2 = new K.PauliOp("Z");
    const composed = p1.compose(p2);
    expect(composed).toBeDefined();
  });

  test("ListOp tensor", () => {
    const p1 = new K.PauliOp("X");
    const p2 = new K.PauliOp("Z");
    const tensored = p1.tensor(p2);
    const m = tensored.toMatrix();
    expect(m.rows).toBe(4);
  });

  test("ListOp pow", () => {
    const p = new K.PauliOp("X");
    const pow = p.pow(2);
    expect(pow).toBeDefined();
  });

  test("ListOp scale", () => {
    const p = new K.PauliOp("X");
    const scaled = p.scale(2);
    expect(scaled).toBeDefined();
  });

  test("ListOp adjoint", () => {
    const p = new K.PauliOp("X");
    const adj = p.adjoint();
    expect(adj).toBeDefined();
  });

  test("ListOp reduce", () => {
    const p1 = new K.PauliOp("X");
    const p2 = new K.PauliOp("X");
    const sum = p1.add(p2);
    const reduced = sum.reduce();
    expect(reduced).toBeDefined();
  });

  test("StateFn adjoint", () => {
    const sf = new K.StateFn("0");
    const adj = sf.adjoint();
    expect(adj.isMeasurement).toBe(true);
  });

  test("StateFn compose with operator", () => {
    const sf = new K.StateFn("0");
    const op = new K.PauliOp("X");
    const result = sf.compose(op);
    expect(result).toBeDefined();
  });

  test("EvolvedOp toMatrix", () => {
    const p = new K.PauliOp("Z");
    const e = new K.EvolvedOp(p, K.Complex.ONE, 1.0);
    const m = e.toMatrix();
    expect(m.rows).toBe(2);
  });

  test("EvolvedOp adjoint", () => {
    const p = new K.PauliOp("Z");
    const e = new K.EvolvedOp(p, K.Complex.ONE, 1.0);
    const adj = e.adjoint();
    expect(adj.time).toBe(-1.0);
  });
});

describe("CircuitOp and CircuitStateFn", () => {
  test("CircuitOp toCircuit", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const c = new K.CircuitOp(qc);
    expect(c.toCircuit()).toBe(qc);
  });

  test("CircuitOp adjoint", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const c = new K.CircuitOp(qc);
    const adj = c.adjoint();
    expect(adj).toBeDefined();
  });

  test("CircuitStateFn toCircuit", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const csf = new K.CircuitStateFn(qc);
    expect(csf.toCircuit()).toBe(qc);
  });
});

describe("PauliSumOp", () => {
  test("reduce", () => {
    const p = K.PauliSumOp.fromList([["Z", 1.0], ["Z", 2.0]]);
    const reduced = p.reduce();
    expect(reduced).toBeDefined();
  });

  test("toList", () => {
    const p = K.PauliSumOp.fromList([["Z", 1.0]]);
    const list = p.toList();
    expect(list.length).toBe(1);
  });

  test("adjoint", () => {
    const p = K.PauliSumOp.fromList([["Z", 1.0]]);
    const adj = p.adjoint();
    expect(adj).toBeDefined();
  });

  test("isHermitian", () => {
    const p = K.PauliSumOp.fromList([["Z", 1.0]]);
    expect(p.isHermitian()).toBe(true);
  });
});
