import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

function matEq(m1, m2, tol = 1e-9) {
  if (m1.rows !== m2.rows || m1.cols !== m2.cols) return false;
  for (let i = 0; i < m1.rows; i++) {
    for (let j = 0; j < m1.cols; j++) {
      if (!approxEq(m1.get(i, j).re, m2.get(i, j).re, tol)) return false;
      if (!approxEq(m1.get(i, j).im, m2.get(i, j).im, tol)) return false;
    }
  }
  return true;
}

describe("QuantumCircuit construction", () => {
  test("basic circuit", () => {
    const qc = new K.QuantumCircuit(3, 2);
    expect(qc.numQubits).toBe(3);
    expect(qc.numClbits).toBe(2);
    qc.h(0);
    qc.cx(0, 1);
    qc.x(2);
    qc.measure(0, 0);
    qc.measure(1, 1);
    expect(qc.data.length).toBe(5);
    expect(qc.depth()).toBeGreaterThanOrEqual(2);
    const ops = qc.countOps();
    expect(ops["h"]).toBe(1);
    expect(ops["cx"]).toBe(1);
    expect(ops["x"]).toBe(1);
    expect(ops["measure"]).toBe(2);
  });

  test("standard gates", () => {
    const qc = new K.QuantumCircuit(3, 3);
    const gates = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "sx", "sxdg", "id"];
    for (const g of gates) qc[g](0);
    const twoQubit = [["cx", 0, 1], ["cy", 0, 1], ["cz", 0, 1], ["ch", 0, 1], ["csx", 0, 1], ["swap", 0, 1], ["iswap", 0, 1], ["dcx", 0, 1]];
    for (const [g, a, b] of twoQubit) qc[g](a, b);
    const threeQubit = [["ccx", 0, 1, 2], ["cswap", 0, 1, 2]];
    for (const [g, a, b, c] of threeQubit) qc[g](a, b, c);
    expect(qc.data.length).toBeGreaterThan(0);
  });

  test("parameterized gates", () => {
    const qc = new K.QuantumCircuit(2);
    qc.rx(Math.PI / 2, 0);
    qc.ry(Math.PI / 4, 0);
    qc.rz(Math.PI / 8, 0);
    qc.rxx(Math.PI / 2, 0, 1);
    qc.ryy(Math.PI / 4, 0, 1);
    qc.rzz(Math.PI / 8, 0, 1);
    qc.p(Math.PI / 3, 0);
    qc.u(Math.PI / 2, Math.PI / 4, Math.PI / 8, 0);
    qc.u1(Math.PI / 4, 0);
    qc.u2(Math.PI / 4, Math.PI / 8, 0);
    qc.u3(Math.PI / 2, Math.PI / 4, Math.PI / 8, 0);
    expect(qc.data.length).toBe(11);
  });

  test("controlled parameterized gates", () => {
    const qc = new K.QuantumCircuit(2);
    qc.cp(Math.PI / 4, 0, 1);
    qc.crx(Math.PI / 4, 0, 1);
    qc.cry(Math.PI / 4, 0, 1);
    qc.crz(Math.PI / 4, 0, 1);
    qc.cu(Math.PI / 4, Math.PI / 8, Math.PI / 16, Math.PI / 32, 0, 1);
    qc.cu1(Math.PI / 4, 0, 1);
    qc.cu3(Math.PI / 4, Math.PI / 8, Math.PI / 16, 0, 1);
    expect(qc.data.length).toBe(7);
  });

  test("multi-controlled gates", () => {
    const qc = new K.QuantumCircuit(4);
    qc.mcx([0, 1, 2], 3);
    qc.mcy([0, 1, 2], 3);
    qc.mcz([0, 1, 2], 3);
    qc.mcu1(Math.PI / 4, [0, 1, 2], 3);
    qc.mcp(Math.PI / 4, [0, 1, 2], 3);
    qc.mcrx(Math.PI / 4, [0, 1, 2], 3);
    qc.mcry(Math.PI / 4, [0, 1, 2], 3);
    qc.mcrz(Math.PI / 4, [0, 1, 2], 3);
    expect(qc.data.length).toBe(8);
  });
});

describe("Circuit methods", () => {
  test("unitary gate", () => {
    const qc = new K.QuantumCircuit(1);
    qc.unitary(K.PAULI.X, 0);
    const op = K.Operator.fromCircuit(qc);
    expect(matEq(op.data, K.PAULI.X, 1e-9)).toBe(true);
  });

  test("initialize with label", () => {
    const qc = new K.QuantumCircuit(2);
    qc.initialize("11");
    const sv = K.Statevector.fromCircuit(qc);
    expect(sv.equals(K.Statevector.fromLabel("11"), 1e-9)).toBe(true);
  });

  test("initialize with amplitudes", () => {
    const qc = new K.QuantumCircuit(1);
    qc.initialize([1 / Math.sqrt(2), 1 / Math.sqrt(2)], 0);
    const sv = K.Statevector.fromCircuit(qc);
    expect(sv.equals(K.Statevector.fromLabel("+"), 1e-9)).toBe(true);
  });

  test("diagonal gate", () => {
    const qc = new K.QuantumCircuit(1);
    qc.diagonal([0, Math.PI], 0);
    const op = K.Operator.fromCircuit(qc);
    expect(matEq(op.data, K.PAULI.Z, 1e-9)).toBe(true);
  });

  test("permutation gate", () => {
    const qc = new K.QuantumCircuit(2);
    qc.permute([1, 0]);
    const op = K.Operator.fromCircuit(qc);
    expect(matEq(op.data, K.standardGates.SwapGate.toMatrix(), 1e-9)).toBe(true);
  });

  test("hamiltonian gate", () => {
    const qc = new K.QuantumCircuit(1);
    qc.hamiltonian(K.PAULI.X, Math.PI / 2, 0);
    const op = K.Operator.fromCircuit(qc);
    // e^{-i pi X / 2} = iX (up to phase)
    const expected = K.PAULI.X.scale(new K.Complex(0, 1));
    let phase = null, ok = true;
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
      const a = op.data.get(i, j); const b = expected.get(i, j);
      if (a.abs() < 1e-9) { if (b.abs() > 1e-9) ok = false; }
      else { const r = b.div(a); if (phase === null) phase = r; else if (!phase.equals(r, 1e-6)) ok = false; }
    }
    expect(ok).toBe(true);
  });

  test("toGate and toInstruction", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const g = qc.toGate("my_h");
    expect(g.numQubits).toBe(1);
    expect(matEq(g.toMatrix(), K.CONSTANTS.H, 1e-9)).toBe(true);
    const instr = qc.toInstruction();
    expect(matEq(instr.toMatrix(), K.CONSTANTS.H, 1e-9)).toBe(true);
  });

  test("repeat and power", () => {
    const qc = new K.QuantumCircuit(1); qc.x(0);
    const rep = qc.repeat(2);
    expect(rep.data.length).toBe(2);
    expect(matEq(K.Operator.fromCircuit(rep).data, K.ComplexMatrix.identity(2), 1e-9)).toBe(true);
    const pow = qc.power(3);
    expect(pow.data.length).toBe(3);
    expect(matEq(K.Operator.fromCircuit(pow).data, K.PAULI.X, 1e-9)).toBe(true);
  });

  test("decompose", () => {
    const sub = new K.QuantumCircuit(2); sub.h(0); sub.cx(0, 1);
    const qc = new K.QuantumCircuit(2);
    qc.append(sub.toGate("bell"), [0, 1]);
    const dec = qc.decompose();
    expect(dec.data.length).toBe(2);
    expect(dec.data[0].operation.name).toBe("h");
    expect(dec.data[1].operation.name).toBe("cx");
  });

  test("qasm export", () => {
    const qc = new K.QuantumCircuit(2, 2); qc.h(0); qc.cx(0, 1);
    const s = qc.qasm();
    expect(s.includes("OPENQASM 2.0")).toBe(true);
    expect(s.includes("h")).toBe(true);
  });

  test("globalPhase setter", () => {
    const qc = new K.QuantumCircuit(1);
    qc.globalPhase = Math.PI / 2;
    expect(approxEq(qc.globalPhase, Math.PI / 2, 1e-9)).toBe(true);
  });

  test("measureActive on multi-qubit gates", () => {
    const qc = new K.QuantumCircuit(3);
    qc.cx(0, 2);
    const measured = qc.measureActive(false);
    let measureCount = 0;
    for (const ci of measured.data) {
      if (ci.operation.name === "measure") measureCount++;
    }
    expect(measureCount).toBe(2);
  });

  test("parameters with ParameterExpression", () => {
    const qc = new K.QuantumCircuit(2);
    const p = new K.Parameter("gamma");
    const expr = p.mul(2);
    qc.rz(expr, 0);
    const params = qc.parameters;
    expect(params.size).toBe(1);
    expect(Array.from(params)[0].name).toBe("gamma");
  });
});

describe("Circuit analysis", () => {
  test("depth", () => {
    const qc = new K.QuantumCircuit(3);
    qc.h(0);
    qc.cx(0, 1);
    qc.x(2);
    expect(qc.depth()).toBeGreaterThan(0);
  });

  test("numNonlocalGates", () => {
    const qc = new K.QuantumCircuit(3);
    qc.h(0);
    qc.cx(0, 1);
    qc.ccx(0, 1, 2);
    expect(qc.numNonlocalGates()).toBe(2);
  });

  test("numTensorFactors", () => {
    const qc = new K.QuantumCircuit(4);
    qc.h(0);
    qc.cx(0, 1);
    qc.h(2);
    qc.x(3);
    // Two separate groups: {0,1} and {2}, {3}
    expect(qc.numTensorFactors()).toBeGreaterThanOrEqual(2);
  });

  test("inverse", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    const inv = qc.inverse();
    // H * H = I
    const composed = new K.QuantumCircuit(1);
    composed.compose(qc);
    composed.compose(inv);
    const op = K.Operator.fromCircuit(composed);
    expect(matEq(op.data, K.ComplexMatrix.identity(2), 1e-9)).toBe(true);
  });
});
