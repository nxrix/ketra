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

describe("Pauli", () => {
  test("basic construction", () => {
    const p = new K.Pauli("XYZ");
    expect(p.numQubits).toBe(3);
    expect(p.weight()).toBe(3);
    const m = p.toMatrix();
    expect(m.rows).toBe(8);
  });

  test("Y matrix no double phase", () => {
    const y = new K.Pauli("Y");
    expect(matEq(y.toMatrix(), K.PAULI.Y, 1e-9)).toBe(true);
    const x = new K.Pauli("X");
    expect(matEq(x.toMatrix(), K.PAULI.X, 1e-9)).toBe(true);
    const z = new K.Pauli("Z");
    expect(matEq(z.toMatrix(), K.PAULI.Z, 1e-9)).toBe(true);
  });

  test("compose", () => {
    const x = new K.Pauli("X");
    const y = new K.Pauli("Y");
    const r = x.compose(y);
    expect(r.pauli.label).toBe("Z");
    expect(approxEq(r.phase.im, 1)).toBe(true); // +i
  });

  test("anticommutation", () => {
    const z = new K.Pauli("Z");
    const x = new K.Pauli("X");
    expect(z.commutes(x)).toBe(false);
    expect(z.commutes(z)).toBe(true);
  });

  test("evolve", () => {
    const p1 = new K.Pauli("X");
    const p2 = new K.Pauli("Y");
    const e = p1.evolve(p2);
    // X * Y (ignoring phase) = (1,0) XOR (1,1) = (0,1) = Z
    expect(e.label).toBe("Z");
  });
});

describe("SparsePauliOp", () => {
  test("construction", () => {
    const spo = K.SparsePauliOp.fromList([["ZZ", 1.0], ["ZI", 0.5], ["IZ", 0.5]]);
    expect(spo.numQubits).toBe(2);
    expect(spo.paulis.size).toBe(3);
    const m = spo.toMatrix();
    expect(m.rows).toBe(4);
  });

  test("expectation value", () => {
    const spo = K.SparsePauliOp.fromList([["ZZ", 1.0], ["ZI", 0.5], ["IZ", 0.5]]);
    const sv = K.Statevector.fromLabel("00");
    const exp = spo.expectationValue(sv);
    expect(approxEq(exp.re, 2.0, 1e-9)).toBe(true);
    const sv11 = K.Statevector.fromLabel("11");
    const exp11 = spo.expectationValue(sv11);
    expect(approxEq(exp11.re, 0.0, 1e-9)).toBe(true);
  });

  test("dot product", () => {
    const z = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const x = K.SparsePauliOp.fromList([["X", 1.0]]);
    const zx = z.dot(x);
    const expected = K.PAULI.Z.mul(K.PAULI.X);
    expect(matEq(zx.toMatrix(), expected, 1e-9)).toBe(true);
  });

  test("simplify", () => {
    const spo = K.SparsePauliOp.fromList([["Z", 1.0], ["Z", 2.0], ["X", 0.5]]);
    const simplified = spo.simplify();
    expect(simplified.paulis.size).toBe(2);
  });

  test("noncommutationGroups", () => {
    const spo = K.SparsePauliOp.fromList([["X", 1.0], ["Y", 1.0], ["Z", 1.0], ["I", 1.0]]);
    const groups = spo.noncommutationGroups();
    expect(groups.length).toBeGreaterThanOrEqual(2);
  });

  test("sort", () => {
    const spo = K.SparsePauliOp.fromList([["II", 1.0], ["IX", 1.0], ["XX", 1.0]]);
    const sorted = spo.sort();
    expect(sorted.paulis.get(0).weight()).toBeGreaterThanOrEqual(
      sorted.paulis.get(sorted.paulis.size - 1).weight()
    );
  });
});

describe("Statevector", () => {
  test("fromLabel", () => {
    const sv = K.Statevector.fromLabel("00");
    expect(sv.numQubits).toBe(2);
    expect(approxEq(sv.probabilities()[0], 1)).toBe(true);
  });

  test("evolve with circuit", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.cx(0, 1);
    const sv = K.Statevector.fromCircuit(qc);
    // Bell state
    expect(approxEq(sv.probabilities()[0], 0.5, 1e-9)).toBe(true);
    expect(approxEq(sv.probabilities()[3], 0.5, 1e-9)).toBe(true);
  });

  test("measure", () => {
    const plus = K.Statevector.fromLabel("+");
    const meas = plus.measure(0, () => 0);
    expect(meas.bit).toBe(0);
    expect(meas.statevector.equals(K.Statevector.fromLabel("0"), 1e-9)).toBe(true);
  });

  test("reset", () => {
    const one = K.Statevector.fromLabel("1");
    const reset = one.reset(0);
    expect(reset.equals(K.Statevector.fromLabel("0"), 1e-9)).toBe(true);
  });

  test("partialTrace", () => {
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    const pt = bell.partialTrace([0]);
    expect(pt.numQubits).toBe(1);
    expect(approxEq(K.purity(pt), 0.5, 1e-9)).toBe(true);
  });

  test("expandDims", () => {
    const sv = K.Statevector.fromLabel("0");
    const expanded = sv.expandDims(1);
    expect(expanded.numQubits).toBe(2);
    expect(expanded.equals(K.Statevector.fromLabel("00"), 1e-9)).toBe(true);
  });

  test("sample", () => {
    const sv = K.Statevector.fromLabel("0");
    const samples = sv.sample(100);
    expect(samples["0"]).toBe(100);
  });
});

describe("Operator", () => {
  test("fromCircuit", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const op = K.Operator.fromCircuit(qc);
    expect(matEq(op.data, K.CONSTANTS.H, 1e-9)).toBe(true);
  });

  test("sum and subtract", () => {
    const id = K.Operator.identity(1);
    const x = K.Operator.fromGate(K.standardGates.XGate);
    const sum = id.sum(x);
    expect(sum.data.rows).toBe(2);
    const diff = id.subtract(x);
    const prod = diff.data.mul(x.data);
    const expected = x.data.sub(id.data);
    expect(matEq(prod, expected, 1e-9)).toBe(true);
  });

  test("partialTrace", () => {
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    const op = bell.toOperator();
    const reduced = op.partialTrace([1]);
    const expected = K.ComplexMatrix.identity(2).scale(new K.Complex(0.5, 0));
    expect(matEq(reduced.data, expected, 1e-9)).toBe(true);
  });
});

describe("Clifford", () => {
  test("fromCircuit H matches Operator", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const cliff = K.Clifford.fromCircuit(qc);
    const op = K.Operator.fromCircuit(qc);
    // Compare up to global phase
    let phase = null, ok = true;
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
      const a = cliff.toMatrix().get(i, j); const b = op.data.get(i, j);
      if (a.abs() < 1e-9) { if (b.abs() > 1e-9) ok = false; }
      else { const r = b.div(a); if (phase === null) phase = r; else if (!phase.equals(r, 1e-6)) ok = false; }
    }
    expect(ok).toBe(true);
  });

  test("compose: S · S = Z", () => {
    const qcS = new K.QuantumCircuit(1); qcS.s(0);
    const cliffS = K.Clifford.fromCircuit(qcS);
    const ss = cliffS.compose(cliffS);
    // Compare to Z up to phase
    let phase = null, ok = true;
    for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) {
      const a = ss.toMatrix().get(i, j); const b = K.PAULI.Z.get(i, j);
      if (a.abs() < 1e-9) { if (b.abs() > 1e-9) ok = false; }
      else { const r = b.div(a); if (phase === null) phase = r; else if (!phase.equals(r, 1e-6)) ok = false; }
    }
    expect(ok).toBe(true);
  });

  test("adjoint: H · H^adj = I", () => {
    const qcH = new K.QuantumCircuit(1); qcH.h(0);
    const cliffH = K.Clifford.fromCircuit(qcH);
    const id = K.Clifford.fromLabel("0");
    expect(cliffH.adjoint().compose(cliffH).equals(id)).toBe(true);
  });

  test("power: S^4 = I", () => {
    const qcS = new K.QuantumCircuit(1); qcS.s(0);
    const cliffS = K.Clifford.fromCircuit(qcS);
    const id = K.Clifford.fromLabel("0");
    expect(cliffS.power(4).equals(id)).toBe(true);
  });

  test("StabilizerState toStatevector", () => {
    expect(K.StabilizerState.fromLabel("0").toStatevector().equals(
      K.Statevector.fromLabel("0"), 1e-9)).toBe(true);
    expect(K.StabilizerState.fromLabel("1").toStatevector().equals(
      K.Statevector.fromLabel("1"), 1e-9)).toBe(true);
    const qcBell = new K.QuantumCircuit(2); qcBell.h(0); qcBell.cx(0, 1);
    const cliffBell = K.Clifford.fromCircuit(qcBell);
    const svBell = new K.StabilizerState(cliffBell).toStatevector();
    const expected = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    expect(svBell.equals(expected, 1e-9)).toBe(true);
  });

  test("StabilizerState expectationValue", () => {
    expect(K.StabilizerState.fromLabel("1").expectationValue(new K.Pauli("Z")).re).toBe(-1);
    expect(K.StabilizerState.fromLabel("0").expectationValue(new K.Pauli("Z")).re).toBe(1);
    expect(K.StabilizerState.fromLabel("0").expectationValue(new K.Pauli("X")).re).toBe(0);
  });
});

describe("DensityMatrix", () => {
  test("from statevector", () => {
    const sv = K.Statevector.fromLabel("0");
    const dm = new K.DensityMatrix(sv.toOperator().data);
    expect(dm.numQubits).toBe(1);
  });

  test("purity of pure state", () => {
    const sv = K.Statevector.fromLabel("00");
    expect(approxEq(K.purity(sv.toOperator()), 1, 1e-9)).toBe(true);
  });
});
