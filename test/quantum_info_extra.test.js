import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("ScalarOp", () => {
  test("construction", () => {
    const s = new K.ScalarOp(2, 3);
    expect(s.numQubits).toBe(2);
    expect(s.coeff.re).toBe(3);
  });

  test("toMatrix", () => {
    const s = new K.ScalarOp(1, 2);
    const m = s.toMatrix();
    expect(approxEq(m.get(0, 0).re, 2)).toBe(true);
    expect(approxEq(m.get(1, 1).re, 2)).toBe(true);
  });

  test("compose", () => {
    const s1 = new K.ScalarOp(1, 2);
    const s2 = new K.ScalarOp(1, 3);
    const s3 = s1.compose(s2);
    expect(s3.coeff.re).toBe(6);
  });

  test("tensor", () => {
    const s1 = new K.ScalarOp(1, 2);
    const s2 = new K.ScalarOp(1, 3);
    const s3 = s1.tensor(s2);
    expect(s3.numQubits).toBe(2);
    expect(s3.coeff.re).toBe(6);
  });

  test("isUnitary", () => {
    const s = new K.ScalarOp(1, 1);
    expect(s.isUnitary()).toBe(true);
    const s2 = new K.ScalarOp(1, 2);
    expect(s2.isUnitary()).toBe(false);
  });

  test("trace and det", () => {
    const s = new K.ScalarOp(1, 2);
    expect(s.trace().re).toBe(4);
    expect(s.det().re).toBe(4);
  });

  test("fromOperator", () => {
    const op = K.Operator.identity(1);
    const s = K.ScalarOp.fromOperator(op);
    expect(s.coeff.re).toBe(1);
  });

  test("adjoint, conjugate, transpose", () => {
    const s = new K.ScalarOp(1, new K.Complex(1, 2));
    expect(s.adjoint().coeff.im).toBe(-2);
    expect(s.conjugate().coeff.im).toBe(-2);
    expect(s.transpose().coeff.im).toBe(2);
  });

  test("applyToVector", () => {
    const s = new K.ScalarOp(1, 2);
    const v = new K.ComplexVector([new K.Complex(1, 0), new K.Complex(0, 0)]);
    const result = s.applyToVector(v);
    expect(approxEq(result.data[0].re, 2)).toBe(true);
  });

  test("expectationValue", () => {
    const s = new K.ScalarOp(1, 2);
    const sv = K.Statevector.fromLabel("0");
    const ev = s.expectationValue(sv);
    expect(ev.re).toBe(2);
  });

  test("equals", () => {
    const s1 = new K.ScalarOp(1, 2);
    const s2 = new K.ScalarOp(1, 2);
    const s3 = new K.ScalarOp(1, 3);
    expect(s1.equals(s2)).toBe(true);
    expect(s1.equals(s3)).toBe(false);
  });

  test("toString", () => {
    const s = new K.ScalarOp(1, 2);
    expect(typeof s.toString()).toBe("string");
  });
});

describe("SchmidtDecomposition", () => {
  test("separable state has rank 1", () => {
    const sv = K.Statevector.fromLabel("00");
    const schmidt = new K.SchmidtDecomposition(sv);
    expect(schmidt.schmidtRank).toBe(1);
  });

  test("Bell state has rank 2", () => {
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    const schmidt = new K.SchmidtDecomposition(bell);
    expect(schmidt.schmidtRank).toBe(2);
  });

  test("entropy of Bell state = 1", () => {
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    const schmidt = new K.SchmidtDecomposition(bell);
    expect(approxEq(schmidt.entropy(), 1, 1e-9)).toBe(true);
  });

  test("isEntangled", () => {
    const sv = K.Statevector.fromLabel("00");
    const s1 = new K.SchmidtDecomposition(sv);
    expect(s1.isEntangled()).toBe(false);
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    const s2 = new K.SchmidtDecomposition(bell);
    expect(s2.isEntangled()).toBe(true);
  });
});

describe("DensityMatrix", () => {
  test("from statevector", () => {
    const sv = K.Statevector.fromLabel("0");
    const dm = K.DensityMatrix.fromStatevector(sv);
    expect(dm.numQubits).toBe(1);
  });

  test("from label", () => {
    const sv = K.Statevector.fromLabel("00");
    const dm = K.DensityMatrix.fromStatevector(sv);
    expect(dm.numQubits).toBe(2);
  });

  test("evolve under unitary", () => {
    const sv = K.Statevector.fromLabel("0");
    const dm = K.DensityMatrix.fromStatevector(sv);
    const xOp = new K.Operator(K.PAULI.X);
    const evolved = dm.evolve(xOp);
    expect(evolved).toBeDefined();
  });

  test("probabilities", () => {
    const sv = K.Statevector.fromLabel("0");
    const dm = K.DensityMatrix.fromStatevector(sv);
    const probs = dm.probabilities();
    expect(approxEq(probs[0], 1)).toBe(true);
  });

  test("toStatevector for pure state", () => {
    const sv = K.Statevector.fromLabel("0");
    const dm = K.DensityMatrix.fromStatevector(sv);
    const out = dm.toStatevector();
    expect(out).toBeDefined();
  });

  test("tensor", () => {
    const dm1 = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    const dm2 = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("1"));
    const dm = dm1.tensor(dm2);
    expect(dm.numQubits).toBe(2);
  });

  test("expand", () => {
    const dm1 = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    const dm2 = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("1"));
    const dm = dm1.expand(dm2);
    expect(dm.numQubits).toBe(2);
  });

  test("partialTrace", () => {
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    const dm = K.DensityMatrix.fromStatevector(bell);
    const reduced = dm.partialTrace([1]);
    expect(reduced.numQubits).toBe(1);
  });

  test("expectationValue", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    const ev = dm.expectationValue(K.PAULI.Z);
    expect(approxEq(ev.re, 1)).toBe(true);
  });

  test("purity", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    expect(approxEq(dm.purity(), 1, 1e-9)).toBe(true);
  });

  test("trace", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    expect(approxEq(dm.trace().re, 1, 1e-9)).toBe(true);
  });
});
