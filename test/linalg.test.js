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

describe("Complex", () => {
  test("basic arithmetic", () => {
    const a = new K.Complex(3, 4);
    const b = new K.Complex(1, 2);
    expect(approxEq(a.abs(), 5)).toBe(true);
    expect(approxEq(a.add(b).re, 4) && approxEq(a.add(b).im, 6)).toBe(true);
    expect(approxEq(a.mul(b).re, -5) && approxEq(a.mul(b).im, 10)).toBe(true);
    expect(approxEq(a.conjugate().re, 3) && approxEq(a.conjugate().im, -4)).toBe(true);
    expect(approxEq(a.exp().re, Math.exp(3) * Math.cos(4))).toBe(true);
  });

  test("constants", () => {
    expect(K.Complex.ZERO.re).toBe(0);
    expect(K.Complex.ZERO.im).toBe(0);
    expect(approxEq(K.Complex.I.mul(K.Complex.I).re, -1)).toBe(true);
    expect(approxEq(K.Complex.I.mul(K.Complex.I).im, 0)).toBe(true);
  });

  test("division", () => {
    const a = new K.Complex(6, 8);
    const b = new K.Complex(2, 0);
    const q = a.div(b);
    expect(approxEq(q.re, 3) && approxEq(q.im, 4)).toBe(true);
  });
});

describe("ComplexVector", () => {
  test("norm and probabilities", () => {
    const v = new K.ComplexVector([
      new K.Complex(1, 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(0, 0),
    ]);
    expect(approxEq(v.norm(), 1)).toBe(true);
    const v2 = new K.ComplexVector([new K.Complex(1, 0), new K.Complex(1, 0)]);
    expect(approxEq(v2.norm(), Math.sqrt(2))).toBe(true);
    const probs = v2.probabilities();
    expect(approxEq(probs[0], 0.5) && approxEq(probs[1], 0.5)).toBe(true);
  });

  test("inner product", () => {
    const v1 = new K.ComplexVector([new K.Complex(1, 0), new K.Complex(0, 0)]);
    const v2 = new K.ComplexVector([new K.Complex(1, 0), new K.Complex(0, 0)]);
    const inner = v1.inner(v2);
    expect(approxEq(inner.re, 1)).toBe(true);
  });
});

describe("ComplexMatrix", () => {
  test("Pauli matrices", () => {
    expect(K.PAULI.X.rows).toBe(2);
    expect(approxEq(K.PAULI.X.get(0, 1).re, 1)).toBe(true);
    expect(approxEq(K.PAULI.X.get(1, 0).re, 1)).toBe(true);
    expect(K.PAULI.Y.isHermitian()).toBe(true);
    expect(K.PAULI.X.isUnitary()).toBe(true);
  });

  test("Hadamard", () => {
    const H = K.CONSTANTS.H;
    expect(H.isUnitary()).toBe(true);
    expect(H.isHermitian()).toBe(true);
  });

  test("determinant, trace, inverse", () => {
    const m = K.ComplexMatrix.fromRows([
      [new K.Complex(1), new K.Complex(2)],
      [new K.Complex(3), new K.Complex(4)],
    ]);
    expect(approxEq(m.det().re, -2)).toBe(true);
    expect(approxEq(m.trace().re, 5)).toBe(true);
    const inv = m.inverse();
    const prod = m.mul(inv);
    expect(matEq(prod, K.ComplexMatrix.identity(2), 1e-9)).toBe(true);
  });

  test("matrix exponential", () => {
    const iPiX2 = K.PAULI.X.scale(new K.Complex(0, Math.PI / 2));
    const expm = iPiX2.expm();
    const expected = K.PAULI.X.scale(new K.Complex(0, 1));
    expect(matEq(expm, expected, 1e-6)).toBe(true);
  });

  test("partial trace", () => {
    const id4 = K.ComplexMatrix.identity(4);
    const traced = id4.partialTrace(2, 0);
    expect(traced.rows).toBe(2);
    expect(approxEq(traced.get(0, 0).re, 2)).toBe(true);
  });

  test("SVD", () => {
    const m = K.ComplexMatrix.fromRows([
      [new K.Complex(1), new K.Complex(0)],
      [new K.Complex(0), new K.Complex(1)],
    ]);
    const { U, S, Vh } = m.svd();
    expect(approxEq(S[0], 1) && approxEq(S[1], 1)).toBe(true);
  });

  test("tensor product", () => {
    const X = K.PAULI.X;
    const XX = X.tensor(X);
    expect(XX.rows).toBe(4);
    expect(XX.cols).toBe(4);
  });

  test("conjugate and dagger", () => {
    const Y = K.PAULI.Y;
    const Ydag = Y.dagger();
    expect(matEq(Y, Ydag, 1e-9)).toBe(true);
  });
});
