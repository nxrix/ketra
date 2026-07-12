/* Copyright nxrix, 2026 */
var __defProp = Object.defineProperty;
var __export = (target, all) => {
  for (var name in all)
    __defProp(target, name, { get: all[name], enumerable: true });
};

// src/math/linalg.js
var Complex = class _Complex {
  constructor(real, imag) {
    if (typeof real !== "number") {
      throw new TypeError("Complex real part must be a number");
    }
    this.re = real;
    this.im = typeof imag === "number" ? imag : 0;
  }
  static from(v) {
    if (v instanceof _Complex) return v;
    if (typeof v === "number") return new _Complex(v, 0);
    throw new TypeError("Cannot coerce value to Complex");
  }
  static get ZERO() {
    return new _Complex(0, 0);
  }
  static get ONE() {
    return new _Complex(1, 0);
  }
  static get I() {
    return new _Complex(0, 1);
  }
  add(other) {
    const o = _Complex.from(other);
    return new _Complex(this.re + o.re, this.im + o.im);
  }
  sub(other) {
    const o = _Complex.from(other);
    return new _Complex(this.re - o.re, this.im - o.im);
  }
  mul(other) {
    const o = _Complex.from(other);
    return new _Complex(
      this.re * o.re - this.im * o.im,
      this.re * o.im + this.im * o.re
    );
  }
  div(other) {
    const o = _Complex.from(other);
    const denom = o.re * o.re + o.im * o.im;
    if (denom === 0) throw new Error("Complex division by zero");
    return new _Complex(
      (this.re * o.re + this.im * o.im) / denom,
      (this.im * o.re - this.re * o.im) / denom
    );
  }
  conjugate() {
    return new _Complex(this.re, -this.im);
  }
  negate() {
    return new _Complex(-this.re, -this.im);
  }
  // Multiply by a real scalar (convenience method)
  scale(s) {
    return new _Complex(this.re * s, this.im * s);
  }
  abs() {
    return Math.hypot(this.re, this.im);
  }
  abs2() {
    return this.re * this.re + this.im * this.im;
  }
  arg() {
    return Math.atan2(this.im, this.re);
  }
  exp() {
    const r4 = Math.exp(this.re);
    return new _Complex(r4 * Math.cos(this.im), r4 * Math.sin(this.im));
  }
  pow(n) {
    const r4 = Math.pow(this.abs(), n);
    const theta = this.arg() * n;
    return new _Complex(r4 * Math.cos(theta), r4 * Math.sin(theta));
  }
  sqrt() {
    const r4 = this.abs();
    const theta = this.arg() / 2;
    const sr = Math.sqrt(r4);
    return new _Complex(sr * Math.cos(theta), sr * Math.sin(theta));
  }
  log() {
    return new _Complex(Math.log(this.abs()), this.arg());
  }
  equals(other, tol) {
    const o = _Complex.from(other);
    const eps = typeof tol === "number" ? tol : 1e-12;
    return Math.abs(this.re - o.re) < eps && Math.abs(this.im - o.im) < eps;
  }
  toString() {
    const sign = this.im < 0 ? "-" : "+";
    return `${this.re}${sign}${Math.abs(this.im)}i`;
  }
  toJSON() {
    return { re: this.re, im: this.im };
  }
};
var ComplexVector = class _ComplexVector {
  constructor(values) {
    if (Array.isArray(values) && values.length > 0 && values[0] instanceof Complex) {
      this.data = values.slice();
    } else if (Array.isArray(values)) {
      this.data = values.map(Complex.from);
    } else if (values && typeof values.length === "number") {
      const out = new Array(values.length);
      for (let i = 0; i < values.length; i++) out[i] = Complex.from(values[i]);
      this.data = out;
    } else {
      throw new TypeError("ComplexVector requires an array-like input");
    }
    this.size = this.data.length;
  }
  static zeros(n) {
    const data = new Array(n);
    for (let i = 0; i < n; i++) data[i] = Complex.ZERO;
    return new _ComplexVector(data);
  }
  static basis(n, index) {
    const v = _ComplexVector.zeros(n);
    v.data[index] = Complex.ONE;
    return v;
  }
  static fromArray(arr) {
    return new _ComplexVector(arr);
  }
  get(i) {
    return this.data[i];
  }
  set(i, v) {
    this.data[i] = Complex.from(v);
  }
  add(other) {
    this._checkSameSize(other);
    const out = new Array(this.size);
    for (let i = 0; i < this.size; i++) out[i] = this.data[i].add(other.data[i]);
    return new _ComplexVector(out);
  }
  sub(other) {
    this._checkSameSize(other);
    const out = new Array(this.size);
    for (let i = 0; i < this.size; i++) out[i] = this.data[i].sub(other.data[i]);
    return new _ComplexVector(out);
  }
  scale(c) {
    const f = Complex.from(c);
    const out = new Array(this.size);
    for (let i = 0; i < this.size; i++) out[i] = this.data[i].mul(f);
    return new _ComplexVector(out);
  }
  inner(other) {
    this._checkSameSize(other);
    let acc = Complex.ZERO;
    for (let i = 0; i < this.size; i++) {
      acc = acc.add(this.data[i].conjugate().mul(other.data[i]));
    }
    return acc;
  }
  dot(other) {
    return this.inner(other);
  }
  norm() {
    return Math.sqrt(this.inner(this).re);
  }
  normalize() {
    const n = this.norm();
    if (n === 0) throw new Error("Cannot normalize zero vector");
    return this.scale(1 / n);
  }
  tensor(other) {
    const out = new Array(this.size * other.size);
    for (let i = 0; i < this.size; i++) {
      for (let j = 0; j < other.size; j++) {
        out[i * other.size + j] = this.data[i].mul(other.data[j]);
      }
    }
    return new _ComplexVector(out);
  }
  conjugate() {
    return new _ComplexVector(this.data.map((c) => c.conjugate()));
  }
  probabilities() {
    const out = new Array(this.size);
    let total = 0;
    for (let i = 0; i < this.size; i++) {
      const p = this.data[i].abs2();
      out[i] = p;
      total += p;
    }
    if (total > 0) for (let i = 0; i < this.size; i++) out[i] /= total;
    return out;
  }
  equals(other, tol) {
    if (!(other instanceof _ComplexVector) || other.size !== this.size) return false;
    for (let i = 0; i < this.size; i++) {
      if (!this.data[i].equals(other.data[i], tol)) return false;
    }
    return true;
  }
  toArray() {
    return this.data.slice();
  }
  toRealArray() {
    return this.data.map((c) => c.re);
  }
  _checkSameSize(other) {
    if (!(other instanceof _ComplexVector) || other.size !== this.size) {
      throw new Error(`Size mismatch: ${this.size} vs ${other && other.size}`);
    }
  }
};
var ComplexMatrix2 = class _ComplexMatrix {
  constructor(rows, cols, data) {
    this.rows = rows;
    this.cols = cols;
    if (data) {
      if (data.length !== rows * cols) {
        throw new Error("ComplexMatrix data length mismatch");
      }
      this.data = data.map(Complex.from);
    } else {
      this.data = new Array(rows * cols);
      for (let i = 0; i < rows * cols; i++) this.data[i] = Complex.ZERO;
    }
  }
  static fromRows(rows) {
    const r4 = rows.length;
    const c = rows[0].length;
    const data = new Array(r4 * c);
    for (let i = 0; i < r4; i++) {
      for (let j = 0; j < c; j++) data[i * c + j] = Complex.from(rows[i][j]);
    }
    return new _ComplexMatrix(r4, c, data);
  }
  static identity(n) {
    const m = new _ComplexMatrix(n, n);
    for (let i = 0; i < n; i++) m.set(i, i, Complex.ONE);
    return m;
  }
  static zeros(n, m) {
    return new _ComplexMatrix(n, m || n);
  }
  static fromDiagonal(diag) {
    const n = diag.length;
    const m = new _ComplexMatrix(n, n);
    for (let i = 0; i < n; i++) m.set(i, i, diag[i]);
    return m;
  }
  get(i, j) {
    return this.data[i * this.cols + j];
  }
  set(i, j, v) {
    this.data[i * this.cols + j] = Complex.from(v);
  }
  getRow(i) {
    const out = new Array(this.cols);
    for (let j = 0; j < this.cols; j++) out[j] = this.get(i, j);
    return out;
  }
  getCol(j) {
    const out = new Array(this.rows);
    for (let i = 0; i < this.rows; i++) out[i] = this.get(i, j);
    return out;
  }
  add(other) {
    this._checkSameShape(other);
    const out = new Array(this.rows * this.cols);
    for (let k = 0; k < out.length; k++) out[k] = this.data[k].add(other.data[k]);
    return new _ComplexMatrix(this.rows, this.cols, out);
  }
  // In-place addition (avoids allocating a new matrix for hot paths like
  // Kraus composition and superoperator assembly).
  add_inplace(other) {
    this._checkSameShape(other);
    for (let k = 0; k < this.data.length; k++) {
      this.data[k] = this.data[k].add(other.data[k]);
    }
    return this;
  }
  _checkSameShape(other) {
    if (!(other instanceof _ComplexMatrix) || other.rows !== this.rows || other.cols !== this.cols) {
      throw new Error(`Matrix shape mismatch: ${this.rows}x${this.cols} vs ${other && other.rows}x${other && other.cols}`);
    }
  }
  sub(other) {
    this._checkSameShape(other);
    const out = new Array(this.rows * this.cols);
    for (let k = 0; k < out.length; k++) out[k] = this.data[k].sub(other.data[k]);
    return new _ComplexMatrix(this.rows, this.cols, out);
  }
  scale(c) {
    const f = Complex.from(c);
    const out = new Array(this.rows * this.cols);
    for (let k = 0; k < out.length; k++) out[k] = this.data[k].mul(f);
    return new _ComplexMatrix(this.rows, this.cols, out);
  }
  mul(other) {
    if (other instanceof _ComplexMatrix) {
      if (this.cols !== other.rows) {
        throw new Error(`Matrix multiply shape mismatch: ${this.rows}x${this.cols} * ${other.rows}x${other.cols}`);
      }
      const out = new _ComplexMatrix(this.rows, other.cols);
      for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < other.cols; j++) {
          let acc = Complex.ZERO;
          for (let k = 0; k < this.cols; k++) {
            acc = acc.add(this.get(i, k).mul(other.get(k, j)));
          }
          out.set(i, j, acc);
        }
      }
      return out;
    }
    return this.scale(other);
  }
  matvec(v) {
    if (v.size !== this.cols) throw new Error("Matrix-vector size mismatch");
    const out = new Array(this.rows);
    for (let i = 0; i < this.rows; i++) {
      let acc = Complex.ZERO;
      for (let j = 0; j < this.cols; j++) {
        acc = acc.add(this.get(i, j).mul(v.data[j]));
      }
      out[i] = acc;
    }
    return new ComplexVector(out);
  }
  tensor(other) {
    const out = new _ComplexMatrix(this.rows * other.rows, this.cols * other.cols);
    for (let i1 = 0; i1 < this.rows; i1++) {
      for (let j1 = 0; j1 < this.cols; j1++) {
        const a = this.get(i1, j1);
        for (let i2 = 0; i2 < other.rows; i2++) {
          for (let j2 = 0; j2 < other.cols; j2++) {
            out.set(i1 * other.rows + i2, j1 * other.cols + j2, a.mul(other.get(i2, j2)));
          }
        }
      }
    }
    return out;
  }
  conjugate() {
    const out = new Array(this.data.length);
    for (let k = 0; k < out.length; k++) out[k] = this.data[k].conjugate();
    return new _ComplexMatrix(this.rows, this.cols, out);
  }
  transpose() {
    const out = new _ComplexMatrix(this.cols, this.rows);
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        out.set(j, i, this.get(i, j));
      }
    }
    return out;
  }
  dagger() {
    const out = new _ComplexMatrix(this.cols, this.rows);
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        out.set(j, i, this.get(i, j).conjugate());
      }
    }
    return out;
  }
  trace() {
    if (this.rows !== this.cols) throw new Error("Trace requires square matrix");
    let acc = Complex.ZERO;
    for (let i = 0; i < this.rows; i++) acc = acc.add(this.get(i, i));
    return acc;
  }
  det() {
    if (this.rows !== this.cols) throw new Error("det requires square matrix");
    const n = this.rows;
    const A = this._clone();
    let det = Complex.ONE;
    for (let i = 0; i < n; i++) {
      let maxRow = i;
      let maxVal = A.get(i, i).abs();
      for (let k = i + 1; k < n; k++) {
        const v = A.get(k, i).abs();
        if (v > maxVal) {
          maxVal = v;
          maxRow = k;
        }
      }
      if (maxVal < 1e-15) return Complex.ZERO;
      if (maxRow !== i) {
        for (let j = 0; j < n; j++) {
          const tmp = A.get(i, j);
          A.set(i, j, A.get(maxRow, j));
          A.set(maxRow, j, tmp);
        }
        det = det.negate();
      }
      det = det.mul(A.get(i, i));
      const pivot = A.get(i, i);
      for (let k = i + 1; k < n; k++) {
        const factor = A.get(k, i).div(pivot);
        for (let j = i; j < n; j++) {
          A.set(k, j, A.get(k, j).sub(factor.mul(A.get(i, j))));
        }
      }
    }
    return det;
  }
  isUnitary(tol) {
    const eps = typeof tol === "number" ? tol : 1e-9;
    if (this.rows !== this.cols) return false;
    const prod = this.dagger().mul(this);
    const id = _ComplexMatrix.identity(this.rows);
    return prod.equals(id, eps);
  }
  isHermitian(tol) {
    const eps = typeof tol === "number" ? tol : 1e-9;
    return this.equals(this.dagger(), eps);
  }
  equals(other, tol) {
    const eps = typeof tol === "number" ? tol : 1e-12;
    if (!(other instanceof _ComplexMatrix)) return false;
    if (this.rows !== other.rows || this.cols !== other.cols) return false;
    for (let k = 0; k < this.data.length; k++) {
      if (!this.data[k].equals(other.data[k], eps)) return false;
    }
    return true;
  }
  // Partial trace over a subsystem. qubits: number of qubits the matrix
  // represents (rows must equal 2^qubits). subsystem: index of qubit to
  // trace out (0 = least significant).
  partialTrace(qubits, subsystem) {
    if (this.rows !== this.cols || this.rows !== 1 << qubits) {
      throw new Error("partialTrace: matrix must be 2^n x 2^n for n qubits");
    }
    const kept = qubits - 1;
    const dim = 1 << kept;
    const out = _ComplexMatrix.zeros(dim, dim);
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        let acc = Complex.ZERO;
        for (let s = 0; s < 2; s++) {
          const bi = _insertBit(i, subsystem, s);
          const bj = _insertBit(j, subsystem, s);
          acc = acc.add(this.get(bi, bj));
        }
        out.set(i, j, acc);
      }
    }
    return out;
  }
  // Eigenvalue decomposition for Hermitian matrices via Jacobi rotations.
  // Returns { eigenvalues: number[], eigenvectors: ComplexMatrix } where
  // eigenvectors columns are the eigenvectors.
  eigh() {
    if (!this.isHermitian(1e-9)) {
      throw new Error("eigh only supports Hermitian matrices");
    }
    const n = this.rows;
    const A = this._clone();
    const V = _ComplexMatrix.identity(n);
    const maxSweeps = 100;
    for (let sweep = 0; sweep < maxSweeps; sweep++) {
      let off = 0;
      for (let p = 0; p < n; p++) {
        for (let q = p + 1; q < n; q++) {
          const apq = A.get(p, q);
          off += apq.abs2();
        }
      }
      if (off < 1e-20) break;
      for (let p = 0; p < n; p++) {
        for (let q = p + 1; q < n; q++) {
          const apq = A.get(p, q);
          if (apq.abs() < 1e-15) continue;
          const app = A.get(p, p).re;
          const aqq = A.get(q, q).re;
          const phi = Math.atan2(apq.im, apq.re);
          const offdiag = Math.hypot(apq.re, apq.im);
          const theta = Math.atan2(2 * offdiag, aqq - app) / 2;
          const c = Math.cos(theta);
          const s = Math.sin(theta);
          const phase = new Complex(Math.cos(phi / 2), Math.sin(phi / 2));
          for (let i = 0; i < n; i++) {
            const aip = A.get(i, p);
            const aiq = A.get(i, q);
            A.set(i, p, aip.scale(c).add(aiq.mul(phase.conjugate()).scale(s)));
            A.set(i, q, aiq.scale(c).sub(aip.mul(phase).scale(s)));
          }
          for (let j = 0; j < n; j++) {
            const apj = A.get(p, j);
            const aqj = A.get(q, j);
            A.set(p, j, apj.scale(c).add(aqj.mul(phase).scale(s)));
            A.set(q, j, aqj.scale(c).sub(apj.mul(phase.conjugate()).scale(s)));
          }
          for (let i = 0; i < n; i++) {
            const vip = V.get(i, p);
            const viq = V.get(i, q);
            V.set(i, p, vip.scale(c).add(viq.mul(phase).scale(s)));
            V.set(i, q, viq.scale(c).sub(vip.mul(phase.conjugate()).scale(s)));
          }
        }
      }
    }
    const eigenvalues = new Array(n);
    for (let i = 0; i < n; i++) eigenvalues[i] = A.get(i, i).re;
    return { eigenvalues, eigenvectors: V };
  }
  // Singular value decomposition for any matrix. Returns { U, S, Vh } where
  // A = U * diag(S) * Vh. Uses eigendecomposition of A^dagger A.
  svd() {
    const m = this.rows;
    const n = this.cols;
    if (m >= n) {
      const AtA = this.dagger().mul(this);
      const { eigenvalues, eigenvectors } = AtA.eigh();
      const indices = eigenvalues.map((v, i) => i).sort((a, b) => eigenvalues[b] - eigenvalues[a]);
      const S = indices.map((i) => Math.sqrt(Math.max(0, eigenvalues[i])));
      const V = eigenvectors;
      const U = _ComplexMatrix.zeros(m, Math.min(m, n));
      const k = Math.min(m, n);
      for (let j = 0; j < k; j++) {
        const s = S[j] > 1e-12 ? S[j] : 1;
        for (let i = 0; i < m; i++) {
          let acc = Complex.ZERO;
          for (let l = 0; l < n; l++) {
            const vIdx = indices[j];
            acc = acc.add(this.get(i, l).mul(V.get(l, vIdx)));
          }
          U.set(i, j, acc.scale(1 / s));
        }
      }
      const Vh = _ComplexMatrix.zeros(k, n);
      for (let i = 0; i < k; i++) {
        for (let j = 0; j < n; j++) {
          Vh.set(i, j, V.get(j, indices[i]).conjugate());
        }
      }
      return { U, S, Vh };
    } else {
      const t = this.transpose();
      const r4 = t.svd();
      return { U: r4.Vh.transpose(), S: r4.S, Vh: r4.U.transpose() };
    }
  }
  // Matrix inverse via Gauss-Jordan elimination
  inverse() {
    if (this.rows !== this.cols) throw new Error("inverse requires square matrix");
    const n = this.rows;
    const A = this._clone();
    const I = _ComplexMatrix.identity(n);
    for (let i = 0; i < n; i++) {
      let maxRow = i;
      let maxVal = A.get(i, i).abs();
      for (let k = i + 1; k < n; k++) {
        const v = A.get(k, i).abs();
        if (v > maxVal) {
          maxVal = v;
          maxRow = k;
        }
      }
      if (maxVal < 1e-15) throw new Error("Matrix is singular");
      if (maxRow !== i) {
        for (let j = 0; j < n; j++) {
          const tmp = A.get(i, j);
          A.set(i, j, A.get(maxRow, j));
          A.set(maxRow, j, tmp);
          const tmp2 = I.get(i, j);
          I.set(i, j, I.get(maxRow, j));
          I.set(maxRow, j, tmp2);
        }
      }
      const pivot = A.get(i, i);
      for (let j = 0; j < n; j++) {
        A.set(i, j, A.get(i, j).div(pivot));
        I.set(i, j, I.get(i, j).div(pivot));
      }
      for (let k = 0; k < n; k++) {
        if (k === i) continue;
        const factor = A.get(k, i);
        for (let j = 0; j < n; j++) {
          A.set(k, j, A.get(k, j).sub(factor.mul(A.get(i, j))));
          I.set(k, j, I.get(k, j).sub(factor.mul(I.get(i, j))));
        }
      }
    }
    return I;
  }
  // Matrix exponentiation via Taylor series with scaling and squaring
  expm() {
    if (this.rows !== this.cols) throw new Error("expm requires square matrix");
    const n = this.rows;
    let norm = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) norm = Math.max(norm, this.get(i, j).abs());
    let s = 0;
    while (norm > 0.5) {
      norm /= 2;
      s++;
    }
    const scaled = this.scale(1 / Math.pow(2, s));
    let result = _ComplexMatrix.identity(n);
    let term = _ComplexMatrix.identity(n);
    for (let k = 1; k <= 20; k++) {
      term = term.mul(scaled).scale(1 / k);
      result = result.add(term);
      if (term.data.every((v) => v.abs() < 1e-18)) break;
    }
    for (let i = 0; i < s; i++) result = result.mul(result);
    return result;
  }
  // Kronecker product alias
  kron(other) {
    return this.tensor(other);
  }
  _clone() {
    return new _ComplexMatrix(this.rows, this.cols, this.data.slice());
  }
  toRows() {
    const out = new Array(this.rows);
    for (let i = 0; i < this.rows; i++) out[i] = this.getRow(i);
    return out;
  }
  toString() {
    const rows = this.toRows().map((r4) => r4.map((c) => c.toString()).join("  "));
    return rows.join("\n");
  }
};
function _insertBit(bits, pos, value) {
  const mask = (1 << pos) - 1;
  const low = bits & mask;
  const high = (bits & ~mask) << 1;
  return high | low | value << pos;
}
function kronMatrices(matrices) {
  if (matrices.length === 0) return ComplexMatrix2.identity(1);
  let acc = matrices[0];
  for (let i = 1; i < matrices.length; i++) acc = acc.tensor(matrices[i]);
  return acc;
}
function kronVectors(vectors) {
  if (vectors.length === 0) return new ComplexVector([Complex.ONE]);
  let acc = vectors[0];
  for (let i = 1; i < vectors.length; i++) acc = acc.tensor(vectors[i]);
  return acc;
}
function randomUniform() {
  const c = typeof globalThis !== "undefined" && globalThis.crypto && globalThis.crypto.getRandomValues;
  if (c) {
    const buf = new Uint32Array(1);
    globalThis.crypto.getRandomValues(buf);
    return buf[0] / 4294967296;
  }
  return Math.random();
}
function sampleDistribution(probs, rng) {
  const r4 = typeof rng === "function" ? rng() : randomUniform();
  let cum = 0;
  for (let i = 0; i < probs.length; i++) {
    cum += probs[i];
    if (r4 < cum) return i;
  }
  return probs.length - 1;
}
var PAULI = {
  I: ComplexMatrix2.identity(2),
  X: ComplexMatrix2.fromRows([[new Complex(0), new Complex(1)], [new Complex(1), new Complex(0)]]),
  Y: ComplexMatrix2.fromRows([[new Complex(0), new Complex(0, -1)], [new Complex(0, 1), new Complex(0)]]),
  Z: ComplexMatrix2.fromRows([[new Complex(1), new Complex(0)], [new Complex(0), new Complex(-1)]])
};
var CONSTANTS = {
  ZERO: Complex.ZERO,
  ONE: Complex.ONE,
  I: Complex.I,
  H: ComplexMatrix2.fromRows([
    [new Complex(1 / Math.SQRT2), new Complex(1 / Math.SQRT2)],
    [new Complex(1 / Math.SQRT2), new Complex(-1 / Math.SQRT2)]
  ]),
  S: ComplexMatrix2.fromRows([
    [new Complex(1), new Complex(0)],
    [new Complex(0), new Complex(0, 1)]
  ]),
  T: ComplexMatrix2.fromRows([
    [new Complex(1), new Complex(0)],
    [new Complex(0), new Complex(Math.cos(Math.PI / 4), Math.sin(Math.PI / 4))]
  ])
};

// src/core/bit.js
var _bitCounter = 0;
var Bit = class _Bit {
  constructor(register, index) {
    this.register = register;
    this.index = index;
    this._uuid = ++_bitCounter;
  }
  equals(other) {
    return other instanceof _Bit && this.register === other.register && this.index === other.index;
  }
  toString() {
    const name = this.register && this.register.name || "?";
    return `${name}[${this.index}]`;
  }
};
var Qubit = class extends Bit {
};
var Clbit = class extends Bit {
};
var Register = class {
  constructor(size, name, bitClass) {
    if (typeof size !== "number" || size <= 0) {
      throw new TypeError("Register size must be a positive integer");
    }
    this._size = size;
    this.name = name;
    this._bits = new Array(size);
    for (let i = 0; i < size; i++) {
      this._bits[i] = new bitClass(this, i);
    }
  }
  get size() {
    return this._size;
  }
  set size(v) {
    this._size = v;
  }
  get nbits() {
    return this._size;
  }
  [Symbol.iterator]() {
    return this._bits[Symbol.iterator]();
  }
  get(i) {
    if (i < 0 || i >= this._size) throw new RangeError(`Register index ${i} out of range`);
    return this._bits[i];
  }
  slice(start, end) {
    return this._bits.slice(start, end);
  }
  len() {
    return this._size;
  }
};
var _qrCounter = 0;
var _crCounter = 0;
var QuantumRegister = class extends Register {
  constructor(size, name) {
    super(size, name || null, Qubit);
    if (!name) this.name = `q${_qrCounter++}`;
  }
  get nqubits() {
    return this._size;
  }
};
var ClassicalRegister = class extends Register {
  constructor(size, name) {
    super(size, name || null, Clbit);
    if (!name) this.name = `c${_crCounter++}`;
  }
  get nclbits() {
    return this._size;
  }
};

// src/core/parameter.js
var _paramCounter = 0;
var Parameter = class _Parameter {
  constructor(name) {
    this.name = name || `param_${_paramCounter++}`;
    this._uuid = ++_paramCounter;
  }
  add(other) {
    return new ParameterExpression("+", this, _coerce(other));
  }
  sub(other) {
    return new ParameterExpression("-", this, _coerce(other));
  }
  mul(other) {
    return new ParameterExpression("*", this, _coerce(other));
  }
  div(other) {
    return new ParameterExpression("/", this, _coerce(other));
  }
  neg() {
    return new ParameterExpression("*", this, _coerce(-1));
  }
  pow(n) {
    return new ParameterExpression("**", this, _coerce(n));
  }
  bind(values) {
    if (typeof values === "number") return values;
    if (!(this.name in values)) {
      throw new Error(`Missing value for parameter ${this.name}`);
    }
    return Number(values[this.name]);
  }
  get parameters() {
    return /* @__PURE__ */ new Set([this]);
  }
  toString() {
    return this.name;
  }
  equals(other) {
    return other instanceof _Parameter && other._uuid === this._uuid;
  }
  // Reject implicit numeric coercion so that bugs like `2 * param.mul(c)`
  // fail loudly instead of silently producing NaN.
  valueOf() {
    throw new TypeError(
      `Parameter "${this.name}" cannot be implicitly coerced to a number. Call .bind(values) first, or chain ParameterExpression methods: use param.mul(2) instead of 2 * param.`
    );
  }
};
var ParameterExpression = class _ParameterExpression {
  constructor(op, left, right) {
    this.op = op;
    this.left = left;
    this.right = right;
  }
  add(other) {
    return new _ParameterExpression("+", this, _coerce(other));
  }
  sub(other) {
    return new _ParameterExpression("-", this, _coerce(other));
  }
  mul(other) {
    return new _ParameterExpression("*", this, _coerce(other));
  }
  div(other) {
    return new _ParameterExpression("/", this, _coerce(other));
  }
  pow(n) {
    return new _ParameterExpression("**", this, _coerce(n));
  }
  bind(values) {
    const l = _bindNode(this.left, values);
    const r4 = _bindNode(this.right, values);
    switch (this.op) {
      case "+":
        return l + r4;
      case "-":
        return l - r4;
      case "*":
        return l * r4;
      case "/":
        return l / r4;
      case "**":
        return Math.pow(l, r4);
      default:
        throw new Error(`Unknown op ${this.op}`);
    }
  }
  get parameters() {
    const set = /* @__PURE__ */ new Set();
    _collectParams(this, set);
    return set;
  }
  toString() {
    return `(${_nodeToString(this.left)} ${this.op} ${_nodeToString(this.right)})`;
  }
  // Same rationale as Parameter.valueOf: catch implicit coercion bugs.
  valueOf() {
    throw new TypeError(
      `ParameterExpression "${this.toString()}" cannot be implicitly coerced to a number. Call .bind(values) first.`
    );
  }
};
function _coerce(v) {
  if (v instanceof Parameter || v instanceof ParameterExpression) return v;
  if (typeof v === "number") return v;
  throw new TypeError("Cannot coerce value to ParameterExpression");
}
function _bindNode(node, values) {
  if (typeof node === "number") return node;
  if (node instanceof Parameter) return node.bind(values);
  if (node instanceof ParameterExpression) return node.bind(values);
  throw new TypeError("Unbindable parameter node");
}
function _collectParams(node, set) {
  if (node instanceof Parameter) set.add(node);
  else if (node instanceof ParameterExpression) {
    _collectParams(node.left, set);
    _collectParams(node.right, set);
  }
}
function _nodeToString(node) {
  if (typeof node === "number") return String(node);
  return node.toString();
}
var ParameterVector = class {
  constructor(name, length) {
    this.name = name || `pv_${_paramCounter++}`;
    this._params = new Array(length || 0);
    for (let i = 0; i < this._params.length; i++) {
      this._params[i] = new Parameter(`${this.name}[${i}]`);
    }
  }
  get length() {
    return this._params.length;
  }
  get(i) {
    return this._params[i];
  }
  [Symbol.iterator]() {
    return this._params[Symbol.iterator]();
  }
  slice(start, end) {
    return this._params.slice(start, end);
  }
  push(name) {
    const p = new Parameter(name || `${this.name}[${this._params.length}]`);
    this._params.push(p);
    return p;
  }
};

// src/core/gate.js
var Instruction = class _Instruction {
  constructor(name, num_qubits, num_clbits, params) {
    this.name = name;
    this.num_qubits = num_qubits;
    this.num_clbits = num_clbits;
    this.params = params || [];
    this.label = null;
    this.condition = null;
    this._definition = null;
    this._matrixBuilder = null;
  }
  get num_params() {
    return this.params.length;
  }
  get definition() {
    if (this._definition) return this._definition(this);
    return null;
  }
  set definition(fn) {
    this._definition = fn;
  }
  to_matrix() {
    if (this._matrixBuilder) {
      return this._matrixBuilder(this.params);
    }
    throw new Error(`Instruction ${this.name} does not define a matrix`);
  }
  has_matrix() {
    return this._matrixBuilder !== null;
  }
  copy() {
    const c = new _Instruction(this.name, this.num_qubits, this.num_clbits, this.params.slice());
    c.label = this.label;
    c.condition = this.condition;
    c._definition = this._definition;
    c._matrixBuilder = this._matrixBuilder;
    return c;
  }
  control(num_ctrl_qubits = 1, label = null, ctrl_state = null) {
    return new ControlledGate(this, num_ctrl_qubits, label, ctrl_state);
  }
  inverse() {
    const selfInverse = [
      "h",
      "x",
      "y",
      "z",
      "cx",
      "cy",
      "cz",
      "swap",
      "ccx",
      "cswap",
      "id",
      "ch"
    ];
    const inverseMap = {
      "s": "sdg",
      "sdg": "s",
      "t": "tdg",
      "tdg": "t",
      "sx": "sxdg",
      "sxdg": "sx"
    };
    const negateParam = [
      "rx",
      "ry",
      "rz",
      "p",
      "u1",
      "rxx",
      "ryy",
      "rzz",
      "rzx",
      "cp",
      "crx",
      "cry",
      "crz",
      "cu1"
    ];
    const daggerGates = ["iswap", "csx", "dcx", "u2", "u3", "u", "cu3", "cu"];
    const name = this.name.toLowerCase();
    if (selfInverse.includes(name)) {
      return this.copy();
    }
    if (inverseMap[name]) {
      const inv2 = this.copy();
      inv2.name = inverseMap[name];
      if (this._matrixBuilder) {
        const originalBuilder = this._matrixBuilder;
        inv2._matrixBuilder = (params) => originalBuilder(params).dagger();
      }
      return inv2;
    }
    if (negateParam.includes(name)) {
      const inv2 = this.copy();
      inv2.params = this.params.map((p) => {
        if (typeof p === "number") return -p;
        if (p && typeof p.negate === "function") return p.negate();
        if (p && typeof p.mul === "function") return p.mul(-1);
        return p;
      });
      return inv2;
    }
    if (daggerGates.includes(name)) {
      const inv2 = this.copy();
      if (this._matrixBuilder) {
        const originalBuilder = this._matrixBuilder;
        inv2._matrixBuilder = (params) => originalBuilder(params).dagger();
      }
      return inv2;
    }
    if (this instanceof ControlledGate) {
      const invBase = this.base_gate.inverse();
      const inv2 = new ControlledGate(invBase, this.num_ctrl_qubits, this.label, this.ctrl_state);
      inv2.name = this.name;
      return inv2;
    }
    const inv = this.copy();
    inv.name = this.name + "_dg";
    if (this._matrixBuilder) {
      const originalBuilder = this._matrixBuilder;
      inv._matrixBuilder = (params) => originalBuilder(params).dagger();
    }
    if (this._definition) {
      const originalDef = this._definition;
      inv._definition = (instr) => {
        const def = originalDef(instr);
        return def.map(([i, q, c]) => [i.inverse(), q.slice().reverse(), c.slice().reverse()]);
      };
    }
    return inv;
  }
  toString() {
    const p = this.params.length ? `(${this.params.map(_p2s).join(",")})` : "";
    return `${this.name}${p}`;
  }
};
function _p2s(p) {
  if (p && typeof p.toString === "function") return p.toString();
  return String(p);
}
var Gate = class extends Instruction {
  constructor(name, num_qubits, params) {
    super(name, num_qubits, 0, params);
  }
};
var ControlledGate = class _ControlledGate extends Gate {
  constructor(baseGate, num_ctrl_qubits, label, ctrl_state) {
    super(
      `c${num_ctrl_qubits}_${baseGate.name}`,
      baseGate.num_qubits + num_ctrl_qubits,
      baseGate.params.slice()
    );
    this.base_gate = baseGate;
    this.num_ctrl_qubits = num_ctrl_qubits;
    this.label = label;
    this.ctrl_state = ctrl_state !== null && ctrl_state !== void 0 ? ctrl_state : (1 << num_ctrl_qubits) - 1;
    this._matrixBuilder = this._controlledMatrixBuilder.bind(this);
  }
  _controlledMatrixBuilder(params) {
    const baseMatrix = this.base_gate.to_matrix();
    const dim = baseMatrix.rows;
    const numCtrl = 1 << this.num_ctrl_qubits;
    const total = dim * numCtrl;
    const out = ComplexMatrix2.zeros(total, total);
    const ctrlOn = this.ctrl_state;
    for (let c = 0; c < numCtrl; c++) {
      if (c === ctrlOn) {
        for (let bi = 0; bi < dim; bi++) {
          for (let bj = 0; bj < dim; bj++) {
            const row = c + (bi << this.num_ctrl_qubits);
            const col = c + (bj << this.num_ctrl_qubits);
            out.set(row, col, baseMatrix.get(bi, bj));
          }
        }
      } else {
        for (let bi = 0; bi < dim; bi++) {
          const idx = c + (bi << this.num_ctrl_qubits);
          out.set(idx, idx, Complex.ONE);
        }
      }
    }
    return out;
  }
  copy() {
    const inv = new _ControlledGate(this.base_gate.copy(), this.num_ctrl_qubits, this.label, this.ctrl_state);
    inv.name = this.name;
    inv.params = this.params.slice();
    return inv;
  }
  inverse() {
    const inv = new _ControlledGate(this.base_gate.inverse(), this.num_ctrl_qubits, this.label, this.ctrl_state);
    inv.name = this.name + "_dg";
    return inv;
  }
};

// src/core/circuit.js
var _extraGateClasses = {};
function _registerExtraGateClass(name, cls) {
  _extraGateClasses[name] = cls;
}
function _extraGateClass(name) {
  if (!_extraGateClasses[name]) {
    throw new Error(
      `Extra gate class ${name} not registered. Import '@nxrix/ketra/library/extra_gates.js' before using circuit.unitary() / .initialize() / .diagonal() / .permute() / .hamiltonian().`
    );
  }
  return _extraGateClasses[name];
}
var _OperatorClass = null;
var _StatevectorClass = null;
var _QASMExporterClass = null;
function _registerOperatorClass(cls) {
  _OperatorClass = cls;
}
function _registerStatevectorClass(cls) {
  _StatevectorClass = cls;
}
function _registerQASMExporterClass(cls) {
  _QASMExporterClass = cls;
}
var CircuitInstruction = class _CircuitInstruction {
  constructor(operation, qubits, clbits) {
    this.operation = operation;
    this.qubits = qubits || [];
    this.clbits = clbits || [];
  }
  copy() {
    return new _CircuitInstruction(this.operation.copy(), this.qubits.slice(), this.clbits.slice());
  }
};
var _stdCache = {};
var _paramBuilders = {};
function _registerStd(name, gate) {
  _stdCache[name] = gate;
}
function _registerParamBuilder(name, builder) {
  _paramBuilders[name] = builder;
}
function _getStdGate(name) {
  if (!_stdCache[name]) {
    throw new Error(`Standard gate ${name} not registered. Did you import library/standard_gates?`);
  }
  return _stdCache[name].copy();
}
var QuantumCircuit = class _QuantumCircuit {
  constructor(...regs) {
    this.qregs = [];
    this.cregs = [];
    this.qubits = [];
    this.clbits = [];
    this.data = [];
    this._qubit_index = /* @__PURE__ */ new Map();
    this._clbit_index = /* @__PURE__ */ new Map();
    this._global_phase = 0;
    this.name = "circuit";
    this.metadata = null;
    if (regs.length > 0 && regs.every((r4) => typeof r4 === "number")) {
      const nQ = regs[0];
      const nC = regs.length > 1 ? regs[1] : 0;
      if (nQ > 0) this.add_register(new QuantumRegister(nQ, "q"));
      if (nC > 0) this.add_register(new ClassicalRegister(nC, "c"));
      return;
    }
    for (const r4 of regs) {
      if (r4 instanceof QuantumRegister) this.add_register(r4);
      else if (r4 instanceof ClassicalRegister) this.add_register(r4);
      else if (typeof r4 === "number") this.add_register(new QuantumRegister(r4));
      else if (Array.isArray(r4)) {
        for (const sub of r4) {
          if (sub instanceof QuantumRegister) this.add_register(sub);
          else if (sub instanceof ClassicalRegister) this.add_register(sub);
          else throw new TypeError("Invalid register in list");
        }
      } else throw new TypeError("Invalid register argument to QuantumCircuit");
    }
  }
  add_register(register) {
    if (register instanceof QuantumRegister) {
      this.qregs.push(register);
      for (const b of register._bits) {
        this._qubit_index.set(b, this.qubits.length);
        this.qubits.push(b);
      }
    } else if (register instanceof ClassicalRegister) {
      this.cregs.push(register);
      for (const b of register._bits) {
        this._clbit_index.set(b, this.clbits.length);
        this.clbits.push(b);
      }
    } else {
      throw new TypeError("add_register requires QuantumRegister or ClassicalRegister");
    }
    return register;
  }
  get num_qubits() {
    return this.qubits.length;
  }
  get num_clbits() {
    return this.clbits.length;
  }
  get nqubits() {
    return this.qubits.length;
  }
  get nclbits() {
    return this.clbits.length;
  }
  qubit_indices(qubit) {
    if (qubit instanceof Qubit) {
      const i = this._qubit_index.get(qubit);
      if (i === void 0) throw new Error("Qubit not in circuit");
      return i;
    }
    throw new TypeError("qubit_indices expects a Qubit");
  }
  clbit_indices(clbit) {
    if (clbit instanceof Clbit) {
      const i = this._clbit_index.get(clbit);
      if (i === void 0) throw new Error("Clbit not in circuit");
      return i;
    }
    throw new TypeError("clbit_indices expects a Clbit");
  }
  _resolveQubits(qargs) {
    if (!Array.isArray(qargs)) qargs = [qargs];
    return qargs.map((q) => {
      if (q instanceof Qubit) return q;
      if (typeof q === "number") {
        if (q < 0 || q >= this.qubits.length) throw new RangeError(`Qubit index ${q} out of range`);
        return this.qubits[q];
      }
      if (q && q.register && typeof q.index === "number") return q;
      throw new TypeError(`Cannot resolve qubit: ${q}`);
    });
  }
  _resolveClbits(cargs) {
    if (!Array.isArray(cargs)) cargs = [cargs];
    return cargs.map((c) => {
      if (c instanceof Clbit) return c;
      if (typeof c === "number") {
        if (c < 0 || c >= this.clbits.length) throw new RangeError(`Clbit index ${c} out of range`);
        return this.clbits[c];
      }
      if (c && c.register && typeof c.index === "number") return c;
      throw new TypeError(`Cannot resolve clbit: ${c}`);
    });
  }
  append(instruction, qargs, cargs) {
    if (!(instruction instanceof Instruction)) {
      throw new TypeError("append expects an Instruction / Gate");
    }
    const qubits = this._resolveQubits(this._normalizeArgs(qargs));
    const clbits = this._resolveClbits(this._normalizeArgs(cargs));
    if (qubits.length !== instruction.num_qubits) {
      throw new Error(`Qubit count mismatch: gate ${instruction.name} expects ${instruction.num_qubits} but got ${qubits.length}`);
    }
    if (clbits.length !== instruction.num_clbits) {
      throw new Error(`Clbit count mismatch: gate ${instruction.name} expects ${instruction.num_clbits} but got ${clbits.length}`);
    }
    const ci = new CircuitInstruction(instruction, qubits, clbits);
    this.data.push(ci);
    return ci;
  }
  _normalizeArgs(args) {
    if (args == null) return [];
    if (Array.isArray(args)) return args;
    return [args];
  }
  compose(other, qubits, clbits, front = false, inplace = false) {
    const target = inplace ? this : this.copy();
    const qMap = other._composeQubitMap(target, qubits);
    const cMap = other._composeClbitMap(target, clbits);
    const newData = other.data.map((ci) => {
      return new CircuitInstruction(
        ci.operation.copy(),
        ci.qubits.map((q) => qMap.get(q)),
        ci.clbits.map((c) => cMap.get(c))
      );
    });
    if (front) target.data = newData.concat(target.data);
    else target.data = target.data.concat(newData);
    return target;
  }
  _composeQubitMap(target, qubits) {
    const map = /* @__PURE__ */ new Map();
    if (qubits == null) {
      for (let i = 0; i < this.qubits.length && i < target.qubits.length; i++) {
        map.set(this.qubits[i], target.qubits[i]);
      }
      return map;
    }
    if (Array.isArray(qubits)) {
      qubits.forEach((q, i) => {
        const tq = target._resolveQubits([q])[0];
        map.set(this.qubits[i], tq);
      });
    }
    return map;
  }
  _composeClbitMap(target, clbits) {
    const map = /* @__PURE__ */ new Map();
    if (clbits == null) {
      for (let i = 0; i < this.clbits.length && i < target.clbits.length; i++) {
        map.set(this.clbits[i], target.clbits[i]);
      }
      return map;
    }
    if (Array.isArray(clbits)) {
      clbits.forEach((c, i) => {
        const tc = target._resolveClbits([c])[0];
        map.set(this.clbits[i], tc);
      });
    }
    return map;
  }
  copy() {
    const c = new _QuantumCircuit();
    const qubitMap = /* @__PURE__ */ new Map();
    const clbitMap = /* @__PURE__ */ new Map();
    for (const r4 of this.qregs) {
      const newReg = new QuantumRegister(r4.size, r4.name);
      c.add_register(newReg);
      for (let i = 0; i < r4.size; i++) qubitMap.set(r4._bits[i], newReg._bits[i]);
    }
    for (const r4 of this.cregs) {
      const newReg = new ClassicalRegister(r4.size, r4.name);
      c.add_register(newReg);
      for (let i = 0; i < r4.size; i++) clbitMap.set(r4._bits[i], newReg._bits[i]);
    }
    c.data = this.data.map((ci) => {
      const newOp = ci.operation.copy();
      const newQubits = ci.qubits.map((q) => qubitMap.get(q) || q);
      const newClbits = ci.clbits.map((cl) => clbitMap.get(cl) || cl);
      return new CircuitInstruction(newOp, newQubits, newClbits);
    });
    c.global_phase = this.global_phase;
    c.name = this.name;
    c.metadata = this.metadata;
    return c;
  }
  // Standard gates
  h(q) {
    return this.append(_getStdGate("H"), q);
  }
  x(q) {
    return this.append(_getStdGate("X"), q);
  }
  y(q) {
    return this.append(_getStdGate("Y"), q);
  }
  z(q) {
    return this.append(_getStdGate("Z"), q);
  }
  s(q) {
    return this.append(_getStdGate("S"), q);
  }
  sdg(q) {
    return this.append(_getStdGate("SDG"), q);
  }
  t(q) {
    return this.append(_getStdGate("T"), q);
  }
  tdg(q) {
    return this.append(_getStdGate("TDG"), q);
  }
  sx(q) {
    return this.append(_getStdGate("SX"), q);
  }
  sxdg(q) {
    return this.append(_getStdGate("SXDG"), q);
  }
  id(q) {
    return this.append(_getStdGate("I"), q);
  }
  cx(control, target) {
    return this.append(_getStdGate("CX"), [control, target]);
  }
  cy(control, target) {
    return this.append(_getStdGate("CY"), [control, target]);
  }
  cz(control, target) {
    return this.append(_getStdGate("CZ"), [control, target]);
  }
  ch(control, target) {
    return this.append(_getStdGate("CH"), [control, target]);
  }
  csx(control, target) {
    return this.append(_getStdGate("CSX"), [control, target]);
  }
  swap(q1, q2) {
    return this.append(_getStdGate("SWAP"), [q1, q2]);
  }
  iswap(q1, q2) {
    return this.append(_getStdGate("ISWAP"), [q1, q2]);
  }
  dcx(q1, q2) {
    return this.append(_getStdGate("DCX"), [q1, q2]);
  }
  ccx(c1, c2, target) {
    return this.append(_getStdGate("CCX"), [c1, c2, target]);
  }
  cswap(c, t1, t2) {
    return this.append(_getStdGate("CSWAP"), [c, t1, t2]);
  }
  rx(theta, q) {
    return this.append(_mkParamGate("RX", 1, [theta], 1), q);
  }
  ry(theta, q) {
    return this.append(_mkParamGate("RY", 1, [theta], 1), q);
  }
  rz(phi, q) {
    return this.append(_mkParamGate("RZ", 1, [phi], 1), q);
  }
  rxx(theta, q1, q2) {
    return this.append(_mkParamGate("RXX", 2, [theta], 1), [q1, q2]);
  }
  ryy(theta, q1, q2) {
    return this.append(_mkParamGate("RYY", 2, [theta], 1), [q1, q2]);
  }
  rzz(theta, q1, q2) {
    return this.append(_mkParamGate("RZZ", 2, [theta], 1), [q1, q2]);
  }
  rzx(theta, q1, q2) {
    return this.append(_mkParamGate("RZX", 2, [theta], 1), [q1, q2]);
  }
  p(theta, q) {
    return this.append(_mkParamGate("P", 1, [theta], 1), q);
  }
  u(theta, phi, lam, q) {
    return this.append(_mkParamGate("U", 1, [theta, phi, lam], 1), q);
  }
  u1(lam, q) {
    return this.append(_mkParamGate("U1", 1, [lam], 1), q);
  }
  u2(phi, lam, q) {
    return this.append(_mkParamGate("U2", 1, [phi, lam], 1), q);
  }
  u3(theta, phi, lam, q) {
    return this.append(_mkParamGate("U3", 1, [theta, phi, lam], 1), q);
  }
  // Controlled parameterized gates
  cp(theta, control, target) {
    const g = _paramBuilders.P ? _paramBuilders.P(theta) : null;
    if (!g) throw new Error("P gate not registered");
    const cg = new ControlledGate(g, 1);
    cg.name = "cp";
    return this.append(cg, [control, target]);
  }
  crx(theta, control, target) {
    const g = _paramBuilders.RX(theta);
    const cg = new ControlledGate(g, 1);
    cg.name = "crx";
    return this.append(cg, [control, target]);
  }
  cry(theta, control, target) {
    const g = _paramBuilders.RY(theta);
    const cg = new ControlledGate(g, 1);
    cg.name = "cry";
    return this.append(cg, [control, target]);
  }
  crz(theta, control, target) {
    const g = _paramBuilders.RZ(theta);
    const cg = new ControlledGate(g, 1);
    cg.name = "crz";
    return this.append(cg, [control, target]);
  }
  cu(theta, phi, lam, gamma, control, target) {
    if (!_paramBuilders.CU) {
      throw new Error("CU gate builder not registered. Did you import library/generalized_gates?");
    }
    const g = _paramBuilders.CU(theta, phi, lam, gamma);
    return this.append(g, [control, target]);
  }
  cu1(lam, control, target) {
    const g = _paramBuilders.U1(lam);
    const cg = new ControlledGate(g, 1);
    cg.name = "cu1";
    return this.append(cg, [control, target]);
  }
  cu3(theta, phi, lam, control, target) {
    const g = _paramBuilders.U3(theta, phi, lam);
    const cg = new ControlledGate(g, 1);
    cg.name = "cu3";
    return this.append(cg, [control, target]);
  }
  // Multi-controlled gates
  mcx(controls, target, ancilla_qubits = null, mode = "noancilla") {
    const ctrls = Array.isArray(controls) ? controls : [controls];
    const base = _stdCache["X"];
    if (!base) throw new Error("X gate not registered");
    const cg = new ControlledGate(base, ctrls.length);
    cg.name = `mcx`;
    return this.append(cg, ctrls.concat([target]));
  }
  mcy(controls, target) {
    const ctrls = Array.isArray(controls) ? controls : [controls];
    const base = _stdCache["Y"];
    const cg = new ControlledGate(base, ctrls.length);
    cg.name = `mcy`;
    return this.append(cg, ctrls.concat([target]));
  }
  mcz(controls, target) {
    const ctrls = Array.isArray(controls) ? controls : [controls];
    const base = _stdCache["Z"];
    const cg = new ControlledGate(base, ctrls.length);
    cg.name = `mcz`;
    return this.append(cg, ctrls.concat([target]));
  }
  mcu1(lam, controls, target) {
    const ctrls = Array.isArray(controls) ? controls : [controls];
    const g = _paramBuilders.U1(lam);
    const cg = new ControlledGate(g, ctrls.length);
    cg.name = `mcu1`;
    return this.append(cg, ctrls.concat([target]));
  }
  // Multi-controlled phase gate: applies e^{i*lam} when all controls are |1>.
  mcp(lam, controls, target) {
    const ctrls = Array.isArray(controls) ? controls : [controls];
    const g = _paramBuilders.P(lam);
    const cg = new ControlledGate(g, ctrls.length);
    cg.name = "mcp";
    return this.append(cg, ctrls.concat([target]));
  }
  // Multi-controlled rotations.
  mcrx(theta, controls, target) {
    const ctrls = Array.isArray(controls) ? controls : [controls];
    const g = _paramBuilders.RX(theta);
    const cg = new ControlledGate(g, ctrls.length);
    cg.name = "mcrx";
    return this.append(cg, ctrls.concat([target]));
  }
  mcry(theta, controls, target) {
    const ctrls = Array.isArray(controls) ? controls : [controls];
    const g = _paramBuilders.RY(theta);
    const cg = new ControlledGate(g, ctrls.length);
    cg.name = "mcry";
    return this.append(cg, ctrls.concat([target]));
  }
  mcrz(theta, controls, target) {
    const ctrls = Array.isArray(controls) ? controls : [controls];
    const g = _paramBuilders.RZ(theta);
    const cg = new ControlledGate(g, ctrls.length);
    cg.name = "mcrz";
    return this.append(cg, ctrls.concat([target]));
  }
  // Append an arbitrary unitary matrix as a gate.
  unitary(matrix, qubits, label = null) {
    const UnitaryGate2 = _extraGateClass("UnitaryGate");
    const g = new UnitaryGate2(matrix, label);
    const qs = Array.isArray(qubits) ? qubits : [qubits];
    return this.append(g, qs);
  }
  // Append a state-preparation instruction. Amplitudes can be an array of
  // Complex/numbers or a state label like "01" or "+-".
  initialize(amplitudes, qubits) {
    const Initialize2 = _extraGateClass("Initialize");
    if (qubits == null) {
      let n2;
      if (typeof amplitudes === "string") n2 = amplitudes.length;
      else if (Array.isArray(amplitudes)) n2 = Math.log2(amplitudes.length);
      else throw new TypeError("initialize requires an array or label string");
      qubits = [];
      for (let i = 0; i < n2; i++) qubits.push(i);
    }
    const qs = Array.isArray(qubits) ? qubits : [qubits];
    let amps;
    if (typeof amplitudes === "string") {
      if (!_StatevectorClass) {
        throw new Error("Statevector class not registered. Import quantum_info/statevector.js first.");
      }
      amps = _StatevectorClass.fromLabel(amplitudes)._data.data;
    } else if (Array.isArray(amplitudes)) {
      amps = amplitudes;
    } else {
      throw new TypeError("initialize requires an array or label string");
    }
    const n = Math.log2(amps.length);
    if (!Number.isInteger(n)) {
      throw new Error(`initialize: amplitudes length ${amps.length} is not 2^n`);
    }
    const instr = new Initialize2(amps, n);
    return this.append(instr, qs);
  }
  // Append a diagonal unitary, given as a list of phases (radians).
  diagonal(diag, qubits) {
    const DiagonalGate2 = _extraGateClass("DiagonalGate");
    const g = new DiagonalGate2(diag);
    const qs = Array.isArray(qubits) ? qubits : [qubits];
    return this.append(g, qs);
  }
  // Append a qubit-permutation gate. `pattern[i]` is the qubit index that
  // ends up at position i.
  permute(pattern, qubits) {
    const PermutationGate2 = _extraGateClass("PermutationGate");
    const g = new PermutationGate2(pattern);
    const qs = Array.isArray(qubits) ? qubits : qubits !== void 0 ? [qubits] : null;
    if (qs === null) {
      const qs2 = [];
      for (let i = 0; i < pattern.length; i++) qs2.push(i);
      return this.append(g, qs2);
    }
    return this.append(g, qs);
  }
  // Append a Hamiltonian evolution gate e^{-i*t*H}.
  hamiltonian(operator, time, qubits, label = null) {
    const HamiltonianGate2 = _extraGateClass("HamiltonianGate");
    const g = new HamiltonianGate2(operator, time, label);
    const qs = Array.isArray(qubits) ? qubits : [qubits];
    return this.append(g, qs);
  }
  // Convert this circuit to a reusable Gate (with the same qubit count).
  to_gate(label = null) {
    const g = new Gate(this.name || "circuit", this.num_qubits, []);
    g.label = label;
    g._matrixBuilder = () => {
      if (!_OperatorClass) {
        throw new Error("Operator class not registered. Import quantum_info/operator.js first.");
      }
      return _OperatorClass.fromCircuit(this)._data;
    };
    g._definition = () => {
      return this.data.map((ci) => [ci.operation, ci.qubits, ci.clbits]);
    };
    return g;
  }
  // Convert this circuit to an Instruction (can have clbits).
  to_instruction(label = null) {
    const instr = new Instruction(this.name || "circuit", this.num_qubits, this.num_clbits, []);
    instr.label = label;
    instr._matrixBuilder = () => {
      if (!_OperatorClass) {
        throw new Error("Operator class not registered. Import quantum_info/operator.js first.");
      }
      return _OperatorClass.fromCircuit(this)._data;
    };
    instr._definition = () => {
      return this.data.map((ci) => [ci.operation, ci.qubits, ci.clbits]);
    };
    return instr;
  }
  // Repeat the circuit `n` times on the same qubits (in series).
  repeat(n) {
    if (n < 0) throw new Error("repeat: n must be non-negative");
    const c = new _QuantumCircuit();
    for (const r4 of this.qregs) c.add_register(new QuantumRegister(r4.size, r4.name));
    for (const r4 of this.cregs) c.add_register(new ClassicalRegister(r4.size, r4.name));
    c.global_phase = this.global_phase;
    c.name = this.name;
    const qMap = /* @__PURE__ */ new Map();
    const cMap = /* @__PURE__ */ new Map();
    for (let i = 0; i < this.qubits.length; i++) qMap.set(this.qubits[i], c.qubits[i]);
    for (let i = 0; i < this.clbits.length; i++) cMap.set(this.clbits[i], c.clbits[i]);
    for (let rep = 0; rep < n; rep++) {
      for (const ci of this.data) {
        c.append(
          ci.operation.copy(),
          ci.qubits.map((q) => qMap.get(q)),
          ci.clbits.map((cl) => cMap.get(cl))
        );
      }
    }
    return c;
  }
  // Like repeat, but composes the circuit with itself n times. For unitary
  // circuits this is equivalent to repeat(n); for parameterized circuits
  // power(n) leaves the parameters intact (matching qiskit's behavior).
  power(n) {
    if (n < 0) {
      const inv = this.inverse();
      return inv.repeat(-n);
    }
    return this.repeat(n);
  }
  // Decompose one level: replace any instruction that has a definition
  // with the instructions in its definition.
  decompose(times = 1) {
    let c = this;
    for (let i = 0; i < times; i++) {
      c = _decomposeOnce(c);
    }
    return c;
  }
  // Emit OpenQASM 2.0 source for this circuit.
  qasm() {
    if (!_QASMExporterClass) {
      throw new Error("QASMExporter not registered. Import qasm/qasm_exporter.js first.");
    }
    return new _QASMExporterClass().export(this);
  }
  // Add a global phase (radians).
  set global_phase(value) {
    this._global_phase = value;
  }
  get global_phase() {
    return this._global_phase || 0;
  }
  // Measurement and friends
  measure(qubit, clbit) {
    const qubits = this._resolveQubits(Array.isArray(qubit) ? qubit : [qubit]);
    const clbits = this._resolveClbits(Array.isArray(clbit) ? clbit : [clbit]);
    if (qubits.length !== clbits.length) {
      throw new Error("measure: number of qubits and clbits must match");
    }
    const instructions = [];
    for (let i = 0; i < qubits.length; i++) {
      const m = new Instruction("measure", 1, 1, []);
      instructions.push(this.append(m, qubits[i], clbits[i]));
    }
    return instructions.length === 1 ? instructions[0] : instructions;
  }
  measure_all(inplace = true) {
    const target = inplace ? this : this.copy();
    if (target.clbits.length === 0) {
      target.add_register(new ClassicalRegister(target.qubits.length));
    }
    for (let i = 0; i < target.qubits.length; i++) {
      target.measure(target.qubits[i], target.clbits[i]);
    }
    return target;
  }
  measure_active(inplace = true) {
    const target = inplace ? this : this.copy();
    const active = /* @__PURE__ */ new Set();
    for (const ci of target.data) {
      if (ci.operation.name === "barrier") continue;
      for (const q of ci.qubits) active.add(q);
    }
    const activeBits = target.qubits.filter((q) => active.has(q));
    if (activeBits.length === 0) return target;
    const startClbit = target.clbits.length;
    target.add_register(new ClassicalRegister(activeBits.length));
    activeBits.forEach((q, i) => target.measure(q, target.clbits[startClbit + i]));
    return target;
  }
  reset(qubit) {
    const qubits = this._resolveQubits(Array.isArray(qubit) ? qubit : [qubit]);
    const out = [];
    for (const q of qubits) {
      const r4 = new Instruction("reset", 1, 0, []);
      out.push(this.append(r4, q));
    }
    return out.length === 1 ? out[0] : out;
  }
  barrier(qubits) {
    let qs;
    if (qubits == null) qs = this.qubits.slice();
    else qs = this._resolveQubits(Array.isArray(qubits) ? qubits : [qubits]);
    const b = new Instruction("barrier", qs.length, 0, []);
    return this.append(b, qs);
  }
  delay(duration, qubit = null, unit = "dt") {
    const d = new Instruction("delay", qubit == null ? this.qubits.length : 1, 0, [duration, unit]);
    if (qubit == null) return this.append(d, this.qubits.slice());
    return this.append(d, qubit);
  }
  if_test(condition, true_predicate, false_predicate = null) {
    const instr = new Instruction("if_else", 0, 0, [condition, true_predicate, false_predicate]);
    this.data.push(new CircuitInstruction(instr, [], []));
    return instr;
  }
  while_loop(condition, body, qubits, clbits) {
    const instr = new Instruction("while_loop", 0, 0, [condition, body]);
    this.data.push(new CircuitInstruction(instr, [], []));
    return instr;
  }
  // Parameter binding
  bind_parameters(values) {
    const c = this.copy();
    c.data = c.data.map((ci) => {
      const newOp = ci.operation.copy();
      if (newOp.params && newOp.params.length) {
        newOp.params = newOp.params.map((p) => {
          if (p && typeof p.bind === "function") return p.bind(values);
          return p;
        });
      }
      return new CircuitInstruction(newOp, ci.qubits.slice(), ci.clbits.slice());
    });
    return c;
  }
  get parameters() {
    const set = /* @__PURE__ */ new Set();
    const seen = [];
    for (const ci of this.data) {
      for (const p of ci.operation.params || []) {
        if (p == null) continue;
        if (typeof p === "number") continue;
        if (p instanceof Parameter) {
          _addToSet(set, seen, p);
          continue;
        }
        let subParams = null;
        if (p && p.parameters !== void 0) {
          subParams = typeof p.parameters === "function" ? p.parameters() : p.parameters;
        }
        if (subParams) {
          for (const sp of subParams) _addToSet(set, seen, sp);
        }
      }
    }
    return set;
  }
  // Analysis methods
  depth() {
    const lastOnQubit = new Array(this.qubits.length).fill(-1);
    const lastOnClbit = new Array(this.clbits.length).fill(-1);
    let maxDepth = 0;
    for (const ci of this.data) {
      if (ci.operation.name === "barrier") continue;
      let d = 0;
      for (const q of ci.qubits) {
        const qi = this._qubit_index.get(q);
        d = Math.max(d, lastOnQubit[qi] + 1);
      }
      for (const c of ci.clbits) {
        const ci2 = this._clbit_index.get(c);
        d = Math.max(d, lastOnClbit[ci2] + 1);
      }
      for (const q of ci.qubits) {
        const qi = this._qubit_index.get(q);
        lastOnQubit[qi] = d;
      }
      for (const c of ci.clbits) {
        const ci2 = this._clbit_index.get(c);
        lastOnClbit[ci2] = d;
      }
      if (d + 1 > maxDepth) maxDepth = d + 1;
    }
    return maxDepth;
  }
  count_ops() {
    const counts = {};
    for (const ci of this.data) {
      counts[ci.operation.name] = (counts[ci.operation.name] || 0) + 1;
    }
    return counts;
  }
  num_nonlocal_gates() {
    let n = 0;
    for (const ci of this.data) {
      if (ci.operation.num_qubits > 1 && ci.operation.name !== "barrier") n++;
    }
    return n;
  }
  num_tensor_factors() {
    const parent = new Array(this.qubits.length).fill(0).map((_, i) => i);
    const find = (x) => parent[x] === x ? x : parent[x] = find(parent[x]);
    const union = (a, b) => {
      parent[find(a)] = find(b);
    };
    for (const ci of this.data) {
      if (ci.operation.num_qubits < 2) continue;
      const idx = ci.qubits.map((q) => this._qubit_index.get(q));
      for (let i = 1; i < idx.length; i++) union(idx[0], idx[i]);
    }
    const roots = /* @__PURE__ */ new Set();
    for (let i = 0; i < this.qubits.length; i++) roots.add(find(i));
    return roots.size;
  }
  size() {
    let n = 0;
    for (const ci of this.data) {
      if (ci.operation.name === "barrier") continue;
      n++;
    }
    return n;
  }
  width() {
    return this.qubits.length + this.clbits.length;
  }
  // Transforms
  reverse_ops() {
    const c = this.copy();
    c.data.reverse();
    return c;
  }
  inverse() {
    const c = this.copy();
    c.data = c.data.filter((ci) => !["measure", "reset", "barrier"].includes(ci.operation.name)).reverse().map((ci) => new CircuitInstruction(ci.operation.inverse(), ci.qubits.slice(), ci.clbits.slice()));
    return c;
  }
  remove_final_measurements(inplace = false) {
    const target = inplace ? this : this.copy();
    let lastNonMeasure = -1;
    for (let i = 0; i < target.data.length; i++) {
      if (target.data[i].operation.name !== "measure") lastNonMeasure = i;
    }
    target.data = target.data.slice(0, lastNonMeasure + 1);
    return target;
  }
  // Drawing - lazy-loaded to avoid circular imports
  draw(output = "text", kwargs) {
    if (_drawHook) {
      return _drawHook(this, output, kwargs);
    }
    return _fallbackDraw(this);
  }
  toString() {
    return _fallbackDraw(this);
  }
};
var _drawHook = null;
function _setDrawHook(fn) {
  _drawHook = fn;
}
function _mkParamGate(name, numQubits, params, numClbits) {
  if (!_paramBuilders[name]) {
    throw new Error(`Parameterized gate ${name} not registered. Did you import library/generalized_gates?`);
  }
  return _paramBuilders[name].apply(null, params);
}
function _fallbackDraw(circuit) {
  const lines = [];
  lines.push(`QuantumCircuit "${circuit.name}" (qubits=${circuit.num_qubits}, clbits=${circuit.num_clbits})`);
  for (const ci of circuit.data) {
    const qs = ci.qubits.map((q) => q.toString()).join(",");
    const cs = ci.clbits.map((c) => c.toString()).join(",");
    const cstr = cs ? ` -> ${cs}` : "";
    const pstr = ci.operation.params.length ? `(${ci.operation.params.map(String).join(",")})` : "";
    lines.push(`  ${ci.operation.name}${pstr} ${qs}${cstr}`);
  }
  return lines.join("\n");
}
function _addToSet(set, seen, item) {
  for (const s of seen) {
    if (s === item) return;
    if (s && item && typeof s.equals === "function" && typeof item.equals === "function" && s.equals(item)) return;
  }
  seen.push(item);
  set.add(item);
}
function _decomposeOnce(circuit) {
  const c = circuit.copy();
  const newData = [];
  for (const ci of c.data) {
    const op = ci.operation;
    const defFn = op._definition;
    if (typeof defFn === "function") {
      try {
        const def = defFn(op);
        if (def && def.length > 0) {
          const defUniqueQubits = [];
          const seen = /* @__PURE__ */ new Set();
          for (const [_, sqs, scs] of def) {
            for (const sq of sqs) {
              if (!seen.has(sq)) {
                seen.add(sq);
                defUniqueQubits.push(sq);
              }
            }
          }
          const defUniqueClbits = [];
          const seenC = /* @__PURE__ */ new Set();
          for (const [_, sqs, scs] of def) {
            for (const sc of scs) {
              if (!seenC.has(sc)) {
                seenC.add(sc);
                defUniqueClbits.push(sc);
              }
            }
          }
          for (const [subOp, subQubits, subClbits] of def) {
            newData.push(new CircuitInstruction(
              subOp.copy(),
              subQubits.map((sq) => {
                const idx = defUniqueQubits.indexOf(sq);
                return idx >= 0 && idx < ci.qubits.length ? ci.qubits[idx] : sq;
              }),
              subClbits.map((sc) => {
                const idx = defUniqueClbits.indexOf(sc);
                return idx >= 0 && idx < ci.clbits.length ? ci.clbits[idx] : sc;
              })
            ));
          }
          continue;
        }
      } catch (e) {
      }
    }
    newData.push(ci);
  }
  c.data = newData;
  return c;
}

// src/library/standard_gates.js
var standard_gates_exports = {};
__export(standard_gates_exports, {
  CCXGate: () => CCXGate,
  CHGate: () => CHGate,
  CSXGate: () => CSXGate,
  CSwapGate: () => CSwapGate,
  CXGate: () => CXGate,
  CYGate: () => CYGate,
  CZGate: () => CZGate,
  DCXGate: () => DCXGate,
  HGate: () => HGate,
  IGate: () => IGate,
  ISwapGate: () => ISwapGate,
  SGate: () => SGate,
  SdgGate: () => SdgGate,
  SwapGate: () => SwapGate,
  SxGate: () => SxGate,
  SxdgGate: () => SxdgGate,
  TGate: () => TGate,
  TdgGate: () => TdgGate,
  XGate: () => XGate,
  YGate: () => YGate,
  ZGate: () => ZGate,
  makeCCXGate: () => makeCCXGate,
  makeCHGate: () => makeCHGate,
  makeCSXGate: () => makeCSXGate,
  makeCSwapGate: () => makeCSwapGate,
  makeCXGate: () => makeCXGate,
  makeCYGate: () => makeCYGate,
  makeCZGate: () => makeCZGate,
  makeDCXGate: () => makeDCXGate,
  makeHGate: () => makeHGate,
  makeIGate: () => makeIGate,
  makeISwapGate: () => makeISwapGate,
  makeSGate: () => makeSGate,
  makeSdgGate: () => makeSdgGate,
  makeSwapGate: () => makeSwapGate,
  makeSxGate: () => makeSxGate,
  makeSxdgGate: () => makeSxdgGate,
  makeTGate: () => makeTGate,
  makeTdgGate: () => makeTdgGate,
  makeXGate: () => makeXGate,
  makeYGate: () => makeYGate,
  makeZGate: () => makeZGate
});
var PI = Math.PI;
var SQRT2 = Math.SQRT2;
var r = (re, im = 0) => new Complex(re, im);
function makeHGate() {
  const g = new Gate("h", 1, []);
  const v = 1 / SQRT2;
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(v), r(v)],
    [r(v), r(-v)]
  ]);
  return g;
}
function makeXGate() {
  const g = new Gate("x", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(0), r(1)],
    [r(1), r(0)]
  ]);
  return g;
}
function makeYGate() {
  const g = new Gate("y", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(0), r(0, -1)],
    [r(0, 1), r(0)]
  ]);
  return g;
}
function makeZGate() {
  const g = new Gate("z", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0)],
    [r(0), r(-1)]
  ]);
  return g;
}
function makeSGate() {
  const g = new Gate("s", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0)],
    [r(0), r(0, 1)]
  ]);
  return g;
}
function makeSdgGate() {
  const g = new Gate("sdg", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0)],
    [r(0), r(0, -1)]
  ]);
  return g;
}
function makeTGate() {
  const g = new Gate("t", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0)],
    [r(0), r(Math.cos(PI / 4), Math.sin(PI / 4))]
  ]);
  return g;
}
function makeTdgGate() {
  const g = new Gate("tdg", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0)],
    [r(0), r(Math.cos(PI / 4), -Math.sin(PI / 4))]
  ]);
  return g;
}
function makeSxGate() {
  const g = new Gate("sx", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(0.5, 0.5), r(0.5, -0.5)],
    [r(0.5, -0.5), r(0.5, 0.5)]
  ]);
  return g;
}
function makeSxdgGate() {
  const g = new Gate("sxdg", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(0.5, -0.5), r(0.5, 0.5)],
    [r(0.5, 0.5), r(0.5, -0.5)]
  ]);
  return g;
}
function makeIGate() {
  const g = new Gate("id", 1, []);
  g._matrixBuilder = () => ComplexMatrix2.identity(2);
  return g;
}
function makeCXGate() {
  const g = new Gate("cx", 2, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0), r(0), r(1)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(1), r(0), r(0)]
  ]);
  return g;
}
function makeCYGate() {
  const g = new Gate("cy", 2, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0), r(0), r(0, -1)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(0, 1), r(0), r(0)]
  ]);
  return g;
}
function makeCZGate() {
  const g = new Gate("cz", 2, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(1), r(0), r(0)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(0), r(0), r(-1)]
  ]);
  return g;
}
function makeCHGate() {
  const g = new Gate("ch", 2, []);
  const v = 1 / SQRT2;
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(v), r(0), r(v)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(v), r(0), r(-v)]
  ]);
  return g;
}
function makeCSXGate() {
  const g = new Gate("csx", 2, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0.5, 0.5), r(0), r(0.5, -0.5)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(0.5, -0.5), r(0), r(0.5, 0.5)]
  ]);
  return g;
}
function makeSwapGate() {
  const g = new Gate("swap", 2, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(1), r(0), r(0)],
    [r(0), r(0), r(0), r(1)]
  ]);
  return g;
}
function makeISwapGate() {
  const g = new Gate("iswap", 2, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0), r(0, 1), r(0)],
    [r(0), r(0, 1), r(0), r(0)],
    [r(0), r(0), r(0), r(1)]
  ]);
  return g;
}
function makeDCXGate() {
  const g = new Gate("dcx", 2, []);
  g._matrixBuilder = () => ComplexMatrix2.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0), r(0), r(1)],
    [r(0), r(1), r(0), r(0)],
    [r(0), r(0), r(1), r(0)]
  ]);
  return g;
}
function makeCCXGate() {
  const g = new Gate("ccx", 3, []);
  g._matrixBuilder = () => {
    const m = ComplexMatrix2.identity(8);
    m.set(3, 3, r(0));
    m.set(7, 7, r(0));
    m.set(3, 7, r(1));
    m.set(7, 3, r(1));
    return m;
  };
  return g;
}
function makeCSwapGate() {
  const g = new Gate("cswap", 3, []);
  g._matrixBuilder = () => {
    const m = ComplexMatrix2.identity(8);
    m.set(3, 3, r(0));
    m.set(5, 5, r(0));
    m.set(3, 5, r(1));
    m.set(5, 3, r(1));
    return m;
  };
  return g;
}
var HGate = makeHGate();
var XGate = makeXGate();
var YGate = makeYGate();
var ZGate = makeZGate();
var SGate = makeSGate();
var SdgGate = makeSdgGate();
var TGate = makeTGate();
var TdgGate = makeTdgGate();
var SxGate = makeSxGate();
var SxdgGate = makeSxdgGate();
var IGate = makeIGate();
var CXGate = makeCXGate();
var CYGate = makeCYGate();
var CZGate = makeCZGate();
var CHGate = makeCHGate();
var CSXGate = makeCSXGate();
var SwapGate = makeSwapGate();
var ISwapGate = makeISwapGate();
var DCXGate = makeDCXGate();
var CCXGate = makeCCXGate();
var CSwapGate = makeCSwapGate();
_registerStd("H", HGate);
_registerStd("X", XGate);
_registerStd("Y", YGate);
_registerStd("Z", ZGate);
_registerStd("S", SGate);
_registerStd("SDG", SdgGate);
_registerStd("T", TGate);
_registerStd("TDG", TdgGate);
_registerStd("SX", SxGate);
_registerStd("SXDG", SxdgGate);
_registerStd("I", IGate);
_registerStd("CX", CXGate);
_registerStd("CY", CYGate);
_registerStd("CZ", CZGate);
_registerStd("CH", CHGate);
_registerStd("CSX", CSXGate);
_registerStd("SWAP", SwapGate);
_registerStd("ISWAP", ISwapGate);
_registerStd("DCX", DCXGate);
_registerStd("CCX", CCXGate);
_registerStd("CSWAP", CSwapGate);

// src/library/generalized_gates.js
var generalized_gates_exports = {};
__export(generalized_gates_exports, {
  makeCUGate: () => makeCUGate,
  makePGate: () => makePGate,
  makeRXGate: () => makeRXGate,
  makeRXXGate: () => makeRXXGate,
  makeRYGate: () => makeRYGate,
  makeRYYGate: () => makeRYYGate,
  makeRZGate: () => makeRZGate,
  makeRZXGate: () => makeRZXGate,
  makeRZZGate: () => makeRZZGate,
  makeU1Gate: () => makeU1Gate,
  makeU2Gate: () => makeU2Gate,
  makeU3Gate: () => makeU3Gate,
  makeUGate: () => makeUGate
});
var r2 = (re, im = 0) => new Complex(re, im);
function _val(p) {
  if (typeof p === "number") return p;
  if (p && typeof p.bind === "function") {
    throw new Error("Cannot build matrix with unbound parameter");
  }
  throw new TypeError(`Unsupported parameter type: ${typeof p}`);
}
function makeRXGate(theta) {
  const g = new Gate("rx", 1, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix2.fromRows([
      [r2(ct), r2(0, -st)],
      [r2(0, -st), r2(ct)]
    ]);
  };
  return g;
}
function makeRYGate(theta) {
  const g = new Gate("ry", 1, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix2.fromRows([
      [r2(ct), r2(st)],
      [r2(-st), r2(ct)]
    ]);
  };
  return g;
}
function makeRZGate(phi) {
  const g = new Gate("rz", 1, [phi]);
  g._matrixBuilder = (params) => {
    const p = _val(params[0]);
    return ComplexMatrix2.fromRows([
      [r2(Math.cos(p / 2), -Math.sin(p / 2)), r2(0)],
      [r2(0), r2(Math.cos(p / 2), Math.sin(p / 2))]
    ]);
  };
  return g;
}
function makePGate(theta) {
  const g = new Gate("p", 1, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    return ComplexMatrix2.fromRows([
      [r2(1), r2(0)],
      [r2(0), r2(Math.cos(t), Math.sin(t))]
    ]);
  };
  return g;
}
function makeU1Gate(lambda) {
  const g = new Gate("u1", 1, [lambda]);
  g._matrixBuilder = (params) => {
    const l = _val(params[0]);
    return ComplexMatrix2.fromRows([
      [r2(1), r2(0)],
      [r2(0), r2(Math.cos(l), Math.sin(l))]
    ]);
  };
  return g;
}
function makeU2Gate(phi, lambda) {
  const g = new Gate("u2", 1, [phi, lambda]);
  g._matrixBuilder = (params) => {
    const phi2 = _val(params[0]);
    const lam = _val(params[1]);
    const v = 1 / Math.SQRT2;
    return ComplexMatrix2.fromRows([
      [r2(v), r2(-v * Math.cos(lam), -v * Math.sin(lam))],
      [r2(v * Math.cos(phi2), v * Math.sin(phi2)), r2(v * Math.cos(lam + phi2), v * Math.sin(lam + phi2))]
    ]);
  };
  return g;
}
function makeU3Gate(theta, phi, lambda) {
  const g = new Gate("u3", 1, [theta, phi, lambda]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const p = _val(params[1]);
    const l = _val(params[2]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix2.fromRows([
      [r2(ct), r2(-st * Math.cos(l), -st * Math.sin(l))],
      [r2(st * Math.cos(p), st * Math.sin(p)), r2(ct * Math.cos(l + p), ct * Math.sin(l + p))]
    ]);
  };
  return g;
}
function makeUGate(theta, phi, lambda) {
  const g = makeU3Gate(theta, phi, lambda);
  g.name = "u";
  return g;
}
function makeCUGate(theta, phi, lambda, gamma) {
  const g = new Gate("cu", 2, [theta, phi, lambda, gamma]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const p = _val(params[1]);
    const l = _val(params[2]);
    const gm = _val(params[3]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    const gRe = Math.cos(gm);
    const gIm = Math.sin(gm);
    const cmul = (a, b) => [gRe * a - gIm * b, gRe * b + gIm * a];
    const m00 = cmul(ct, 0);
    const e_il = [Math.cos(l), Math.sin(l)];
    const tr_prod = cmul(e_il[0] * st, e_il[1] * st);
    const m01 = [-tr_prod[0], -tr_prod[1]];
    const e_ip = [Math.cos(p), Math.sin(p)];
    const bl_prod = cmul(e_ip[0] * st, e_ip[1] * st);
    const m10 = bl_prod;
    const e_ipl = [Math.cos(p + l), Math.sin(p + l)];
    const br_prod = cmul(e_ipl[0] * ct, e_ipl[1] * ct);
    const m11 = br_prod;
    return ComplexMatrix2.fromRows([
      [r2(1), r2(0), r2(0), r2(0)],
      [r2(0), r2(1), r2(0), r2(0)],
      [r2(0), r2(0), r2(m00[0], m00[1]), r2(m01[0], m01[1])],
      [r2(0), r2(0), r2(m10[0], m10[1]), r2(m11[0], m11[1])]
    ]);
  };
  return g;
}
function makeRXXGate(theta) {
  const g = new Gate("rxx", 2, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix2.fromRows([
      [r2(ct), r2(0), r2(0), r2(0, -st)],
      [r2(0), r2(ct), r2(0, -st), r2(0)],
      [r2(0), r2(0, -st), r2(ct), r2(0)],
      [r2(0, -st), r2(0), r2(0), r2(ct)]
    ]);
  };
  return g;
}
function makeRYYGate(theta) {
  const g = new Gate("ryy", 2, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix2.fromRows([
      [r2(ct), r2(0), r2(0), r2(0, st)],
      [r2(0), r2(ct), r2(0, -st), r2(0)],
      [r2(0), r2(0, -st), r2(ct), r2(0)],
      [r2(0, st), r2(0), r2(0), r2(ct)]
    ]);
  };
  return g;
}
function makeRZZGate(theta) {
  const g = new Gate("rzz", 2, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix2.fromRows([
      [r2(ct, -st), r2(0), r2(0), r2(0)],
      [r2(0), r2(ct, st), r2(0), r2(0)],
      [r2(0), r2(0), r2(ct, st), r2(0)],
      [r2(0), r2(0), r2(0), r2(ct, -st)]
    ]);
  };
  return g;
}
function makeRZXGate(theta) {
  const g = new Gate("rzx", 2, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix2.fromRows([
      [r2(ct), r2(0, -st), r2(0), r2(0)],
      [r2(0, -st), r2(ct), r2(0), r2(0)],
      [r2(0), r2(0), r2(ct), r2(0, st)],
      [r2(0), r2(0), r2(0, st), r2(ct)]
    ]);
  };
  return g;
}
_registerParamBuilder("RX", makeRXGate);
_registerParamBuilder("RY", makeRYGate);
_registerParamBuilder("RZ", makeRZGate);
_registerParamBuilder("P", makePGate);
_registerParamBuilder("U1", makeU1Gate);
_registerParamBuilder("U2", makeU2Gate);
_registerParamBuilder("U3", makeU3Gate);
_registerParamBuilder("U", makeUGate);
_registerParamBuilder("CU", makeCUGate);
_registerParamBuilder("RXX", makeRXXGate);
_registerParamBuilder("RYY", makeRYYGate);
_registerParamBuilder("RZZ", makeRZZGate);
_registerParamBuilder("RZX", makeRZXGate);

// src/quantum_info/operator.js
var Operator = class _Operator {
  constructor(data) {
    if (data instanceof ComplexMatrix2) {
      this._data = data;
    } else if (Array.isArray(data) && Array.isArray(data[0])) {
      this._data = ComplexMatrix2.fromRows(data);
    } else {
      throw new TypeError("Operator requires a ComplexMatrix or 2D array");
    }
    if (this._data.rows !== this._data.cols) {
      throw new Error("Operator must be square");
    }
    this._numQubits = Math.log2(this._data.rows);
    if (!Number.isInteger(this._numQubits)) {
      throw new Error("Operator dimension must be 2^n");
    }
  }
  static fromMatrix(matrix) {
    return new _Operator(matrix);
  }
  static fromLabel(label) {
    const c = (re, im = 0) => new Complex(re, im);
    const gates = {
      "I": ComplexMatrix2.identity(2),
      "X": PAULI.X,
      "Y": PAULI.Y,
      "Z": PAULI.Z,
      "0": ComplexMatrix2.fromRows([[c(1), c(0)], [c(0), c(0)]]),
      "1": ComplexMatrix2.fromRows([[c(0), c(0)], [c(0), c(1)]]),
      "+": ComplexMatrix2.fromRows([[c(0.5), c(0.5)], [c(0.5), c(0.5)]]),
      "-": ComplexMatrix2.fromRows([[c(0.5), c(-0.5)], [c(-0.5), c(0.5)]]),
      "r": ComplexMatrix2.fromRows([[c(0.5, 0.5), c(0)], [c(0.5, -0.5), c(0)]]),
      "l": ComplexMatrix2.fromRows([[c(0.5, -0.5), c(0)], [c(0.5, 0.5), c(0)]])
    };
    const mats = [];
    for (let i = label.length - 1; i >= 0; i--) {
      const ch = label[i].toUpperCase();
      if (!gates[ch]) throw new Error(`Unknown label character: ${ch}`);
      mats.push(gates[ch]);
    }
    return new _Operator(kronMatrices(mats));
  }
  static fromCircuit(circuit) {
    const op = _Operator.identity(circuit.num_qubits);
    for (const ci of circuit.data) {
      if (ci.operation.name === "barrier" || ci.operation.name === "measure" || ci.operation.name === "reset") continue;
      if (ci.operation.num_qubits === 0) continue;
      const gateOp = _Operator.fromGate(ci.operation);
      const fullOp = gateOp.embedIntoCircuit(circuit.num_qubits, ci.qubits.map((q) => circuit._qubit_index.get(q)));
      op._data = fullOp._data.mul(op._data);
    }
    return op;
  }
  static fromGate(gate) {
    if (typeof gate.to_matrix === "function") {
      return new _Operator(gate.to_matrix());
    }
    throw new TypeError("Cannot construct Operator from gate without matrix");
  }
  static identity(numQubits) {
    return new _Operator(ComplexMatrix2.identity(1 << numQubits));
  }
  static zero(numQubits) {
    return new _Operator(ComplexMatrix2.zeros(1 << numQubits));
  }
  static fromPhase(phase, numQubits = 0) {
    return new _Operator(ComplexMatrix2.identity(1 << numQubits).scale(phase));
  }
  get data() {
    return this._data;
  }
  get dim() {
    return this._data.rows;
  }
  get num_qubits() {
    return this._numQubits;
  }
  compose(other) {
    if (this._numQubits !== other._numQubits) {
      throw new Error("compose: qubit count mismatch");
    }
    return new _Operator(this._data.mul(other._data));
  }
  tensor(other) {
    return new _Operator(this._data.tensor(other._data));
  }
  expand(other) {
    return new _Operator(other._data.tensor(this._data));
  }
  adjoint() {
    return new _Operator(this._data.dagger());
  }
  conjugate() {
    return new _Operator(this._data.conjugate());
  }
  transpose() {
    return new _Operator(this._data.transpose());
  }
  trace() {
    return this._data.trace();
  }
  det() {
    return this._data.det();
  }
  is_unitary(tol) {
    return this._data.isUnitary(tol);
  }
  is_hermitian(tol) {
    return this._data.isHermitian(tol);
  }
  pow(n) {
    if (n === 0) return _Operator.identity(this._numQubits);
    let result = this;
    for (let i = 1; i < n; i++) result = result.compose(this);
    return result;
  }
  exp() {
    return new _Operator(this._data.expm());
  }
  equals(other, tol) {
    return other instanceof _Operator && this._data.equals(other._data, tol);
  }
  equals_up_to_phase(other, tol) {
    const eps = typeof tol === "number" ? tol : 1e-9;
    if (!(other instanceof _Operator) || this._numQubits !== other._numQubits) return false;
    let phase = null;
    for (let i = 0; i < this._data.rows; i++) {
      for (let j = 0; j < this._data.cols; j++) {
        const a = this._data.get(i, j);
        const b = other._data.get(i, j);
        if (a.abs() > eps) {
          if (b.abs() < eps) return false;
          const r4 = b.div(a);
          if (phase === null) phase = r4;
          else if (!phase.equals(r4, eps)) return false;
        } else if (b.abs() > eps) {
          return false;
        }
      }
    }
    return true;
  }
  embedIntoCircuit(numQubits, qubitIndices) {
    if (qubitIndices.length !== this._numQubits) {
      throw new Error("embedIntoCircuit: qubit count mismatch");
    }
    const dim = 1 << numQubits;
    const out = ComplexMatrix2.zeros(dim, dim);
    const indices = qubitIndices;
    for (let row = 0; row < dim; row++) {
      for (let col = 0; col < dim; col++) {
        let subRow = 0;
        let subCol = 0;
        let restMatches = true;
        for (let qi = 0; qi < indices.length; qi++) {
          const bit = row >> indices[qi] & 1;
          subRow |= bit << qi;
        }
        for (let qi = 0; qi < indices.length; qi++) {
          const bit = col >> indices[qi] & 1;
          subCol |= bit << qi;
        }
        for (let q = 0; q < numQubits; q++) {
          if (indices.indexOf(q) === -1) {
            const rb = row >> q & 1;
            const cb = col >> q & 1;
            if (rb !== cb) {
              restMatches = false;
              break;
            }
          }
        }
        if (!restMatches) continue;
        const elem = this._data.get(subRow, subCol);
        out.set(row, col, elem);
      }
    }
    return new _Operator(out);
  }
  // Apply to a statevector
  apply_to_vector(vector) {
    return this._data.matvec(vector);
  }
  // Partial trace: trace out the qubits in `qubitsToTraceOut`, returning an
  // Operator on the remaining qubits.
  partial_trace(qubitsToTraceOut) {
    const traceOut = (Array.isArray(qubitsToTraceOut) ? qubitsToTraceOut : [qubitsToTraceOut]).slice();
    if (traceOut.length === 0) return this;
    const n = this._numQubits;
    let current = this._data;
    let currentN = n;
    traceOut.sort((a, b) => b - a);
    for (const q of traceOut) {
      current = current.partialTrace(currentN, q);
      currentN -= 1;
    }
    return new _Operator(current);
  }
  // Sum of two operators: returns this + other.
  sum(other) {
    if (this._numQubits !== other._numQubits) {
      throw new Error("sum: qubit count mismatch");
    }
    const result = ComplexMatrix2.zeros(this._data.rows, this._data.cols);
    for (let i = 0; i < this._data.rows; i++) {
      for (let j = 0; j < this._data.cols; j++) {
        result.set(i, j, this._data.get(i, j).add(other._data.get(i, j)));
      }
    }
    return new _Operator(result);
  }
  // Subtract: returns this - other.
  subtract(other) {
    if (this._numQubits !== other._numQubits) {
      throw new Error("subtract: qubit count mismatch");
    }
    const result = ComplexMatrix2.zeros(this._data.rows, this._data.cols);
    for (let i = 0; i < this._data.rows; i++) {
      for (let j = 0; j < this._data.cols; j++) {
        result.set(i, j, this._data.get(i, j).sub(other._data.get(i, j)));
      }
    }
    return new _Operator(result);
  }
  // Convert to a measurement channel: project the input state onto the operator
  to_matrix() {
    return this._data;
  }
  // Eigenvalues (for Hermitian operators)
  eigvals() {
    if (!this.is_hermitian()) {
      throw new Error("eigvals only supports Hermitian operators (use svd for general)");
    }
    return this._data.eigh().eigenvalues;
  }
  // SVD decomposition
  svd() {
    return this._data.svd();
  }
  toDict() {
    const rows = this._data.toRows();
    return {
      rows: this._data.rows,
      cols: this._data.cols,
      data: rows.map((r4) => r4.map((c) => ({ re: c.re, im: c.im })))
    };
  }
  toString() {
    return `Operator(num_qubits=${this._numQubits})
${this._data.toString()}`;
  }
};
_registerOperatorClass(Operator);

// src/quantum_info/statevector.js
var _DensityMatrixClass = null;
function _registerDensityMatrixClass(cls) {
  _DensityMatrixClass = cls;
}
var ZERO = new ComplexVector([new Complex(1, 0), new Complex(0, 0)]);
var ONE = new ComplexVector([new Complex(0, 0), new Complex(1, 0)]);
var LABEL_VECTORS = {
  "0": ZERO,
  "1": ONE,
  "+": new ComplexVector([new Complex(1 / Math.SQRT2, 0), new Complex(1 / Math.SQRT2, 0)]),
  "-": new ComplexVector([new Complex(1 / Math.SQRT2, 0), new Complex(-1 / Math.SQRT2, 0)]),
  "r": new ComplexVector([new Complex(1 / Math.SQRT2, 0), new Complex(0, 1 / Math.SQRT2)]),
  "l": new ComplexVector([new Complex(1 / Math.SQRT2, 0), new Complex(0, -1 / Math.SQRT2)])
};
var Statevector = class _Statevector {
  constructor(data, numQubits = null) {
    if (data instanceof ComplexVector) {
      this._data = data;
    } else if (Array.isArray(data)) {
      this._data = new ComplexVector(data);
    } else {
      throw new TypeError("Statevector requires a ComplexVector or array");
    }
    const dim = this._data.size;
    const nq = Math.log2(dim);
    if (!Number.isInteger(nq)) {
      throw new Error(`Statevector dimension must be 2^n, got ${dim}`);
    }
    this._numQubits = numQubits !== null ? numQubits : nq;
  }
  static fromLabel(label) {
    const parts = [];
    for (let i = 0; i < label.length; i++) {
      const ch = label[i].toLowerCase();
      if (!LABEL_VECTORS[ch]) throw new Error(`Unknown state label: ${ch}`);
      parts.push(LABEL_VECTORS[ch]);
    }
    if (parts.length === 0) return new _Statevector(new ComplexVector([Complex.ONE]), 0);
    return new _Statevector(kronVectors(parts), parts.length);
  }
  static zero(numQubits) {
    const data = ComplexVector.zeros(1 << numQubits);
    data.data[0] = Complex.ONE;
    return new _Statevector(data, numQubits);
  }
  static one(numQubits) {
    const data = ComplexVector.zeros(1 << numQubits);
    data.data[(1 << numQubits) - 1] = Complex.ONE;
    return new _Statevector(data, numQubits);
  }
  static fromInstruction(instruction) {
    const nq = instruction.num_qubits;
    const sv = _Statevector.zero(nq);
    if (typeof instruction.to_matrix === "function") {
      const m = instruction.to_matrix();
      return new _Statevector(m.matvec(sv._data), nq);
    }
    throw new Error("Instruction has no matrix");
  }
  static fromCircuit(circuit, initState = null) {
    let sv = initState ? initState.copy() : _Statevector.zero(circuit.num_qubits);
    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier") continue;
      if (op.name === "measure" || op.name === "reset") continue;
      if (op.num_qubits === 0) continue;
      if (typeof op.to_matrix !== "function") {
        throw new Error(`Cannot evolve Statevector under non-unitary ${op.name}`);
      }
      const fullOp = _embedGate(op, circuit.num_qubits, ci.qubits.map((q) => circuit._qubit_index.get(q)));
      sv = new _Statevector(fullOp.matvec(sv._data), circuit.num_qubits);
    }
    return sv;
  }
  static from_int(index, numQubits) {
    const data = ComplexVector.zeros(1 << numQubits);
    data.data[index] = Complex.ONE;
    return new _Statevector(data, numQubits);
  }
  get data() {
    return this._data;
  }
  get dim() {
    return this._data.size;
  }
  get num_qubits() {
    return this._numQubits;
  }
  copy() {
    return new _Statevector(new ComplexVector(this._data.data.slice()), this._numQubits);
  }
  conjugate() {
    return new _Statevector(this._data.conjugate(), this._numQubits);
  }
  inner(other) {
    return this._data.inner(other._data);
  }
  dot(other) {
    return this._data.inner(other._data);
  }
  norm() {
    return this._data.norm();
  }
  is_unitary() {
    return false;
  }
  probabilities(qargs = null) {
    const probs = this._data.probabilities();
    if (qargs == null) return probs;
    const indices = Array.isArray(qargs) ? qargs : [qargs];
    return _marginalProbabilities(probs, this._numQubits, indices);
  }
  evolve(other) {
    if (other instanceof _Statevector) {
      return new _Statevector(this._data.tensor(other._data), this._numQubits + other._numQubits);
    }
    if (other && other._data && other._data instanceof ComplexMatrix2) {
      if (other.num_qubits !== this._numQubits) {
        throw new Error("evolve: operator qubit count mismatch");
      }
      return new _Statevector(other._data.matvec(this._data), this._numQubits);
    }
    if (other && typeof other.to_matrix === "function") {
      if (other.num_qubits !== this._numQubits) {
        throw new Error("evolve: gate qubit count mismatch");
      }
      return new _Statevector(other.to_matrix().matvec(this._data), this._numQubits);
    }
    if (other && other.data && other.data instanceof Array) {
      if (typeof _Statevector.fromCircuit !== "function") {
        throw new Error("Cannot evolve by circuit without circuit module loaded");
      }
      return _Statevector.fromCircuit(other, this);
    }
    throw new TypeError("evolve: unsupported operand type");
  }
  tensor(other) {
    return new _Statevector(this._data.tensor(other._data), this._numQubits + other._numQubits);
  }
  expand(other) {
    return new _Statevector(other._data.tensor(this._data), this._numQubits + other._numQubits);
  }
  sample(counts = 1024, qargs = null, rng) {
    const probs = this.probabilities(qargs);
    const results = {};
    const nq = qargs == null ? this._numQubits : Array.isArray(qargs) ? qargs.length : 1;
    for (let i = 0; i < counts; i++) {
      const idx = sampleDistribution(probs, rng);
      const bits = new Array(nq).fill("0");
      for (let q = 0; q < nq; q++) bits[q] = (idx >> q & 1).toString();
      const key = bits.reverse().join("");
      results[key] = (results[key] || 0) + 1;
    }
    return results;
  }
  sample_memory(counts = 1024, qargs = null, rng) {
    const probs = this.probabilities(qargs);
    const nq = qargs == null ? this._numQubits : Array.isArray(qargs) ? qargs.length : 1;
    const results = new Array(counts);
    for (let i = 0; i < counts; i++) {
      const idx = sampleDistribution(probs, rng);
      const bits = new Array(nq).fill("0");
      for (let q = 0; q < nq; q++) bits[q] = (idx >> q & 1).toString();
      results[i] = bits.reverse().join("");
    }
    return results;
  }
  expectation_value(operator) {
    const m = operator instanceof ComplexMatrix2 ? operator : operator._data;
    const op = m.matvec(this._data);
    return this._data.inner(op);
  }
  // Compute the density matrix rho = |psi><psi|
  to_operator() {
    const dim = this._data.size;
    const m = ComplexMatrix2.zeros(dim, dim);
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        m.set(i, j, this._data.get(i).mul(this._data.get(j).conjugate()));
      }
    }
    return new Operator(m);
  }
  equals(other, tol) {
    if (!(other instanceof _Statevector)) return false;
    if (other._numQubits !== this._numQubits) return false;
    return this._data.equals(other._data, tol);
  }
  // Check equality up to global phase
  equiv(other, tol = 1e-9) {
    if (!(other instanceof _Statevector) || other._numQubits !== this._numQubits) return false;
    let phase = null;
    for (let i = 0; i < this._data.size; i++) {
      const a = this._data.get(i);
      const b = other._data.get(i);
      if (a.abs() > tol) {
        if (b.abs() < tol) return false;
        const r4 = b.div(a);
        if (phase === null) phase = r4;
        else if (!phase.equals(r4, tol)) return false;
      } else if (b.abs() > tol) {
        return false;
      }
    }
    return true;
  }
  to_dict() {
    const out = {};
    for (let i = 0; i < this._data.size; i++) {
      const bits = i.toString(2).padStart(this._numQubits, "0");
      out[bits] = this._data.get(i);
    }
    return out;
  }
  // Measure a single qubit (or all qubits if no argument) and collapse
  // the state. Returns a { bit, statevector } pair.
  measure(qubit = null, rng) {
    if (qubit == null) {
      const probs2 = this.probabilities();
      const idx = sampleDistribution(probs2, rng);
      const bits = new Array(this._numQubits).fill("0");
      for (let q = 0; q < this._numQubits; q++) bits[q] = (idx >> q & 1).toString();
      const bitstring = bits.reverse().join("");
      const collapsed = new Array(this._data.size).fill(Complex.ZERO);
      collapsed[idx] = Complex.ONE;
      return { bits: bitstring, statevector: new _Statevector(new ComplexVector(collapsed), this._numQubits) };
    }
    const probs = this.probabilities([qubit]);
    const prob0 = probs[0];
    const r4 = rng ? rng() : Math.random();
    const outcome = r4 < prob0 ? 0 : 1;
    const dim = this._data.size;
    const newData = new Array(dim);
    let normSq = 0;
    for (let i = 0; i < dim; i++) {
      if ((i >> qubit & 1) === outcome) {
        newData[i] = this._data.get(i);
        normSq += newData[i].re * newData[i].re + newData[i].im * newData[i].im;
      } else {
        newData[i] = Complex.ZERO;
      }
    }
    const norm = Math.sqrt(normSq);
    for (let i = 0; i < dim; i++) {
      newData[i] = new Complex(newData[i].re / norm, newData[i].im / norm);
    }
    return { bit: outcome, statevector: new _Statevector(new ComplexVector(newData), this._numQubits) };
  }
  // Reset a qubit to |0>: measure it, and if outcome is |1>, apply X.
  reset(qubit) {
    const { bit, statevector } = this.measure(qubit);
    if (bit === 0) return statevector;
    const X = PAULI.X;
    const dim = statevector._data.size;
    const newData = new Array(dim);
    for (let i = 0; i < dim; i++) {
      const j = i ^ 1 << qubit;
      newData[i] = statevector._data.get(j);
    }
    return new _Statevector(new ComplexVector(newData), this._numQubits);
  }
  // Partial trace: trace out the qubits NOT in `keep`, returning a
  // DensityMatrix on the kept qubits.
  partial_trace(qargsToKeep) {
    const keep = Array.isArray(qargsToKeep) ? qargsToKeep : [qargsToKeep];
    const rho = this.to_operator();
    const traceOut = _complement(keep, this._numQubits);
    const reducedOp = rho.partial_trace(traceOut);
    if (!_DensityMatrixClass) {
      throw new Error("DensityMatrix class not registered. Import quantum_info/density_matrix.js first.");
    }
    return new _DensityMatrixClass(reducedOp._data);
  }
  // Expand dimensions: tensor with |0...0> on `numQubits` additional qubits
  // at the front (high-order bits).
  expand_dims(numQubits) {
    const newDim = this._data.size << numQubits;
    const newData = new Array(newDim).fill(Complex.ZERO);
    for (let i = 0; i < this._data.size; i++) {
      newData[i << numQubits] = this._data.get(i);
    }
    return new _Statevector(new ComplexVector(newData), this._numQubits + numQubits);
  }
  toString() {
    const terms = [];
    for (let i = 0; i < this._data.size; i++) {
      const amp = this._data.get(i);
      if (amp.abs() < 1e-10) continue;
      const bits = i.toString(2).padStart(this._numQubits, "0");
      terms.push(`${amp.toString()}|${bits}>`);
    }
    return terms.join(" + ");
  }
};
function _embedGate(gate, numQubits, qubitIndices) {
  if (typeof gate.to_matrix !== "function") {
    throw new Error("Cannot embed gate without matrix");
  }
  const gateMatrix = gate.to_matrix();
  if (gate.num_qubits !== qubitIndices.length) {
    throw new Error("embedGate: qubit count mismatch");
  }
  const dim = 1 << numQubits;
  const out = ComplexMatrix2.zeros(dim, dim);
  const k = qubitIndices.length;
  for (let row = 0; row < dim; row++) {
    for (let col = 0; col < dim; col++) {
      let subRow = 0, subCol = 0, restMatches = true;
      for (let qi = 0; qi < k; qi++) {
        const rb = row >> qubitIndices[qi] & 1;
        const cb = col >> qubitIndices[qi] & 1;
        subRow |= rb << qi;
        subCol |= cb << qi;
      }
      for (let q = 0; q < numQubits; q++) {
        if (qubitIndices.indexOf(q) === -1) {
          const rb = row >> q & 1;
          const cb = col >> q & 1;
          if (rb !== cb) {
            restMatches = false;
            break;
          }
        }
      }
      if (!restMatches) continue;
      out.set(row, col, gateMatrix.get(subRow, subCol));
    }
  }
  return out;
}
function _marginalProbabilities(probs, numQubits, indices) {
  const sortedIndices = indices.slice().sort((a, b) => a - b);
  const outDim = 1 << indices.length;
  const out = new Array(outDim).fill(0);
  for (let i = 0; i < probs.length; i++) {
    let subIdx = 0;
    for (let qi = 0; qi < indices.length; qi++) {
      const bit = i >> indices[qi] & 1;
      subIdx |= bit << qi;
    }
    out[subIdx] += probs[i];
  }
  return out;
}
_registerStatevectorClass(Statevector);
function _complement(subset, n) {
  const set = new Set(subset);
  const out = [];
  for (let i = 0; i < n; i++) {
    if (!set.has(i)) out.push(i);
  }
  return out;
}

// src/quantum_info/pauli.js
var PAULI_MATRICES = {
  "I": ComplexMatrix2.identity(2),
  "X": PAULI.X,
  "Y": PAULI.Y,
  "Z": PAULI.Z
};
var Pauli = class _Pauli {
  constructor(label) {
    if (typeof label !== "string") {
      throw new TypeError("Pauli label must be a string");
    }
    this.label = label.toUpperCase();
    const n = label.length;
    this.x = new Array(n).fill(0);
    this.z = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
      const ch = this.label[n - 1 - i];
      if (ch === "X") {
        this.x[i] = 1;
      } else if (ch === "Z") {
        this.z[i] = 1;
      } else if (ch === "Y") {
        this.x[i] = 1;
        this.z[i] = 1;
      } else if (ch === "I") {
      } else throw new Error(`Invalid Pauli label char: ${ch}`);
    }
    this.num_qubits = n;
  }
  static fromSymplectic(x, z) {
    if (x.length !== z.length) throw new Error("x and z must have same length");
    const n = x.length;
    let label = "";
    for (let i = n - 1; i >= 0; i--) {
      if (x[i] && z[i]) label += "Y";
      else if (x[i]) label += "X";
      else if (z[i]) label += "Z";
      else label += "I";
    }
    return new _Pauli(label);
  }
  // Build the matrix for this Pauli operator.
  //
  // Symplectic convention used throughout Ketra (matching the
  // SparsePauliOp expectation_value fast path):
  //
  //   P(x, z) = i^(x·z) · X^x · Z^z    (per qubit, then tensor product)
  //
  // For one qubit this gives:
  //   P(0,0) = I
  //   P(1,0) = X
  //   P(0,1) = Z
  //   P(1,1) = i · X · Z = i · (XZ) = i · (−iY) = Y     ← the i^(xz)
  //            factor is exactly what turns XZ into Y, so the matrix for
  //            the (1,1) entry is the standard Pauli-Y matrix and NO
  //            extra global phase is applied on top. (Applying i^#Y again
  //            would double-count the i and produce a wrong sign — that
  //            was the previous bug.)
  to_matrix() {
    const mats = [];
    for (let i = 0; i < this.num_qubits; i++) {
      if (this.x[i] && this.z[i]) mats.push(PAULI_MATRICES["Y"]);
      else if (this.x[i]) mats.push(PAULI_MATRICES["X"]);
      else if (this.z[i]) mats.push(PAULI_MATRICES["Z"]);
      else mats.push(PAULI_MATRICES["I"]);
    }
    let result = mats[0];
    for (let i = 1; i < mats.length; i++) result = result.tensor(mats[i]);
    return result;
  }
  // Multiply two Paulis in the symplectic representation.
  //
  //   P(x1,z1) · P(x2,z2) = phase · P(x1⊕x2, z1⊕z2)
  //
  // Derivation (per qubit, all arithmetic mod 4 for the i-exponent):
  //   P(x1,z1) = i^(x1 z1) X^x1 Z^z1
  //   P(x2,z2) = i^(x2 z2) X^x2 Z^z2
  //   Product  = i^(x1 z1 + x2 z2) · X^x1 Z^z1 X^x2 Z^z2
  // X and Z anticommute: Z·X = −X·Z, so moving Z^z1 past X^x2 contributes
  // a factor (−1)^(z1 x2) = i^(2 z1 x2). The remaining X and Z parts add
  // mod 2, and the new (x·z) factor for the result Pauli is (x1⊕x2)(z1⊕z2),
  // which we have to remove (subtract) since the convention encodes it
  // implicitly via Y = i·XZ. Putting it together:
  //
  //   phase_exp = x1 z1 + x2 z2 − (x1⊕x2)(z1⊕z2) + 2 z1 x2     (mod 4)
  //
  // Verified by exhaustive enumeration over all 16 single-qubit cases
  // (X·Y = iZ, Y·X = −iZ, Z·X = iY, Y·Z = iX, Z·Y = −iX, X·Z = −iY,
  // Y·Y = I, X·X = I, Z·Z = I, etc.).
  compose(other) {
    if (this.num_qubits !== other.num_qubits) {
      throw new Error("Pauli composition: qubit count mismatch");
    }
    const newX = new Array(this.num_qubits).fill(0);
    const newZ = new Array(this.num_qubits).fill(0);
    let phase = 0;
    for (let i = 0; i < this.num_qubits; i++) {
      const x1 = this.x[i], z1 = this.z[i];
      const x2 = other.x[i], z2 = other.z[i];
      newX[i] = x1 ^ x2;
      newZ[i] = z1 ^ z2;
      const yzProd = (x1 ^ x2) & (z1 ^ z2);
      phase += x1 * z1 + x2 * z2 - yzProd + 2 * (z1 & x2);
    }
    phase = (phase % 4 + 4) % 4;
    const result = _Pauli.fromSymplectic(newX, newZ);
    const phaseFactor = [new Complex(1, 0), new Complex(0, 1), new Complex(-1, 0), new Complex(0, -1)][phase];
    return { pauli: result, phase: phaseFactor };
  }
  multiply(other) {
    const r4 = this.compose(other);
    return r4.pauli.to_matrix().scale(r4.phase);
  }
  to_label() {
    return this.label;
  }
  equals(other) {
    if (!(other instanceof _Pauli)) return false;
    if (this.num_qubits !== other.num_qubits) return false;
    for (let i = 0; i < this.num_qubits; i++) {
      if (this.x[i] !== other.x[i] || this.z[i] !== other.z[i]) return false;
    }
    return true;
  }
  weight() {
    let w = 0;
    for (let i = 0; i < this.num_qubits; i++) {
      if (this.x[i] || this.z[i]) w++;
    }
    return w;
  }
  commutes(other) {
    if (this.num_qubits !== other.num_qubits) {
      throw new Error("Pauli commutes: qubit count mismatch");
    }
    let anticommute = 0;
    for (let i = 0; i < this.num_qubits; i++) {
      anticommute += (this.x[i] * other.z[i] + this.z[i] * other.x[i]) % 2;
    }
    return anticommute % 2 === 0;
  }
  anticommutes(other) {
    return !this.commutes(other);
  }
  // Apply X/Y/Z to specific qubits (single-qubit Pauli evolution)
  evolve(other) {
    if (this.num_qubits !== other.num_qubits) {
      throw new Error("evolve: qubit count mismatch");
    }
    const newX = new Array(this.num_qubits).fill(0);
    const newZ = new Array(this.num_qubits).fill(0);
    for (let i = 0; i < this.num_qubits; i++) {
      newX[i] = this.x[i] ^ other.x[i];
      newZ[i] = this.z[i] ^ other.z[i];
    }
    return _Pauli.fromSymplectic(newX, newZ);
  }
  toString() {
    return `Pauli(${this.label})`;
  }
};
var PauliList = class _PauliList {
  constructor(labels) {
    if (Array.isArray(labels)) {
      this.paulis = labels.map((l) => l instanceof Pauli ? l : new Pauli(l));
    } else {
      this.paulis = [labels instanceof Pauli ? labels : new Pauli(labels)];
    }
    if (this.paulis.length > 0) {
      this.num_qubits = this.paulis[0].num_qubits;
      for (const p of this.paulis) {
        if (p.num_qubits !== this.num_qubits) throw new Error("All Paulis must have same length");
      }
    } else {
      this.num_qubits = 0;
    }
  }
  get size() {
    return this.paulis.length;
  }
  get length() {
    return this.paulis.length;
  }
  get(i) {
    return this.paulis[i];
  }
  [Symbol.iterator]() {
    return this.paulis[Symbol.iterator]();
  }
  to_labels() {
    return this.paulis.map((p) => p.label);
  }
  to_matrices() {
    return this.paulis.map((p) => p.to_matrix());
  }
  weights() {
    return this.paulis.map((p) => p.weight());
  }
  append(other) {
    if (other instanceof Pauli) this.paulis.push(other);
    else if (other instanceof _PauliList) this.paulis.push(...other.paulis);
    else throw new TypeError("append: expects Pauli or PauliList");
  }
};
var SparsePauliOp = class _SparsePauliOp {
  constructor(paulis, coeffs = null) {
    if (paulis instanceof PauliList) this.paulis = paulis;
    else this.paulis = new PauliList(paulis);
    const n = this.paulis.size;
    if (coeffs == null) {
      this.coeffs = new Array(n).fill(0).map(() => new Complex(1, 0));
    } else {
      this.coeffs = coeffs.map((c) => c instanceof Complex ? c : new Complex(c, 0));
    }
    if (this.coeffs.length !== n) {
      throw new Error("SparsePauliOp: paulis and coeffs length mismatch");
    }
    this.num_qubits = this.paulis.num_qubits;
  }
  static from_list(list) {
    const paulis = list.map(([p]) => p instanceof Pauli ? p : new Pauli(p));
    const coeffs = list.map(([_, c]) => c instanceof Complex ? c : new Complex(c, 0));
    return new _SparsePauliOp(new PauliList(paulis), coeffs);
  }
  static from_sparse_list(zlist, numQubits) {
    const paulis = zlist.map(([label]) => {
      let padded = label;
      while (padded.length < numQubits) padded = "I" + padded;
      return new Pauli(padded);
    });
    const coeffs = zlist.map(([_, c]) => c instanceof Complex ? c : new Complex(c, 0));
    return new _SparsePauliOp(new PauliList(paulis), coeffs);
  }
  static from_operator(operator) {
    const n = operator.num_qubits;
    const dim = 1 << n;
    const pauliLabels = [];
    function genLabels(prefix, depth) {
      if (depth === 0) {
        pauliLabels.push(prefix);
        return;
      }
      for (const c of ["I", "X", "Y", "Z"]) genLabels(prefix + c, depth - 1);
    }
    genLabels("", n);
    const results = [];
    for (const label of pauliLabels) {
      const p = new Pauli(label);
      const m = p.to_matrix();
      let coeff = new Complex(0, 0);
      for (let i = 0; i < dim; i++) {
        for (let j = 0; j < dim; j++) {
          coeff = coeff.add(m.get(i, j).conjugate().mul(operator._data.get(i, j)));
        }
      }
      coeff = coeff.scale(1 / dim);
      if (coeff.abs() > 1e-10) {
        results.push([p, coeff]);
      }
    }
    return _SparsePauliOp.from_list(results);
  }
  to_matrix() {
    const dim = 1 << this.num_qubits;
    const result = ComplexMatrix2.zeros(dim, dim);
    for (let i = 0; i < this.paulis.size; i++) {
      const m = this.paulis.get(i).to_matrix();
      const c = this.coeffs[i];
      for (let k = 0; k < m.data.length; k++) {
        result.data[k] = result.data[k].add(m.data[k].mul(c));
      }
    }
    return result;
  }
  to_list() {
    const out = [];
    for (let i = 0; i < this.paulis.size; i++) {
      out.push([this.paulis.get(i).label, this.coeffs[i]]);
    }
    return out;
  }
  add(other) {
    if (other instanceof _SparsePauliOp) {
      return new _SparsePauliOp(
        new PauliList(this.paulis.paulis.concat(other.paulis.paulis)),
        this.coeffs.concat(other.coeffs)
      );
    }
    throw new TypeError("add: unsupported type");
  }
  // Matrix product: this @ other. For SparsePauliOp this corresponds to
  // multiplying the two Pauli sums together, which produces up to N*M terms
  // (one for each pair). The result is then simplified.
  dot(other) {
    if (other instanceof _SparsePauliOp) {
      const newPaulis = [];
      const newCoeffs = [];
      for (let i = 0; i < this.paulis.size; i++) {
        const pi = this.paulis.get(i);
        const ci = this.coeffs[i];
        for (let j = 0; j < other.paulis.size; j++) {
          const pj = other.paulis.get(j);
          const cj = other.coeffs[j];
          const r4 = pi.compose(pj);
          newPaulis.push(r4.pauli);
          newCoeffs.push(ci.mul(r4.phase).mul(cj));
        }
      }
      return new _SparsePauliOp(new PauliList(newPaulis), newCoeffs).simplify();
    }
    throw new TypeError("dot: unsupported type");
  }
  // Iterate over non-zero (row, col, value) entries of the operator's
  // matrix representation. Yields [row, col, Complex] tuples. Useful for
  // sparse matrix construction without materializing the full matrix.
  *matrix_iter() {
    const dim = 1 << this.num_qubits;
    for (let ti = 0; ti < this.paulis.size; ti++) {
      const pauli = this.paulis.get(ti);
      const coeff = this.coeffs[ti];
      const x = pauli.x, z = pauli.z;
      let mask = 0, yMask = 0, zOnlyMask = 0, yCount = 0;
      for (let q = 0; q < this.num_qubits; q++) {
        if (x[q]) {
          mask |= 1 << q;
          if (z[q]) {
            yMask |= 1 << q;
            yCount++;
          }
        } else if (z[q]) {
          zOnlyMask |= 1 << q;
        }
      }
      const iPow = [new Complex(1, 0), new Complex(0, 1), new Complex(-1, 0), new Complex(0, -1)][yCount % 4];
      for (let i = 0; i < dim; i++) {
        const j = i ^ mask;
        let parity = 0;
        for (let q = 0; q < this.num_qubits; q++) {
          if (zOnlyMask & 1 << q & i) parity ^= 1;
          if (yMask & 1 << q & j) parity ^= 1;
        }
        const sign = parity ? -1 : 1;
        const vRe = sign * iPow.re * coeff.re - sign * iPow.im * coeff.im;
        const vIm = sign * iPow.re * coeff.im + sign * iPow.im * coeff.re;
        if (Math.abs(vRe) > 1e-15 || Math.abs(vIm) > 1e-15) {
          yield [i, j, new Complex(vRe, vIm)];
        }
      }
    }
  }
  // Group the Paulis into mutually-commuting subsets (for measurement
  // optimization). Returns an array of SparsePauliOp instances.
  noncommutation_groups() {
    const groups = [];
    for (let i = 0; i < this.paulis.size; i++) {
      const pi = this.paulis.get(i);
      const ci = this.coeffs[i];
      let placed = false;
      for (const g of groups) {
        let commutes = true;
        for (const pj of g.paulis) {
          if (!pi.commutes(pj)) {
            commutes = false;
            break;
          }
        }
        if (commutes) {
          g.paulis.push(pi);
          g.coeffs.push(ci);
          placed = true;
          break;
        }
      }
      if (!placed) {
        groups.push({ paulis: [pi], coeffs: [ci] });
      }
    }
    return groups.map((g) => new _SparsePauliOp(new PauliList(g.paulis), g.coeffs));
  }
  // Sort terms by weight (descending by default) for measurement efficiency.
  sort(reverse = true) {
    const indexed = this.paulis.paulis.map((p, i) => ({ p, c: this.coeffs[i], w: p.weight() }));
    indexed.sort((a, b) => reverse ? b.w - a.w : a.w - b.w);
    return new _SparsePauliOp(
      new PauliList(indexed.map((x) => x.p)),
      indexed.map((x) => x.c)
    );
  }
  // Split into chunks of mutually-commuting terms (alias for noncommutation_groups
  // with a size cap).
  chunk(maxSize = Infinity) {
    const groups = this.noncommutation_groups();
    if (maxSize === Infinity) return groups;
    const result = [];
    for (const g of groups) {
      if (g.paulis.size <= maxSize) {
        result.push(g);
      } else {
        for (let i = 0; i < g.paulis.size; i += maxSize) {
          const slice = g.paulis.paulis.slice(i, i + maxSize);
          const coeffs = g.coeffs.slice(i, i + maxSize);
          result.push(new _SparsePauliOp(new PauliList(slice), coeffs));
        }
      }
    }
    return result;
  }
  conjugate() {
    const newCoeffs = this.coeffs.map((c) => c.conjugate());
    return new _SparsePauliOp(new PauliList(this.paulis.paulis), newCoeffs);
  }
  transpose() {
    const newCoeffs = this.coeffs.map((c, i) => {
      let yCount = 0;
      const p = this.paulis.get(i);
      for (let j = 0; j < p.num_qubits; j++) if (p.x[j] && p.z[j]) yCount++;
      const sign = yCount % 2 === 0 ? 1 : -1;
      return new Complex(c.re * sign, c.im * sign);
    });
    return new _SparsePauliOp(new PauliList(this.paulis.paulis), newCoeffs);
  }
  adjoint() {
    const newCoeffs = this.coeffs.map((c) => c.conjugate());
    const adjustedPaulis = this.paulis.paulis.map((p) => {
      let yCount = 0;
      for (let i = 0; i < p.num_qubits; i++) if (p.x[i] && p.z[i]) yCount++;
      return { pauli: p, sign: yCount % 2 === 0 ? 1 : -1 };
    });
    const finalCoeffs = newCoeffs.map((c, i) => new Complex(c.re * adjustedPaulis[i].sign, c.im * adjustedPaulis[i].sign));
    return new _SparsePauliOp(new PauliList(this.paulis.paulis), finalCoeffs);
  }
  copy() {
    return new _SparsePauliOp(
      new PauliList(this.paulis.paulis.map((p) => new Pauli(p.label))),
      this.coeffs.map((c) => new Complex(c.re, c.im))
    );
  }
  simplify(tol = 1e-12) {
    const map = /* @__PURE__ */ new Map();
    for (let i = 0; i < this.paulis.size; i++) {
      const label = this.paulis.get(i).label;
      if (map.has(label)) {
        const existing = map.get(label);
        existing.coeff = existing.coeff.add(this.coeffs[i]);
      } else {
        map.set(label, { pauli: this.paulis.get(i), coeff: this.coeffs[i] });
      }
    }
    const newPaulis = [];
    const newCoeffs = [];
    for (const v of map.values()) {
      if (v.coeff.abs() > tol) {
        newPaulis.push(v.pauli);
        newCoeffs.push(v.coeff);
      }
    }
    return new _SparsePauliOp(new PauliList(newPaulis), newCoeffs);
  }
  // Apply to a statevector
  apply_to_vector(vector) {
    return this.to_matrix().matvec(vector);
  }
  // Expectation value <psi|H|psi>.
  // Uses a fast O(2^n * m) algorithm that applies each Pauli term directly
  // to the statevector amplitudes without materializing the full 2^n x 2^n
  // matrix. For n=10 qubits and m=100 terms this is ~100x faster than the
  // matrix-based approach.
  expectation_value(statevector) {
    const sv = statevector._data || statevector.data;
    const n = this.num_qubits;
    const dim = 1 << n;
    let totalRe = 0, totalIm = 0;
    for (let ti = 0; ti < this.paulis.size; ti++) {
      const pauli = this.paulis.get(ti);
      const coeff = this.coeffs[ti];
      const x = pauli.x, z = pauli.z;
      let mask = 0, yMask = 0, zOnlyMask = 0, yCount = 0;
      for (let q = 0; q < n; q++) {
        if (x[q]) {
          mask |= 1 << q;
          if (z[q]) {
            yMask |= 1 << q;
            yCount++;
          }
        } else if (z[q]) {
          zOnlyMask |= 1 << q;
        }
      }
      const iPow = [new Complex(1, 0), new Complex(0, 1), new Complex(-1, 0), new Complex(0, -1)][yCount % 4];
      const iPowRe = iPow.re, iPowIm = iPow.im;
      const cr = coeff.re, ci = coeff.im;
      for (let i = 0; i < dim; i++) {
        const j = i ^ mask;
        const zPhase = _popcountParity(i & zOnlyMask);
        const yPhase = _popcountParity(j & yMask);
        const phase = zPhase * yPhase;
        const psiI = sv.data[i];
        const psiJ = sv.data[j];
        const pRe = iPowRe * phase;
        const pIm = iPowIm * phase;
        const cRe = psiI.re * pRe - -psiI.im * pIm;
        const cIm = psiI.re * pIm + -psiI.im * pRe;
        const tRe = cRe * psiJ.re - cIm * psiJ.im;
        const tIm = cRe * psiJ.im + cIm * psiJ.re;
        totalRe += tRe * cr - tIm * ci;
        totalIm += tRe * ci + tIm * cr;
      }
    }
    return new Complex(totalRe, totalIm);
  }
  // Legacy matrix-based expectation value (kept for testing / verification).
  expectation_value_matrix(statevector) {
    const m = this.to_matrix();
    const opResult = m.matvec(statevector._data);
    return statevector._data.inner(opResult);
  }
  toString() {
    return this.to_list().map(([l, c]) => `${c}*${l}`).join(" + ");
  }
};
function _popcountParity(x) {
  x ^= x >> 16;
  x ^= x >> 8;
  x ^= x >> 4;
  x &= 15;
  return 27030 >> x & 1 ? -1 : 1;
}

// src/library/extra_gates.js
var r3 = (re, im = 0) => new Complex(re, im);
var UnitaryGate = class _UnitaryGate extends Gate {
  constructor(matrix, label = null) {
    let mat;
    if (matrix instanceof ComplexMatrix2) {
      mat = matrix;
    } else if (matrix && matrix._data instanceof ComplexMatrix2) {
      mat = matrix._data;
    } else if (Array.isArray(matrix)) {
      mat = ComplexMatrix2.fromRows(matrix.map(
        (row) => row.map((v) => v instanceof Complex ? v : r3(v))
      ));
    } else {
      throw new TypeError("UnitaryGate requires a ComplexMatrix, Operator, or 2D array");
    }
    const n = Math.log2(mat.rows);
    if (!Number.isInteger(n)) {
      throw new Error(`UnitaryGate: matrix dimension ${mat.rows} is not 2^n`);
    }
    super("unitary", n, []);
    this._matrix = mat;
    this.label = label;
    this._matrixBuilder = () => mat;
  }
  copy() {
    const g = new _UnitaryGate(this._matrix, this.label);
    return g;
  }
  // Allow a unitary gate to be parameter-bound (it has no parameters, so
  // this is a no-op, but bind_parameters calls it).
  bind() {
    return this;
  }
  control(numCtrl = 1) {
    return new ControlledGate(this, numCtrl);
  }
};
var DiagonalGate = class _DiagonalGate extends Gate {
  constructor(diag) {
    const n = Math.log2(diag.length);
    if (!Number.isInteger(n)) {
      throw new Error(`DiagonalGate: diag length ${diag.length} is not 2^n`);
    }
    super("diagonal", n, diag);
    this._diag = diag.slice();
    this._matrixBuilder = (params) => {
      const d = params && params.length ? params : this._diag;
      const dim = 1 << n;
      const m = ComplexMatrix2.zeros(dim, dim);
      for (let i = 0; i < dim; i++) {
        m.set(i, i, r3(Math.cos(d[i]), Math.sin(d[i])));
      }
      return m;
    };
  }
  copy() {
    return new _DiagonalGate(this._diag);
  }
};
var PermutationGate = class _PermutationGate extends Gate {
  constructor(pattern) {
    const n = pattern.length;
    super("permutation", n, []);
    this._pattern = pattern.slice();
    this._matrixBuilder = () => {
      const dim = 1 << n;
      const m = ComplexMatrix2.zeros(dim, dim);
      for (let i = 0; i < dim; i++) {
        let j = 0;
        for (let q = 0; q < n; q++) {
          if (i >> q & 1) {
            j |= 1 << pattern.indexOf(q);
          }
        }
        m.set(j, i, r3(1));
      }
      return m;
    };
  }
  copy() {
    return new _PermutationGate(this._pattern);
  }
};
var HamiltonianGate = class _HamiltonianGate extends Gate {
  constructor(operator, time, label = null) {
    let mat;
    if (operator instanceof ComplexMatrix2) {
      mat = operator;
    } else if (operator && operator._data instanceof ComplexMatrix2) {
      mat = operator._data;
    } else if (operator instanceof SparsePauliOp) {
      mat = operator.to_matrix();
    } else if (Array.isArray(operator)) {
      mat = ComplexMatrix2.fromRows(operator.map(
        (row) => row.map((v) => v instanceof Complex ? v : r3(v))
      ));
    } else {
      throw new TypeError("HamiltonianGate requires a matrix, Operator, or SparsePauliOp");
    }
    const n = Math.log2(mat.rows);
    if (!Number.isInteger(n)) {
      throw new Error(`HamiltonianGate: dimension ${mat.rows} is not 2^n`);
    }
    super("hamiltonian", n, [time]);
    this._hamiltonian = mat;
    this._time = time;
    this.label = label;
    this._matrixBuilder = (params) => {
      const t = params && params.length ? params[0] : this._time;
      const scaledH = mat.scale(new Complex(0, -t));
      return scaledH.expm();
    };
  }
  copy() {
    return new _HamiltonianGate(this._hamiltonian, this._time, this.label);
  }
};
var Initialize = class _Initialize extends Instruction {
  constructor(amplitudes, numQubits = null) {
    let amps;
    let n;
    if (typeof amplitudes === "string") {
      const sv = Statevector.fromLabel(amplitudes);
      amps = sv._data.data;
      n = amplitudes.length;
    } else if (Array.isArray(amplitudes)) {
      amps = amplitudes.map((v) => v instanceof Complex ? v : r3(v));
      n = Math.log2(amps.length);
      if (!Number.isInteger(n)) {
        throw new Error(`Initialize: amplitudes length ${amps.length} is not 2^n`);
      }
    } else {
      throw new TypeError("Initialize requires an array of amplitudes or a label string");
    }
    if (numQubits !== null) n = numQubits;
    super("initialize", n, 0, amps);
    this._amplitudes = amps;
    this._matrixBuilder = () => {
      const dim = amps.length;
      const m = ComplexMatrix2.identity(dim);
      for (let i = 0; i < dim; i++) m.set(i, 0, amps[i]);
      for (let j = 1; j < dim; j++) m.set(0, j, r3(0));
      return m;
    };
  }
  copy() {
    return new _Initialize(this._amplitudes);
  }
};
var MCPhaseGate = class _MCPhaseGate extends Gate {
  constructor(lambda, numCtrlQubits) {
    super("mcphase", numCtrlQubits + 1, [lambda]);
    this._lambda = lambda;
    this._numCtrl = numCtrlQubits;
    this._matrixBuilder = (params) => {
      const l = params && params.length ? params[0] : this._lambda;
      const numCtrl = this._numCtrl;
      const dim = 1 << numCtrl + 1;
      const m = ComplexMatrix2.identity(dim);
      m.set(dim - 1, dim - 1, r3(Math.cos(l), Math.sin(l)));
      return m;
    };
  }
  copy() {
    return new _MCPhaseGate(this._lambda, this._numCtrl);
  }
};
var MCRXGate = class _MCRXGate extends Gate {
  constructor(theta, numCtrlQubits) {
    super("mcrx", numCtrlQubits + 1, [theta]);
    this._theta = theta;
    this._numCtrl = numCtrlQubits;
    this._matrixBuilder = (params) => {
      const t = params && params.length ? params[0] : this._theta;
      const numCtrl = this._numCtrl;
      const dim = 1 << numCtrl + 1;
      const m = ComplexMatrix2.identity(dim);
      const ct = Math.cos(t / 2);
      const st = Math.sin(t / 2);
      const idx00 = dim - 2;
      const idx11 = dim - 1;
      m.set(idx00, idx00, r3(ct));
      m.set(idx00, idx11, r3(0, -st));
      m.set(idx11, idx00, r3(0, -st));
      m.set(idx11, idx11, r3(ct));
      return m;
    };
  }
  copy() {
    return new _MCRXGate(this._theta, this._numCtrl);
  }
};
var MCRYGate = class _MCRYGate extends Gate {
  constructor(theta, numCtrlQubits) {
    super("mcry", numCtrlQubits + 1, [theta]);
    this._theta = theta;
    this._numCtrl = numCtrlQubits;
    this._matrixBuilder = (params) => {
      const t = params && params.length ? params[0] : this._theta;
      const numCtrl = this._numCtrl;
      const dim = 1 << numCtrl + 1;
      const m = ComplexMatrix2.identity(dim);
      const ct = Math.cos(t / 2);
      const st = Math.sin(t / 2);
      const idx00 = dim - 2;
      const idx11 = dim - 1;
      m.set(idx00, idx00, r3(ct));
      m.set(idx00, idx11, r3(-st));
      m.set(idx11, idx00, r3(st));
      m.set(idx11, idx11, r3(ct));
      return m;
    };
  }
  copy() {
    return new _MCRYGate(this._theta, this._numCtrl);
  }
};
var MCRZGate = class _MCRZGate extends Gate {
  constructor(theta, numCtrlQubits) {
    super("mcrz", numCtrlQubits + 1, [theta]);
    this._theta = theta;
    this._numCtrl = numCtrlQubits;
    this._matrixBuilder = (params) => {
      const t = params && params.length ? params[0] : this._theta;
      const numCtrl = this._numCtrl;
      const dim = 1 << numCtrl + 1;
      const m = ComplexMatrix2.identity(dim);
      const ct = Math.cos(t / 2);
      const st = Math.sin(t / 2);
      const idx00 = dim - 2;
      const idx11 = dim - 1;
      m.set(idx00, idx00, r3(ct, -st));
      m.set(idx11, idx11, r3(ct, st));
      return m;
    };
  }
  copy() {
    return new _MCRZGate(this._theta, this._numCtrl);
  }
};
var MCMTGate = class _MCMTGate extends Gate {
  constructor(baseGate, numCtrlQubits, numTargetQubits) {
    super("mcmt", numCtrlQubits + numTargetQubits, baseGate.params.slice());
    this._baseGate = baseGate;
    this._numCtrl = numCtrlQubits;
    this._numTarget = numTargetQubits;
    this._matrixBuilder = () => {
      const baseMat = baseGate.to_matrix();
      const baseDim = baseMat.rows;
      if (numTargetQubits === 1 && baseDim === 2) {
        const cg = new ControlledGate(baseGate, numCtrlQubits);
        return cg.to_matrix();
      }
      const totalDim = 1 << numCtrlQubits + numTargetQubits;
      const m = ComplexMatrix2.identity(totalDim);
      const ctrlOn = (1 << numCtrlQubits) - 1;
      const targetBlockDim = 1 << numTargetQubits;
      let targetU = baseMat;
      for (let i = 1; i < numTargetQubits; i++) {
        targetU = targetU.tensor(baseMat);
      }
      for (let bi = 0; bi < targetBlockDim; bi++) {
        for (let bj = 0; bj < targetBlockDim; bj++) {
          const row = ctrlOn + (bi << numCtrlQubits);
          const col = ctrlOn + (bj << numCtrlQubits);
          m.set(row, col, targetU.get(bi, bj));
        }
      }
      return m;
    };
  }
  copy() {
    return new _MCMTGate(this._baseGate.copy(), this._numCtrl, this._numTarget);
  }
};
_registerParamBuilder("UNITARY", (matrix) => new UnitaryGate(matrix));
_registerParamBuilder("DIAGONAL", (diag) => new DiagonalGate(diag));
_registerParamBuilder("PERMUTATION", (pattern) => new PermutationGate(pattern));
_registerParamBuilder("HAMILTONIAN", (op, t) => new HamiltonianGate(op, t));
_registerExtraGateClass("UnitaryGate", UnitaryGate);
_registerExtraGateClass("DiagonalGate", DiagonalGate);
_registerExtraGateClass("PermutationGate", PermutationGate);
_registerExtraGateClass("HamiltonianGate", HamiltonianGate);
_registerExtraGateClass("Initialize", Initialize);
_registerExtraGateClass("MCPhaseGate", MCPhaseGate);
_registerExtraGateClass("MCRXGate", MCRXGate);
_registerExtraGateClass("MCRYGate", MCRYGate);
_registerExtraGateClass("MCRZGate", MCRZGate);
_registerExtraGateClass("MCMTGate", MCMTGate);

// src/library/circuits.js
var circuits_exports = {};
__export(circuits_exports, {
  bellState: () => bellState,
  efficientSU2: () => efficientSU2,
  ghzState: () => ghzState,
  graphState: () => graphState,
  hiddenLinearFunction: () => hiddenLinearFunction,
  iqft: () => iqft,
  pauliEvolution: () => pauliEvolution,
  pauliTwoDesign: () => pauliTwoDesign,
  qft: () => qft,
  qftInverse: () => qftInverse,
  quantumVolume: () => quantumVolume,
  realAmplitudes: () => realAmplitudes,
  twoLocal: () => twoLocal
});
function bellState() {
  const qc = new QuantumCircuit(2);
  qc.h(0);
  qc.cx(0, 1);
  return qc;
}
function ghzState(numQubits) {
  if (numQubits < 2) throw new Error("GHZ requires at least 2 qubits");
  const qc = new QuantumCircuit(numQubits);
  qc.h(0);
  for (let i = 0; i < numQubits - 1; i++) {
    qc.cx(i, i + 1);
  }
  return qc;
}
function qft(numQubits) {
  const qc = new QuantumCircuit(numQubits);
  for (let i = numQubits - 1; i >= 0; i--) {
    qc.h(i);
    for (let j = i - 1; j >= 0; j--) {
      const lam = Math.PI / Math.pow(2, i - j);
      qc.cp(lam, j, i);
    }
  }
  for (let i = 0; i < Math.floor(numQubits / 2); i++) {
    qc.swap(i, numQubits - 1 - i);
  }
  return qc;
}
function qftInverse(numQubits) {
  return qft(numQubits).inverse();
}
function quantumVolume(numQubits, depth = null, seed = null) {
  const d = depth || numQubits;
  const qc = new QuantumCircuit(numQubits);
  const rng = seed != null ? _makeRng(seed) : Math.random;
  for (let layer = 0; layer < d; layer++) {
    const perm = _randomPermutation(numQubits, rng);
    for (let i = 0; i + 1 < numQubits; i += 2) {
      const q1 = perm[i];
      const q2 = perm[i + 1];
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q1);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q2);
      qc.cx(q1, q2);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q1);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q2);
      qc.cx(q1, q2);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q1);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q2);
      qc.cx(q1, q2);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q1);
      qc.u(rng() * 2 * Math.PI, rng() * 2 * Math.PI, rng() * 2 * Math.PI, q2);
    }
    qc.barrier();
  }
  return qc;
}
function realAmplitudes(numQubits, reps = 3, entanglement = "full") {
  const qc = new QuantumCircuit(numQubits);
  for (let q = 0; q < numQubits; q++) {
    qc.ry(new Parameter(`theta[0][${q}]`), q);
  }
  qc.barrier();
  for (let r4 = 1; r4 <= reps; r4++) {
    _applyEntanglement(qc, numQubits, entanglement, "cx");
    qc.barrier();
    for (let q = 0; q < numQubits; q++) {
      qc.ry(new Parameter(`theta[${r4}][${q}]`), q);
    }
    qc.barrier();
  }
  return qc;
}
function efficientSU2(numQubits, reps = 3, entanglement = "full") {
  const qc = new QuantumCircuit(numQubits);
  for (let q = 0; q < numQubits; q++) {
    qc.ry(new Parameter(`theta[0][${q}]`), q);
    qc.rz(new Parameter(`phi[0][${q}]`), q);
  }
  qc.barrier();
  for (let r4 = 1; r4 <= reps; r4++) {
    _applyEntanglement(qc, numQubits, entanglement, "cx");
    qc.barrier();
    for (let q = 0; q < numQubits; q++) {
      qc.ry(new Parameter(`theta[${r4}][${q}]`), q);
      qc.rz(new Parameter(`phi[${r4}][${q}]`), q);
    }
    qc.barrier();
  }
  return qc;
}
function pauliTwoDesign(numQubits, reps = 3, seed = 12345) {
  const qc = new QuantumCircuit(numQubits);
  const rng = _makeRng(seed);
  for (let q = 0; q < numQubits; q++) {
    qc.ry(new Parameter(`p[0][${q}]`), q);
  }
  const paulis = ["x", "y", "z", "xx", "yy", "zz"];
  for (let r4 = 0; r4 < reps; r4++) {
    const layer = paulis[Math.floor(rng() * paulis.length)];
    if (layer === "x") qc.x(r4 % numQubits);
    else if (layer === "y") qc.y(r4 % numQubits);
    else if (layer === "z") qc.z(r4 % numQubits);
    else if (layer === "xx") qc.cx(r4 % numQubits, (r4 + 1) % numQubits);
    else if (layer === "yy") qc.cy(r4 % numQubits, (r4 + 1) % numQubits);
    else if (layer === "zz") qc.cz(r4 % numQubits, (r4 + 1) % numQubits);
    for (let q = 0; q < numQubits; q++) {
      qc.rz(new Parameter(`p[${r4 + 1}][${q}]`), q);
    }
    qc.barrier();
  }
  return qc;
}
function twoLocal(numQubits, rotationBlocks = ["ry"], entanglementBlocks = ["cx"], reps = 3, entanglement = "full") {
  const qc = new QuantumCircuit(numQubits);
  const rotName = (block) => block.toLowerCase();
  const paramFor = (layer, q, block) => {
    if (rotationBlocks.length === 1) {
      return new Parameter(`theta[${layer}][${q}]`);
    }
    return new Parameter(`${rotName(block)}[${layer}][${q}]`);
  };
  const applyRotation = (q, layer, block) => {
    const p = paramFor(layer, q, block);
    if (block === "rx") qc.rx(p, q);
    else if (block === "ry") qc.ry(p, q);
    else if (block === "rz") qc.rz(p, q);
    else throw new Error(`twoLocal: unknown rotation block ${block}`);
  };
  for (let q = 0; q < numQubits; q++) {
    for (const block of rotationBlocks) applyRotation(q, 0, block);
  }
  for (let r4 = 1; r4 <= reps; r4++) {
    for (const block of entanglementBlocks) {
      _applyEntanglement(qc, numQubits, entanglement, block.toLowerCase());
    }
    for (let q = 0; q < numQubits; q++) {
      for (const block of rotationBlocks) applyRotation(q, r4, block);
    }
  }
  return qc;
}
function graphState(numQubits, edges) {
  const qc = new QuantumCircuit(numQubits);
  for (let q = 0; q < numQubits; q++) qc.h(q);
  for (const [i, j] of edges) qc.cz(i, j);
  return qc;
}
function pauliEvolution(pauliLabel, time) {
  const n = pauliLabel.length;
  const qc = new QuantumCircuit(n);
  const chars = pauliLabel.toUpperCase().split("");
  for (let i = 0; i < n; i++) {
    const ch = chars[n - 1 - i];
    if (ch === "X") qc.h(i);
    else if (ch === "Y") {
      qc.sdg(i);
      qc.h(i);
    }
  }
  const nonI = [];
  for (let i = 0; i < n; i++) {
    if (chars[n - 1 - i] !== "I") nonI.push(i);
  }
  for (let k = 0; k + 1 < nonI.length; k++) {
    qc.cx(nonI[k], nonI[k + 1]);
  }
  if (nonI.length > 0) {
    const angle = time && typeof time.mul === "function" ? time.mul(2) : 2 * time;
    qc.rz(angle, nonI[nonI.length - 1]);
  }
  for (let k = nonI.length - 2; k >= 0; k--) {
    qc.cx(nonI[k], nonI[k + 1]);
  }
  for (let i = 0; i < n; i++) {
    const ch = chars[n - 1 - i];
    if (ch === "X") qc.h(i);
    else if (ch === "Y") {
      qc.h(i);
      qc.s(i);
    }
  }
  return qc;
}
function hiddenLinearFunction(numQubits, matrix) {
  const qc = new QuantumCircuit(numQubits);
  for (let q = 0; q < numQubits; q++) qc.h(q);
  for (let i = 0; i < numQubits; i++) {
    for (let j = i + 1; j < numQubits; j++) {
      if (matrix[i][j]) qc.cz(i, j);
    }
  }
  for (let q = 0; q < numQubits; q++) qc.h(q);
  return qc;
}
var iqft = qftInverse;
function _applyEntanglement(qc, n, entanglement, gate) {
  const apply2q = (i, j) => {
    switch (gate) {
      case "cx":
        qc.cx(i, j);
        break;
      case "cz":
        qc.cz(i, j);
        break;
      case "cy":
        qc.cy(i, j);
        break;
      case "ch":
        qc.ch(i, j);
        break;
      case "swap":
        qc.swap(i, j);
        break;
      case "iswap":
        qc.iswap(i, j);
        break;
      case "dcx":
        qc.dcx(i, j);
        break;
      case "csx":
        qc.csx(i, j);
        break;
      default:
        throw new Error(`Unknown entanglement block: ${gate}`);
    }
  };
  if (entanglement === "full") {
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) apply2q(i, j);
    }
  } else if (entanglement === "linear") {
    for (let i = 0; i + 1 < n; i++) apply2q(i, i + 1);
  } else if (entanglement === "circular") {
    for (let i = 0; i + 1 < n; i++) apply2q(i, i + 1);
    if (n > 1) apply2q(n - 1, 0);
  } else if (entanglement === "sca") {
    for (let i = 0; i + 1 < n; i++) apply2q(i, i + 1);
    if (n > 1) apply2q(n - 1, 0);
  } else if (entanglement === "pairwise") {
    for (let i = 0; i + 1 < n; i += 2) apply2q(i, i + 1);
  } else {
    throw new Error(`Unknown entanglement type: ${entanglement}`);
  }
}
function _randomPermutation(n, rng) {
  const arr = Array.from({ length: n }, (_, i) => i);
  for (let i = n - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    const tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
  }
  return arr;
}
function _makeRng(seed) {
  let s = seed >>> 0;
  return function() {
    s |= 0;
    s = s + 1831565813 | 0;
    let t = Math.imul(s ^ s >>> 15, 1 | s);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

// src/quantum_info/density_matrix.js
var DensityMatrix = class _DensityMatrix {
  constructor(data, numQubits = null) {
    if (data instanceof ComplexMatrix2) {
      this._data = data;
    } else if (Array.isArray(data) && Array.isArray(data[0])) {
      this._data = ComplexMatrix2.fromRows(data);
    } else {
      throw new TypeError("DensityMatrix requires a ComplexMatrix or 2D array");
    }
    if (this._data.rows !== this._data.cols) {
      throw new Error("DensityMatrix must be square");
    }
    const nq = Math.log2(this._data.rows);
    if (!Number.isInteger(nq)) {
      throw new Error(`DensityMatrix dimension must be 2^n, got ${this._data.rows}`);
    }
    this._numQubits = numQubits !== null ? numQubits : nq;
  }
  static fromStatevector(statevector) {
    const dim = statevector.dim;
    const data = ComplexMatrix2.zeros(dim, dim);
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        const a = statevector.data.get(i);
        const b = statevector.data.get(j);
        data.set(i, j, a.mul(b.conjugate()));
      }
    }
    return new _DensityMatrix(data, statevector.num_qubits);
  }
  static fromOperator(operator) {
    return new _DensityMatrix(operator.data, operator.num_qubits);
  }
  static fromCircuit(circuit, initState = null) {
    const sv = initState ? initState : Statevector.zero(circuit.num_qubits);
    let dm = _DensityMatrix.fromStatevector(sv);
    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (op.num_qubits === 0) continue;
      if (typeof op.to_matrix !== "function") continue;
      const gateMatrix = op.to_matrix();
      const qubitIndices = ci.qubits.map((q) => circuit._qubit_index.get(q));
      const fullOp = _embedGate2(gateMatrix, circuit.num_qubits, qubitIndices);
      const newRho = fullOp.mul(dm._data).mul(fullOp.dagger());
      dm = new _DensityMatrix(newRho, circuit.num_qubits);
    }
    return dm;
  }
  static fromLabel(label) {
    return _DensityMatrix.fromStatevector(Statevector.fromLabel(label));
  }
  static zero(numQubits) {
    return new _DensityMatrix(ComplexMatrix2.zeros(1 << numQubits), numQubits);
  }
  static identity(numQubits) {
    return new _DensityMatrix(ComplexMatrix2.identity(1 << numQubits).scale(1 / (1 << numQubits)), numQubits);
  }
  get data() {
    return this._data;
  }
  get dim() {
    return this._data.rows;
  }
  get num_qubits() {
    return this._numQubits;
  }
  copy() {
    return new _DensityMatrix(new ComplexMatrix2(this._data.rows, this._data.cols, this._data.data.slice()), this._numQubits);
  }
  // Trace of density matrix (should be 1 for valid quantum states)
  trace() {
    return this._data.trace();
  }
  // Purity: Tr(rho^2). Pure states have purity 1, mixed < 1.
  purity() {
    const sq = this._data.mul(this._data);
    const tr = sq.trace();
    return tr.re;
  }
  is_valid(tol = 1e-9) {
    if (!this._data.isHermitian(tol)) return false;
    const tr = this.trace();
    if (Math.abs(tr.re - 1) > tol || Math.abs(tr.im) > tol) return false;
    try {
      const { eigenvalues } = this._data.eigh();
      return eigenvalues.every((ev) => ev >= -tol);
    } catch (e) {
      return false;
    }
  }
  is_pure(tol = 1e-9) {
    return Math.abs(this.purity() - 1) < tol;
  }
  // Evolve under a unitary: rho' = U rho U^dagger
  evolve(other) {
    if (other && other._data && other._data instanceof ComplexMatrix2) {
      if (other.num_qubits !== this._numQubits) {
        throw new Error("evolve: operator qubit count mismatch");
      }
      return new _DensityMatrix(other._data.mul(this._data).mul(other._data.dagger()), this._numQubits);
    }
    if (other && typeof other.to_matrix === "function") {
      if (other.num_qubits !== this._numQubits) {
        throw new Error("evolve: gate qubit count mismatch");
      }
      const m = other.to_matrix();
      return new _DensityMatrix(m.mul(this._data).mul(m.dagger()), this._numQubits);
    }
    throw new TypeError("evolve: unsupported operand type");
  }
  // Tensor product with another DensityMatrix
  tensor(other) {
    return new _DensityMatrix(this._data.tensor(other._data), this._numQubits + other._numQubits);
  }
  expand(other) {
    return new _DensityMatrix(other._data.tensor(this._data), this._numQubits + other._numQubits);
  }
  // Add another density matrix (for mixture)
  add(other) {
    if (this._numQubits !== other._numQubits) {
      throw new Error("add: qubit count mismatch");
    }
    return new _DensityMatrix(this._data.add(other._data), this._numQubits);
  }
  scale(c) {
    return new _DensityMatrix(this._data.scale(c), this._numQubits);
  }
  // Partial trace over a qubit
  partial_trace(qubits) {
    const qList = Array.isArray(qubits) ? qubits : [qubits];
    let result = this;
    const sorted = qList.slice().sort((a, b) => b - a);
    let currentQubits = this._numQubits;
    let currentData = this._data;
    for (const q of sorted) {
      currentData = currentData.partialTrace(currentQubits, q);
      currentQubits--;
    }
    return new _DensityMatrix(currentData, currentQubits);
  }
  // Probabilities of measuring each computational basis state
  probabilities(qargs = null) {
    const diag = new Array(this._data.rows);
    for (let i = 0; i < this._data.rows; i++) {
      diag[i] = this._data.get(i, i).re;
    }
    if (qargs == null) return diag;
    const indices = Array.isArray(qargs) ? qargs : [qargs];
    return _marginalProbabilitiesDM(diag, this._numQubits, indices);
  }
  // Sample measurement outcomes
  sample(counts = 1024, qargs = null, rng) {
    const probs = this.probabilities(qargs);
    const results = {};
    const nq = qargs == null ? this._numQubits : Array.isArray(qargs) ? qargs.length : 1;
    for (let i = 0; i < counts; i++) {
      const idx = sampleDistribution(probs, rng);
      const bits = new Array(nq).fill("0");
      for (let q = 0; q < nq; q++) bits[q] = (idx >> q & 1).toString();
      const key = bits.reverse().join("");
      results[key] = (results[key] || 0) + 1;
    }
    return results;
  }
  // Expectation value <O> = Tr(O rho)
  expectation_value(operator) {
    const m = operator instanceof ComplexMatrix2 ? operator : operator._data;
    const prod = m.mul(this._data);
    const tr = prod.trace();
    return tr;
  }
  // Convert to a Statevector if pure (returns null if not pure)
  to_statevector() {
    if (!this.is_pure()) return null;
    const { eigenvalues, eigenvectors } = this._data.eigh();
    let maxIdx = 0;
    let maxVal = -Infinity;
    for (let i = 0; i < eigenvalues.length; i++) {
      if (eigenvalues[i] > maxVal) {
        maxVal = eigenvalues[i];
        maxIdx = i;
      }
    }
    const dim = this._data.rows;
    const data = new Array(dim);
    for (let i = 0; i < dim; i++) {
      data[i] = eigenvectors.get(i, maxIdx);
    }
    return new Statevector(new ComplexVector(data), this._numQubits);
  }
  equals(other, tol) {
    return other instanceof _DensityMatrix && this._data.equals(other._data, tol);
  }
  to_matrix() {
    return this._data;
  }
  to_dict() {
    const rows = this._data.toRows();
    return rows.map((r4) => r4.map((c) => ({ re: c.re, im: c.im })));
  }
  toString() {
    return `DensityMatrix(num_qubits=${this._numQubits})
${this._data.toString()}`;
  }
};
function _embedGate2(gateMatrix, numQubits, qubitIndices) {
  const dim = 1 << numQubits;
  const out = ComplexMatrix2.zeros(dim, dim);
  const k = qubitIndices.length;
  for (let row = 0; row < dim; row++) {
    for (let col = 0; col < dim; col++) {
      let subRow = 0, subCol = 0, restMatches = true;
      for (let qi = 0; qi < k; qi++) {
        const rb = row >> qubitIndices[qi] & 1;
        const cb = col >> qubitIndices[qi] & 1;
        subRow |= rb << qi;
        subCol |= cb << qi;
      }
      for (let q = 0; q < numQubits; q++) {
        if (qubitIndices.indexOf(q) === -1) {
          const rb = row >> q & 1;
          const cb = col >> q & 1;
          if (rb !== cb) {
            restMatches = false;
            break;
          }
        }
      }
      if (!restMatches) continue;
      out.set(row, col, gateMatrix.get(subRow, subCol));
    }
  }
  return out;
}
function _marginalProbabilitiesDM(diag, numQubits, indices) {
  const outDim = 1 << indices.length;
  const out = new Array(outDim).fill(0);
  for (let i = 0; i < diag.length; i++) {
    let subIdx = 0;
    for (let qi = 0; qi < indices.length; qi++) {
      const bit = i >> indices[qi] & 1;
      subIdx |= bit << qi;
    }
    out[subIdx] += diag[i];
  }
  return out;
}
_registerDensityMatrixClass(DensityMatrix);

// src/quantum_info/clifford.js
var Clifford = class _Clifford {
  constructor(data) {
    if (data && data.x && data.z) {
      this.x = data.x.map((row) => row.slice());
      this.z = data.z.map((row) => row.slice());
      this.signs = (data.signs || new Array(data.x.length).fill(0)).slice();
      this.num_qubits = data.x.length / 2;
      if (!Number.isInteger(this.num_qubits)) {
        throw new Error(
          `Clifford: tableau must have 2n rows (got ${data.x.length}); num_qubits would be non-integer.`
        );
      }
    } else if (typeof data === "string") {
      this.num_qubits = data.length;
      this.x = [];
      this.z = [];
      this.signs = [];
      for (let i = 0; i < this.num_qubits; i++) {
        const xRow = new Array(this.num_qubits).fill(0);
        const zRow = new Array(this.num_qubits).fill(0);
        xRow[i] = 1;
        this.x.push(xRow);
        this.z.push(zRow);
        this.signs.push(0);
      }
      for (let i = 0; i < this.num_qubits; i++) {
        const xRow = new Array(this.num_qubits).fill(0);
        const zRow = new Array(this.num_qubits).fill(0);
        const bit = parseInt(data[data.length - 1 - i], 10);
        if (bit === 1) {
          zRow[i] = 1;
          this.signs.push(1);
        } else {
          zRow[i] = 1;
          this.signs.push(0);
        }
        this.x.push(xRow);
        this.z.push(zRow);
      }
    } else if (Array.isArray(data) && data.length > 0 && typeof data[0] === "string") {
      const paulis = data.map((l) => new Pauli(l));
      this.num_qubits = paulis[0].num_qubits;
      this.x = [];
      this.z = [];
      this.signs = [];
      for (let i = 0; i < this.num_qubits; i++) {
        const xRow = new Array(this.num_qubits).fill(0);
        const zRow = new Array(this.num_qubits).fill(0);
        xRow[i] = 1;
        this.x.push(xRow);
        this.z.push(zRow);
        this.signs.push(0);
      }
      for (const p of paulis) {
        this.x.push(p.x.slice());
        this.z.push(p.z.slice());
        this.signs.push(0);
      }
    } else {
      throw new TypeError("Clifford: invalid constructor argument");
    }
  }
  static fromCircuit(circuit) {
    const n = circuit.num_qubits;
    const cliff = _Clifford.fromLabel("0".repeat(n));
    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (op.num_qubits === 0) continue;
      const targets = ci.qubits.map((q) => circuit._qubit_index.get(q));
      cliff._applyGate(op.name, targets, op.params || []);
    }
    return cliff;
  }
  static fromLabel(label) {
    return new _Clifford(label);
  }
  // Build a Clifford from a unitary matrix. The matrix must actually be a
  // Clifford unitary (i.e. map Paulis to Paulis); otherwise an error is
  // thrown. The algorithm reads off the destabilizer and stabilizer
  // generators directly by computing U X_i U^† and U Z_i U^† for each
  // qubit i and identifying the resulting Paulis.
  static fromMatrix(matrix) {
    const n = Math.log2(matrix.rows);
    if (!Number.isInteger(n)) {
      throw new Error(`Clifford.fromMatrix: matrix dimension ${matrix.rows} is not 2^n`);
    }
    const dim = 1 << n;
    const matrixDag = matrix.dagger();
    const paulis = [
      { name: "X", mat: PAULI.X },
      { name: "Z", mat: PAULI.Z }
    ];
    const destabs = [];
    const stabs = [];
    for (let q = 0; q < n; q++) {
      for (const { name, mat } of paulis) {
        const fullP = _embedPauli(mat, n, q);
        const conj = matrix.mul(fullP).mul(matrixDag);
        const result = _identifyPauli(conj, n);
        if (result === null) {
          throw new Error(
            `Clifford.fromMatrix: U ${name}_${q} U^\u2020 is not a Pauli (matrix is not a Clifford)`
          );
        }
        if (name === "X") destabs.push(result);
        else stabs.push(result);
      }
    }
    const x = destabs.map((d) => d.x).concat(stabs.map((s) => s.x));
    const z = destabs.map((d) => d.z).concat(stabs.map((s) => s.z));
    const signs = destabs.map((d) => d.sign).concat(stabs.map((s) => s.sign));
    return new _Clifford({ x, z, signs });
  }
  static random(numQubits, seed = null) {
    let rng = Math.random;
    if (seed != null) {
      let s = seed >>> 0;
      rng = () => {
        s = s + 1831565813 | 0;
        let t = Math.imul(s ^ s >>> 15, 1 | s);
        t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
        return ((t ^ t >>> 14) >>> 0) / 4294967296;
      };
    }
    const c = _Clifford.fromLabel("0".repeat(numQubits));
    for (let q = 0; q < numQubits; q++) {
      const r4 = rng();
      if (r4 < 0.25) {
      } else if (r4 < 0.5) c._h(q);
      else if (r4 < 0.75) {
        c._h(q);
        c._s(q);
      } else c._s(q);
    }
    for (let i = 0; i < numQubits; i++) {
      for (let j = i + 1; j < numQubits; j++) {
        if (rng() < 0.5) c._cx(i, j);
      }
    }
    return c;
  }
  // Apply gates directly to the stabilizer table.
  // Decompositions of X, Y, Z into H and S are written so that the
  // conjugation action matches the standard Pauli definitions exactly
  // (the previous implementation used sequences that produced Z when X
  // was requested, etc.). Verified against the canonical Clifford-action
  // table for single-qubit Paulis.
  _applyGate(name, targets, params) {
    switch (name) {
      case "h":
      case "H":
        this._h(targets[0]);
        break;
      case "s":
      case "S":
        this._s(targets[0]);
        break;
      case "sdg":
      case "SDG":
        this._s(targets[0]);
        this._s(targets[0]);
        this._s(targets[0]);
        break;
      // X = H · Z · H  (H swaps X and Z, so conjugating Z by H gives X).
      case "x":
      case "X":
        this._h(targets[0]);
        this._s(targets[0]);
        this._s(targets[0]);
        this._h(targets[0]);
        break;
      // Z = Z (identity on the tableau for Z itself is just applying Z,
      // which is S^2 — but to keep the convention consistent we apply
      // S twice, which is exactly the Z action in the H/S basis).
      case "z":
      case "Z":
        this._s(targets[0]);
        this._s(targets[0]);
        break;
      // Y = i·X·Z = X·Z·(phase). On the tableau, Y acts as (X·Z) up to a
      // phase, which we can realize as Z then X (or equivalently S^2 then
      // HSH). Using X·Z means: apply Z first (S^2), then X (H S^2 H).
      case "y":
      case "Y":
        this._s(targets[0]);
        this._s(targets[0]);
        this._h(targets[0]);
        this._s(targets[0]);
        this._s(targets[0]);
        this._h(targets[0]);
        break;
      case "cx":
      case "CX":
        this._cx(targets[0], targets[1]);
        break;
      case "cz":
      case "CZ":
        this._h(targets[1]);
        this._cx(targets[0], targets[1]);
        this._h(targets[1]);
        break;
      // CY = (I⊗Sdg) · CX · (I⊗S). Decompose Sdg as S^3.
      case "cy":
      case "CY":
        this._s(targets[1]);
        this._s(targets[1]);
        this._s(targets[1]);
        this._cx(targets[0], targets[1]);
        this._s(targets[1]);
        break;
      case "swap":
      case "SWAP":
        this._cx(targets[0], targets[1]);
        this._cx(targets[1], targets[0]);
        this._cx(targets[0], targets[1]);
        break;
      case "sx":
      case "SX":
        this._h(targets[0]);
        this._s(targets[0]);
        this._h(targets[0]);
        break;
      case "sxdg":
      case "SXDG":
        this._h(targets[0]);
        this._s(targets[0]);
        this._s(targets[0]);
        this._s(targets[0]);
        this._h(targets[0]);
        break;
      case "id":
      case "I":
        break;
      default:
        throw new Error(`Clifford: cannot apply gate ${name} (not Clifford)`);
    }
  }
  // Apply H to qubit q
  _h(q) {
    for (let i = 0; i < this.num_qubits * 2; i++) {
      const xBit = this.x[i][q];
      const zBit = this.z[i][q];
      this.x[i][q] = zBit;
      this.z[i][q] = xBit;
      if (xBit === 1 && zBit === 1) {
        this.signs[i] ^= 1;
      }
    }
  }
  // Apply S to qubit q
  _s(q) {
    for (let i = 0; i < this.num_qubits * 2; i++) {
      const xBit = this.x[i][q];
      const zBit = this.z[i][q];
      this.z[i][q] = xBit ^ zBit;
      if (xBit === 1 && zBit === 1) {
        this.signs[i] ^= 1;
      }
    }
  }
  // Apply CNOT(control, target)
  _cx(control, target) {
    for (let i = 0; i < this.num_qubits * 2; i++) {
      const xc = this.x[i][control];
      const zc = this.z[i][control];
      const xt = this.x[i][target];
      const zt = this.z[i][target];
      if (xc === 1 && zc === 0 && xt === 0 && zt === 1) {
        this.signs[i] ^= 1;
      }
      this.x[i][target] = xc ^ xt;
      this.z[i][control] = zc ^ zt;
    }
  }
  // Build the unitary matrix representing this Clifford.
  //
  // Strategy: U|i> = U X^i |0> = (U X^i U^†) (U |0>) = P_i · (U|0>), where
  // P_i = prod_k (destab[k])^{i_k} is a Pauli determined by the destabilizer
  // rows of the tableau. We compute U|0> once (via StabilizerState), then
  // apply each P_i to get the i-th column. This keeps a consistent global
  // phase across columns, so the resulting matrix is unitary (not just
  // unitary up to per-column phases).
  //
  // Cost: O(2^n · (n + 2^n)) — exponential, but unavoidable since the
  // output matrix has 4^n entries. Practical for n ≤ 8 or so.
  to_matrix() {
    const n = this.num_qubits;
    const dim = 1 << n;
    const uZero = new StabilizerState(this).to_statevector();
    const uZeroVec = uZero._data.data;
    const destabs = [];
    for (let k = 0; k < n; k++) {
      destabs.push({
        x: this.x[k].slice(),
        z: this.z[k].slice(),
        sign: this.signs[k]
      });
    }
    const m = ComplexMatrix2.zeros(dim, dim);
    for (let i = 0; i < dim; i++) {
      let pX = new Array(n).fill(0);
      let pZ = new Array(n).fill(0);
      let pPhase = new Complex(1, 0);
      for (let k = 0; k < n; k++) {
        if (i >> k & 1) {
          const d = destabs[k];
          let dFullPhase = d.sign ? new Complex(-1, 0) : new Complex(1, 0);
          let dIPhase = 0;
          for (let q = 0; q < n; q++) {
            if (d.x[q] && d.z[q]) dIPhase = (dIPhase + 1) % 4;
          }
          const iPhaseFactor = [
            new Complex(1, 0),
            new Complex(0, 1),
            new Complex(-1, 0),
            new Complex(0, -1)
          ][dIPhase];
          dFullPhase = dFullPhase.mul(iPhaseFactor);
          let antiSign = 0;
          for (let q = 0; q < n; q++) {
            if (pZ[q] && d.x[q]) antiSign ^= 1;
          }
          if (antiSign) pPhase = pPhase.mul(new Complex(-1, 0));
          pPhase = pPhase.mul(dFullPhase);
          for (let q = 0; q < n; q++) {
            pX[q] ^= d.x[q];
            pZ[q] ^= d.z[q];
          }
        }
      }
      for (let j = 0; j < dim; j++) {
        const srcIdx = j ^ _bitsToNum(pX, n);
        let parity = 0;
        for (let q = 0; q < n; q++) {
          if (pZ[q] && srcIdx >> q & 1) parity ^= 1;
        }
        const sign = parity ? -1 : 1;
        const amp = uZeroVec[srcIdx];
        const re = sign * pPhase.re * amp.re - sign * pPhase.im * amp.im;
        const im = sign * pPhase.re * amp.im + sign * pPhase.im * amp.re;
        m.set(j, i, new Complex(re, im));
      }
    }
    return m;
  }
  // Compose with another Clifford: returns self * other (apply other first,
  // then self), matching qiskit's Clifford.compose convention.
  //
  // The resulting tableau is computed by conjugating each of other's
  // destabilizers and stabilizers by self's Clifford action. Conjugation
  // preserves the Pauli group structure, so the result is again a valid
  // Clifford tableau.
  compose(other) {
    if (this.num_qubits !== other.num_qubits) {
      throw new Error("compose: qubit count mismatch");
    }
    const n = this.num_qubits;
    const newX = [];
    const newZ = [];
    const newSigns = [];
    for (let k = 0; k < n; k++) {
      const r4 = this._conjugatePauli(other.x[k], other.z[k], other.signs[k]);
      newX.push(r4.x);
      newZ.push(r4.z);
      newSigns.push(r4.sign);
    }
    for (let k = 0; k < n; k++) {
      const r4 = this._conjugatePauli(other.x[n + k], other.z[n + k], other.signs[n + k]);
      newX.push(r4.x);
      newZ.push(r4.z);
      newSigns.push(r4.sign);
    }
    return new _Clifford({
      x: newX,
      z: newZ,
      signs: newSigns
    });
  }
  // Compute U * P * U^† where U is this Clifford and P is the Pauli
  // (-1)^inSign * P(inX, inZ) (using the convention P(x,z) = i^(x·z) X^x Z^z).
  //
  // Returns { x, z, sign } where the result is (-1)^sign * P(x, z).
  // Conjugation by a Clifford maps Paulis to Paulis (up to phase ±1, never
  // ±i), so the result is always expressible in the standard tableau form.
  _conjugatePauli(inX, inZ, inSign) {
    const n = this.num_qubits;
    let currentPauli = new Pauli("I".repeat(n));
    let currentPhase = new Complex(1, 0);
    if (inSign) currentPhase = currentPhase.mul(new Complex(-1, 0));
    let inIPhase = 0;
    for (let q = 0; q < n; q++) {
      if (inX[q] && inZ[q]) inIPhase = (inIPhase + 1) % 4;
    }
    const iPhaseFactor = [
      new Complex(1, 0),
      new Complex(0, 1),
      new Complex(-1, 0),
      new Complex(0, -1)
    ][inIPhase];
    currentPhase = currentPhase.mul(iPhaseFactor);
    for (let k = 0; k < n; k++) {
      if (inX[k]) {
        const destab = Pauli.fromSymplectic(this.x[k], this.z[k]);
        const r4 = currentPauli.compose(destab);
        currentPauli = r4.pauli;
        currentPhase = currentPhase.mul(r4.phase);
        if (this.signs[k]) currentPhase = currentPhase.mul(new Complex(-1, 0));
      }
    }
    for (let k = 0; k < n; k++) {
      if (inZ[k]) {
        const stab = Pauli.fromSymplectic(this.x[n + k], this.z[n + k]);
        const r4 = currentPauli.compose(stab);
        currentPauli = r4.pauli;
        currentPhase = currentPhase.mul(r4.phase);
        if (this.signs[n + k]) currentPhase = currentPhase.mul(new Complex(-1, 0));
      }
    }
    if (Math.abs(currentPhase.im) > 1e-9) {
      throw new Error(
        `Clifford._conjugatePauli: result has imaginary phase ${currentPhase}; expected \xB11 (input was a Pauli from another Clifford tableau).`
      );
    }
    const sign = currentPhase.re < 0 ? 1 : 0;
    return { x: currentPauli.x, z: currentPauli.z, sign };
  }
  to_dict() {
    return {
      stabilizer: this.x.slice(0, this.num_qubits).map((xRow, i) => ({
        x: xRow.slice(),
        z: this.z[i].slice(),
        sign: this.signs[i]
      })),
      destabilizer: this.x.slice(this.num_qubits).map((xRow, i) => ({
        x: xRow.slice(),
        z: this.z[this.num_qubits + i].slice(),
        sign: this.signs[this.num_qubits + i]
      }))
    };
  }
  // Adjoint (inverse) of this Clifford: returns the Clifford U^† such that
  // U · U^† = I. Computed by inverting the tableau's symplectic part and
  // recomputing the signs.
  adjoint() {
    const n = this.num_qubits;
    const m = this.to_matrix();
    const mDag = m.dagger();
    return _Clifford.fromMatrix(mDag);
  }
  // Alias for adjoint.
  inverse() {
    return this.adjoint();
  }
  // Compose this Clifford with itself n times. For n=0 returns identity.
  power(n) {
    if (n === 0) return _Clifford.fromLabel("0".repeat(this.num_qubits));
    if (n < 0) return this.adjoint().power(-n);
    let result = this;
    for (let i = 1; i < n; i++) result = result.compose(this);
    return result;
  }
  // Append another Clifford's action to this one (in-place semantic, but
  // returns a new Clifford). Equivalent to self.compose(other).
  append(other) {
    return this.compose(other);
  }
  equals(other) {
    if (!(other instanceof _Clifford)) return false;
    if (this.num_qubits !== other.num_qubits) return false;
    for (let i = 0; i < this.x.length; i++) {
      for (let q = 0; q < this.num_qubits; q++) {
        if (this.x[i][q] !== other.x[i][q]) return false;
        if (this.z[i][q] !== other.z[i][q]) return false;
      }
      if (this.signs[i] !== other.signs[i]) return false;
    }
    return true;
  }
  toString() {
    return `Clifford(num_qubits=${this.num_qubits})`;
  }
};
Clifford.fromStabilizers = function(labels) {
  const paulis = labels.map((l) => new Pauli(l));
  const n = paulis[0].num_qubits;
  const x = [];
  const z = [];
  const signs = [];
  for (let i = 0; i < n; i++) {
    const xRow = new Array(n).fill(0);
    const zRow = new Array(n).fill(0);
    xRow[i] = 1;
    x.push(xRow);
    z.push(zRow);
    signs.push(0);
  }
  for (const p of paulis) {
    x.push(p.x.slice());
    z.push(p.z.slice());
    signs.push(0);
  }
  return new Clifford({ x, z, signs });
};
var StabilizerState = class _StabilizerState {
  constructor(clifford) {
    if (clifford instanceof Clifford) {
      this.clifford = clifford;
    } else {
      this.clifford = new Clifford(clifford);
    }
    this.num_qubits = this.clifford.num_qubits;
  }
  static fromLabel(label) {
    return new _StabilizerState(Clifford.fromLabel(label));
  }
  static zero(numQubits) {
    return _StabilizerState.fromLabel("0".repeat(numQubits));
  }
  // Probabilities of measuring each basis state
  probabilities(qargs = null) {
    const sv = this.to_statevector();
    if (!sv) return [];
    return sv.probabilities(qargs);
  }
  // Convert to a statevector.
  //
  // Algorithm: a stabilizer state |ψ> is the +1 eigenstate of its n
  // stabilizer generators. Applying the projector Π = ∏_k (I + S_k) / 2 to
  // any reference state |φ> with ⟨ψ|φ⟩ ≠ 0 gives a vector proportional
  // to |ψ>:
  //
  //   Π |φ⟩ = (1/2ⁿ) · Σ_{P ∈ stab group} P |φ⟩ ∝ |ψ⟩
  //
  // We try computational basis states |0⟩, |1⟩, … in turn until we find
  // one with non-zero overlap. (For a generic stabilizer state, |0⟩ works
  // most of the time, but states like |1⟩, |−⟩, etc. have ⟨ψ|0⟩ = 0.)
  //
  // For P = fullPhase · X^x · Z^z and reference |i*⟩:
  //   P |i*⟩ = fullPhase · (-1)^{popcount(i* AND z)} · |i* XOR x⟩
  // so the amplitude of |j⟩ in Π|i*⟩ is:
  //   amp_j = (1/2ⁿ) · Σ_{P with x_P = j XOR i*} fullPhase_P · (-1)^{popcount(i* AND z_P)}
  //
  // Cost: O(4ⁿ) — fine for n ≤ 10.
  to_statevector() {
    const n = this.num_qubits;
    const dim = 1 << n;
    const groupElements = [];
    for (let subset = 0; subset < dim; subset++) {
      let x = new Array(n).fill(0);
      let z = new Array(n).fill(0);
      let fullPhase = new Complex(1, 0);
      for (let k = 0; k < n; k++) {
        if (subset >> k & 1) {
          const gx = this.clifford.x[n + k];
          const gz = this.clifford.z[n + k];
          const gSign = this.clifford.signs[n + k];
          let genFullPhase = gSign ? new Complex(-1, 0) : new Complex(1, 0);
          let genIPhase = 0;
          for (let q = 0; q < n; q++) {
            if (gx[q] && gz[q]) genIPhase = (genIPhase + 1) % 4;
          }
          const iPhaseFactor = [
            new Complex(1, 0),
            new Complex(0, 1),
            new Complex(-1, 0),
            new Complex(0, -1)
          ][genIPhase];
          genFullPhase = genFullPhase.mul(iPhaseFactor);
          let antiSign = 0;
          for (let q = 0; q < n; q++) {
            if (z[q] && gx[q]) antiSign ^= 1;
          }
          if (antiSign) fullPhase = fullPhase.mul(new Complex(-1, 0));
          fullPhase = fullPhase.mul(genFullPhase);
          for (let q = 0; q < n; q++) {
            x[q] ^= gx[q];
            z[q] ^= gz[q];
          }
        }
      }
      groupElements.push({ x, z, fullPhase });
    }
    let refState = -1;
    let refNormSq = 0;
    for (let iStar = 0; iStar < dim; iStar++) {
      let sumRe = 0;
      for (const g of groupElements) {
        let xZero = true;
        for (let q = 0; q < n; q++) {
          if (g.x[q]) {
            xZero = false;
            break;
          }
        }
        if (!xZero) continue;
        let parity = 0;
        for (let q = 0; q < n; q++) {
          if (g.z[q] && iStar >> q & 1) parity ^= 1;
        }
        sumRe += parity ? -g.fullPhase.re : g.fullPhase.re;
      }
      if (sumRe > 1e-9) {
        refState = iStar;
        refNormSq = sumRe / dim;
        break;
      }
    }
    if (refState < 0) {
      throw new Error(
        "StabilizerState.to_statevector: could not find a reference state with non-zero overlap (invalid stabilizer state)"
      );
    }
    const refNorm = Math.sqrt(refNormSq);
    const amps = new Array(dim);
    for (let j = 0; j < dim; j++) {
      const targetX = j ^ refState;
      let sumRe = 0, sumIm = 0;
      for (const g of groupElements) {
        let xMatches = true;
        for (let q = 0; q < n; q++) {
          if ((targetX >> q & 1) !== g.x[q]) {
            xMatches = false;
            break;
          }
        }
        if (!xMatches) continue;
        let parity = 0;
        for (let q = 0; q < n; q++) {
          if (g.z[q] && refState >> q & 1) parity ^= 1;
        }
        const sign = parity ? -1 : 1;
        sumRe += sign * g.fullPhase.re;
        sumIm += sign * g.fullPhase.im;
      }
      amps[j] = new Complex(sumRe / (dim * refNorm), sumIm / (dim * refNorm));
    }
    return new Statevector(new ComplexVector(amps), n);
  }
  // Expectation value ⟨ψ|P|ψ⟩ of a Pauli operator P on this stabilizer state.
  //
  // For a stabilizer state, the expectation is:
  //   - 0 if P is NOT in the stabilizer group (which includes both Paulis
  //     that anticommute with some generator AND Paulis that commute with
  //     every generator but aren't group elements),
  //   - ±1 if P is in the stabilizer group (the sign is the eigenvalue).
  //
  // The previous implementation returned +1 for every Pauli that merely
  // commuted with all stabilizer generators, which gave wrong answers for
  // states like |1⟩ (where ⟨1|Z|1⟩ = −1, not +1) and for Bell states
  // (where ⟨Bell|ZI|Bell⟩ = 0 since ZI commutes with ZZ and XX but is not
  // itself a group element).
  //
  // We brute-force enumerate the 2^n products of the n stabilizer
  // generators and check if any product equals P (up to phase). The
  // accumulated phase of the matching product is the eigenvalue.
  // Cost: O(2^n · n) — fine for n ≤ 16 or so.
  expectation_value(pauli) {
    if (!(pauli instanceof Pauli)) pauli = new Pauli(pauli);
    const n = this.num_qubits;
    if (pauli.label === "I".repeat(n)) {
      return new Complex(1, 0);
    }
    const dim = 1 << n;
    const targetX = pauli.x;
    const targetZ = pauli.z;
    for (let subset = 0; subset < dim; subset++) {
      let x = new Array(n).fill(0);
      let z = new Array(n).fill(0);
      let fullPhase = new Complex(1, 0);
      for (let k = 0; k < n; k++) {
        if (subset >> k & 1) {
          const gx = this.clifford.x[n + k];
          const gz = this.clifford.z[n + k];
          const gSign = this.clifford.signs[n + k];
          let genFullPhase = gSign ? new Complex(-1, 0) : new Complex(1, 0);
          let genIPhase = 0;
          for (let q = 0; q < n; q++) {
            if (gx[q] && gz[q]) genIPhase = (genIPhase + 1) % 4;
          }
          const iPhaseFactor = [
            new Complex(1, 0),
            new Complex(0, 1),
            new Complex(-1, 0),
            new Complex(0, -1)
          ][genIPhase];
          genFullPhase = genFullPhase.mul(iPhaseFactor);
          let antiSign = 0;
          for (let q = 0; q < n; q++) {
            if (z[q] && gx[q]) antiSign ^= 1;
          }
          if (antiSign) fullPhase = fullPhase.mul(new Complex(-1, 0));
          fullPhase = fullPhase.mul(genFullPhase);
          for (let q = 0; q < n; q++) {
            x[q] ^= gx[q];
            z[q] ^= gz[q];
          }
        }
      }
      let matches = true;
      for (let q = 0; q < n; q++) {
        if (x[q] !== targetX[q] || z[q] !== targetZ[q]) {
          matches = false;
          break;
        }
      }
      if (matches) {
        let targetIPhase = 0;
        for (let q = 0; q < n; q++) {
          if (targetX[q] && targetZ[q]) targetIPhase = (targetIPhase + 1) % 4;
        }
        const invIPhase = (4 - targetIPhase) % 4;
        const invIPhaseFactor = [
          new Complex(1, 0),
          new Complex(0, -1),
          new Complex(-1, 0),
          new Complex(0, 1)
        ][invIPhase];
        return fullPhase.mul(invIPhaseFactor);
      }
    }
    return Complex.ZERO;
  }
  equals(other, tol) {
    if (!(other instanceof _StabilizerState)) return false;
    if (this.num_qubits !== other.num_qubits) return false;
    for (let i = 0; i < this.clifford.x.length; i++) {
      for (let q = 0; q < this.num_qubits; q++) {
        if (this.clifford.x[i][q] !== other.clifford.x[i][q]) return false;
        if (this.clifford.z[i][q] !== other.clifford.z[i][q]) return false;
      }
      if (this.clifford.signs[i] !== other.clifford.signs[i]) return false;
    }
    return true;
  }
};
function _bitsToNum(bits, n) {
  let num = 0;
  for (let q = 0; q < n; q++) {
    if (bits[q]) num |= 1 << q;
  }
  return num;
}
function _embedPauli(pauli2x2, n, q) {
  const dim = 1 << n;
  const out = ComplexMatrix2.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      const ib = i >> q & 1;
      const jb = j >> q & 1;
      let ok = true;
      for (let k = 0; k < n; k++) {
        if (k === q) continue;
        if ((i >> k & 1) !== (j >> k & 1)) {
          ok = false;
          break;
        }
      }
      if (ok) out.set(i, j, pauli2x2.get(ib, jb));
    }
  }
  return out;
}
function _identifyPauli(m, n) {
  const dim = 1 << n;
  let firstNonZero = null;
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      if (m.get(i, j).abs() > 1e-9) {
        firstNonZero = { i, j, val: m.get(i, j) };
        break;
      }
    }
    if (firstNonZero) break;
  }
  if (!firstNonZero) return null;
  const phase = firstNonZero.val;
  const mag = phase.abs();
  if (mag < 1e-9) return null;
  const x = new Array(n).fill(0);
  const z = new Array(n).fill(0);
  for (let q = 0; q < n; q++) {
    let hasFlip = false;
    let hasSame = false;
    for (let i = 0; i < dim; i++) {
      const j = i ^ 1 << q;
      const v = m.get(i, j);
      if (v.abs() > 1e-9) hasFlip = true;
      const v2 = m.get(i, i);
      if (v2.abs() > 1e-9) hasSame = true;
    }
    if (hasFlip) x[q] = 1;
  }
  let j0 = 0;
  for (let q = 0; q < n; q++) if (x[q]) j0 |= 1 << q;
  const v00 = m.get(0, j0);
  if (v00.abs() < 1e-9) return null;
  const rel = v00.div(phase);
  let yCount = 0;
  if (Math.abs(rel.re - 1) < 1e-9 && Math.abs(rel.im) < 1e-9) yCount = 0;
  else if (Math.abs(rel.re) < 1e-9 && Math.abs(rel.im - 1) < 1e-9) yCount = 1;
  else if (Math.abs(rel.re + 1) < 1e-9 && Math.abs(rel.im) < 1e-9) yCount = 2;
  else if (Math.abs(rel.re) < 1e-9 && Math.abs(rel.im + 1) < 1e-9) yCount = 3;
  else return null;
  if (yCount === 0) {
    for (let q = 0; q < n; q++) {
      if (x[q]) continue;
      let hasNeg = false, hasPos = false;
      for (let i = 0; i < dim; i++) {
        if ((i >> q & 1) !== 1) continue;
        const v = m.get(i, i);
        if (v.abs() < 1e-9) continue;
        const r4 = v.div(phase);
        if (r4.re < -0.5) hasNeg = true;
        else if (r4.re > 0.5) hasPos = true;
      }
      if (hasNeg && !hasPos) z[q] = 1;
      else if (hasNeg && hasPos) return null;
    }
  } else {
    for (let q = 0; q < n; q++) {
      if (!x[q]) {
        let hasNeg = false, hasPos = false;
        for (let i = 0; i < dim; i++) {
          if ((i >> q & 1) !== 1) continue;
          const v = m.get(i, i);
          if (v.abs() < 1e-9) continue;
          const r4 = v.div(phase);
          if (r4.re < -0.5 || r4.im < -0.5) hasNeg = true;
          else if (r4.re > 0.5 || r4.im > 0.5) hasPos = true;
        }
        if (hasNeg && !hasPos) z[q] = 1;
        continue;
      }
      let found = false;
      for (let i = 0; i < dim && !found; i++) {
        const j = i ^ j0;
        const i2 = i ^ 1 << q;
        const j2 = i2 ^ j0;
        const v1 = m.get(i, j);
        const v2 = m.get(i2, j2);
        if (v1.abs() < 1e-9 || v2.abs() < 1e-9) continue;
        const r4 = v2.div(v1);
        if (Math.abs(r4.re + 1) < 1e-9 && Math.abs(r4.im) < 1e-9) {
          z[q] = 1;
          found = true;
        } else if (Math.abs(r4.re - 1) < 1e-9 && Math.abs(r4.im) < 1e-9) {
          found = true;
        } else {
          return null;
        }
      }
    }
  }
  const expectedPhase = [new Complex(1, 0), new Complex(0, 1), new Complex(-1, 0), new Complex(0, -1)][yCount];
  const signPhase = phase.div(expectedPhase);
  if (Math.abs(signPhase.im) > 1e-9) return null;
  const sign = signPhase.re < 0 ? 1 : 0;
  let actualYCount = 0;
  for (let q = 0; q < n; q++) if (x[q] && z[q]) actualYCount++;
  if (actualYCount !== yCount && actualYCount % 4 !== yCount % 4) {
    return null;
  }
  return { x, z, sign };
}

// src/quantum_info/scalar_op.js
var ScalarOp = class _ScalarOp {
  constructor(numQubits, coeff = null) {
    this._numQubits = numQubits;
    if (coeff instanceof Complex) {
      this._coeff = coeff;
    } else if (typeof coeff === "number") {
      this._coeff = new Complex(coeff, 0);
    } else if (coeff == null) {
      this._coeff = Complex.ONE;
    } else {
      throw new TypeError("ScalarOp coeff must be a number or Complex");
    }
  }
  static fromOperator(operator) {
    const dim = operator.dim;
    const c0 = operator.data.get(0, 0);
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        if (i === j) {
          if (!operator.data.get(i, i).equals(c0, 1e-9)) {
            throw new Error("Operator is not a scalar multiple of identity");
          }
        } else if (operator.data.get(i, j).abs() > 1e-9) {
          throw new Error("Operator is not a scalar multiple of identity");
        }
      }
    }
    return new _ScalarOp(operator.num_qubits, c0);
  }
  get num_qubits() {
    return this._numQubits;
  }
  get coeff() {
    return this._coeff;
  }
  get dim() {
    return 1 << this._numQubits;
  }
  to_matrix() {
    return ComplexMatrix2.identity(this.dim).scale(this._coeff);
  }
  to_operator() {
    return new Operator(this.to_matrix());
  }
  compose(other) {
    if (other instanceof _ScalarOp) {
      return new _ScalarOp(this._numQubits, this._coeff.mul(other._coeff));
    }
    return other.scale(this._coeff);
  }
  tensor(other) {
    if (other instanceof _ScalarOp) {
      return new _ScalarOp(this._numQubits + other._numQubits, this._coeff.mul(other._coeff));
    }
    return other.scale(this._coeff);
  }
  expand(other) {
    return this.tensor(other);
  }
  adjoint() {
    return new _ScalarOp(this._numQubits, this._coeff.conjugate());
  }
  conjugate() {
    return new _ScalarOp(this._numQubits, this._coeff.conjugate());
  }
  transpose() {
    return this;
  }
  trace() {
    return this._coeff.scale(this.dim);
  }
  det() {
    return this._coeff.pow(this.dim);
  }
  is_unitary(tol = 1e-9) {
    return Math.abs(this._coeff.abs() - 1) < tol;
  }
  equals(other, tol) {
    if (!(other instanceof _ScalarOp)) return false;
    if (this._numQubits !== other._numQubits) return false;
    return this._coeff.equals(other._coeff, tol);
  }
  apply_to_vector(vector) {
    return vector.scale(this._coeff);
  }
  expectation_value(statevector) {
    const norm2 = statevector.norm() ** 2;
    return this._coeff.scale(norm2);
  }
  toString() {
    return `ScalarOp(${this._coeff.toString()}, num_qubits=${this._numQubits})`;
  }
};

// src/quantum_info/schmidt.js
var SchmidtDecomposition = class {
  constructor(statevector, num_qubits_a = null) {
    let nA;
    if (num_qubits_a !== null) {
      nA = num_qubits_a;
    } else {
      nA = Math.floor(statevector.num_qubits / 2);
    }
    this.num_qubits_a = nA;
    this.num_qubits_b = statevector.num_qubits - nA;
    this.dim_a = 1 << nA;
    this.dim_b = 1 << this.num_qubits_b;
    const m = ComplexMatrix2.zeros(this.dim_a, this.dim_b);
    for (let i = 0; i < this.dim_a; i++) {
      for (let j = 0; j < this.dim_b; j++) {
        const idx = i << this.num_qubits_b | j;
        m.set(i, j, statevector.data.get(idx));
      }
    }
    const { U, S, Vh } = m.svd();
    this.U = U;
    this.S = S;
    this.Vh = Vh;
    this._statevector = statevector;
    this.schmidt_rank = S.filter((s) => s > 1e-12).length;
  }
  // Get the Schmidt coefficients
  schmidt_coefficients() {
    return this.S.slice(0, this.schmidt_rank);
  }
  // Get the k-th Schmidt state of subsystem A
  schmidt_state_a(k) {
    if (k >= this.schmidt_rank) return null;
    const data = new Array(this.dim_a);
    for (let i = 0; i < this.dim_a; i++) {
      data[i] = this.U.get(i, k);
    }
    return new Statevector(new ComplexVector(data), this.num_qubits_a);
  }
  // Get the k-th Schmidt state of subsystem B
  schmidt_state_b(k) {
    if (k >= this.schmidt_rank) return null;
    const data = new Array(this.dim_b);
    for (let j = 0; j < this.dim_b; j++) {
      data[j] = this.Vh.get(k, j).conjugate();
    }
    return new Statevector(new ComplexVector(data), this.num_qubits_b);
  }
  // Entanglement entropy (von Neumann) = -sum p_k log2(p_k)
  entropy() {
    let h = 0;
    for (const s of this.schmidt_coefficients()) {
      if (s > 1e-12) {
        const p = s * s;
        h -= p * Math.log2(p);
      }
    }
    return h;
  }
  // Schmidt number (rank)
  rank() {
    return this.schmidt_rank;
  }
  // Check if the state is entangled (schmidt_rank > 1)
  is_entangled(tol = 1e-9) {
    return this.schmidt_rank > 1 || this.schmidt_rank === 1 && Math.abs(this.S[0] - 1) > tol;
  }
  // Reconstruct the original statevector from the decomposition
  to_statevector() {
    return this._statevector;
  }
  toString() {
    const coeffs = this.schmidt_coefficients();
    return `SchmidtDecomposition(rank=${this.schmidt_rank}, entropy=${this.entropy().toFixed(4)}, coeffs=[${coeffs.map((c) => c.toFixed(4)).join(", ")}])`;
  }
};

// src/quantum_info/channels.js
var QuantumChannel = class {
  constructor(numQubits) {
    this.num_qubits = numQubits;
    this.dim = 1 << numQubits;
  }
  // Convert to superoperator (dim^2 x dim^2 matrix)
  to_superop() {
    throw new Error("to_superop not implemented");
  }
  // Convert to Kraus list
  to_kraus() {
    throw new Error("to_kraus not implemented");
  }
  // Convert to Choi matrix
  to_choi() {
    throw new Error("to_choi not implemented");
  }
  // Convert to PTM
  to_ptm() {
    throw new Error("to_ptm not implemented");
  }
  // Apply channel to a statevector -> density matrix
  apply(statevector) {
    const rho = _stateToDensity(statevector);
    return this.apply_density(rho);
  }
  // Apply channel to a density matrix
  apply_density(rho) {
    const superop = this.to_superop();
    const vec = _matrixToVector(rho);
    const result = superop.matvec(vec);
    return _vectorToMatrix(result, this.dim);
  }
  // Compose two channels (this after other)
  compose(other) {
    const s1 = this.to_superop();
    const s2 = other.to_superop();
    return SuperOp.from_matrix(s1.mul(s2), this.num_qubits);
  }
  // Tensor product
  tensor(other) {
    const s1 = this.to_superop();
    const s2 = other.to_superop();
    return SuperOp.from_matrix(s1.tensor(s2), this.num_qubits + other.num_qubits);
  }
};
var Kraus = class _Kraus extends QuantumChannel {
  constructor(data, numQubits = null) {
    const n = numQubits !== null ? numQubits : Math.log2(data[0].rows);
    super(n);
    this.data = data;
  }
  static from_operator(operator) {
    return new _Kraus([operator.to_matrix ? operator.to_matrix() : operator], operator.num_qubits || Math.log2(operator.rows));
  }
  to_kraus() {
    return this.data;
  }
  to_superop() {
    const dim = this.dim;
    const superDim = dim * dim;
    const result = ComplexMatrix2.zeros(superDim, superDim);
    for (const K of this.data) {
      const Kconj = K.conjugate();
      const Kkron = Kconj.tensor(K);
      result.add_inplace(Kkron);
    }
    return result;
  }
  to_choi() {
    const dim = this.dim;
    const choiDim = dim * dim;
    const result = ComplexMatrix2.zeros(choiDim, choiDim);
    for (const K of this.data) {
      const vec = _matrixToVector(K);
      for (let i = 0; i < vec.size; i++) {
        for (let j = 0; j < vec.size; j++) {
          const val = result.get(i, j).add(vec.data[i].mul(vec.data[j].conjugate()));
          result.set(i, j, val);
        }
      }
    }
    return result;
  }
  to_ptm() {
    const n = this.num_qubits;
    const dim = this.dim;
    const paulis = _generatePauliOps(n);
    const numPaulis = paulis.length;
    const ptm = ComplexMatrix2.zeros(numPaulis, numPaulis);
    for (let i = 0; i < numPaulis; i++) {
      for (let j = 0; j < numPaulis; j++) {
        const ej = this.apply_density(paulis[j]);
        const trace = _traceOfProduct(paulis[i], ej);
        ptm.set(i, j, trace.scale(1 / dim));
      }
    }
    return ptm;
  }
};
var SuperOp = class _SuperOp extends QuantumChannel {
  constructor(data, numQubits) {
    super(numQubits);
    this.data = data;
  }
  static from_matrix(matrix, numQubits) {
    return new _SuperOp(matrix, numQubits);
  }
  to_superop() {
    return this.data;
  }
  to_kraus() {
    const choi = this.to_choi();
    const { eigenvalues, eigenvectors } = choi.eigh();
    const krausOps = [];
    for (let i = 0; i < eigenvalues.length; i++) {
      if (eigenvalues[i] > 1e-10) {
        const sqrtEv = Math.sqrt(eigenvalues[i]);
        const vec = new ComplexVector(new Array(this.dim * this.dim));
        for (let j = 0; j < this.dim * this.dim; j++) {
          vec.data[j] = eigenvectors.get(j, i);
        }
        const K = _vectorToMatrix(vec, this.dim).scale(sqrtEv);
        krausOps.push(K);
      }
    }
    return krausOps;
  }
  to_choi() {
    const dim = this.dim;
    const choi = ComplexMatrix2.zeros(dim * dim, dim * dim);
    for (let i = 0; i < dim; i++) {
      for (let j = 0; j < dim; j++) {
        for (let k = 0; k < dim; k++) {
          for (let l = 0; l < dim; l++) {
            const sRow = i * dim + j;
            const sCol = k * dim + l;
            const cRow = i * dim + k;
            const cCol = j * dim + l;
            choi.set(cRow, cCol, this.data.get(sRow, sCol));
          }
        }
      }
    }
    return choi;
  }
};
var Chi = class extends QuantumChannel {
  constructor(data, numQubits) {
    super(numQubits);
    this.data = data;
  }
  to_superop() {
    const n = this.num_qubits;
    const paulis = _generatePauliOps(n);
    const dim = this.dim;
    const superDim = dim * dim;
    const result = ComplexMatrix2.zeros(superDim, superDim);
    for (let i = 0; i < paulis.length; i++) {
      for (let j = 0; j < paulis.length; j++) {
        const chiVal = this.data.get(i, j);
        if (chiVal.abs() < 1e-15) continue;
        const term = paulis[i].conjugate().tensor(paulis[j]).scale(chiVal);
        result.add_inplace(term);
      }
    }
    return result;
  }
  to_kraus() {
    return new Kraus(this.to_kraus(), this.num_qubits);
  }
  to_choi() {
    const superop = this.to_superop();
    return new SuperOp(superop, this.num_qubits).to_choi();
  }
};
var PTM = class extends QuantumChannel {
  constructor(data, numQubits) {
    super(numQubits);
    this.data = data;
  }
  to_superop() {
    const n = this.num_qubits;
    const paulis = _generatePauliOps(n);
    const dim = this.dim;
    const superDim = dim * dim;
    const result = ComplexMatrix2.zeros(superDim, superDim);
    for (let i = 0; i < paulis.length; i++) {
      for (let j = 0; j < paulis.length; j++) {
        const rVal = this.data.get(i, j);
        if (rVal.re < 1e-15 && rVal.im < 1e-15) continue;
        const term = paulis[i].conjugate().tensor(paulis[j]).scale(rVal.scale(1 / dim));
        result.add_inplace(term);
      }
    }
    return result;
  }
};
function state_fidelity(state1, state2) {
  const isPure1 = state1._data && state1._data.size !== void 0;
  const isPure2 = state2._data && state2._data.size !== void 0;
  if (isPure1 && isPure2) {
    const inner = state1.data.inner(state2.data);
    return inner.abs2();
  }
  const rho1 = isPure1 ? _stateToDensity(state1) : state1;
  const rho2 = isPure2 ? _stateToDensity(state2) : state2;
  const product = rho1.mul(rho2);
  return product.trace().re;
}
function process_fidelity(channel, operator = null) {
  if (operator === null) {
    const choi = channel.to_choi();
    const dim = channel.dim;
    const identity = ComplexMatrix2.identity(dim * dim).scale(1 / dim);
    return choi.mul(identity).trace().re;
  }
  const sE = channel.to_superop();
  const sU = _unitaryToSuperop(operator.to_matrix ? operator.to_matrix() : operator);
  const product = sE.mul(sU.dagger());
  return product.trace().re / (channel.dim * channel.dim);
}
function average_gate_fidelity(channel, operator = null) {
  const d = channel.dim;
  const fproc = process_fidelity(channel, operator);
  return (d * fproc + 1) / (d + 1);
}
function diamond_norm(channel) {
  const choi = channel.to_choi();
  const { eigenvalues } = choi.eigh();
  return eigenvalues.reduce((s, v) => s + Math.abs(v), 0);
}
function _stateToDensity(statevector) {
  const dim = statevector.size || statevector.dim;
  const data = statevector.data || statevector._data;
  const rho = ComplexMatrix2.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      rho.set(i, j, data.get(i).mul(data.get(j).conjugate()));
    }
  }
  return rho;
}
function _matrixToVector(matrix) {
  const dim = matrix.rows;
  const vec = new ComplexVector(new Array(dim * dim));
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      vec.data[i * dim + j] = matrix.get(i, j);
    }
  }
  return vec;
}
function _vectorToMatrix(vec, dim) {
  const m = ComplexMatrix2.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      m.set(i, j, vec.data[i * dim + j]);
    }
  }
  return m;
}
function _traceOfProduct(m1, m2) {
  const product = m1.mul(m2);
  return product.trace();
}
function _unitaryToSuperop(U) {
  return U.conjugate().tensor(U);
}
function _generatePauliOps(numQubits) {
  const paulis1 = [ComplexMatrix2.identity(2), PAULI.X, PAULI.Y, PAULI.Z];
  if (numQubits === 1) return paulis1;
  let result = paulis1;
  for (let i = 1; i < numQubits; i++) {
    const newResult = [];
    for (const p1 of result) {
      for (const p2 of paulis1) {
        newResult.push(p1.tensor(p2));
      }
    }
    result = newResult;
  }
  return result;
}

// src/quantum_info/quantum_info_extra.js
function random_unitary(numQubits, seed = null) {
  const dim = 1 << numQubits;
  const rng = _makeRng2(seed);
  const m = ComplexMatrix2.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      m.set(i, j, new Complex(_gaussian(rng), _gaussian(rng)));
    }
  }
  const Q = ComplexMatrix2.zeros(dim, dim);
  const R = ComplexMatrix2.zeros(dim, dim);
  for (let j = 0; j < dim; j++) {
    const v = new Array(dim);
    for (let i = 0; i < dim; i++) v[i] = m.get(i, j);
    for (let k = 0; k < j; k++) {
      let rkj = Complex.ZERO;
      for (let i = 0; i < dim; i++) {
        rkj = rkj.add(Q.get(i, k).conjugate().mul(v[i]));
      }
      R.set(k, j, rkj);
      for (let i = 0; i < dim; i++) {
        v[i] = v[i].sub(rkj.mul(Q.get(i, k)));
      }
    }
    let normSq = Complex.ZERO;
    for (let i = 0; i < dim; i++) normSq = normSq.add(v[i].mul(v[i].conjugate()));
    const norm = Math.sqrt(normSq.re);
    R.set(j, j, new Complex(norm, 0));
    if (norm < 1e-12) {
      for (let i = 0; i < dim; i++) v[i] = new Complex(_gaussian(rng), _gaussian(rng));
      let ns = Complex.ZERO;
      for (let i = 0; i < dim; i++) ns = ns.add(v[i].mul(v[i].conjugate()));
      const nn = Math.sqrt(ns.re);
      for (let i = 0; i < dim; i++) Q.set(i, j, v[i].scale(1 / nn));
    } else {
      for (let i = 0; i < dim; i++) Q.set(i, j, v[i].scale(1 / norm));
    }
  }
  for (let j = 0; j < dim; j++) {
    const rjj = R.get(j, j);
    if (rjj.abs() > 1e-12) {
      const phase = rjj.scale(1 / rjj.abs());
      for (let i = 0; i < dim; i++) {
        Q.set(i, j, Q.get(i, j).mul(phase.conjugate()));
      }
    }
  }
  return new Operator(Q);
}
function random_statevector(numQubits, seed = null) {
  const dim = 1 << numQubits;
  const rng = _makeRng2(seed);
  const data = new Array(dim);
  let normSq = 0;
  for (let i = 0; i < dim; i++) {
    const re = _gaussian(rng);
    const im = _gaussian(rng);
    data[i] = new Complex(re, im);
    normSq += re * re + im * im;
  }
  const norm = Math.sqrt(normSq);
  const normalized = data.map((c) => new Complex(c.re / norm, c.im / norm));
  return new Statevector(new ComplexVector(normalized), numQubits);
}
function random_pauli(numQubits, seed = null) {
  const rng = _makeRng2(seed);
  let label = "";
  const chars = ["I", "X", "Y", "Z"];
  for (let i = 0; i < numQubits; i++) {
    label += chars[Math.floor(rng() * 4)];
  }
  return new Pauli(label);
}
function random_clifford(numQubits, seed = null) {
  return Clifford.random(numQubits, seed);
}
function random_density_matrix(numQubits, seed = null, mixedWeight = null) {
  const dim = 1 << numQubits;
  const rng = _makeRng2(seed);
  const psi = random_statevector(numQubits, seed);
  const pureRho = psi.to_operator()._data;
  const w = mixedWeight !== null ? mixedWeight : rng();
  const mixedRho = ComplexMatrix2.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      mixedRho.set(i, j, pureRho.get(i, j).scale(1 - w).add(
        i === j ? new Complex(w / dim, 0) : Complex.ZERO
      ));
    }
  }
  return new DensityMatrix(mixedRho);
}
function purity(densityMatrix) {
  const rho = densityMatrix._data;
  const dim = rho.rows;
  let tr = 0;
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      tr += rho.get(i, j).mul(rho.get(j, i)).re;
    }
  }
  return tr;
}
function concurrence(densityMatrix) {
  if (densityMatrix.num_qubits !== 2) {
    throw new Error("concurrence is only defined for 2-qubit states");
  }
  const rho = densityMatrix._data;
  const YY = PAULI.Y.tensor(PAULI.Y);
  const rhoStar = rho.conjugate();
  const R = rho.mul(YY).mul(rhoStar).mul(YY);
  const { eigenvalues, eigenvectors } = R.eigh();
  const sqrtEigs = eigenvalues.map((v) => Math.sqrt(Math.abs(v)));
  sqrtEigs.sort((a, b) => b - a);
  const C = Math.max(0, sqrtEigs[0] - sqrtEigs[1] - sqrtEigs[2] - sqrtEigs[3]);
  return C;
}
function entanglement_of_formation(densityMatrix) {
  const C = concurrence(densityMatrix);
  if (C < 1e-12) return 0;
  const x = (1 + Math.sqrt(1 - C * C)) / 2;
  const h = -x * Math.log2(x) - (1 - x) * Math.log2(1 - x);
  return h;
}
function mutual_information(densityMatrix, subsystemA, subsystemB = null) {
  const n = densityMatrix.num_qubits;
  if (subsystemB == null) {
    subsystemB = [];
    for (let q = 0; q < n; q++) {
      if (!subsystemA.includes(q)) subsystemB.push(q);
    }
  }
  const rhoAB = densityMatrix;
  const rhoA = rhoAB.partial_trace(_complement2(subsystemA, n));
  const rhoB = rhoAB.partial_trace(_complement2(subsystemB, n));
  const sAB = _vonNeumannEntropy(rhoAB);
  const sA = _vonNeumannEntropy(rhoA);
  const sB = _vonNeumannEntropy(rhoB);
  return sA + sB - sAB;
}
function gate_fidelity(unitary1, unitary2) {
  const u1 = unitary1._data || unitary1;
  const u2 = unitary2._data || unitary2;
  if (u1.rows !== u2.rows) {
    throw new Error("gate_fidelity: unitaries must have the same dimension");
  }
  const d = u1.rows;
  const u1dag = u1.dagger();
  const product = u1dag.mul(u2);
  const tr = product.trace();
  return (tr.re * tr.re + tr.im * tr.im) / (d * d);
}
function unitarity(channel) {
  const d = channel.dim;
  const choi = channel.to_choi();
  const { eigenvalues } = choi.eigh();
  let sumSq = 0;
  for (const e of eigenvalues) sumSq += e * e;
  return Math.sqrt(sumSq) / d;
}
function _makeRng2(seed) {
  if (seed == null) return Math.random;
  let s = seed >>> 0;
  return () => {
    s = s + 1831565813 | 0;
    let t = Math.imul(s ^ s >>> 15, 1 | s);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}
function _gaussian(rng) {
  let u = 0, v = 0;
  while (u === 0) u = rng();
  while (v === 0) v = rng();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
function _complement2(subset, n) {
  const set = new Set(subset);
  const out = [];
  for (let i = 0; i < n; i++) {
    if (!set.has(i)) out.push(i);
  }
  return out;
}
function _vonNeumannEntropy(densityMatrix) {
  const rho = densityMatrix._data;
  const { eigenvalues } = rho.eigh();
  let s = 0;
  for (const e of eigenvalues) {
    if (e > 1e-12) s -= e * Math.log2(e);
  }
  return s;
}

// src/result/result.js
var Counts = class _Counts {
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
  get(key) {
    return this._data[key] || 0;
  }
  set(key, value) {
    this._data[key] = value;
  }
  keys() {
    return Object.keys(this._data);
  }
  values() {
    return Object.values(this._data);
  }
  items() {
    return Object.entries(this._data);
  }
  most_frequent() {
    let bestKey = null, bestCount = -1;
    for (const k in this._data) {
      if (this._data[k] > bestCount) {
        bestCount = this._data[k];
        bestKey = k;
      }
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
  to_dict() {
    return Object.assign({}, this._data);
  }
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
        const bit = c >> idxs[i] & 1;
        new_c |= bit << i;
      }
      let newKey = "";
      for (let i = idxs.length - 1; i >= 0; i--) {
        newKey += (new_c >> i & 1).toString();
      }
      out[newKey] = (out[newKey] || 0) + this._data[k];
    }
    return new _Counts(out, this.hexKeys);
  }
  clear() {
    this._data = {};
  }
};
var Result = class {
  constructor(kwargs) {
    this.backend_name = kwargs.backend_name || "unknown";
    this.backend_version = kwargs.backend_version || "0.0.0";
    this.qobj_id = kwargs.qobj_id || "qobj";
    this.job_id = kwargs.job_id || `job_${Date.now()}`;
    this.success = kwargs.success !== void 0 ? kwargs.success : true;
    this.results = kwargs.results || [];
    this.date = kwargs.date || (/* @__PURE__ */ new Date()).toISOString();
    this.status = kwargs.status || "COMPLETED";
  }
  get_counts(experiment_id = 0) {
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
  to_dict() {
    return {
      backend_name: this.backend_name,
      backend_version: this.backend_version,
      qobj_id: this.qobj_id,
      job_id: this.job_id,
      success: this.success,
      results: this.results,
      date: this.date,
      status: this.status
    };
  }
};

// src/simulator/statevector_simulator.js
function _applyGateInPlace(state, gateMatrix, k, qubitIndices) {
  const n = state.size;
  const dimK = 1 << k;
  const processed = new Array(n).fill(false);
  const newState = new Array(n);
  for (let i = 0; i < n; i++) newState[i] = state.data[i];
  for (let i = 0; i < n; i++) {
    if (processed[i]) continue;
    const subs = new Array(dimK);
    const idxs = new Array(dimK);
    for (let s = 0; s < dimK; s++) {
      let idx = i;
      for (let qi = 0; qi < k; qi++) {
        const mask = 1 << qubitIndices[qi];
        idx &= ~mask;
      }
      for (let qi = 0; qi < k; qi++) {
        if (s >> qi & 1) idx |= 1 << qubitIndices[qi];
      }
      idxs[s] = idx;
      subs[s] = state.data[idx];
      processed[idx] = true;
    }
    const newSubs = new Array(dimK);
    for (let rr = 0; rr < dimK; rr++) {
      let acc = Complex.ZERO;
      for (let cc = 0; cc < dimK; cc++) {
        acc = acc.add(gateMatrix.get(rr, cc).mul(subs[cc]));
      }
      newSubs[rr] = acc;
    }
    for (let s = 0; s < dimK; s++) newState[idxs[s]] = newSubs[s];
  }
  state.data = newState;
}
var StatevectorSimulator = class {
  constructor(options = {}) {
    this.method = options.method || "statevector";
    this.precision = options.precision || "double";
    this.max_memory_mb = options.max_memory_mb || 8192;
    this.initial_statevector = options.initial_statevector || null;
  }
  run(circuit, shots = 1024, options = {}) {
    const experiments = Array.isArray(circuit) ? circuit : [circuit];
    const results = experiments.map((c) => this._runOne(c, shots, options));
    return new Result({
      backend_name: "statevector_simulator",
      backend_version: "1.0.0",
      qobj_id: options.qobj_id || "qobj",
      job_id: options.job_id || `job_${Date.now()}`,
      success: true,
      results
    });
  }
  _runOne(circuit, shots, options) {
    const n = circuit.num_qubits;
    const state = this.initial_statevector ? new ComplexVector(this.initial_statevector.data.slice()) : ComplexVector.zeros(1 << n);
    if (!this.initial_statevector) state.data[0] = Complex.ONE;
    const finalMeasures = [];
    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier") continue;
      if (op.name === "delay") continue;
      if (op.name === "measure") {
        const qIdx = circuit._qubit_index.get(ci.qubits[0]);
        const cIdx = circuit._clbit_index.get(ci.clbits[0]);
        finalMeasures.push({ qIdx, cIdx });
        continue;
      }
      if (op.name === "reset") {
        const qIdx = circuit._qubit_index.get(ci.qubits[0]);
        const prob0 = _probOfZero(state, qIdx);
        if (randomUniform() < prob0) {
          _projectToZero(state, qIdx, Math.sqrt(prob0));
        } else {
          _projectToOne(state, qIdx, Math.sqrt(1 - prob0));
        }
        continue;
      }
      if (op.name === "if_else" || op.name === "while_loop") continue;
      if (typeof op.to_matrix !== "function") {
        throw new Error(`Cannot simulate ${op.name} (no matrix)`);
      }
      const gateMatrix = op.to_matrix();
      const qubitIndices = ci.qubits.map((q) => circuit._qubit_index.get(q));
      _applyGateInPlace(state, gateMatrix, op.num_qubits, qubitIndices);
    }
    let counts = null;
    let memory = null;
    if (shots > 0) {
      const probs = state.probabilities();
      counts = {};
      memory = new Array(shots);
      const numClbits = circuit.num_clbits;
      const numQubits = circuit.num_qubits;
      for (let s = 0; s < shots; s++) {
        const idx = sampleDistribution(probs, options.rng);
        const bits = new Array(numClbits).fill("0");
        if (finalMeasures.length === 0) {
          for (let q = 0; q < numQubits && q < numClbits; q++) {
            bits[q] = (idx >> q & 1).toString();
          }
        } else {
          for (const m of finalMeasures) {
            const b = idx >> m.qIdx & 1;
            bits[m.cIdx] = b.toString();
          }
        }
        const key = bits.reverse().join("");
        counts[key] = (counts[key] || 0) + 1;
        memory[s] = key;
      }
    }
    return {
      success: true,
      shots,
      data: {
        statevector: new ComplexVector(state.data.slice()),
        counts: counts ? new Counts(counts) : null,
        memory
      },
      status: "DONE",
      time_taken: 0
    };
  }
};
function _probOfZero(state, qIdx) {
  let p0 = 0;
  for (let i = 0; i < state.size; i++) {
    if ((i >> qIdx & 1) === 0) p0 += state.data[i].abs2();
  }
  return p0;
}
function _projectToZero(state, qIdx, norm) {
  for (let i = 0; i < state.size; i++) {
    if ((i >> qIdx & 1) === 1) state.data[i] = Complex.ZERO;
    else state.data[i] = state.data[i].scale(1 / norm);
  }
}
function _projectToOne(state, qIdx, norm) {
  for (let i = 0; i < state.size; i++) {
    if ((i >> qIdx & 1) === 0) state.data[i] = Complex.ZERO;
    else state.data[i] = state.data[i].scale(1 / norm);
  }
}
function simulate(circuit, shots = 1024, options = {}) {
  const sim = new StatevectorSimulator(options);
  return sim.run(circuit, shots, options);
}

// src/simulator/qasm_simulator.js
var QasmSimulator = class {
  constructor(options = {}) {
    this.method = options.method || "automatic";
    this.noise_model = options.noise_model || null;
    this.basis_gates = options.basis_gates || null;
    this.max_memory_mb = options.max_memory_mb || 8192;
    this.seed = options.seed || null;
  }
  // Run one or more circuits
  run(circuits2, shots = 1024, options = {}) {
    const circuitList = Array.isArray(circuits2) ? circuits2 : [circuits2];
    const noiseModel = options.noise_model || this.noise_model;
    const results = circuitList.map((c) => this._runOne(c, shots, noiseModel, options));
    return new Result({
      backend_name: "qasm_simulator",
      backend_version: "1.0.0",
      qobj_id: options.qobj_id || "qobj",
      job_id: options.job_id || `job_${Date.now()}`,
      success: true,
      results
    });
  }
  _runOne(circuit, shots, noiseModel, options) {
    const n = circuit.num_qubits;
    const hasNoise = noiseModel && !noiseModel.is_empty();
    if (hasNoise) {
      return this._runNoisy(circuit, shots, noiseModel, options);
    }
    return this._runIdeal(circuit, shots, options);
  }
  // Ideal simulation: statevector + sampling
  _runIdeal(circuit, shots, options) {
    const n = circuit.num_qubits;
    const state = ComplexVector.zeros(1 << n);
    state.data[0] = Complex.ONE;
    const finalMeasures = [];
    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "delay") continue;
      if (op.name === "measure") {
        const qIdx = circuit._qubit_index.get(ci.qubits[0]);
        const cIdx = circuit._clbit_index.get(ci.clbits[0]);
        finalMeasures.push({ qIdx, cIdx });
        continue;
      }
      if (op.name === "reset") {
        const qIdx = circuit._qubit_index.get(ci.qubits[0]);
        const prob0 = _probOfZero2(state, qIdx);
        if (randomUniform() < prob0) {
          _projectToZero2(state, qIdx, Math.sqrt(prob0));
        } else {
          _projectToOne2(state, qIdx, Math.sqrt(1 - prob0));
        }
        continue;
      }
      if (op.name === "if_else" || op.name === "while_loop") continue;
      if (typeof op.to_matrix !== "function") continue;
      const gateMatrix = op.to_matrix();
      const qubitIndices = ci.qubits.map((q) => circuit._qubit_index.get(q));
      _applyGateInPlace2(state, gateMatrix, op.num_qubits, qubitIndices);
    }
    const probs = state.probabilities();
    const counts = {};
    const memory = new Array(shots);
    const numClbits = circuit.num_clbits;
    const numQubits = circuit.num_qubits;
    for (let s = 0; s < shots; s++) {
      const idx = sampleDistribution(probs, options.rng);
      const bits = new Array(numClbits).fill("0");
      if (finalMeasures.length === 0) {
        for (let q = 0; q < numQubits && q < numClbits; q++) {
          bits[q] = (idx >> q & 1).toString();
        }
      } else {
        for (const m of finalMeasures) {
          bits[m.cIdx] = (idx >> m.qIdx & 1).toString();
        }
      }
      const key = bits.reverse().join("");
      counts[key] = (counts[key] || 0) + 1;
      memory[s] = key;
    }
    return {
      success: true,
      shots,
      data: { counts: new Counts(counts), memory },
      status: "DONE"
    };
  }
  // Noisy simulation: density matrix evolution with Kraus operators
  _runNoisy(circuit, shots, noiseModel, options) {
    const n = circuit.num_qubits;
    const counts = {};
    const memory = new Array(shots);
    const finalMeasures = [];
    for (const ci of circuit.data) {
      if (ci.operation.name === "measure") {
        finalMeasures.push({
          op: ci.operation,
          qIdx: circuit._qubit_index.get(ci.qubits[0]),
          cIdx: circuit._clbit_index.get(ci.clbits[0]),
          qubits: ci.qubits.map((q) => circuit._qubit_index.get(q)),
          clbits: ci.clbits.map((c) => circuit._clbit_index.get(c))
        });
      }
    }
    for (let shot = 0; shot < shots; shot++) {
      const state = ComplexVector.zeros(1 << n);
      state.data[0] = Complex.ONE;
      let measuredBits = new Array(circuit.num_clbits).fill(0);
      for (const ci of circuit.data) {
        const op = ci.operation;
        if (op.name === "barrier" || op.name === "delay") continue;
        if (op.name === "measure") {
          const qIdx = circuit._qubit_index.get(ci.qubits[0]);
          const cIdx = circuit._clbit_index.get(ci.clbits[0]);
          const prob1 = _probOfOne(state, qIdx);
          if (randomUniform() < prob1) {
            _projectToOne2(state, qIdx, Math.sqrt(prob1));
            measuredBits[cIdx] = 1;
          } else {
            _projectToZero2(state, qIdx, Math.sqrt(1 - prob1));
            measuredBits[cIdx] = 0;
          }
          continue;
        }
        if (op.name === "reset") {
          const qIdx = circuit._qubit_index.get(ci.qubits[0]);
          const prob0 = _probOfZero2(state, qIdx);
          if (randomUniform() < prob0) {
            _projectToZero2(state, qIdx, Math.sqrt(prob0));
          } else {
            _projectToOne2(state, qIdx, Math.sqrt(1 - prob0));
          }
          continue;
        }
        if (op.name === "if_else" || op.name === "while_loop") continue;
        if (typeof op.to_matrix !== "function") continue;
        const gateMatrix = op.to_matrix();
        const qubitIndices = ci.qubits.map((q) => circuit._qubit_index.get(q));
        _applyGateInPlace2(state, gateMatrix, op.num_qubits, qubitIndices);
        const error = noiseModel.get_quantum_error(op.name, qubitIndices);
        if (error) {
          _applyQuantumError(state, error, qubitIndices, n);
        }
      }
      for (const re of noiseModel.get_readout_errors()) {
        const targetQubits = re.qubits;
        for (const m of finalMeasures) {
          let appliesThisQubit = false;
          if (targetQubits === "all") appliesThisQubit = true;
          else if (Array.isArray(targetQubits) && targetQubits.includes(m.qIdx)) {
            appliesThisQubit = true;
          }
          if (!appliesThisQubit) continue;
          const r4 = randomUniform();
          const errProb = re.error.probabilities[measuredBits[m.cIdx]][1 - measuredBits[m.cIdx]];
          if (r4 < errProb) measuredBits[m.cIdx] = 1 - measuredBits[m.cIdx];
        }
      }
      const key = measuredBits.reverse().join("");
      counts[key] = (counts[key] || 0) + 1;
      memory[shot] = key;
    }
    return {
      success: true,
      shots,
      data: { counts: new Counts(counts), memory },
      status: "DONE"
    };
  }
};
function _applyGateInPlace2(state, gateMatrix, k, qubitIndices) {
  const n = state.size;
  const dimK = 1 << k;
  const processed = new Array(n).fill(false);
  const newState = new Array(n);
  for (let i = 0; i < n; i++) newState[i] = state.data[i];
  for (let i = 0; i < n; i++) {
    if (processed[i]) continue;
    const subs = new Array(dimK);
    const idxs = new Array(dimK);
    for (let s = 0; s < dimK; s++) {
      let idx = i;
      for (let qi = 0; qi < k; qi++) {
        const mask = 1 << qubitIndices[qi];
        idx &= ~mask;
      }
      for (let qi = 0; qi < k; qi++) {
        if (s >> qi & 1) idx |= 1 << qubitIndices[qi];
      }
      idxs[s] = idx;
      subs[s] = state.data[idx];
      processed[idx] = true;
    }
    const newSubs = new Array(dimK);
    for (let rr = 0; rr < dimK; rr++) {
      let acc = Complex.ZERO;
      for (let cc = 0; cc < dimK; cc++) {
        acc = acc.add(gateMatrix.get(rr, cc).mul(subs[cc]));
      }
      newSubs[rr] = acc;
    }
    for (let s = 0; s < dimK; s++) newState[idxs[s]] = newSubs[s];
  }
  state.data = newState;
}
function _probOfZero2(state, qIdx) {
  let p0 = 0;
  for (let i = 0; i < state.size; i++) {
    if ((i >> qIdx & 1) === 0) p0 += state.data[i].abs2();
  }
  return p0;
}
function _probOfOne(state, qIdx) {
  let p1 = 0;
  for (let i = 0; i < state.size; i++) {
    if ((i >> qIdx & 1) === 1) p1 += state.data[i].abs2();
  }
  return p1;
}
function _projectToZero2(state, qIdx, norm) {
  for (let i = 0; i < state.size; i++) {
    if ((i >> qIdx & 1) === 1) state.data[i] = Complex.ZERO;
    else state.data[i] = state.data[i].scale(1 / norm);
  }
}
function _projectToOne2(state, qIdx, norm) {
  for (let i = 0; i < state.size; i++) {
    if ((i >> qIdx & 1) === 0) state.data[i] = Complex.ZERO;
    else state.data[i] = state.data[i].scale(1 / norm);
  }
}
function _applyQuantumError(state, error, qubitIndices, totalQubits) {
  const r4 = randomUniform();
  let cum = 0;
  let chosenTerm = error.terms[0];
  for (const term of error.terms) {
    cum += term.probability;
    if (r4 < cum) {
      chosenTerm = term;
      break;
    }
  }
  if (chosenTerm.operators.length === 0) return;
  const k0 = chosenTerm.operators[0];
  const krausQubits = Math.log2(k0.rows);
  if (!Number.isInteger(krausQubits) || krausQubits < 1) {
    return;
  }
  let krausOp;
  if (chosenTerm.operators.length === 1) {
    krausOp = chosenTerm.operators[0];
  } else {
    const probs = chosenTerm.operators.map((K) => {
      const KdK = K.dagger().mul(K);
      const tempState = { size: state.size, data: state.data.slice() };
      _applyGateInPlace2(tempState, KdK, krausQubits, qubitIndices.slice(0, krausQubits));
      let p = 0;
      for (let i = 0; i < state.size; i++) {
        p += state.data[i].conjugate().mul(tempState.data[i]).re;
      }
      return Math.max(0, p);
    });
    const totalP = probs.reduce((a, b) => a + b, 0);
    if (totalP < 1e-15) return;
    const r22 = randomUniform() * totalP;
    let cum2 = 0;
    let idx = 0;
    for (let i = 0; i < probs.length; i++) {
      cum2 += probs[i];
      if (r22 < cum2) {
        idx = i;
        break;
      }
    }
    krausOp = chosenTerm.operators[idx];
  }
  _applyGateInPlace2(state, krausOp, krausQubits, qubitIndices.slice(0, krausQubits));
  const norm = state.norm();
  if (norm > 1e-15) {
    const inv = 1 / norm;
    for (let i = 0; i < state.size; i++) {
      state.data[i] = state.data[i].scale(inv);
    }
  }
}
function simulate_noisy(circuit, shots = 1024, options = {}) {
  const sim = new QasmSimulator(options);
  return sim.run(circuit, shots, options);
}

// src/noise/noise_models.js
var QuantumError = class _QuantumError {
  constructor(terms = []) {
    this.terms = terms;
    const total = terms.reduce((s, t) => s + t.probability, 0);
    if (terms.length > 0 && Math.abs(total - 1) > 1e-9) {
      for (const t of this.terms) t.probability /= total;
    }
  }
  get size() {
    return this.terms.length;
  }
  get num_qubits() {
    if (this.terms.length === 0) return 0;
    const rows = this.terms[0].operators[0].rows;
    const n = Math.log2(rows);
    if (!Number.isInteger(n)) {
      throw new Error(`QuantumError: Kraus operator has non-power-of-2 dimension ${rows}`);
    }
    return n;
  }
  // Apply this error to a density matrix (Kraus representation)
  apply(densityMatrix) {
    const dim = densityMatrix.rows;
    const result = ComplexMatrix2.zeros(dim, dim);
    for (const term of this.terms) {
      for (const K of term.operators) {
        const Krho = K.mul(densityMatrix);
        const KrhoKd = Krho.mul(K.dagger());
        result.add_inplace(KrhoKd.scale(term.probability));
      }
    }
    return result;
  }
  // Compose two errors (sequential application)
  compose(other) {
    const newTerms = [];
    for (const t1 of this.terms) {
      for (const t2 of other.terms) {
        const newOps = [];
        for (const K1 of t1.operators) {
          for (const K2 of t2.operators) {
            newOps.push(K2.mul(K1));
          }
        }
        newTerms.push({ operators: newOps, probability: t1.probability * t2.probability });
      }
    }
    return new _QuantumError(newTerms);
  }
  // Tensor product of two errors
  tensor(other) {
    const newTerms = [];
    for (const t1 of this.terms) {
      for (const t2 of other.terms) {
        const newOps = [];
        for (const K1 of t1.operators) {
          for (const K2 of t2.operators) {
            newOps.push(K1.tensor(K2));
          }
        }
        newTerms.push({ operators: newOps, probability: t1.probability * t2.probability });
      }
    }
    return new _QuantumError(newTerms);
  }
  copy() {
    return new _QuantumError(this.terms.map((t) => ({
      operators: t.operators.map((K) => new ComplexMatrix2(K.rows, K.cols, K.data.slice())),
      probability: t.probability
    })));
  }
};
function depolarizing_error(prob, numQubits) {
  const dim = 1 << numQubits;
  const identity = ComplexMatrix2.identity(dim);
  if (numQubits === 1) {
    const p = prob;
    const K0 = identity.scale(Math.sqrt(1 - p));
    const K1 = PAULI.X.scale(Math.sqrt(p / 3));
    const K2 = PAULI.Y.scale(Math.sqrt(p / 3));
    const K3 = PAULI.Z.scale(Math.sqrt(p / 3));
    return new QuantumError([
      { operators: [K0], probability: 1 - p },
      { operators: [K1], probability: p / 3 },
      { operators: [K2], probability: p / 3 },
      { operators: [K3], probability: p / 3 }
    ]);
  }
  const paulis = _generatePauliGroup(numQubits);
  const numNonId = paulis.length - 1;
  const terms = [{ operators: [identity], probability: 1 - prob }];
  for (let i = 1; i < paulis.length; i++) {
    terms.push({ operators: [paulis[i]], probability: prob / numNonId });
  }
  return new QuantumError(terms);
}
function bit_flip_error(prob) {
  const K0 = ComplexMatrix2.identity(2).scale(Math.sqrt(1 - prob));
  const K1 = PAULI.X.scale(Math.sqrt(prob));
  return new QuantumError([
    { operators: [K0], probability: 1 - prob },
    { operators: [K1], probability: prob }
  ]);
}
function phase_flip_error(prob) {
  const K0 = ComplexMatrix2.identity(2).scale(Math.sqrt(1 - prob));
  const K1 = PAULI.Z.scale(Math.sqrt(prob));
  return new QuantumError([
    { operators: [K0], probability: 1 - prob },
    { operators: [K1], probability: prob }
  ]);
}
function pauli_x_error(prob) {
  return bit_flip_error(prob);
}
function pauli_y_error(prob) {
  const K0 = ComplexMatrix2.identity(2).scale(Math.sqrt(1 - prob));
  const K1 = PAULI.Y.scale(Math.sqrt(prob));
  return new QuantumError([
    { operators: [K0], probability: 1 - prob },
    { operators: [K1], probability: prob }
  ]);
}
function pauli_z_error(prob) {
  return phase_flip_error(prob);
}
function amplitude_damping_error(gamma) {
  const K0 = ComplexMatrix2.fromRows([
    [new Complex(1, 0), new Complex(0, 0)],
    [new Complex(0, 0), new Complex(Math.sqrt(1 - gamma), 0)]
  ]);
  const K1 = ComplexMatrix2.fromRows([
    [new Complex(0, 0), new Complex(Math.sqrt(gamma), 0)],
    [new Complex(0, 0), new Complex(0, 0)]
  ]);
  return new QuantumError([
    { operators: [K0, K1], probability: 1 }
  ]);
}
function phase_damping_error(gamma) {
  const K0 = ComplexMatrix2.fromRows([
    [new Complex(1, 0), new Complex(0, 0)],
    [new Complex(0, 0), new Complex(Math.sqrt(1 - gamma), 0)]
  ]);
  const K1 = ComplexMatrix2.fromRows([
    [new Complex(0, 0), new Complex(0, 0)],
    [new Complex(0, 0), new Complex(Math.sqrt(gamma), 0)]
  ]);
  return new QuantumError([
    { operators: [K0, K1], probability: 1 }
  ]);
}
function reset_error(prob, numQubits = 1) {
  const dim = 1 << numQubits;
  const proj0 = ComplexMatrix2.zeros(dim, dim);
  proj0.set(0, 0, new Complex(1, 0));
  const identity = ComplexMatrix2.identity(dim);
  const K0 = identity.scale(Math.sqrt(1 - prob));
  const K1 = proj0.scale(Math.sqrt(prob));
  return new QuantumError([
    { operators: [K0], probability: 1 - prob },
    { operators: [K1], probability: prob }
  ]);
}
function kraus_error(krausOps) {
  return new QuantumError([
    { operators: krausOps, probability: 1 }
  ]);
}
function mixed_unitary_error(errors) {
  return new QuantumError(errors.map((e) => ({
    operators: [e.unitary],
    probability: e.probability
  })));
}
function ReadoutError(probabilities) {
  this.probabilities = probabilities;
  this.num_qubits = Math.log2(probabilities.length);
  this.apply = function(classicalBits) {
    const measured = [];
    for (let i = 0; i < classicalBits.length; i++) {
      const actual = classicalBits[i];
      const r4 = Math.random();
      const probActual = this.probabilities[actual];
      measured.push(r4 < probActual[0] ? 0 : 1);
    }
    return measured;
  };
  this.copy = function() {
    return new ReadoutError(this.probabilities.map((row) => row.slice()));
  };
}
var NoiseModel = class _NoiseModel {
  constructor() {
    this._local_quantum_errors = /* @__PURE__ */ new Map();
    this._local_readout_errors = /* @__PURE__ */ new Map();
    this._basis_gate_errors = /* @__PURE__ */ new Map();
    this._readout_errors = [];
  }
  get basis_gates() {
    const gates = /* @__PURE__ */ new Set();
    for (const g of this._local_quantum_errors.keys()) gates.add(g);
    for (const g of this._basis_gate_errors.keys()) gates.add(g);
    return Array.from(gates);
  }
  // Add a quantum error for a specific gate on specific qubits
  add_quantum_error(error, gates, qubits = null) {
    const gateList = Array.isArray(gates) ? gates : [gates];
    for (const gate of gateList) {
      const key = qubits ? `${gate}:${qubits.join(",")}` : gate;
      if (qubits) {
        if (!this._local_quantum_errors.has(gate)) {
          this._local_quantum_errors.set(gate, /* @__PURE__ */ new Map());
        }
        this._local_quantum_errors.get(gate).set(qubits.join(","), error);
      } else {
        this._basis_gate_errors.set(gate, error);
      }
    }
    return this;
  }
  // Add a readout error for specific qubits
  add_readout_error(error, qubits = null) {
    if (qubits) {
      this._readout_errors.push({ error, qubits });
    } else {
      this._readout_errors.push({ error, qubits: "all" });
    }
    return this;
  }
  // Get the quantum error for a gate on specific qubits (or null)
  get_quantum_error(gate, qubits) {
    if (this._local_quantum_errors.has(gate)) {
      const local = this._local_quantum_errors.get(gate);
      const key = qubits.join(",");
      if (local.has(key)) return local.get(key);
    }
    if (this._basis_gate_errors.has(gate)) {
      return this._basis_gate_errors.get(gate);
    }
    return null;
  }
  // Check if a gate has an error
  has_quantum_error(gate, qubits) {
    return this.get_quantum_error(gate, qubits) !== null;
  }
  // Get all readout errors
  get_readout_errors() {
    return this._readout_errors;
  }
  copy() {
    const nm = new _NoiseModel();
    nm._local_quantum_errors = new Map(this._local_quantum_errors);
    nm._basis_gate_errors = new Map(this._basis_gate_errors);
    nm._readout_errors = this._readout_errors.slice();
    return nm;
  }
  is_empty() {
    return this._local_quantum_errors.size === 0 && this._basis_gate_errors.size === 0 && this._readout_errors.length === 0;
  }
};
function _generatePauliGroup(numQubits) {
  const paulis1 = [ComplexMatrix2.identity(2), PAULI.X, PAULI.Y, PAULI.Z];
  if (numQubits === 1) return paulis1;
  let result = paulis1;
  for (let i = 1; i < numQubits; i++) {
    const newResult = [];
    for (const p1 of result) {
      for (const p2 of paulis1) {
        newResult.push(p1.tensor(p2));
      }
    }
    result = newResult;
  }
  return result;
}
function combine_errors(...errors) {
  if (errors.length === 0) return new QuantumError([]);
  let result = errors[0];
  for (let i = 1; i < errors.length; i++) {
    result = result.compose(errors[i]);
  }
  return result;
}

// src/transpiler/transpiler.js
var DEFAULT_BASIS = ["cx", "id", "rz", "sx", "x"];
var DECOMP_RULES = {
  "h": () => [
    ["rz", [0], [Math.PI / 2]],
    ["sx", [0], []],
    ["rz", [0], [Math.PI / 2]]
  ],
  "y": () => [
    ["sx", [0], []],
    ["sx", [0], []],
    ["rz", [0], [Math.PI]]
  ],
  "z": () => [["rz", [0], [Math.PI]]],
  "s": () => [["rz", [0], [Math.PI / 2]]],
  "sdg": () => [["rz", [0], [-Math.PI / 2]]],
  "t": () => [["rz", [0], [Math.PI / 4]]],
  "tdg": () => [["rz", [0], [-Math.PI / 4]]],
  "sx": () => null,
  "sxdg": () => [
    ["sx", [0], []],
    ["x", [0], []]
  ],
  "rx": (params) => [
    ["h", [0], []],
    ["rz", [0], [params[0]]],
    ["h", [0], []]
  ],
  "ry": (params) => [
    // RY(theta) = RZ(-pi/2) H RZ(theta) H RZ(pi/2).
    ["rz", [0], [Math.PI / 2]],
    ["h", [0], []],
    ["rz", [0], [params[0]]],
    ["h", [0], []],
    ["rz", [0], [-Math.PI / 2]]
  ],
  "p": (params) => [["rz", [0], [params[0]]]],
  "u1": (params) => [["rz", [0], [params[0]]]],
  "u2": (params) => [
    // U2(phi, lambda) = RZ(phi) H RZ(lambda) up to global phase.
    ["rz", [0], [params[1]]],
    ["h", [0], []],
    ["rz", [0], [params[0]]]
  ],
  "u3": (params) => [
    // U3(theta, phi, lambda) = RZ(phi) H RZ(theta) H RZ(lambda) up to global phase.
    ["rz", [0], [params[2]]],
    ["h", [0], []],
    ["rz", [0], [params[0]]],
    ["h", [0], []],
    ["rz", [0], [params[1]]]
  ],
  "u": (params) => [
    ["rz", [0], [params[2]]],
    ["h", [0], []],
    ["rz", [0], [params[0]]],
    ["h", [0], []],
    ["rz", [0], [params[1]]]
  ],
  "cy": () => [
    ["sdg", [1], []],
    ["cx", [0, 1], []],
    ["s", [1], []]
  ],
  "cz": () => [
    ["h", [1], []],
    ["cx", [0, 1], []],
    ["h", [1], []]
  ],
  "ch": () => [
    // CH = S(target) H(target) CX(control, target) H(target) Sdg(target).
    ["s", [1], []],
    ["h", [1], []],
    ["cx", [0, 1], []],
    ["h", [1], []],
    ["sdg", [1], []]
  ],
  "swap": () => [
    ["cx", [0, 1], []],
    ["cx", [1, 0], []],
    ["cx", [0, 1], []]
  ],
  "iswap": () => [
    // iSWAP = S(0) S(1) H(0) CX(0,1) CX(1,0) H(1) (verified up to global phase).
    ["s", [0], []],
    ["s", [1], []],
    ["h", [0], []],
    ["cx", [0, 1], []],
    ["cx", [1, 0], []],
    ["h", [1], []]
  ],
  "ccx": () => [
    ["h", [2], []],
    ["cx", [1, 2], []],
    ["tdg", [2], []],
    ["cx", [0, 2], []],
    ["t", [2], []],
    ["cx", [1, 2], []],
    ["tdg", [2], []],
    ["cx", [0, 2], []],
    ["t", [1], []],
    ["t", [2], []],
    ["cx", [0, 1], []],
    ["h", [2], []],
    ["t", [0], []],
    ["tdg", [1], []],
    ["cx", [0, 1], []]
  ],
  "cswap": () => [
    ["cx", [2, 1], []],
    ["ccx", [0, 1, 2], []],
    ["cx", [2, 1], []]
  ],
  "crx": (params) => [
    // CRX(theta) = RZ(-pi/2) RY(theta/2) CX RY(-theta/2) CX RZ(pi/2) on target.
    ["rz", [1], [-Math.PI / 2]],
    ["ry", [1], [params[0] / 2]],
    ["cx", [0, 1], []],
    ["ry", [1], [-params[0] / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [Math.PI / 2]]
  ],
  "cry": (params) => [
    // CRY(theta) = RY(theta/2) CX RY(-theta/2) CX on target.
    ["ry", [1], [params[0] / 2]],
    ["cx", [0, 1], []],
    ["ry", [1], [-params[0] / 2]],
    ["cx", [0, 1], []]
  ],
  "crz": (params) => [
    // CRZ(theta) = RZ(theta/2) CX RZ(-theta/2) CX on target.
    ["rz", [1], [params[0] / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [-params[0] / 2]],
    ["cx", [0, 1], []]
  ],
  "cp": (params) => [
    // CP(theta) = RZ(theta/2) on control, CX, RZ(-theta/2) on target, CX, RZ(theta/2) on target.
    ["rz", [0], [params[0] / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [-params[0] / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [params[0] / 2]]
  ],
  "rzz": (params) => [
    ["cx", [0, 1], []],
    ["rz", [1], [params[0]]],
    ["cx", [0, 1], []]
  ],
  "rxx": (params) => [
    ["h", [0], []],
    ["h", [1], []],
    ["cx", [0, 1], []],
    ["rz", [1], [params[0]]],
    ["cx", [0, 1], []],
    ["h", [1], []],
    ["h", [0], []]
  ],
  "ryy": (params) => [
    // RYY(theta) = RY(-pi/2) on both, CX, RZ(theta) on target, CX, RY(pi/2) on both.
    ["ry", [0], [-Math.PI / 2]],
    ["ry", [1], [-Math.PI / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [params[0]]],
    ["cx", [0, 1], []],
    ["ry", [0], [Math.PI / 2]],
    ["ry", [1], [Math.PI / 2]]
  ]
};
function decomposeGate(gateName, params, basis) {
  basis = basis || DEFAULT_BASIS;
  if (basis.indexOf(gateName) !== -1) {
    return [[gateName, [0], params]];
  }
  const rule = DECOMP_RULES[gateName];
  if (!rule) return null;
  const expansion = rule(params || []);
  if (expansion == null) return null;
  const result = [];
  for (const [name, qubits, p] of expansion) {
    if (basis.indexOf(name) !== -1) {
      result.push([name, qubits, p]);
    } else {
      const sub = decomposeGate(name, p, basis);
      if (!sub) return null;
      for (const [subName, subQubits, subP] of sub) {
        const mappedQubits = subQubits.map((qi) => qubits[qi]);
        result.push([subName, mappedQubits, subP]);
      }
    }
  }
  return result;
}
function transpile(circuits2, options = {}) {
  const basis = options.basis_gates || DEFAULT_BASIS;
  const optimizationLevel = options.optimization_level !== void 0 ? options.optimization_level : 1;
  const single = !Array.isArray(circuits2) || circuits2 instanceof QuantumCircuit;
  const list = single ? [circuits2] : circuits2;
  const transpiled = list.map((c) => _transpileOne(c, basis, optimizationLevel, options));
  return single ? transpiled[0] : transpiled;
}
function _transpileOne(circuit, basis, optLevel, options) {
  const out = new QuantumCircuit();
  const qubitMap = /* @__PURE__ */ new Map();
  const clbitMap = /* @__PURE__ */ new Map();
  if (circuit.qregs.length === 0) {
    const newReg = new QuantumRegister(circuit.num_qubits, "q");
    out.add_register(newReg);
    for (let i = 0; i < circuit.num_qubits; i++) qubitMap.set(circuit.qubits[i], newReg._bits[i]);
  } else {
    for (const r4 of circuit.qregs) {
      const newReg = new QuantumRegister(r4.size, r4.name);
      out.add_register(newReg);
      for (let i = 0; i < r4.size; i++) qubitMap.set(r4._bits[i], newReg._bits[i]);
    }
  }
  if (circuit.cregs.length > 0) {
    for (const r4 of circuit.cregs) {
      const newReg = new ClassicalRegister(r4.size, r4.name);
      out.add_register(newReg);
      for (let i = 0; i < r4.size; i++) clbitMap.set(r4._bits[i], newReg._bits[i]);
    }
  } else if (circuit.num_clbits > 0) {
    const newReg = new ClassicalRegister(circuit.num_clbits, "c");
    out.add_register(newReg);
    for (let i = 0; i < circuit.num_clbits; i++) clbitMap.set(circuit.clbits[i], newReg._bits[i]);
  }
  out.name = circuit.name + "_transpiled";
  out.global_phase = circuit.global_phase;
  const remapQubits = (qs) => qs.map((q) => qubitMap.get(q) || q);
  const remapClbits = (cs) => cs.map((c) => clbitMap.get(c) || c);
  for (const ci of circuit.data) {
    const opName = ci.operation.name;
    const newQubits = remapQubits(ci.qubits);
    const newClbits = remapClbits(ci.clbits);
    if (opName === "barrier") {
      out.barrier(newQubits);
      continue;
    }
    if (opName === "measure") {
      out.append(ci.operation.copy(), newQubits, newClbits);
      continue;
    }
    if (opName === "reset") {
      out.append(ci.operation.copy(), newQubits);
      continue;
    }
    if (opName === "delay") {
      out.append(ci.operation.copy(), newQubits);
      continue;
    }
    if (opName === "if_else" || opName === "while_loop") {
      out.data.push(new CircuitInstruction(ci.operation.copy(), [], []));
      continue;
    }
    const params = ci.operation.params || [];
    if (basis.indexOf(opName) !== -1) {
      out.append(ci.operation.copy(), newQubits);
      continue;
    }
    const expansion = decomposeGate(opName, params, basis);
    if (!expansion) {
      out.append(ci.operation.copy(), newQubits);
      continue;
    }
    for (const [name, qubitOffsets, p] of expansion) {
      const mappedQubits = qubitOffsets.map((qi) => newQubits[qi]);
      const gate = _makeBasicGate(name, p);
      out.append(gate, mappedQubits);
    }
  }
  if (optLevel >= 1) _removeAdjacentInverseGates(out);
  if (optLevel >= 2) _commuteCZ(out);
  return out;
}
function _makeBasicGate(name, params) {
  params = params || [];
  if (standard_gates_exports) {
    const factories = {
      "h": makeHGate,
      "x": makeXGate,
      "y": makeYGate,
      "z": makeZGate,
      "s": makeSGate,
      "sdg": makeSdgGate,
      "t": makeTGate,
      "tdg": makeTdgGate,
      "sx": makeSxGate,
      "sxdg": makeSxdgGate,
      "id": makeIGate,
      "cx": makeCXGate,
      "cy": makeCYGate,
      "cz": makeCZGate,
      "ch": makeCHGate,
      "csx": makeCSXGate,
      "swap": makeSwapGate,
      "iswap": makeISwapGate,
      "dcx": makeDCXGate,
      "ccx": makeCCXGate,
      "cswap": makeCSwapGate
    };
    const factory = factories[name];
    if (factory) return factory();
  }
  if (generalized_gates_exports) {
    const paramFactories = {
      "rx": makeRXGate,
      "ry": makeRYGate,
      "rz": makeRZGate,
      "p": makePGate,
      "u1": makeU1Gate,
      "u2": makeU2Gate,
      "u3": makeU3Gate,
      "u": makeUGate,
      "rxx": makeRXXGate,
      "ryy": makeRYYGate,
      "rzz": makeRZZGate,
      "rzx": makeRZXGate
    };
    const factory = paramFactories[name];
    if (factory) return factory.apply(null, params);
  }
  if (generalized_gates_exports) {
    const controlledBases = {
      "cp": makePGate,
      "crx": makeRXGate,
      "cry": makeRYGate,
      "crz": makeRZGate,
      "cu1": makeU1Gate,
      "cu3": makeU3Gate,
      "cu": makeUGate
    };
    const baseFactory = controlledBases[name];
    if (baseFactory) {
      const base = baseFactory.apply(null, params);
      const cg = new ControlledGate(base, 1);
      cg.name = name;
      return cg;
    }
  }
  const numQubits = name === "cx" || name === "cy" || name === "cz" || name === "ch" || name === "swap" || name === "iswap" || name === "cp" || name === "crx" || name === "cry" || name === "crz" || name === "rxx" || name === "ryy" || name === "rzz" || name === "rzx" || name === "csx" || name === "cu" || name === "cu1" || name === "cu3" || name === "dcx" ? 2 : name === "ccx" || name === "cswap" ? 3 : 1;
  return new Gate(name, numQubits, params);
}
function _removeAdjacentInverseGates(circuit) {
  const newData = [];
  for (const ci of circuit.data) {
    if (newData.length === 0) {
      newData.push(ci);
      continue;
    }
    const prev = newData[newData.length - 1];
    if (prev.qubits.length === ci.qubits.length && prev.qubits.every((q, i) => q === ci.qubits[i]) && prev.operation.name === ci.operation.name + "_dg") {
      newData.pop();
      continue;
    }
    if (prev.qubits.length === ci.qubits.length && prev.qubits.every((q, i) => q === ci.qubits[i]) && ci.operation.name === prev.operation.name + "_dg") {
      newData.pop();
      continue;
    }
    const selfInverse = ["h", "x", "y", "z", "cx", "cz", "swap"];
    if (selfInverse.indexOf(prev.operation.name) !== -1 && prev.operation.name === ci.operation.name && prev.qubits.length === ci.qubits.length && prev.qubits.every((q, i) => q === ci.qubits[i])) {
      newData.pop();
      continue;
    }
    newData.push(ci);
  }
  circuit.data = newData;
}
function _commuteCZ(circuit) {
  return circuit;
}
var PassManager = class {
  constructor() {
    this.passes = [];
  }
  append(pass) {
    this.passes.push(pass);
  }
  run(circuit) {
    let c = circuit;
    for (const pass of this.passes) c = pass(c);
    return c;
  }
};
var PassManagerConfig = class {
  constructor(kwargs) {
    this.basis_gates = kwargs.basis_gates || DEFAULT_BASIS;
    this.optimization_level = kwargs.optimization_level || 0;
    this.initial_layout = kwargs.initial_layout || null;
    this.routing_method = kwargs.routing_method || "basic";
  }
};
function presetPassManager(optimizationLevel = 1, backend = null, basisGates = null) {
  const pm = new PassManager();
  const basis = basisGates || backend && backend.basis_gates || DEFAULT_BASIS;
  pm.append((c) => transpile(c, { basis_gates: basis, optimization_level: optimizationLevel }));
  return pm;
}

// src/transpiler/layout.js
var Layout = class _Layout {
  constructor(mapping = null) {
    this._p2v = /* @__PURE__ */ new Map();
    this._v2p = /* @__PURE__ */ new Map();
    if (mapping) {
      for (const [k, v] of Object.entries(mapping)) {
        const phys = parseInt(k, 10);
        if (!isNaN(phys)) {
          this.setPhysical(phys, v);
        } else {
          this.setVirtual(k, parseInt(v, 10));
        }
      }
    }
  }
  static trivial(numQubits) {
    const layout = new _Layout();
    for (let i = 0; i < numQubits; i++) {
      layout.setPhysical(i, i);
    }
    return layout;
  }
  static from_dict(d) {
    const layout = new _Layout();
    for (const [k, v] of Object.entries(d)) {
      layout.setPhysical(parseInt(k, 10), v);
    }
    return layout;
  }
  static generate_trivial(numQubits, qubits) {
    return _Layout.trivial(numQubits);
  }
  static generate_from_integers(intList) {
    const layout = new _Layout();
    intList.forEach((v, p) => layout.setPhysical(p, v));
    return layout;
  }
  setPhysical(phys, virtual) {
    if (this._p2v.has(phys)) {
      const oldVirt = this._p2v.get(phys);
      this._v2p.delete(oldVirt);
    }
    if (this._v2p.has(virtual)) {
      const oldPhys = this._v2p.get(virtual);
      this._p2v.delete(oldPhys);
    }
    this._p2v.set(phys, virtual);
    this._v2p.set(virtual, phys);
  }
  setVirtual(virtual, phys) {
    this.setPhysical(phys, virtual);
  }
  getPhysical(phys) {
    return this._p2v.get(phys);
  }
  getVirtual(virtual) {
    return this._v2p.get(virtual);
  }
  // physical -> virtual mapping (object form)
  get_physical_bits() {
    const out = {};
    for (const [k, v] of this._p2v.entries()) out[k] = v;
    return out;
  }
  get_virtual_bits() {
    const out = {};
    for (const [k, v] of this._v2p.entries()) out[k] = v;
    return out;
  }
  // Get all physical qubits in use
  get_physical_bits_list() {
    return Array.from(this._p2v.keys()).sort((a, b) => a - b);
  }
  // Add an ancilla qubit
  add(phys, virtual) {
    this.setPhysical(phys, virtual);
  }
  // Swap two physical qubits' mappings
  swap(phys1, phys2) {
    const v1 = this._p2v.get(phys1);
    const v2 = this._p2v.get(phys2);
    this.setPhysical(phys1, v2);
    this.setPhysical(phys2, v1);
  }
  // Combine two layouts
  combine(other) {
    const result = new _Layout();
    for (const [p, v] of this._p2v.entries()) result.setPhysical(p, v);
    for (const [p, v] of other._p2v.entries()) result.setPhysical(p, v);
    return result;
  }
  copy() {
    const result = new _Layout();
    for (const [p, v] of this._p2v.entries()) result.setPhysical(p, v);
    return result;
  }
  size() {
    return this._p2v.size;
  }
  to_dict() {
    return this.get_physical_bits();
  }
  toString() {
    const pairs = [];
    for (const [p, v] of this._p2v.entries()) pairs.push(`${p}: ${v}`);
    return `Layout({${pairs.join(", ")}})`;
  }
};
var CouplingMap = class _CouplingMap {
  constructor(connections = null, numQubits = null) {
    this._edges = [];
    this._graph = /* @__PURE__ */ new Map();
    this._numQubits = numQubits;
    if (connections) {
      if (typeof connections === "string") {
        const parts = connections.split(",");
        for (const part of parts) {
          const [a, b] = part.trim().split("-").map((s) => parseInt(s.trim(), 10));
          this.addEdge(a, b);
        }
      } else if (Array.isArray(connections)) {
        for (const [a, b] of connections) this.addEdge(a, b);
      } else if (connections instanceof Map) {
        for (const [k, neighbors] of connections.entries()) {
          for (const n of neighbors) this.addEdge(k, n);
        }
      } else {
        throw new TypeError("CouplingMap: invalid connections argument");
      }
    }
  }
  static fromLine(numQubits) {
    const c = new _CouplingMap(null, numQubits);
    for (let i = 0; i + 1 < numQubits; i++) {
      c.addEdge(i, i + 1);
      c.addEdge(i + 1, i);
    }
    return c;
  }
  static fromRing(numQubits) {
    const c = new _CouplingMap(null, numQubits);
    for (let i = 0; i < numQubits; i++) {
      c.addEdge(i, (i + 1) % numQubits);
      c.addEdge((i + 1) % numQubits, i);
    }
    return c;
  }
  static fromGrid(rows, cols) {
    const n = rows * cols;
    const c = new _CouplingMap(null, n);
    for (let r4 = 0; r4 < rows; r4++) {
      for (let col = 0; col < cols; col++) {
        const i = r4 * cols + col;
        if (col + 1 < cols) {
          c.addEdge(i, i + 1);
          c.addEdge(i + 1, i);
        }
        if (r4 + 1 < rows) {
          c.addEdge(i, i + cols);
          c.addEdge(i + cols, i);
        }
      }
    }
    return c;
  }
  static fullyConnected(numQubits) {
    const c = new _CouplingMap(null, numQubits);
    for (let i = 0; i < numQubits; i++) {
      for (let j = 0; j < numQubits; j++) {
        if (i !== j) c.addEdge(i, j);
      }
    }
    return c;
  }
  addEdge(a, b) {
    this._edges.push([a, b]);
    if (!this._graph.has(a)) this._graph.set(a, /* @__PURE__ */ new Set());
    if (!this._graph.has(b)) this._graph.set(b, /* @__PURE__ */ new Set());
    this._graph.get(a).add(b);
    if (this._numQubits == null) {
      this._numQubits = Math.max(a, b) + 1;
    } else {
      this._numQubits = Math.max(this._numQubits, a + 1, b + 1);
    }
  }
  get edges() {
    return this._edges.slice();
  }
  get size() {
    return this._numQubits;
  }
  get num_qubits() {
    return this._numQubits;
  }
  neighbors(node) {
    return Array.from(this._graph.get(node) || []);
  }
  hasEdge(a, b) {
    return this._graph.has(a) && this._graph.get(a).has(b);
  }
  is_connected() {
    if (this._numQubits <= 1) return true;
    const visited = /* @__PURE__ */ new Set([0]);
    const queue = [0];
    while (queue.length > 0) {
      const n = queue.shift();
      for (const m of this.neighbors(n)) {
        if (!visited.has(m)) {
          visited.add(m);
          queue.push(m);
        }
      }
    }
    return visited.size === this._numQubits;
  }
  // Shortest path between two physical qubits (BFS)
  shortestPath(a, b) {
    if (a === b) return [a];
    const visited = /* @__PURE__ */ new Set([a]);
    const queue = [[a]];
    while (queue.length > 0) {
      const path = queue.shift();
      const node = path[path.length - 1];
      for (const next of this.neighbors(node)) {
        if (visited.has(next)) continue;
        if (next === b) return path.concat([next]);
        visited.add(next);
        queue.push(path.concat([next]));
      }
    }
    return null;
  }
  // Distance between two physical qubits
  distance(a, b) {
    const path = this.shortestPath(a, b);
    return path ? path.length - 1 : Infinity;
  }
  // Get the undirected graph (each edge once)
  get_undirected_edges() {
    const seen = /* @__PURE__ */ new Set();
    const result = [];
    for (const [a, b] of this._edges) {
      const key = a < b ? `${a}-${b}` : `${b}-${a}`;
      if (!seen.has(key)) {
        seen.add(key);
        result.push([Math.min(a, b), Math.max(a, b)]);
      }
    }
    return result;
  }
  copy() {
    const c = new _CouplingMap(null, this._numQubits);
    for (const [a, b] of this._edges) c.addEdge(a, b);
    return c;
  }
  to_dict() {
    const out = {};
    for (const [k, neighbors] of this._graph.entries()) {
      out[k] = Array.from(neighbors);
    }
    return out;
  }
};
var AnalysisPass = class {
  constructor() {
    this.requires = [];
    this.preserves = [];
    this.analysis_name = this.constructor.name;
    this.property_set = {};
  }
  run(dag) {
    throw new Error("AnalysisPass.run not implemented");
  }
};
var TransformationPass = class {
  constructor() {
    this.requires = [];
    this.preserves = [];
    this.transformation_name = this.constructor.name;
  }
  run(dag) {
    throw new Error("TransformationPass.run not implemented");
  }
};

// src/dagcircuit/dagcircuit.js
var _nodeCounter = 0;
var DAGNode = class {
  constructor() {
    this._id = ++_nodeCounter;
    this._node_id = this._id;
    this.name = "";
    this.wires = [];
  }
  is_op_node() {
    return false;
  }
  is_in_node() {
    return false;
  }
  is_out_node() {
    return false;
  }
};
var DAGOpNode = class _DAGOpNode extends DAGNode {
  constructor(op, qargs = [], cargs = []) {
    super();
    this.op = op;
    this.qargs = qargs;
    this.cargs = cargs;
    this.name = op.name;
    this.wires = qargs.concat(cargs);
  }
  is_op_node() {
    return true;
  }
  semantic_eq(other) {
    return other instanceof _DAGOpNode && this.op.name === other.op.name && this.qargs.length === other.qargs.length && this.qargs.every((q, i) => q.equals(other.qargs[i]));
  }
};
var DAGInNode = class extends DAGNode {
  constructor(wire) {
    super();
    this.wire = wire;
    this.name = wire.toString();
    this.wires = [wire];
  }
  is_in_node() {
    return true;
  }
};
var DAGOutNode = class extends DAGNode {
  constructor(wire) {
    super();
    this.wire = wire;
    this.name = wire.toString();
    this.wires = [wire];
  }
  is_out_node() {
    return true;
  }
};
var DAGRegister = class {
  constructor(name, bits, bitClass) {
    this.name = name;
    this.bits = bits;
    this.bitClass = bitClass;
  }
};
var DAGCircuit = class _DAGCircuit {
  constructor() {
    this.qregs = /* @__PURE__ */ new Map();
    this.cregs = /* @__PURE__ */ new Map();
    this.qubits = [];
    this.clbits = [];
    this._input_nodes = /* @__PURE__ */ new Map();
    this._output_nodes = /* @__PURE__ */ new Map();
    this._multi_graph = /* @__PURE__ */ new Map();
    this._wire_to_in_edges = /* @__PURE__ */ new Map();
    this._wire_to_out_edges = /* @__PURE__ */ new Map();
    this.global_phase = 0;
    this.name = "dag";
    this.metadata = null;
    this.calibrations = [];
    this.duration = null;
    this.unit = "dt";
  }
  add_qreg(register) {
    if (!(register instanceof QuantumRegister)) {
      throw new TypeError("add_qreg requires a QuantumRegister");
    }
    if (this.qregs.has(register.name)) {
      throw new Error(`Quantum register ${register.name} already exists`);
    }
    const dagReg = new DAGRegister(register.name, register._bits.slice(), Qubit);
    this.qregs.set(register.name, dagReg);
    for (const q of register._bits) {
      this.qubits.push(q);
      const inNode = new DAGInNode(q);
      const outNode = new DAGOutNode(q);
      this._input_nodes.set(q, inNode);
      this._output_nodes.set(q, outNode);
      this._multi_graph.set(inNode, { successors: /* @__PURE__ */ new Set([outNode]), predecessors: /* @__PURE__ */ new Set() });
      this._multi_graph.set(outNode, { successors: /* @__PURE__ */ new Set(), predecessors: /* @__PURE__ */ new Set([inNode]) });
      this._wire_to_in_edges.set(q, [{ from: inNode, to: outNode }]);
      this._wire_to_out_edges.set(q, []);
    }
    return dagReg;
  }
  add_creg(register) {
    if (!(register instanceof ClassicalRegister)) {
      throw new TypeError("add_creg requires a ClassicalRegister");
    }
    if (this.cregs.has(register.name)) {
      throw new Error(`Classical register ${register.name} already exists`);
    }
    const dagReg = new DAGRegister(register.name, register._bits.slice(), Clbit);
    this.cregs.set(register.name, dagReg);
    for (const c of register._bits) {
      this.clbits.push(c);
      const inNode = new DAGInNode(c);
      const outNode = new DAGOutNode(c);
      this._input_nodes.set(c, inNode);
      this._output_nodes.set(c, outNode);
      this._multi_graph.set(inNode, { successors: /* @__PURE__ */ new Set([outNode]), predecessors: /* @__PURE__ */ new Set() });
      this._multi_graph.set(outNode, { successors: /* @__PURE__ */ new Set(), predecessors: /* @__PURE__ */ new Set([inNode]) });
      this._wire_to_in_edges.set(c, [{ from: inNode, to: outNode }]);
      this._wire_to_out_edges.set(c, []);
    }
    return dagReg;
  }
  get num_qubits() {
    return this.qubits.length;
  }
  get num_clbits() {
    return this.clbits.length;
  }
  // Apply an operation to qubits and clbits
  apply_operation(op, qargs = [], cargs = []) {
    if (!Array.isArray(qargs)) qargs = [qargs];
    if (!Array.isArray(cargs)) cargs = [cargs];
    const node = new DAGOpNode(op, qargs, cargs);
    this._multi_graph.set(node, { successors: /* @__PURE__ */ new Set(), predecessors: /* @__PURE__ */ new Set() });
    for (const wire of qargs.concat(cargs)) {
      const inEdges = this._wire_to_in_edges.get(wire);
      if (!inEdges) {
        throw new Error(`Wire ${wire} not in DAG`);
      }
      const lastEdge = inEdges[inEdges.length - 1];
      const fromNode = lastEdge.from;
      const toNode = lastEdge.to;
      this._removeEdge(fromNode, toNode);
      this._addEdge(fromNode, node);
      this._addEdge(node, toNode);
      inEdges[inEdges.length - 1] = { from: node, to: toNode };
      this._wire_to_out_edges.get(wire).push({ from: fromNode, to: node });
    }
    return node;
  }
  _addEdge(from, to) {
    this._multi_graph.get(from).successors.add(to);
    this._multi_graph.get(to).predecessors.add(from);
  }
  _removeEdge(from, to) {
    this._multi_graph.get(from).successors.delete(to);
    this._multi_graph.get(to).predecessors.delete(from);
  }
  // Topologically iterate over operation nodes
  topological_op_nodes() {
    const inDegree = /* @__PURE__ */ new Map();
    const queue = [];
    for (const [node, links] of this._multi_graph.entries()) {
      let deg = 0;
      for (const pred of links.predecessors) {
        if (pred.is_op_node()) deg++;
      }
      inDegree.set(node, deg);
      if (deg === 0 && node.is_op_node()) queue.push(node);
    }
    const result = [];
    while (queue.length > 0) {
      const node = queue.shift();
      if (node.is_op_node()) result.push(node);
      const links = this._multi_graph.get(node);
      for (const succ of links.successors) {
        const newDeg = inDegree.get(succ) - 1;
        inDegree.set(succ, newDeg);
        if (newDeg === 0 && succ.is_op_node()) queue.push(succ);
      }
    }
    return result;
  }
  // Get all operation nodes (unsorted)
  op_nodes() {
    return Array.from(this._multi_graph.keys()).filter((n) => n.is_op_node());
  }
  // Get all nodes (input, output, ops)
  all_nodes() {
    return Array.from(this._multi_graph.keys());
  }
  // Count operations by name
  count_ops() {
    const counts = {};
    for (const node of this.op_nodes()) {
      counts[node.op.name] = (counts[node.op.name] || 0) + 1;
    }
    return counts;
  }
  // Depth (longest path through op nodes)
  depth() {
    const opNodes = this.topological_op_nodes();
    const longestPath = /* @__PURE__ */ new Map();
    for (const node of opNodes) {
      let maxPred = 0;
      const links = this._multi_graph.get(node);
      for (const pred of links.predecessors) {
        if (pred.is_op_node()) {
          const predLen = longestPath.get(pred) || 0;
          if (predLen > maxPred) maxPred = predLen;
        }
      }
      longestPath.set(node, maxPred + 1);
    }
    let max = 0;
    for (const v of longestPath.values()) if (v > max) max = v;
    return max;
  }
  // Get edges on a wire
  edges_on_wire(wire) {
    const result = [];
    const inEdges = this._wire_to_in_edges.get(wire);
    if (!inEdges) return result;
    let current = this._input_nodes.get(wire);
    const outNode = this._output_nodes.get(wire);
    while (current !== outNode) {
      const links = this._multi_graph.get(current);
      let next = null;
      for (const succ of links.successors) {
        if (succ.wires.includes(wire) || succ === outNode) {
          next = succ;
          break;
        }
      }
      if (!next) break;
      result.push({ from: current, to: next });
      current = next;
    }
    return result;
  }
  // Remove an operation node
  remove_op_node(node) {
    if (!node.is_op_node()) throw new Error("remove_op_node requires an op node");
    const links = this._multi_graph.get(node);
    const preds = Array.from(links.predecessors);
    const succs = Array.from(links.successors);
    for (const wire of node.wires) {
      let pred = null, succ = null;
      for (const p of preds) if (p.wires.includes(wire)) {
        pred = p;
        break;
      }
      for (const s of succs) if (s.wires.includes(wire)) {
        succ = s;
        break;
      }
      if (pred && succ) {
        this._removeEdge(pred, node);
        this._removeEdge(node, succ);
        this._addEdge(pred, succ);
        const inEdges = this._wire_to_in_edges.get(wire);
        for (let i = 0; i < inEdges.length; i++) {
          if (inEdges[i].from === node && inEdges[i].to === succ) {
            inEdges[i] = { from: pred, to: succ };
            break;
          }
        }
      }
    }
    this._multi_graph.delete(node);
  }
  // Substitute one node with another (preserves edges)
  substitute_node(old_node, new_node) {
    if (!old_node.is_op_node()) throw new Error("substitute_node: old must be op node");
    const links = this._multi_graph.get(old_node);
    this._multi_graph.set(new_node, {
      successors: new Set(links.successors),
      predecessors: new Set(links.predecessors)
    });
    for (const pred of links.predecessors) {
      const predLinks = this._multi_graph.get(pred);
      predLinks.successors.delete(old_node);
      predLinks.successors.add(new_node);
    }
    for (const succ of links.successors) {
      const succLinks = this._multi_graph.get(succ);
      succLinks.predecessors.delete(old_node);
      succLinks.predecessors.add(new_node);
    }
    this._multi_graph.delete(old_node);
  }
  // Replace a node with a subcircuit (DAGCircuit)
  substitute_node_with_dag(node, sub_dag, wires_map = null) {
    if (!node.is_op_node()) throw new Error("substitute_node_with_dag: old must be op node");
    const subQubits = sub_dag.qubits;
    const subClbits = sub_dag.clbits;
    const defaultMap = /* @__PURE__ */ new Map();
    if (!wires_map) {
      for (let i = 0; i < subQubits.length; i++) defaultMap.set(subQubits[i], node.qargs[i]);
      for (let i = 0; i < subClbits.length; i++) defaultMap.set(subClbits[i], node.cargs[i]);
    } else {
      for (const [k, v] of wires_map.entries()) defaultMap.set(k, v);
    }
    const predLinks = this._multi_graph.get(node).predecessors;
    const succLinks = this._multi_graph.get(node).successors;
    const wirePred = /* @__PURE__ */ new Map();
    const wireSucc = /* @__PURE__ */ new Map();
    for (const pred of predLinks) for (const w of pred.wires) wirePred.set(w, pred);
    for (const succ of succLinks) for (const w of succ.wires) wireSucc.set(w, succ);
    for (const w of node.wires) {
      const p = wirePred.get(w);
      const s = wireSucc.get(w);
      if (p && s) {
        this._removeEdge(p, node);
        this._removeEdge(node, s);
      }
    }
    this._multi_graph.delete(node);
    let lastNodeOnWire = new Map(wirePred);
    for (const opNode of sub_dag.topological_op_nodes()) {
      const mappedQargs = opNode.qargs.map((q) => defaultMap.get(q) || q);
      const mappedCargs = opNode.cargs.map((c) => defaultMap.get(c) || c);
      const newNode = new DAGOpNode(opNode.op, mappedQargs, mappedCargs);
      this._multi_graph.set(newNode, { successors: /* @__PURE__ */ new Set(), predecessors: /* @__PURE__ */ new Set() });
      for (const w of mappedQargs.concat(mappedCargs)) {
        const pred = lastNodeOnWire.get(w) || this._input_nodes.get(w);
        this._addEdge(pred, newNode);
        lastNodeOnWire.set(w, newNode);
      }
    }
    for (const w of node.wires) {
      const last = lastNodeOnWire.get(w) || wirePred.get(w);
      const succ = wireSucc.get(w);
      if (last && succ) this._addEdge(last, succ);
    }
  }
  // Get node predecessors
  predecessors(node) {
    const links = this._multi_graph.get(node);
    return links ? Array.from(links.predecessors) : [];
  }
  // Get node successors
  successors(node) {
    const links = this._multi_graph.get(node);
    return links ? Array.from(links.successors) : [];
  }
  // Get quantum register by name
  qreg(name) {
    return this.qregs.get(name);
  }
  creg(name) {
    return this.cregs.get(name);
  }
  // Copy the DAG
  copy() {
    const newDag = new _DAGCircuit();
    for (const [name, reg] of this.qregs.entries()) {
      newDag.add_qreg(new QuantumRegister(reg.bits.length, name));
    }
    for (const [name, reg] of this.cregs.entries()) {
      newDag.add_creg(new ClassicalRegister(reg.bits.length, name));
    }
    newDag.name = this.name;
    newDag.global_phase = this.global_phase;
    newDag.metadata = this.metadata;
    for (const node of this.topological_op_nodes()) {
      newDag.apply_operation(node.op.copy(), node.qargs.slice(), node.cargs.slice());
    }
    return newDag;
  }
  // Convert DAG to a list of CircuitInstructions (in topological order)
  to_instructions() {
    return this.topological_op_nodes().map((node) => new CircuitInstruction(
      node.op.copy(),
      node.qargs.slice(),
      node.cargs.slice()
    ));
  }
};

// src/transpiler/sabre.js
var SabreSwap = class {
  constructor(couplingMap, heuristic = "lookahead", seed = null) {
    this.coupling_map = couplingMap;
    this.heuristic = heuristic;
    this.seed = seed;
    this._rng = seed != null ? _makeRng3(seed) : Math.random;
    this.swap_weight = 0.5;
    this.decrement = 1e-3;
    this.max_distance = couplingMap ? couplingMap.size + 1 : 100;
  }
  // Run Sabre routing on a DAGCircuit.
  // Returns a routed QuantumCircuit with SWAP gates inserted where needed.
  run(dag) {
    if (!this.coupling_map) {
      return dag;
    }
    const layout = new Layout();
    for (let i = 0; i < dag.qubits.length; i++) {
      layout.setPhysical(i, dag.qubits[i]);
    }
    const gates = dag.topological_op_nodes();
    const frontLayer = [];
    const executedGates = /* @__PURE__ */ new Set();
    const remainingGates = new Set(gates);
    const successors = /* @__PURE__ */ new Map();
    const predecessors = /* @__PURE__ */ new Map();
    for (const gate of gates) {
      successors.set(gate, []);
      predecessors.set(gate, []);
    }
    for (const gate of gates) {
      for (const succ of dag.successors(gate)) {
        if (succ.is_op_node()) {
          successors.get(gate).push(succ);
          predecessors.get(succ).push(gate);
        }
      }
    }
    this._successorsMap = successors;
    for (const gate of gates) {
      if (predecessors.get(gate).length === 0) {
        frontLayer.push(gate);
      }
    }
    const decay = new Array(this.coupling_map.size).fill(1);
    const distMatrix = this._computeDistanceMatrix();
    const schedule = [];
    let iterations = 0;
    const maxIterations = gates.length * 10;
    while (frontLayer.length > 0 && iterations < maxIterations) {
      iterations++;
      const executeList = [];
      const remainingFront = [];
      for (const gate of frontLayer) {
        if (this._canExecute(gate, layout)) {
          executeList.push(gate);
        } else {
          remainingFront.push(gate);
        }
      }
      if (executeList.length > 0) {
        for (const gate of executeList) {
          schedule.push({ kind: "gate", gate, layout: layout.copy() });
          executedGates.add(gate);
          remainingGates.delete(gate);
          for (const succ of successors.get(gate)) {
            const allDone = predecessors.get(succ).every((p) => executedGates.has(p));
            if (allDone && !frontLayer.includes(succ) && !remainingFront.includes(succ)) {
              remainingFront.push(succ);
            }
          }
        }
        frontLayer.length = 0;
        frontLayer.push(...remainingFront);
        decay.fill(1);
        continue;
      }
      const swapCandidates = this._generateSwapCandidates(frontLayer, layout);
      if (swapCandidates.length === 0) {
        break;
      }
      let bestSwap = null;
      let bestCost = Infinity;
      for (const swap of swapCandidates) {
        const cost = this._swapCost(swap, frontLayer, layout, distMatrix, decay);
        if (cost < bestCost) {
          bestCost = cost;
          bestSwap = swap;
        }
      }
      if (bestSwap) {
        schedule.push({ kind: "swap", physical: [bestSwap[0], bestSwap[1]], layout: layout.copy() });
        layout.swap(bestSwap[0], bestSwap[1]);
        decay[bestSwap[0]] += this.decrement;
        decay[bestSwap[1]] += this.decrement;
      }
    }
    return this._buildRoutedCircuit(dag, schedule, layout);
  }
  // Check if a gate can execute on the current layout
  _canExecute(gate, layout) {
    if (gate.op.num_qubits < 2) return true;
    const physicalQubits = gate.qargs.map((q) => layout.getVirtual(q));
    if (physicalQubits.some((p) => p === void 0)) return false;
    if (gate.op.name === "cx" || gate.op.name === "cy" || gate.op.name === "cz" || gate.op.name === "ch") {
      return this.coupling_map.hasEdge(physicalQubits[0], physicalQubits[1]);
    }
    return true;
  }
  // Generate candidate SWAPs from the front layer
  _generateSwapCandidates(frontLayer, layout) {
    const candidates = [];
    const seen = /* @__PURE__ */ new Set();
    for (const gate of frontLayer) {
      if (gate.op.num_qubits < 2) continue;
      const q1 = layout.getVirtual(gate.qargs[0]);
      const q2 = layout.getVirtual(gate.qargs[1]);
      if (q1 === void 0 || q2 === void 0) continue;
      for (const neighbor of this.coupling_map.neighbors(q1)) {
        const key = [Math.min(q1, neighbor), Math.max(q1, neighbor)].join(",");
        if (!seen.has(key)) {
          seen.add(key);
          candidates.push([q1, neighbor]);
        }
      }
      for (const neighbor of this.coupling_map.neighbors(q2)) {
        const key = [Math.min(q2, neighbor), Math.max(q2, neighbor)].join(",");
        if (!seen.has(key)) {
          seen.add(key);
          candidates.push([q2, neighbor]);
        }
      }
    }
    return candidates;
  }
  // Compute the cost of a SWAP using the Sabre heuristic
  _swapCost(swap, frontLayer, layout, distMatrix, decay) {
    const tempLayout = layout.copy();
    tempLayout.swap(swap[0], swap[1]);
    let totalCost = 0;
    for (const gate of frontLayer) {
      if (gate.op.num_qubits < 2) continue;
      const q1 = tempLayout.getVirtual(gate.qargs[0]);
      const q2 = tempLayout.getVirtual(gate.qargs[1]);
      if (q1 === void 0 || q2 === void 0) continue;
      const dist = distMatrix[q1][q2];
      const weight = Math.max(decay[q1], decay[q2]);
      totalCost += dist * weight;
    }
    if (this.heuristic === "lookahead") {
      const nextLayer = this._getNextLayer(frontLayer);
      let nextCost = 0;
      let count = 0;
      for (const gate of nextLayer) {
        if (gate.op.num_qubits < 2) continue;
        const q1 = tempLayout.getVirtual(gate.qargs[0]);
        const q2 = tempLayout.getVirtual(gate.qargs[1]);
        if (q1 === void 0 || q2 === void 0) continue;
        nextCost += distMatrix[q1][q2];
        count++;
      }
      if (count > 0) totalCost += this.swap_weight * nextCost / count;
    }
    return totalCost;
  }
  // Get the next layer (immediate successors of the current front layer).
  // Used by the lookahead term of the cost function.
  _getNextLayer(frontLayer) {
    const next = /* @__PURE__ */ new Set();
    const succMap = this._successorsMap;
    if (!succMap) return [];
    for (const gate of frontLayer) {
      const succs = succMap.get(gate) || [];
      for (const succ of succs) next.add(succ);
    }
    return Array.from(next);
  }
  // Precompute shortest-path distances between all pairs
  _computeDistanceMatrix() {
    const n = this.coupling_map.size;
    const dist = new Array(n);
    for (let i = 0; i < n; i++) {
      dist[i] = new Array(n);
      for (let j = 0; j < n; j++) {
        dist[i][j] = i === j ? 0 : this.coupling_map.hasEdge(i, j) ? 1 : Infinity;
      }
    }
    for (let k = 0; k < n; k++) {
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          if (dist[i][k] + dist[k][j] < dist[i][j]) {
            dist[i][j] = dist[i][k] + dist[k][j];
          }
        }
      }
    }
    return dist;
  }
  // Build the routed circuit from the recorded schedule.
  // `schedule` is an array of {kind: 'gate'|'swap', ...} entries in execution
  // order. `finalLayout` is the layout at the end of routing (used only for
  // diagnostic info on the returned circuit).
  _buildRoutedCircuit(dag, schedule, finalLayout) {
    const n = dag.qubits.length;
    const routed = new QuantumCircuit();
    routed.add_register(new QuantumRegister(n, "q"));
    if (dag.clbits.length > 0) {
      routed.add_register(new ClassicalRegister(dag.clbits.length, "c"));
    }
    const qubitMap = /* @__PURE__ */ new Map();
    for (let i = 0; i < n; i++) {
      qubitMap.set(dag.qubits[i], routed.qubits[i]);
    }
    const clbitMap = /* @__PURE__ */ new Map();
    for (let i = 0; i < dag.clbits.length; i++) {
      clbitMap.set(dag.clbits[i], routed.clbits[i]);
    }
    for (const entry of schedule) {
      if (entry.kind === "gate") {
        const gate = entry.gate;
        const op = gate.op.copy();
        const qargs = gate.qargs.map((q) => {
          const phys = entry.layout.getVirtual(q);
          return phys === void 0 ? qubitMap.get(q) : routed.qubits[phys];
        });
        const cargs = (gate.cargs || []).map((c) => clbitMap.get(c) || c);
        routed.append(op, qargs, cargs);
      } else if (entry.kind === "swap") {
        const [p0, p1] = entry.physical;
        routed.append(_makeStdGate("SWAP"), [routed.qubits[p0], routed.qubits[p1]]);
      }
    }
    return routed;
  }
};
var SabreLayout = class {
  constructor(couplingMap, seed = null, maxIterations = 4) {
    this.coupling_map = couplingMap;
    this.seed = seed;
    this.max_iterations = maxIterations;
  }
  run(dag) {
    if (!this.coupling_map) {
      const layout = new Layout();
      for (let i = 0; i < dag.qubits.length; i++) layout.setPhysical(i, dag.qubits[i]);
      return layout;
    }
    let bestLayout = null;
    let bestScore = Infinity;
    for (let iter = 0; iter < this.max_iterations; iter++) {
      let initialLayout;
      if (iter === 0) {
        initialLayout = new Layout();
        for (let i = 0; i < dag.qubits.length; i++) {
          initialLayout.setPhysical(i, dag.qubits[i]);
        }
      } else {
        initialLayout = this._randomLayout(dag.qubits, this.coupling_map.size);
      }
      const sabre = new SabreSwap(this.coupling_map, "lookahead", this.seed != null ? this.seed + iter : null);
      const routed = sabre.run(dag);
      const score = routed.data.filter((ci) => ci.operation.name === "swap").length;
      if (score < bestScore) {
        bestScore = score;
        bestLayout = initialLayout;
      }
    }
    return bestLayout;
  }
  _randomLayout(virtualQubits, numPhysical) {
    const layout = new Layout();
    const physical = Array.from({ length: numPhysical }, (_, i) => i);
    for (let i = numPhysical - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [physical[i], physical[j]] = [physical[j], physical[i]];
    }
    for (let i = 0; i < virtualQubits.length; i++) {
      layout.setPhysical(physical[i], virtualQubits[i]);
    }
    return layout;
  }
};
function collect_1q_runs(dag) {
  const runs = [];
  const visited = /* @__PURE__ */ new Set();
  for (const node of dag.topological_op_nodes()) {
    if (visited.has(node)) continue;
    if (node.op.num_qubits !== 1) continue;
    const run = [node];
    visited.add(node);
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter((s) => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.num_qubits !== 1) break;
      if (next.qubits[0] !== current.qubits[0]) break;
      run.push(next);
      visited.add(next);
      current = next;
    }
    if (run.length > 1) runs.push(run);
  }
  return runs;
}
function collect_2q_runs(dag) {
  const runs = [];
  const visited = /* @__PURE__ */ new Set();
  for (const node of dag.topological_op_nodes()) {
    if (visited.has(node)) continue;
    if (node.op.num_qubits !== 2) continue;
    const run = [node];
    visited.add(node);
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter((s) => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.num_qubits !== 2) break;
      if (next.qubits[0] !== current.qubits[0] || next.qubits[1] !== current.qubits[1]) break;
      run.push(next);
      visited.add(next);
      current = next;
    }
    if (run.length > 1) runs.push(run);
  }
  return runs;
}
function consolidate_blocks(dag) {
  const runs1q = collect_1q_runs(dag);
  for (const run of runs1q) {
    let combined = null;
    let ok = true;
    for (const node of run) {
      try {
        const m = node.op.to_matrix();
        combined = combined ? m.mul(combined) : m;
      } catch (e) {
        ok = false;
        break;
      }
    }
    if (!ok || !combined) continue;
    const newGate = new Gate("unitary", 1, [combined]);
    newGate._matrixBuilder = () => combined;
    const newNode = new DAGOpNode(newGate, run[0].qargs.slice(), []);
    dag.substitute_node(run[0], newNode);
    for (let i = 1; i < run.length; i++) {
      try {
        dag.remove_op_node(run[i]);
      } catch (e) {
      }
    }
  }
  return dag;
}
function optimize_cliffords(dag) {
  const nodes = dag.topological_op_nodes();
  const selfInverse = /* @__PURE__ */ new Set(["h", "x", "y", "z", "cx", "cz", "swap"]);
  const toRemove = /* @__PURE__ */ new Set();
  const stack = [];
  for (const node of nodes) {
    if (toRemove.has(node)) continue;
    if (stack.length > 0) {
      const prev = stack[stack.length - 1];
      if (selfInverse.has(prev.op.name) && prev.op.name === node.op.name && prev.qubits.length === node.qubits.length && prev.qubits.every((q, i) => q === node.qubits[i])) {
        stack.pop();
        toRemove.add(prev);
        toRemove.add(node);
        continue;
      }
    }
    stack.push(node);
  }
  for (const node of toRemove) {
    try {
      dag.remove_op_node(node);
    } catch (e) {
    }
  }
  return dag;
}
function remove_diagonal_gates_before_measure(dag) {
  const diagonalGates = ["z", "s", "sdg", "t", "tdg", "p", "u1", "rz", "cz", "cp"];
  const toRemove = [];
  for (const node of dag.topological_op_nodes()) {
    if (!diagonalGates.includes(node.op.name)) continue;
    const succs = dag.successors(node).filter((s) => s.is_op_node());
    if (succs.length > 0 && succs.every((s) => s.op.name === "measure")) {
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.remove_op_node(node);
  }
  return dag;
}
function _makeRng3(seed) {
  let s = seed >>> 0;
  return function() {
    s = s + 1831565813 | 0;
    let t = Math.imul(s ^ s >>> 15, 1 | s);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}
function _makeStdGate(name) {
  const upper = name.toUpperCase();
  switch (upper) {
    case "SWAP":
      return makeSwapGate();
    case "CX":
      return makeCXGate();
    case "CZ":
      return makeCZGate();
    case "H":
      return makeHGate();
    case "X":
      return makeXGate();
    case "Y":
      return makeYGate();
    case "Z":
      return makeZGate();
    case "S":
      return makeSGate();
    case "SDG":
      return makeSdgGate();
    case "T":
      return makeTGate();
    case "TDG":
      return makeTdgGate();
    case "SX":
      return makeSxGate();
    case "I":
      return makeIGate();
    default:
      return new Gate(name.toLowerCase(), 1, []);
  }
}

// src/transpiler/passes.js
var TrivialLayout = class extends TransformationPass {
  constructor(couplingMap = null) {
    super();
    this.coupling_map = couplingMap;
  }
  run(dag) {
    const layout = Layout.trivial(dag.qubits.length);
    dag._layout = layout;
    return dag;
  }
};
var DenseLayout = class extends TransformationPass {
  constructor(couplingMap = null) {
    super();
    this.coupling_map = couplingMap;
  }
  run(dag) {
    if (!this.coupling_map) {
      return new TrivialLayout().run(dag);
    }
    const n = dag.qubits.length;
    const degree = /* @__PURE__ */ new Map();
    for (let p = 0; p < this.coupling_map.size; p++) {
      degree.set(p, this.coupling_map.neighbors(p).length);
    }
    const sorted = Array.from(degree.entries()).sort((a, b) => b[1] - a[1]).slice(0, n).map(([p]) => p);
    const layout = new Layout();
    for (let i = 0; i < n; i++) {
      layout.setPhysical(sorted[i], i);
    }
    dag._layout = layout;
    return dag;
  }
};
var SabreLayout2 = class extends TransformationPass {
  constructor(couplingMap = null, seed = null, maxIterations = 3) {
    super();
    this.coupling_map = couplingMap;
    this.seed = seed;
    this.max_iterations = maxIterations;
  }
  run(dag) {
    if (!this.coupling_map) {
      return new TrivialLayout().run(dag);
    }
    const real = new SabreLayout(this.coupling_map, this.seed, this.max_iterations);
    dag._layout = real.run(dag);
    return dag;
  }
};
var ApplyLayout = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    if (!dag._layout) return dag;
    const layout = dag._layout;
    const maxPhys = Math.max(...Array.from(layout._p2v.keys())) + 1;
    if (maxPhys === dag.qubits.length) return dag;
    const newQubitsNeeded = maxPhys - dag.qubits.length;
    if (newQubitsNeeded > 0) {
      const ancReg = new QuantumRegister(newQubitsNeeded, "ancilla");
      dag.add_qreg(ancReg);
    }
    return dag;
  }
};
var BasicSwap = class extends TransformationPass {
  constructor(couplingMap = null) {
    super();
    this.coupling_map = couplingMap;
  }
  run(dag) {
    if (!this.coupling_map) return dag;
    const newDag = new DAGCircuit();
    for (const [name, reg] of dag.qregs.entries()) {
      newDag.add_qreg(new QuantumRegister(reg.bits.length, name));
    }
    for (const [name, reg] of dag.cregs.entries()) {
      newDag.add_creg(new ClassicalRegister(reg.bits.length, name));
    }
    newDag.name = dag.name;
    newDag.global_phase = dag.global_phase;
    newDag.metadata = dag.metadata;
    for (const node of dag.topological_op_nodes()) {
      const op = node.op;
      const qargs = node.qargs;
      const cargs = node.cargs;
      if (op.num_qubits !== 2) {
        newDag.apply_operation(op.copy(), qargs.map((q) => _mapQubit(q, dag, newDag)), cargs.map((c) => _mapClbit(c, dag, newDag)));
        continue;
      }
      const q0 = dag.qubits.indexOf(qargs[0]);
      const q1 = dag.qubits.indexOf(qargs[1]);
      if (this.coupling_map.hasEdge(q0, q1)) {
        newDag.apply_operation(op.copy(), qargs.map((q) => _mapQubit(q, dag, newDag)), cargs.map((c) => _mapClbit(c, dag, newDag)));
        continue;
      }
      const path = this.coupling_map.shortestPath(q0, q1);
      if (!path || path.length < 2) {
        throw new Error(`BasicSwap: no path between ${q0} and ${q1} in coupling map`);
      }
      const forwardSwaps = [];
      for (let i = 0; i + 1 < path.length; i++) {
        const a = path[i];
        const b = path[i + 1];
        forwardSwaps.push([a, b]);
        const swapOp = SwapGate.copy();
        const qa = newDag.qubits[a];
        const qb = newDag.qubits[b];
        newDag.apply_operation(swapOp, [qa, qb], []);
      }
      const lastIdx = path.length - 1;
      const newQargs = [
        newDag.qubits[path[lastIdx]],
        newDag.qubits[path[lastIdx - 1]]
      ];
      newDag.apply_operation(op.copy(), newQargs, cargs.map((c) => _mapClbit(c, dag, newDag)));
      for (let i = forwardSwaps.length - 1; i >= 0; i--) {
        const [a, b] = forwardSwaps[i];
        const swapOp = SwapGate.copy();
        const qa = newDag.qubits[a];
        const qb = newDag.qubits[b];
        newDag.apply_operation(swapOp, [qa, qb], []);
      }
    }
    return newDag;
  }
};
function _mapQubit(q, oldDag, newDag) {
  const idx = oldDag.qubits.indexOf(q);
  return newDag.qubits[idx];
}
function _mapClbit(c, oldDag, newDag) {
  const idx = oldDag.clbits.indexOf(c);
  return newDag.clbits[idx];
}
var LookaheadSwap = class extends TransformationPass {
  constructor(couplingMap = null, searchDepth = 5) {
    super();
    this.coupling_map = couplingMap;
    this.search_depth = searchDepth;
  }
  run(dag) {
    if (!this.coupling_map) return dag;
    const sabre = new SabreSwap(this.coupling_map, "lookahead", null);
    return sabre.run(dag);
  }
};
var StochasticSwap = class extends TransformationPass {
  constructor(couplingMap = null, seed = null, trials = 20) {
    super();
    this.coupling_map = couplingMap;
    this.seed = seed;
    this.trials = trials;
  }
  run(dag) {
    if (!this.coupling_map) return dag;
    let best = null;
    let bestSwaps = Infinity;
    for (let t = 0; t < this.trials; t++) {
      const sabre = new SabreSwap(this.coupling_map, "lookahead", this.seed != null ? this.seed + t : null);
      const routed = sabre.run(dag);
      const numSwaps = routed.data.filter((ci) => ci.operation.name === "swap").length;
      if (numSwaps < bestSwaps) {
        bestSwaps = numSwaps;
        best = routed;
      }
    }
    return best || dag;
  }
};
var SabreSwap2 = class extends TransformationPass {
  constructor(couplingMap = null, heuristic = "lookahead", seed = null) {
    super();
    this.coupling_map = couplingMap;
    this.heuristic = heuristic;
    this.seed = seed;
  }
  run(dag) {
    if (!this.coupling_map) return dag;
    const sabre = new SabreSwap(this.coupling_map, this.heuristic, this.seed);
    return sabre.run(dag);
  }
};
var CXCancellation = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    const toRemove = /* @__PURE__ */ new Set();
    for (const node of opNodes) {
      if (toRemove.has(node)) continue;
      if (node.op.name !== "cx") continue;
      for (const succ of dag.successors(node)) {
        if (succ.is_op_node() && succ.op.name === "cx" && succ.qargs.length === 2 && succ.qargs[0] === node.qargs[0] && succ.qargs[1] === node.qargs[1]) {
          toRemove.add(node);
          toRemove.add(succ);
          break;
        }
      }
    }
    for (const node of toRemove) {
      try {
        dag.remove_op_node(node);
      } catch (e) {
      }
    }
    return dag;
  }
};
var CommutativeCancellation = class extends TransformationPass {
  constructor(basisGates = ["cx", "rz", "sx", "x"]) {
    super();
    this.basis_gates = basisGates;
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    const toRemove = /* @__PURE__ */ new Set();
    const stack = [];
    for (const node of opNodes) {
      if (toRemove.has(node)) continue;
      if (stack.length > 0) {
        const prev = stack[stack.length - 1];
        if (prev.op.name === "rz" && node.op.name === "rz" && prev.qargs[0] === node.qargs[0] && prev.op.params.length > 0 && node.op.params.length > 0) {
          const a = typeof prev.op.params[0] === "number" ? prev.op.params[0] : 0;
          const b = typeof node.op.params[0] === "number" ? node.op.params[0] : 0;
          if (Math.abs(a + b) < 1e-12) {
            stack.pop();
            toRemove.add(prev);
            toRemove.add(node);
            continue;
          }
        }
      }
      stack.push(node);
    }
    for (const node of toRemove) {
      try {
        dag.remove_op_node(node);
      } catch (e) {
      }
    }
    return dag;
  }
};
var Optimize1qGates = class extends TransformationPass {
  constructor(basisGates = ["u1", "u2", "u3", "cx"]) {
    super();
    this.basis_gates = basisGates;
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    for (const node of opNodes) {
      if (node.op.num_qubits !== 1) continue;
      const qarg = node.qargs[0];
      let current = node;
      const chain = [node];
      while (true) {
        const succs = dag.successors(current).filter((s) => s.is_op_node());
        if (succs.length !== 1) break;
        const next = succs[0];
        if (next.op.num_qubits !== 1 || next.qargs[0] !== qarg) break;
        chain.push(next);
        current = next;
      }
      if (chain.length > 1) {
        let combined = null;
        let ok = true;
        for (const n of chain) {
          try {
            const mat = n.op.to_matrix();
            combined = combined ? mat.mul(combined) : mat;
          } catch (e) {
            ok = false;
            break;
          }
        }
        if (!ok || !combined) continue;
        const euler = _decomposeZYZ(combined);
        const newGate = makeU3Gate(euler.theta, euler.phi, euler.lambda);
        const newNode = new DAGOpNode(newGate, [qarg], []);
        dag.substitute_node(chain[0], newNode);
        for (let i = 1; i < chain.length; i++) {
          try {
            dag.remove_op_node(chain[i]);
          } catch (e) {
          }
        }
      }
    }
    return dag;
  }
};
function _decomposeZYZ(m) {
  const a = m.get(0, 0);
  const b = m.get(0, 1);
  const c = m.get(1, 0);
  const d = m.get(1, 1);
  const theta = 2 * Math.atan2(b.abs(), a.abs());
  const halfTheta = theta / 2;
  let phi, lambda;
  if (Math.abs(Math.sin(halfTheta)) < 1e-12) {
    phi = 0;
    lambda = a.arg() - d.arg() >= 0 ? 0 : 0;
    lambda = 2 * a.arg();
  } else if (Math.abs(Math.cos(halfTheta)) < 1e-12) {
    phi = c.arg();
    lambda = b.arg() + Math.PI;
  } else {
    phi = c.arg();
    lambda = b.arg() + Math.PI;
  }
  return { theta, phi, lambda };
}
var OptimizeSwapBeforeMeasure = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    for (const node of opNodes) {
      if (node.op.name !== "swap") continue;
      const succs = dag.successors(node).filter((s) => s.is_op_node());
      if (succs.length === 2 && succs.every((s) => s.op.name === "measure")) {
        const m0 = succs[0], m1 = succs[1];
        const tmp = m0.cargs[0];
        m0.cargs[0] = m1.cargs[0];
        m1.cargs[0] = tmp;
        try {
          dag.remove_op_node(node);
        } catch (e) {
        }
      }
    }
    return dag;
  }
};
var RemoveBarriers = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.op_nodes();
    for (const node of opNodes) {
      if (node.op.name === "barrier") {
        dag.remove_op_node(node);
      }
    }
    return dag;
  }
};
var RemoveResetInZeroState = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    const zeroQubits = /* @__PURE__ */ new Set();
    for (const q of dag.qubits) zeroQubits.add(q);
    for (const node of opNodes) {
      if (node.op.name === "reset") {
        if (zeroQubits.has(node.qargs[0])) {
          dag.remove_op_node(node);
        }
      } else {
        for (const q of node.qargs) zeroQubits.delete(q);
      }
    }
    return dag;
  }
};
var DAGFixedPointPass = class extends AnalysisPass {
  constructor() {
    super();
  }
  run(dag) {
    return dag;
  }
};

// src/transpiler/extra_passes.js
var Decompose = class extends TransformationPass {
  constructor(gatesToDecompose = null) {
    super();
    this.gates = gatesToDecompose ? Array.isArray(gatesToDecompose) ? gatesToDecompose : [gatesToDecompose] : null;
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    for (const node of opNodes) {
      if (this.gates && !this.gates.includes(node.op.name)) continue;
      const defFn = node.op._definition;
      if (typeof defFn !== "function") continue;
      try {
        const def = defFn(node.op);
        if (!def || def.length === 0) continue;
        const subDag = new DAGCircuit();
        const numQubits = node.qargs.length;
        const numClbits = node.cargs.length;
        if (numQubits > 0) subDag.add_qreg(new QuantumRegister(numQubits, "q"));
        if (numClbits > 0) subDag.add_creg(new ClassicalRegister(numClbits, "c"));
        for (const [subOp, subQubits, subClbits] of def) {
          const mappedQargs = subQubits.map((_, i) => {
            return null;
          });
          if (!node._qubitIndexMap) {
            node._qubitIndexMap = /* @__PURE__ */ new Map();
            const seen = /* @__PURE__ */ new Set();
            let posIdx = 0;
            for (const sq of node.qargs) {
              if (!seen.has(sq)) {
                node._qubitIndexMap.set(posIdx, sq);
                seen.add(sq);
                posIdx++;
              }
            }
          }
          if (!node._defUniqueQubits) {
            node._defUniqueQubits = [];
            const seenSet = /* @__PURE__ */ new Set();
            for (const [_, sqs, scs] of def) {
              for (const sq of sqs) {
                if (!seenSet.has(sq)) {
                  seenSet.add(sq);
                  node._defUniqueQubits.push(sq);
                }
              }
            }
          }
          const mappedQ = subQubits.map((sq) => {
            const idx = node._defUniqueQubits.indexOf(sq);
            return idx >= 0 ? subDag.qubits[idx] : sq;
          });
          const mappedC = subClbits.map((sc) => {
            const idx = node.cargs.indexOf(sc);
            return idx >= 0 ? subDag.clbits[idx] : sc;
          });
          subDag.apply_operation(subOp.copy(), mappedQ, mappedC);
        }
        dag.substitute_node_with_dag(node, subDag);
      } catch (e) {
      }
    }
    return dag;
  }
};
var RemoveFinalMeasurements = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    const toRemove = /* @__PURE__ */ new Set();
    const qubitLastOp = /* @__PURE__ */ new Map();
    for (let i = opNodes.length - 1; i >= 0; i--) {
      const node = opNodes[i];
      if (node.op.name !== "measure") continue;
      let isFinal = true;
      for (const q of node.qargs) {
        if (qubitLastOp.has(q)) {
          isFinal = false;
          break;
        }
      }
      if (isFinal) {
        toRemove.add(node);
        for (const q of node.qargs) qubitLastOp.set(q, node);
      }
    }
    for (const node of toRemove) {
      try {
        dag.remove_op_node(node);
      } catch (e) {
      }
    }
    return dag;
  }
};
var BarrierBeforeFinalMeasurements = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    let firstFinalMeasureIdx = -1;
    for (let i = opNodes.length - 1; i >= 0; i--) {
      if (opNodes[i].op.name === "measure") {
        firstFinalMeasureIdx = i;
      } else if (firstFinalMeasureIdx !== -1) {
        break;
      }
    }
    if (firstFinalMeasureIdx === -1) return dag;
    const newDag = new DAGCircuit();
    for (const [name, reg] of dag.qregs.entries()) {
      newDag.add_qreg(new QuantumRegister(reg.bits.length, name));
    }
    for (const [name, reg] of dag.cregs.entries()) {
      newDag.add_creg(new ClassicalRegister(reg.bits.length, name));
    }
    for (let i = 0; i < opNodes.length; i++) {
      if (i === firstFinalMeasureIdx) {
        const barrier = new Instruction("barrier", newDag.qubits.length, 0, []);
        newDag.apply_operation(barrier, newDag.qubits.slice(), []);
      }
      const node = opNodes[i];
      const qargs = node.qargs.map((q) => {
        const idx = dag.qubits.indexOf(q);
        return newDag.qubits[idx];
      });
      const cargs = node.cargs.map((c) => {
        const idx = dag.clbits.indexOf(c);
        return newDag.clbits[idx];
      });
      newDag.apply_operation(node.op.copy(), qargs, cargs);
    }
    return newDag;
  }
};
var Collect2qBlocks = class extends AnalysisPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    const blocks = [];
    let currentBlock = null;
    for (const node of opNodes) {
      if (node.op.num_qubits === 2) {
        const q0 = node.qargs[0];
        const q1 = node.qargs[1];
        if (currentBlock && currentBlock.qubits.includes(q0) && currentBlock.qubits.includes(q1)) {
          currentBlock.nodes.push(node);
        } else {
          if (currentBlock) blocks.push(currentBlock);
          currentBlock = { qubits: [q0, q1], nodes: [node] };
        }
      } else {
        if (currentBlock) blocks.push(currentBlock);
        currentBlock = null;
      }
    }
    if (currentBlock) blocks.push(currentBlock);
    this.property_set["2q_blocks"] = blocks;
    return dag;
  }
};
var ConsolidateBlocks = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const collector = new Collect2qBlocks();
    collector.run(dag);
    const blocks = collector.property_set["2q_blocks"] || [];
    for (const block of blocks) {
      if (block.nodes.length < 2) continue;
      let combined = null;
      for (const node of block.nodes) {
        try {
          const mat = node.op.to_matrix();
          combined = combined ? mat.mul(combined) : mat;
        } catch (e) {
          combined = null;
          break;
        }
      }
      if (!combined) continue;
      const UnitaryGate2 = _lookupExtraGate("UnitaryGate");
      if (!UnitaryGate2) continue;
      const newGate = new UnitaryGate2(combined);
      block.nodes[0].op = newGate;
      for (let i = 1; i < block.nodes.length; i++) {
        try {
          dag.remove_op_node(block.nodes[i]);
        } catch (e) {
        }
      }
    }
    return dag;
  }
};
var Unroll3qOrMore = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    for (const node of opNodes) {
      if (node.op.num_qubits < 3) continue;
      const defFn = node.op._definition;
      if (typeof defFn === "function") {
        try {
          const def = defFn(node.op);
          if (def && def.length > 0) {
            dag.substitute_node_with_dag(node, _buildDefDag(def, node));
            continue;
          }
        } catch (e) {
        }
      }
    }
    return dag;
  }
};
var MergeAdjacentBarriers = class extends TransformationPass {
  constructor() {
    super();
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    const toRemove = /* @__PURE__ */ new Set();
    let prevBarrier = null;
    for (const node of opNodes) {
      if (node.op.name === "barrier") {
        if (prevBarrier && prevBarrier.qargs.length === node.qargs.length && prevBarrier.qargs.every((q) => node.qargs.includes(q))) {
          toRemove.add(node);
        } else {
          prevBarrier = node;
        }
      } else {
        prevBarrier = null;
      }
    }
    for (const node of toRemove) {
      try {
        dag.remove_op_node(node);
      } catch (e) {
      }
    }
    return dag;
  }
};
var CheckMap = class extends AnalysisPass {
  constructor(couplingMap = null) {
    super();
    this.coupling_map = couplingMap;
  }
  run(dag) {
    if (!this.coupling_map) {
      this.property_set["is_mapped"] = true;
      return dag;
    }
    let isMapped = true;
    for (const node of dag.topological_op_nodes()) {
      if (node.op.num_qubits !== 2) continue;
      const q0 = dag.qubits.indexOf(node.qargs[0]);
      const q1 = dag.qubits.indexOf(node.qargs[1]);
      if (!this.coupling_map.hasEdge(q0, q1)) {
        isMapped = false;
        break;
      }
    }
    this.property_set["is_mapped"] = isMapped;
    return dag;
  }
};
var GateDirection = class extends TransformationPass {
  constructor(couplingMap = null) {
    super();
    this.coupling_map = couplingMap;
  }
  run(dag) {
    if (!this.coupling_map) return dag;
    const opNodes = dag.topological_op_nodes();
    for (const node of opNodes) {
      if (node.op.num_qubits !== 2) continue;
      const q0 = dag.qubits.indexOf(node.qargs[0]);
      const q1 = dag.qubits.indexOf(node.qargs[1]);
      if (this.coupling_map.hasEdge(q0, q1)) continue;
      if (this.coupling_map.hasEdge(q1, q0)) {
        const tmp = node.qargs[0];
        node.qargs[0] = node.qargs[1];
        node.qargs[1] = tmp;
      }
    }
    return dag;
  }
};
var CountOps = class extends AnalysisPass {
  constructor() {
    super();
  }
  run(dag) {
    const counts = {};
    for (const node of dag.topological_op_nodes()) {
      counts[node.op.name] = (counts[node.op.name] || 0) + 1;
    }
    this.property_set["count_ops"] = counts;
    return dag;
  }
};
var Depth = class extends AnalysisPass {
  constructor() {
    super();
  }
  run(dag) {
    const lastOnQubit = /* @__PURE__ */ new Map();
    for (const q of dag.qubits) lastOnQubit.set(q, -1);
    for (const c of dag.clbits) lastOnQubit.set(c, -1);
    let maxDepth = 0;
    for (const node of dag.topological_op_nodes()) {
      if (node.op.name === "barrier") continue;
      let d = 0;
      for (const q of node.qargs) d = Math.max(d, lastOnQubit.get(q) + 1);
      for (const c of node.cargs) d = Math.max(d, lastOnQubit.get(c) + 1);
      for (const q of node.qargs) lastOnQubit.set(q, d);
      for (const c of node.cargs) lastOnQubit.set(c, d);
      if (d + 1 > maxDepth) maxDepth = d + 1;
    }
    this.property_set["depth"] = maxDepth;
    return dag;
  }
};
var Size = class extends AnalysisPass {
  constructor() {
    super();
  }
  run(dag) {
    let n = 0;
    for (const node of dag.topological_op_nodes()) {
      if (node.op.name === "barrier") continue;
      n++;
    }
    this.property_set["size"] = n;
    return dag;
  }
};
var Width = class extends AnalysisPass {
  constructor() {
    super();
  }
  run(dag) {
    this.property_set["width"] = dag.qubits.length + dag.clbits.length;
    return dag;
  }
};
var Optimize1qGatesDecomposition = class extends TransformationPass {
  constructor(basisGates = ["u3", "cx"]) {
    super();
    this.basis_gates = basisGates;
  }
  run(dag) {
    const opNodes = dag.topological_op_nodes();
    for (const node of opNodes) {
      if (node.op.num_qubits !== 1) continue;
      const succs = dag.successors(node).filter((s) => s.is_op_node());
      if (succs.length !== 1) continue;
      const next = succs[0];
      if (next.op.num_qubits !== 1 || next.qargs[0] !== node.qargs[0]) continue;
      try {
        const m1 = node.op.to_matrix();
        const m2 = next.op.to_matrix();
        const combined = m2.mul(m1);
        const euler = _decomposeZYZ2(combined);
        const newGate = makeU3Gate(euler.theta, euler.phi, euler.lambda);
        const newNode = new DAGOpNode(newGate, [node.qargs[0]], []);
        dag.substitute_node(node, newNode);
        try {
          dag.remove_op_node(next);
        } catch (e) {
        }
      } catch (e) {
      }
    }
    return dag;
  }
};
var _extraGateRegistry = null;
function _lookupExtraGate(name) {
  if (!_extraGateRegistry) {
    try {
      const mod = globalThis.__ketraExtraGates;
      if (mod && mod[name]) return mod[name];
    } catch (e) {
    }
    return null;
  }
  return _extraGateRegistry[name];
}
function _decomposeZYZ2(m) {
  const a = m.get(0, 0);
  const b = m.get(0, 1);
  const c = m.get(1, 0);
  const d = m.get(1, 1);
  const theta = 2 * Math.atan2(b.abs(), a.abs());
  const halfTheta = theta / 2;
  let phi, lambda;
  if (Math.abs(Math.sin(halfTheta)) < 1e-12) {
    phi = 0;
    lambda = 2 * a.arg();
  } else if (Math.abs(Math.cos(halfTheta)) < 1e-12) {
    phi = c.arg();
    lambda = b.arg() + Math.PI;
  } else {
    phi = c.arg();
    lambda = b.arg() + Math.PI;
  }
  return { theta, phi, lambda };
}

// src/dagcircuit/dag_passes.js
function collect_1q_runs2(dag) {
  const runs = [];
  const visited = /* @__PURE__ */ new Set();
  for (const node of dag.topological_op_nodes()) {
    if (visited.has(node)) continue;
    if (node.op.num_qubits !== 1) continue;
    if (!node.qubits || node.qubits.length === 0) continue;
    if (["measure", "reset", "barrier"].includes(node.op.name)) continue;
    const run = [node];
    visited.add(node);
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter((s) => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.num_qubits !== 1) break;
      if (!next.qubits || next.qubits.length === 0) break;
      if (next.qubits[0] !== current.qubits[0]) break;
      if (["measure", "reset", "barrier"].includes(next.op.name)) break;
      run.push(next);
      visited.add(next);
      current = next;
    }
    if (run.length > 1) runs.push(run);
  }
  return runs;
}
function collect_2q_runs2(dag) {
  const runs = [];
  const visited = /* @__PURE__ */ new Set();
  for (const node of dag.topological_op_nodes()) {
    if (visited.has(node)) continue;
    if (node.op.num_qubits !== 2) continue;
    if (["measure", "reset", "barrier"].includes(node.op.name)) continue;
    const run = [node];
    visited.add(node);
    let current = node;
    while (true) {
      const succs = dag.successors(current).filter((s) => s.is_op_node());
      if (succs.length !== 1) break;
      const next = succs[0];
      if (next.op.num_qubits !== 2) break;
      if (next.qubits[0] !== current.qubits[0] || next.qubits[1] !== current.qubits[1]) break;
      run.push(next);
      visited.add(next);
      current = next;
    }
    if (run.length > 1) runs.push(run);
  }
  return runs;
}
function consolidate_blocks2(dag) {
  const runs = collect_1q_runs2(dag);
  for (const run of runs) {
    let combined = null;
    for (const node of run) {
      try {
        const m = node.op.to_matrix();
        combined = combined ? m.mul(combined) : m;
      } catch (e) {
        break;
      }
    }
    if (combined) {
      const newGate = new Gate("unitary", 1, [combined]);
      newGate._matrixBuilder = () => combined;
      const newNode = new DAGOpNode(newGate, run[0].qargs.slice(), []);
      dag.substitute_node(run[0], newNode);
      for (let i = 1; i < run.length; i++) {
        try {
          dag.remove_op_node(run[i]);
        } catch (e) {
        }
      }
    }
  }
  return dag;
}
function optimize_cliffords2(dag) {
  const selfInverse = /* @__PURE__ */ new Set(["h", "x", "y", "z", "cx", "cz", "swap", "sx", "sxdg"]);
  const nodes = dag.topological_op_nodes();
  const toRemove = /* @__PURE__ */ new Set();
  const sSimplifications = [];
  const stack = [];
  for (const node of nodes) {
    if (toRemove.has(node)) continue;
    if (stack.length > 0) {
      const prev = stack[stack.length - 1];
      if (prev.op.name === node.op.name && prev.qubits.length === node.qubits.length && prev.qubits.every((q, i) => q === node.qubits[i])) {
        if (selfInverse.has(prev.op.name)) {
          stack.pop();
          toRemove.add(prev);
          toRemove.add(node);
          continue;
        }
        if (prev.op.name === "s") {
          stack.pop();
          const newGate = new Gate("sdg", 1, []);
          newGate._matrixBuilder = () => makeSdgGate().to_matrix();
          sSimplifications.push({ replace: prev, with: new DAGOpNode(newGate, node.qargs.slice(), []) });
          toRemove.add(node);
          continue;
        }
      }
    }
    stack.push(node);
  }
  for (const { replace, with: newNode } of sSimplifications) {
    dag.substitute_node(replace, newNode);
  }
  for (const node of toRemove) {
    try {
      dag.remove_op_node(node);
    } catch (e) {
    }
  }
  return dag;
}
function remove_diagonal_gates_before_measure2(dag) {
  const diagonalGates = ["z", "s", "sdg", "t", "tdg", "p", "u1", "rz", "cz", "cp"];
  const toRemove = [];
  for (const node of dag.topological_op_nodes()) {
    if (!diagonalGates.includes(node.op.name)) continue;
    const succs = dag.successors(node).filter((s) => s.is_op_node());
    if (succs.length > 0 && succs.every((s) => s.op.name === "measure")) {
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.remove_op_node(node);
  }
  return dag;
}
function elide_permutations(dag) {
  const toRemove = [];
  for (const node of dag.topological_op_nodes()) {
    if (node.op.name !== "swap") continue;
    const succs = dag.successors(node).filter((s) => s.is_op_node());
    if (succs.length === 2 && succs.every((s) => s.op.name === "measure")) {
      const m0 = succs[0], m1 = succs[1];
      const tmp = m0.cargs[0];
      m0.cargs[0] = m1.cargs[0];
      m1.cargs[0] = tmp;
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.remove_op_node(node);
  }
  return dag;
}
function remove_redundant_gates(dag) {
  const toRemove = [];
  for (const node of dag.topological_op_nodes()) {
    if (node.op.name === "id") toRemove.push(node);
    if (["rz", "rx", "ry", "p", "u1"].includes(node.op.name) && node.op.params.length > 0 && typeof node.op.params[0] === "number" && Math.abs(node.op.params[0]) < 1e-15) {
      toRemove.push(node);
    }
  }
  for (const node of toRemove) {
    dag.remove_op_node(node);
  }
  return dag;
}
function commutative_cancellation(dag) {
  const nodes = dag.topological_op_nodes();
  const selfInverse = /* @__PURE__ */ new Set(["h", "x", "y", "z", "cx", "cz", "swap"]);
  const toRemove = /* @__PURE__ */ new Set();
  const stack = [];
  for (const node of nodes) {
    if (toRemove.has(node)) continue;
    if (stack.length > 0) {
      const prev = stack[stack.length - 1];
      if (selfInverse.has(prev.op.name) && prev.op.name === node.op.name && prev.qubits.length === node.qubits.length && prev.qubits.every((q, i) => q === node.qubits[i])) {
        stack.pop();
        toRemove.add(prev);
        toRemove.add(node);
        continue;
      }
    }
    stack.push(node);
  }
  for (const node of toRemove) {
    try {
      dag.remove_op_node(node);
    } catch (e) {
    }
  }
  return dag;
}
function template_optimization(dag) {
  const nodes = dag.topological_op_nodes();
  const replacements = [];
  const toRemove = [];
  for (let i = 0; i < nodes.length - 2; i++) {
    const n1 = nodes[i], n2 = nodes[i + 1], n3 = nodes[i + 2];
    if (n1.op.name === "h" && n2.op.name === "z" && n3.op.name === "h" && n1.qubits[0] === n2.qubits[0] && n2.qubits[0] === n3.qubits[0]) {
      const newGate = new Gate("x", 1, []);
      newGate._matrixBuilder = () => makeXGate().to_matrix();
      replacements.push({ node: n1, gate: newGate });
      toRemove.push(n2, n3);
    } else if (n1.op.name === "h" && n2.op.name === "x" && n3.op.name === "h" && n1.qubits[0] === n2.qubits[0] && n2.qubits[0] === n3.qubits[0]) {
      const newGate = new Gate("z", 1, []);
      newGate._matrixBuilder = () => makeZGate().to_matrix();
      replacements.push({ node: n1, gate: newGate });
      toRemove.push(n2, n3);
    }
  }
  for (const { node, gate } of replacements) {
    dag.substitute_node(node, new DAGOpNode(gate, node.qargs.slice(), []));
  }
  for (const node of toRemove) {
    try {
      dag.remove_op_node(node);
    } catch (e) {
    }
  }
  return dag;
}

// src/converters/converters.js
function circuit_to_dag(circuit) {
  const dag = new DAGCircuit();
  dag.name = circuit.name;
  dag.global_phase = circuit.global_phase;
  dag.metadata = circuit.metadata;
  for (const r4 of circuit.qregs) {
    dag.add_qreg(new QuantumRegister(r4.size, r4.name));
  }
  for (const r4 of circuit.cregs) {
    dag.add_creg(new ClassicalRegister(r4.size, r4.name));
  }
  if (dag.qubits.length < circuit.qubits.length) {
    const qr = new QuantumRegister(circuit.qubits.length - dag.qubits.length, "q");
    dag.add_qreg(qr);
  }
  if (dag.clbits.length < circuit.clbits.length) {
    const cr = new ClassicalRegister(circuit.clbits.length - dag.clbits.length, "c");
    dag.add_creg(cr);
  }
  const qubitMap = /* @__PURE__ */ new Map();
  for (let i = 0; i < circuit.qubits.length; i++) {
    qubitMap.set(circuit.qubits[i], dag.qubits[i]);
  }
  const clbitMap = /* @__PURE__ */ new Map();
  for (let i = 0; i < circuit.clbits.length; i++) {
    clbitMap.set(circuit.clbits[i], dag.clbits[i]);
  }
  for (const ci of circuit.data) {
    const op = ci.operation.copy();
    const qargs = ci.qubits.map((q) => qubitMap.get(q));
    const cargs = ci.clbits.map((c) => clbitMap.get(c));
    dag.apply_operation(op, qargs, cargs);
  }
  return dag;
}
function dag_to_circuit(dag) {
  const regs = [];
  for (const [name, reg] of dag.qregs.entries()) {
    regs.push(new QuantumRegister(reg.bits.length, name));
  }
  for (const [name, reg] of dag.cregs.entries()) {
    regs.push(new ClassicalRegister(reg.bits.length, name));
  }
  const circuit = new QuantumCircuit(...regs);
  circuit.name = dag.name;
  circuit.global_phase = dag.global_phase;
  circuit.metadata = dag.metadata;
  const qubitMap = /* @__PURE__ */ new Map();
  for (let i = 0; i < dag.qubits.length; i++) {
    qubitMap.set(dag.qubits[i], circuit.qubits[i]);
  }
  const clbitMap = /* @__PURE__ */ new Map();
  for (let i = 0; i < dag.clbits.length; i++) {
    clbitMap.set(dag.clbits[i], circuit.clbits[i]);
  }
  for (const node of dag.topological_op_nodes()) {
    const op = node.op.copy();
    const qargs = node.qargs.map((q) => qubitMap.get(q));
    const cargs = node.cargs.map((c) => clbitMap.get(c));
    circuit.append(op, qargs, cargs);
  }
  return circuit;
}

// src/qasm/qasm_parser.js
var Token = class {
  constructor(type, value, pos) {
    this.type = type;
    this.value = value;
    this.pos = pos;
  }
  toString() {
    return `${this.type}(${this.value})`;
  }
};
var KEYWORDS = /* @__PURE__ */ new Set([
  "OPENQASM",
  "include",
  "qreg",
  "creg",
  "gate",
  "measure",
  "barrier",
  "reset",
  "if",
  "opaque",
  "format"
]);
function tokenize(source) {
  const tokens = [];
  let i = 0;
  let line = 1, col = 1;
  while (i < source.length) {
    const c = source[i];
    if (c === " " || c === "	") {
      i++;
      col++;
      continue;
    }
    if (c === "\n") {
      i++;
      line++;
      col = 1;
      continue;
    }
    if (c === "\r") {
      i++;
      continue;
    }
    if (c === "/" && source[i + 1] === "/") {
      while (i < source.length && source[i] !== "\n") {
        i++;
        col++;
      }
      continue;
    }
    if (c === '"') {
      let str = "";
      i++;
      col++;
      while (i < source.length && source[i] !== '"') {
        str += source[i];
        i++;
        col++;
      }
      i++;
      col++;
      tokens.push(new Token("STRING", str, { line, col }));
      continue;
    }
    if (c === "-" && source[i + 1] !== ">" || c >= "0" && c <= "9") {
      let num = "";
      if (c === "-") {
        num += c;
        i++;
        col++;
      }
      while (i < source.length && (source[i] >= "0" && source[i] <= "9" || source[i] === ".")) {
        num += source[i];
        i++;
        col++;
      }
      if (source[i] === "*") {
        const afterStar = source.slice(i + 1, i + 3);
        if (afterStar === "pi") {
          num = String(parseFloat(num) * Math.PI);
          i += 3;
          col += 3;
        }
      }
      if (source.slice(i, i + 2) === "pi") {
        const piVal = parseFloat(num);
        if (source[i + 2] === "/") {
          const denom = parseInt(source.slice(i + 3), 10);
          num = String(piVal * Math.PI / denom);
          i += 3 + String(denom).length;
          col += 3 + String(denom).length;
        } else {
          num = String(piVal * Math.PI);
          i += 2;
          col += 2;
        }
      }
      tokens.push(new Token("NUMBER", num, { line, col }));
      continue;
    }
    if (c >= "a" && c <= "z" || c >= "A" && c <= "Z" || c === "_") {
      let id = "";
      while (i < source.length && (source[i] >= "a" && source[i] <= "z" || source[i] >= "A" && source[i] <= "Z" || source[i] >= "0" && source[i] <= "9" || source[i] === "_")) {
        id += source[i];
        i++;
        col++;
      }
      if (KEYWORDS.has(id)) {
        tokens.push(new Token("KEYWORD", id, { line, col }));
      } else {
        tokens.push(new Token("IDENT", id, { line, col }));
      }
      if (id === "pi") {
        tokens.pop();
        tokens.push(new Token("NUMBER", String(Math.PI), { line, col }));
      }
      continue;
    }
    if ("()[]{};,->=+*/".indexOf(c) !== -1) {
      if (c === "-" && source[i + 1] === ">") {
        tokens.push(new Token("ARROW", "->", { line, col }));
        i += 2;
        col += 2;
        continue;
      }
      if (c === "=" && source[i + 1] === "=") {
        tokens.push(new Token("EQ", "==", { line, col }));
        i += 2;
        col += 2;
        continue;
      }
      const typeMap = {
        "(": "LPAREN",
        ")": "RPAREN",
        "[": "LBRACKET",
        "]": "RBRACKET",
        "{": "LBRACE",
        "}": "RBRACE",
        ";": "SEMICOLON",
        ",": "COMMA",
        "=": "ASSIGN",
        "+": "PLUS",
        "*": "STAR",
        "/": "SLASH"
      };
      tokens.push(new Token(typeMap[c] || c, c, { line, col }));
      i++;
      col++;
      continue;
    }
    i++;
    col++;
  }
  tokens.push(new Token("EOF", null, { line, col }));
  return tokens;
}
var QASMParser = class {
  constructor(source) {
    this.source = source;
    this.tokens = tokenize(source);
    this.pos = 0;
    this.circuit = null;
    this.gateDefinitions = /* @__PURE__ */ new Map();
    this.qregs = /* @__PURE__ */ new Map();
    this.cregs = /* @__PURE__ */ new Map();
    this.allQubits = [];
    this.allClbits = [];
  }
  peek() {
    return this.tokens[this.pos];
  }
  next() {
    return this.tokens[this.pos++];
  }
  expect(type) {
    const t = this.tokens[this.pos];
    if (t.type !== type) {
      throw new Error(`QASM parse error at line ${t.pos.line}: expected ${type} but got ${t.type} (${t.value})`);
    }
    this.pos++;
    return t;
  }
  accept(type) {
    if (this.tokens[this.pos].type === type) {
      return this.tokens[this.pos++];
    }
    return null;
  }
  parse() {
    const versionTok = this.expect("KEYWORD");
    if (versionTok.value !== "OPENQASM") {
      throw new Error(`Expected OPENQASM, got ${versionTok.value}`);
    }
    this.expect("NUMBER");
    this.expect("SEMICOLON");
    this.circuit = new QuantumCircuit();
    while (this.peek().type !== "EOF") {
      this.parseStatement();
    }
    return this.circuit;
  }
  parseStatement() {
    const t = this.peek();
    if (t.type === "KEYWORD") {
      switch (t.value) {
        case "include":
          this.parseInclude();
          return;
        case "qreg":
          this.parseQregDecl();
          return;
        case "creg":
          this.parseCregDecl();
          return;
        case "gate":
          this.parseGateDef();
          return;
        case "measure":
          this.parseMeasure();
          return;
        case "barrier":
          this.parseBarrier();
          return;
        case "reset":
          this.parseReset();
          return;
        case "if":
          this.parseIf();
          return;
        case "opaque":
          this.parseOpaque();
          return;
      }
    }
    if (t.type === "IDENT") {
      this.parseGateCall();
      return;
    }
    throw new Error(`QASM parse error at line ${t.pos.line}: unexpected token ${t.type} (${t.value})`);
  }
  parseInclude() {
    this.expect("KEYWORD");
    this.expect("STRING");
    this.expect("SEMICOLON");
  }
  parseQregDecl() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const size = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("SEMICOLON");
    const reg = new QuantumRegister(size, name);
    this.circuit.add_register(reg);
    this.qregs.set(name, reg);
    for (const q of reg._bits) this.allQubits.push({ reg: name, index: q.index, bit: q });
  }
  parseCregDecl() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const size = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("SEMICOLON");
    const reg = new ClassicalRegister(size, name);
    this.circuit.add_register(reg);
    this.cregs.set(name, reg);
    for (const c of reg._bits) this.allClbits.push({ reg: name, index: c.index, bit: c });
  }
  parseGateDef() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    let params = [];
    if (this.accept("LBRACKET")) {
      while (this.peek().type !== "RBRACKET") {
        params.push(this.expect("IDENT").value);
        this.accept("COMMA");
      }
      this.expect("RBRACKET");
    }
    let qargs = [];
    if (this.peek().type === "IDENT") {
      while (this.peek().type === "IDENT") {
        qargs.push(this.expect("IDENT").value);
        this.accept("COMMA");
      }
    }
    this.expect("LBRACE");
    const body = [];
    while (this.peek().type !== "RBRACE") {
      const gateName = this.expect("IDENT").value;
      let gateParams = [];
      if (this.accept("LPAREN")) {
        while (this.peek().type !== "RPAREN") {
          gateParams.push(this.parseExpression());
          this.accept("COMMA");
        }
        this.expect("RPAREN");
      }
      const gateQargs = [];
      while (this.peek().type === "IDENT") {
        const qName = this.expect("IDENT").value;
        if (this.accept("LBRACKET")) {
          const qIdx = parseInt(this.expect("NUMBER").value, 10);
          this.expect("RBRACKET");
          gateQargs.push({ name: qName, index: qIdx });
        } else {
          gateQargs.push({ name: qName, index: 0 });
        }
        this.accept("COMMA");
      }
      this.expect("SEMICOLON");
      body.push({ gate: gateName, params: gateParams, qargs: gateQargs });
    }
    this.expect("RBRACE");
    this.gateDefinitions.set(name, { params, qargs, body });
  }
  parseOpaque() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    while (this.peek().type !== "SEMICOLON") this.next();
    this.expect("SEMICOLON");
  }
  parseMeasure() {
    this.expect("KEYWORD");
    const qubit = this.parseBitRef();
    this.expect("ARROW");
    const clbit = this.parseBitRef();
    this.expect("SEMICOLON");
    this.circuit.measure(qubit, clbit);
  }
  parseBarrier() {
    this.expect("KEYWORD");
    const qubits = [];
    if (this.peek().type === "IDENT") {
      do {
        qubits.push(this.parseBitRef());
      } while (this.accept("COMMA"));
    }
    this.expect("SEMICOLON");
    if (qubits.length === 0) this.circuit.barrier();
    else this.circuit.barrier(qubits);
  }
  parseReset() {
    this.expect("KEYWORD");
    const qubit = this.parseBitRef();
    this.expect("SEMICOLON");
    this.circuit.reset(qubit);
  }
  parseIf() {
    this.expect("KEYWORD");
    this.expect("LPAREN");
    const cregName = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const cregIdx = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("EQ");
    const value = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RPAREN");
    const dataLenBefore = this.circuit.data.length;
    this.parseGateCall();
    const creg = this.cregs.get(cregName);
    for (let i = dataLenBefore; i < this.circuit.data.length; i++) {
      const ci = this.circuit.data[i];
      ci.operation.condition = {
        register: cregName,
        index: cregIdx,
        value,
        bit: creg ? creg._bits[cregIdx] : null
      };
    }
  }
  parseGateCall() {
    const name = this.expect("IDENT").value;
    let params = [];
    if (this.accept("LPAREN")) {
      while (this.peek().type !== "RPAREN") {
        params.push(this.parseExpression());
        this.accept("COMMA");
      }
      this.expect("RPAREN");
    }
    const qargs = [];
    const cargs = [];
    while (this.peek().type === "IDENT" || this.peek().type === "KEYWORD") {
      const ref = this.parseBitRef();
      if (ref instanceof Clbit || ref && ref.register instanceof ClassicalRegister) {
        cargs.push(ref);
      } else {
        qargs.push(ref);
      }
      if (!this.accept("COMMA")) break;
    }
    this.expect("SEMICOLON");
    this._applyGateByName(name, params, qargs, cargs);
  }
  _applyGateByName(name, params, qargs, cargs) {
    if (this.gateDefinitions.has(name)) {
      const def = this.gateDefinitions.get(name);
      const argMap = /* @__PURE__ */ new Map();
      for (let i = 0; i < def.qargs.length && i < qargs.length; i++) {
        argMap.set(def.qargs[i], qargs[i]);
      }
      for (const stmt of def.body) {
        const mappedQargs = stmt.qargs.map((q) => argMap.get(q.name) || qargs[0]);
        this._applyGateByName(stmt.gate, stmt.params, mappedQargs, []);
      }
      return;
    }
    const n = name.toLowerCase();
    switch (n) {
      case "h":
        this.circuit.h(qargs[0]);
        break;
      case "x":
        this.circuit.x(qargs[0]);
        break;
      case "y":
        this.circuit.y(qargs[0]);
        break;
      case "z":
        this.circuit.z(qargs[0]);
        break;
      case "s":
        this.circuit.s(qargs[0]);
        break;
      case "sdg":
        this.circuit.sdg(qargs[0]);
        break;
      case "t":
        this.circuit.t(qargs[0]);
        break;
      case "tdg":
        this.circuit.tdg(qargs[0]);
        break;
      case "sx":
        this.circuit.sx(qargs[0]);
        break;
      case "id":
        this.circuit.id(qargs[0]);
        break;
      case "u":
      case "u3":
        this.circuit.u3(params[0], params[1], params[2], qargs[0]);
        break;
      case "u2":
        this.circuit.u2(params[0], params[1], qargs[0]);
        break;
      case "u1":
      case "p":
        this.circuit.u1(params[0], qargs[0]);
        break;
      case "rx":
        this.circuit.rx(params[0], qargs[0]);
        break;
      case "ry":
        this.circuit.ry(params[0], qargs[0]);
        break;
      case "rz":
        this.circuit.rz(params[0], qargs[0]);
        break;
      case "cx":
        this.circuit.cx(qargs[0], qargs[1]);
        break;
      case "cy":
        this.circuit.cy(qargs[0], qargs[1]);
        break;
      case "cz":
        this.circuit.cz(qargs[0], qargs[1]);
        break;
      case "ch":
        this.circuit.ch(qargs[0], qargs[1]);
        break;
      case "csx":
        this.circuit.csx(qargs[0], qargs[1]);
        break;
      case "swap":
        this.circuit.swap(qargs[0], qargs[1]);
        break;
      case "iswap":
        this.circuit.iswap(qargs[0], qargs[1]);
        break;
      case "ccx":
        this.circuit.ccx(qargs[0], qargs[1], qargs[2]);
        break;
      case "cswap":
        this.circuit.cswap(qargs[0], qargs[1], qargs[2]);
        break;
      case "crx":
        this.circuit.crx(params[0], qargs[0], qargs[1]);
        break;
      case "cry":
        this.circuit.cry(params[0], qargs[0], qargs[1]);
        break;
      case "crz":
        this.circuit.crz(params[0], qargs[0], qargs[1]);
        break;
      case "cu1":
      case "cp":
        this.circuit.cp(params[0], qargs[0], qargs[1]);
        break;
      case "cu3":
        this.circuit.cu3(params[0], params[1], params[2], qargs[0], qargs[1]);
        break;
      case "rxx":
        this.circuit.rxx(params[0], qargs[0], qargs[1]);
        break;
      case "ryy":
        this.circuit.ryy(params[0], qargs[0], qargs[1]);
        break;
      case "rzz":
        this.circuit.rzz(params[0], qargs[0], qargs[1]);
        break;
      case "rzx":
        this.circuit.rzx(params[0], qargs[0], qargs[1]);
        break;
      default:
        const numQubits = qargs.length;
        const numClbits = cargs.length;
        const instr = new Instruction(n, numQubits, numClbits, params);
        this.circuit.append(instr, qargs, cargs);
    }
  }
  parseBitRef() {
    const name = this.expect("IDENT").value;
    if (this.accept("LBRACKET")) {
      const idx = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
      if (this.qregs.has(name)) {
        return this.qregs.get(name)._bits[idx];
      }
      if (this.cregs.has(name)) {
        return this.cregs.get(name)._bits[idx];
      }
      throw new Error(`Unknown register: ${name}`);
    }
    if (this.qregs.has(name)) return this.qregs.get(name)._bits[0];
    if (this.cregs.has(name)) return this.cregs.get(name)._bits[0];
    throw new Error(`Unknown register: ${name}`);
  }
  parseExpression() {
    const t = this.peek();
    if (t.type === "NUMBER") {
      this.next();
      let value = parseFloat(t.value);
      while (this.peek().type === "PLUS" || this.peek().type === "STAR" || this.peek().type === "SLASH") {
        const op = this.next().type;
        const right = this.expect("NUMBER").value;
        if (op === "PLUS") value += parseFloat(right);
        else if (op === "STAR") value *= parseFloat(right);
        else if (op === "SLASH") value /= parseFloat(right);
      }
      return value;
    }
    if (t.type === "IDENT") {
      this.next();
      if (t.value === "pi") {
        let value = Math.PI;
        while (this.peek().type === "STAR" || this.peek().type === "SLASH") {
          const op = this.next().type;
          const right = this.expect("NUMBER").value;
          if (op === "STAR") value *= parseFloat(right);
          else if (op === "SLASH") value /= parseFloat(right);
        }
        return value;
      }
      return t.value;
    }
    throw new Error(`QASM parse error: expected expression, got ${t.type}`);
  }
};

// src/qasm/qasm_exporter.js
var QASMExporter = class {
  constructor(options = {}) {
    this.includes = options.includes || ["qelib1.inc"];
    this.basis_gates = options.basis_gates || ["cx", "u3", "u1", "u2"];
  }
  export(circuit) {
    const lines = [];
    lines.push("OPENQASM 2.0;");
    for (const inc of this.includes) {
      lines.push(`include "${inc}";`);
    }
    const qrName = circuit.qregs.length > 0 ? circuit.qregs[0].name : "q";
    const crName = circuit.cregs.length > 0 ? circuit.cregs[0].name : "c";
    if (circuit.num_qubits > 0) {
      lines.push(`qreg ${qrName}[${circuit.num_qubits}];`);
    }
    if (circuit.num_clbits > 0) {
      lines.push(`creg ${crName}[${circuit.num_clbits}];`);
    }
    const qRef = (q) => `${qrName}[${circuit.qubit_indices(q)}]`;
    const cRef = (c) => `${crName}[${circuit.clbit_indices(c)}]`;
    const customGates = /* @__PURE__ */ new Map();
    for (const ci of circuit.data) {
      const op = ci.operation;
      const name = op.name;
      const params = op.params || [];
      const qargs = ci.qubits.map(qRef).join(",");
      const cargs = ci.clbits.map(cRef).join(",");
      if (name === "barrier") {
        const qList = ci.qubits.map(qRef).join(",");
        lines.push(`barrier ${qList};`);
        continue;
      }
      if (name === "measure") {
        lines.push(`measure ${qRef(ci.qubits[0])} -> ${cRef(ci.clbits[0])};`);
        continue;
      }
      if (name === "reset") {
        lines.push(`reset ${qRef(ci.qubits[0])};`);
        continue;
      }
      if (name === "delay") {
        lines.push(`// delay ${params[0]} ${params[1] || "dt"} on ${qargs}`);
        continue;
      }
      if (name === "if_else" || name === "while_loop") {
        continue;
      }
      const cond = op.condition;
      let condPrefix = "";
      if (cond) {
        condPrefix = `if (${cond.register}[${cond.index}] == ${cond.value}) `;
      }
      const paramStr = params.length > 0 ? `(${params.map((p) => _formatParam(p)).join(",")})` : "";
      let qasmName = name;
      if (name === "u") qasmName = "u3";
      if (name === "p") qasmName = "u1";
      if (name === "rx" || name === "ry" || name === "rz" || name === "rxx" || name === "ryy" || name === "rzz" || name === "rzx") {
      }
      if (params.length > 0) {
        lines.push(`${condPrefix}${qasmName}${paramStr} ${qargs};`);
      } else {
        const stdNames = [
          "h",
          "x",
          "y",
          "z",
          "s",
          "sdg",
          "t",
          "tdg",
          "sx",
          "id",
          "cx",
          "cy",
          "cz",
          "ch",
          "swap",
          "ccx",
          "cswap",
          "iswap",
          "crx",
          "cry",
          "crz"
        ];
        if (stdNames.indexOf(name) !== -1) {
          lines.push(`${condPrefix}${qasmName} ${qargs};`);
        } else {
          lines.push(`${condPrefix}${qasmName} ${qargs};`);
        }
      }
    }
    return lines.join("\n");
  }
};
_registerQASMExporterClass(QASMExporter);

// src/qasm/qasm3_parser.js
var Token3 = class {
  constructor(type, value, pos) {
    this.type = type;
    this.value = value;
    this.pos = pos;
  }
};
var KEYWORDS3 = /* @__PURE__ */ new Set([
  "OPENQASM",
  "include",
  "qreg",
  "creg",
  "gate",
  "measure",
  "barrier",
  "reset",
  "if",
  "opaque",
  "format",
  // QASM 3.0 additions
  "let",
  "for",
  "while",
  "switch",
  "case",
  "default",
  "break",
  "continue",
  "def",
  "return",
  "cal",
  "defcal",
  "stretch",
  "durationof",
  "extern",
  "const",
  "int",
  "uint",
  "float",
  "bool",
  "bit",
  "qubit",
  "qubits",
  "in",
  "array",
  "void",
  "input",
  "output"
]);
function tokenize3(source) {
  const tokens = [];
  let i = 0, line = 1, col = 1;
  while (i < source.length) {
    const c = source[i];
    if (c === " " || c === "	") {
      i++;
      col++;
      continue;
    }
    if (c === "\n") {
      i++;
      line++;
      col = 1;
      continue;
    }
    if (c === "\r") {
      i++;
      continue;
    }
    if (c === "/" && source[i + 1] === "/") {
      while (i < source.length && source[i] !== "\n") {
        i++;
        col++;
      }
      continue;
    }
    if (c === "/" && source[i + 1] === "*") {
      while (i < source.length && !(source[i] === "*" && source[i + 1] === "/")) {
        if (source[i] === "\n") {
          line++;
          col = 1;
        } else col++;
        i++;
      }
      i += 2;
      col += 2;
      continue;
    }
    if (c === '"') {
      let str = "";
      i++;
      col++;
      while (i < source.length && source[i] !== '"') {
        str += source[i];
        i++;
        col++;
      }
      i++;
      col++;
      tokens.push(new Token3("STRING", str, { line, col }));
      continue;
    }
    if (c === "-" && source[i + 1] !== ">" || c >= "0" && c <= "9") {
      let num = "";
      if (c === "-") {
        num += c;
        i++;
        col++;
      }
      while (i < source.length && (source[i] >= "0" && source[i] <= "9" || source[i] === ".")) {
        num += source[i];
        i++;
        col++;
      }
      if (source.slice(i, i + 2) === "pi") {
        const val = parseFloat(num);
        if (source[i + 2] === "/") {
          const denom = parseInt(source.slice(i + 3), 10);
          num = String(val * Math.PI / denom);
          i += 3 + String(denom).length;
          col += 3 + String(denom).length;
        } else {
          num = String(val * Math.PI);
          i += 2;
          col += 2;
        }
      }
      tokens.push(new Token3("NUMBER", num, { line, col }));
      continue;
    }
    if (c >= "a" && c <= "z" || c >= "A" && c <= "Z" || c === "_") {
      let id = "";
      while (i < source.length && (source[i] >= "a" && source[i] <= "z" || source[i] >= "A" && source[i] <= "Z" || source[i] >= "0" && source[i] <= "9" || source[i] === "_")) {
        id += source[i];
        i++;
        col++;
      }
      if (KEYWORDS3.has(id)) tokens.push(new Token3("KEYWORD", id, { line, col }));
      else if (id === "pi") tokens.push(new Token3("NUMBER", String(Math.PI), { line, col }));
      else tokens.push(new Token3("IDENT", id, { line, col }));
      continue;
    }
    if (c === "$") {
      let num = "";
      i++;
      col++;
      while (i < source.length && source[i] >= "0" && source[i] <= "9") {
        num += source[i];
        i++;
        col++;
      }
      tokens.push(new Token3("PHYSICAL_QUBIT", num, { line, col }));
      continue;
    }
    if (c === "-" && source[i + 1] === ">") {
      tokens.push(new Token3("ARROW", "->", { line, col }));
      i += 2;
      col += 2;
      continue;
    }
    if (c === "=" && source[i + 1] === "=") {
      tokens.push(new Token3("EQ", "==", { line, col }));
      i += 2;
      col += 2;
      continue;
    }
    if (c === "!" && source[i + 1] === "=") {
      tokens.push(new Token3("NEQ", "!=", { line, col }));
      i += 2;
      col += 2;
      continue;
    }
    if (c === "<" && source[i + 1] === "=") {
      tokens.push(new Token3("LE", "<=", { line, col }));
      i += 2;
      col += 2;
      continue;
    }
    if (c === ">" && source[i + 1] === "=") {
      tokens.push(new Token3("GE", ">=", { line, col }));
      i += 2;
      col += 2;
      continue;
    }
    if (c === "+" && source[i + 1] === "+") {
      tokens.push(new Token3("INC", "++", { line, col }));
      i += 2;
      col += 2;
      continue;
    }
    if (c === "-" && source[i + 1] === "-") {
      tokens.push(new Token3("DEC", "--", { line, col }));
      i += 2;
      col += 2;
      continue;
    }
    const typeMap = {
      "(": "LPAREN",
      ")": "RPAREN",
      "[": "LBRACKET",
      "]": "RBRACKET",
      "{": "LBRACE",
      "}": "RBRACE",
      ";": "SEMICOLON",
      ",": "COMMA",
      "=": "ASSIGN",
      "+": "PLUS",
      "*": "STAR",
      "/": "SLASH",
      "<": "LT",
      ">": "GT",
      "&": "AMP",
      "|": "PIPE",
      "^": "CARET",
      "~": "TILDE",
      "%": "PERCENT",
      "?": "QUESTION",
      ":": "COLON"
    };
    if (typeMap[c]) {
      tokens.push(new Token3(typeMap[c], c, { line, col }));
      i++;
      col++;
      continue;
    }
    i++;
    col++;
  }
  tokens.push(new Token3("EOF", null, { line, col }));
  return tokens;
}
var QASM3Parser = class {
  constructor(source) {
    this.source = source;
    this.tokens = tokenize3(source);
    this.pos = 0;
    this.circuit = null;
    this.qregs = /* @__PURE__ */ new Map();
    this.cregs = /* @__PURE__ */ new Map();
    this.classicalVars = /* @__PURE__ */ new Map();
    this.gateDefinitions = /* @__PURE__ */ new Map();
  }
  peek() {
    return this.tokens[this.pos];
  }
  next() {
    return this.tokens[this.pos++];
  }
  expect(type) {
    const t = this.tokens[this.pos];
    if (t.type !== type) {
      throw new Error(`QASM3 parse error at line ${t.pos.line}: expected ${type} but got ${t.type} (${t.value})`);
    }
    this.pos++;
    return t;
  }
  accept(type) {
    if (this.tokens[this.pos].type === type) return this.tokens[this.pos++];
    return null;
  }
  parse() {
    const versionTok = this.expect("KEYWORD");
    if (versionTok.value !== "OPENQASM") {
      throw new Error(`Expected OPENQASM, got ${versionTok.value}`);
    }
    const verTok = this.next();
    this.expect("SEMICOLON");
    this.circuit = new QuantumCircuit();
    while (this.peek().type !== "EOF") {
      this.parseStatement();
    }
    return this.circuit;
  }
  parseStatement() {
    const t = this.peek();
    if (t.type === "KEYWORD") {
      switch (t.value) {
        case "include":
          this.parseInclude();
          return;
        case "qreg":
        case "qubit":
        case "qubits":
          this.parseQregDecl();
          return;
        case "creg":
        case "bit":
          this.parseCregDecl();
          return;
        case "gate":
          this.parseGateDef();
          return;
        case "def":
          this.parseDef();
          return;
        case "measure":
          this.parseMeasure();
          return;
        case "barrier":
          this.parseBarrier();
          return;
        case "reset":
          this.parseReset();
          return;
        case "if":
          this.parseIf();
          return;
        case "for":
          this.parseFor();
          return;
        case "while":
          this.parseWhile();
          return;
        case "let":
          this.parseLet();
          return;
        case "const":
          this.parseConst();
          return;
        case "opaque":
          this.parseOpaque();
          return;
      }
    }
    if (t.type === "IDENT") {
      const savedPos = this.pos;
      const name = this.next().value;
      if (this.accept("LBRACKET")) {
        const idx = parseInt(this.expect("NUMBER").value, 10);
        this.expect("RBRACKET");
        if (this.accept("ASSIGN") && this.peek().type === "KEYWORD" && this.peek().value === "measure") {
          this.expect("KEYWORD");
          const qubit = this.parseBitRef();
          this.expect("SEMICOLON");
          if (this.cregs.has(name)) {
            this.circuit.measure(qubit, this.cregs.get(name)._bits[idx]);
          }
          return;
        }
      }
      this.pos = savedPos;
      this.parseGateCall();
      return;
    }
    throw new Error(`QASM3 parse error at line ${t.pos.line}: unexpected ${t.type} (${t.value})`);
  }
  parseInclude() {
    this.expect("KEYWORD");
    this.expect("STRING");
    this.expect("SEMICOLON");
  }
  parseQregDecl() {
    this.expect("KEYWORD");
    let name, size;
    if (this.peek().type === "LBRACKET") {
      this.expect("LBRACKET");
      size = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
      name = this.expect("IDENT").value;
    } else {
      name = this.expect("IDENT").value;
      this.expect("LBRACKET");
      size = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
    }
    this.expect("SEMICOLON");
    const reg = new QuantumRegister(size, name);
    this.circuit.add_register(reg);
    this.qregs.set(name, reg);
  }
  parseCregDecl() {
    this.expect("KEYWORD");
    let name, size;
    if (this.peek().type === "LBRACKET") {
      this.expect("LBRACKET");
      size = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
      name = this.expect("IDENT").value;
    } else {
      name = this.expect("IDENT").value;
      this.expect("LBRACKET");
      size = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
    }
    this.expect("SEMICOLON");
    const reg = new ClassicalRegister(size, name);
    this.circuit.add_register(reg);
    this.cregs.set(name, reg);
  }
  parseGateDef() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    let params = [];
    if (this.accept("LBRACKET")) {
      while (this.peek().type !== "RBRACKET") {
        params.push(this.expect("IDENT").value);
        this.accept("COMMA");
      }
      this.expect("RBRACKET");
    }
    let qargs = [];
    while (this.peek().type === "IDENT") {
      qargs.push(this.expect("IDENT").value);
      this.accept("COMMA");
    }
    this.expect("LBRACE");
    const body = [];
    while (this.peek().type !== "RBRACE") {
      const gateName = this.expect("IDENT").value;
      let gateParams = [];
      if (this.accept("LPAREN")) {
        while (this.peek().type !== "RPAREN") {
          gateParams.push(this.parseExpression());
          this.accept("COMMA");
        }
        this.expect("RPAREN");
      }
      const gateQargs = [];
      while (this.peek().type === "IDENT") {
        const qName = this.expect("IDENT").value;
        if (this.accept("LBRACKET")) {
          const qIdx = parseInt(this.expect("NUMBER").value, 10);
          this.expect("RBRACKET");
          gateQargs.push({ name: qName, index: qIdx });
        } else {
          gateQargs.push({ name: qName, index: 0 });
        }
        this.accept("COMMA");
      }
      this.expect("SEMICOLON");
      body.push({ gate: gateName, params: gateParams, qargs: gateQargs });
    }
    this.expect("RBRACE");
    this.gateDefinitions.set(name, { params, qargs, body });
  }
  parseDef() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("LPAREN");
    const args = [];
    while (this.peek().type !== "RPAREN") {
      args.push(this.expect("IDENT").value);
      this.accept("COMMA");
    }
    this.expect("RPAREN");
    if (this.accept("ARROW")) {
      this.next();
    }
    this.expect("LBRACE");
    let depth = 1;
    while (depth > 0) {
      const t = this.next();
      if (t.type === "LBRACE") depth++;
      else if (t.type === "RBRACE") depth--;
      else if (t.type === "EOF") break;
    }
  }
  parseMeasure() {
    this.expect("KEYWORD");
    const qubit = this.parseBitRef();
    this.expect("ARROW");
    const clbit = this.parseBitRef();
    this.expect("SEMICOLON");
    this.circuit.measure(qubit, clbit);
  }
  parseBarrier() {
    this.expect("KEYWORD");
    const qubits = [];
    if (this.peek().type === "IDENT") {
      do {
        qubits.push(this.parseBitRef());
      } while (this.accept("COMMA"));
    }
    this.expect("SEMICOLON");
    if (qubits.length === 0) this.circuit.barrier();
    else this.circuit.barrier(qubits);
  }
  parseReset() {
    this.expect("KEYWORD");
    const qubit = this.parseBitRef();
    this.expect("SEMICOLON");
    this.circuit.reset(qubit);
  }
  parseIf() {
    this.expect("KEYWORD");
    this.expect("LPAREN");
    const cregName = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const cregIdx = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    const opTok = this.next();
    const op = opTok.type === "EQ" ? "==" : opTok.type === "NEQ" ? "!=" : opTok.type === "LT" ? "<" : opTok.type === "GT" ? ">" : opTok.type === "LE" ? "<=" : opTok.type === "GE" ? ">=" : "==";
    const value = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RPAREN");
    const creg = this.cregs.get(cregName);
    const condition = {
      register: cregName,
      index: cregIdx,
      op,
      value,
      bit: creg ? creg._bits[cregIdx] : null
    };
    const dataLenBefore = this.circuit.data.length;
    if (this.accept("LBRACE")) {
      while (this.peek().type !== "RBRACE") {
        this.parseStatement();
      }
      this.expect("RBRACE");
    } else {
      this.parseStatement();
    }
    for (let i = dataLenBefore; i < this.circuit.data.length; i++) {
      this.circuit.data[i].operation.condition = condition;
    }
    if (this.peek().type === "IDENT" && this.peek().value === "else") {
      this.next();
      const elseCondition = { ...condition, invert: true };
      const elseBefore = this.circuit.data.length;
      if (this.accept("LBRACE")) {
        while (this.peek().type !== "RBRACE") this.parseStatement();
        this.expect("RBRACE");
      } else {
        this.parseStatement();
      }
      for (let i = elseBefore; i < this.circuit.data.length; i++) {
        this.circuit.data[i].operation.condition = elseCondition;
      }
    }
  }
  parseFor() {
    this.expect("KEYWORD");
    this.expect("LPAREN");
    this.next();
    const varName = this.expect("IDENT").value;
    this.expect("KEYWORD");
    const values = [];
    if (this.accept("LBRACKET")) {
      const start = parseFloat(this.expect("NUMBER").value);
      let stop, step = 1;
      if (this.accept("COLON")) {
        const next = parseFloat(this.expect("NUMBER").value);
        if (this.accept("COLON")) {
          step = next;
          stop = parseFloat(this.expect("NUMBER").value);
        } else {
          stop = next;
        }
      }
      this.expect("RBRACKET");
      if (Number.isInteger(start) && Number.isInteger(step) && Number.isInteger(stop)) {
        if (step > 0) {
          for (let v = start; v < stop; v += step) values.push(v);
        } else {
          for (let v = start; v > stop; v += step) values.push(v);
        }
      } else {
        const n = Math.ceil(Math.abs((stop - start) / step));
        for (let i = 0; i < n; i++) values.push(start + i * step);
      }
    } else if (this.accept("LBRACE")) {
      while (this.peek().type !== "RBRACE") {
        values.push(this.parseExpression());
        if (!this.accept("COMMA")) break;
      }
      this.expect("RBRACE");
    } else {
      values.push(this.parseExpression());
    }
    this.expect("RPAREN");
    const bodyStartPos = this.pos;
    if (this.accept("LBRACE")) {
      const braceStartPos = this.pos;
      for (const v of values) {
        this.classicalVars.set(varName, v);
        if (v !== values[0]) {
          this.pos = braceStartPos;
        }
        while (this.peek().type !== "RBRACE") {
          this.parseStatement();
        }
      }
      this.expect("RBRACE");
    } else {
      const stmtStartPos = this.pos;
      for (const v of values) {
        this.classicalVars.set(varName, v);
        if (v !== values[0]) {
          this.pos = stmtStartPos;
        }
        this.parseStatement();
      }
    }
    this.classicalVars.delete(varName);
  }
  parseWhile() {
    this.expect("KEYWORD");
    this.expect("LPAREN");
    const cregName = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const cregIdx = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("EQ");
    const value = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RPAREN");
    const MAX_ITERATIONS = 1e3;
    const bodyStartPos = this.pos;
    if (this.accept("LBRACE")) {
      const braceStartPos = this.pos;
      for (let iter = 0; iter < MAX_ITERATIONS; iter++) {
        const dataLenBefore = this.circuit.data.length;
        if (iter > 0) this.pos = braceStartPos;
        while (this.peek().type !== "RBRACE") {
          this.parseStatement();
        }
        const creg = this.cregs.get(cregName);
        const condition = {
          register: cregName,
          index: cregIdx,
          op: "==",
          value,
          bit: creg ? creg._bits[cregIdx] : null,
          isWhile: true
        };
        for (let i = dataLenBefore; i < this.circuit.data.length; i++) {
          this.circuit.data[i].operation.condition = condition;
        }
      }
      this.expect("RBRACE");
    } else {
      const stmtStartPos = this.pos;
      for (let iter = 0; iter < MAX_ITERATIONS; iter++) {
        if (iter > 0) this.pos = stmtStartPos;
        this.parseStatement();
      }
    }
  }
  parseLet() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("ASSIGN");
    const value = this.parseExpression();
    this.expect("SEMICOLON");
    this.classicalVars.set(name, value);
  }
  parseConst() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("ASSIGN");
    const value = this.parseExpression();
    this.expect("SEMICOLON");
    this.classicalVars.set(name, value);
  }
  parseOpaque() {
    this.expect("KEYWORD");
    while (this.peek().type !== "SEMICOLON") this.next();
    this.expect("SEMICOLON");
  }
  parseGateCall() {
    const name = this.expect("IDENT").value;
    let params = [];
    if (this.accept("LPAREN")) {
      while (this.peek().type !== "RPAREN") {
        params.push(this.parseExpression());
        this.accept("COMMA");
      }
      this.expect("RPAREN");
    }
    const qargs = [];
    const cargs = [];
    while (this.peek().type === "IDENT" || this.peek().type === "KEYWORD") {
      const ref = this.parseBitRef();
      if (ref instanceof Clbit || ref && ref.register instanceof ClassicalRegister) {
        cargs.push(ref);
      } else {
        qargs.push(ref);
      }
      if (!this.accept("COMMA")) break;
    }
    this.expect("SEMICOLON");
    this._applyGateByName(name, params, qargs, cargs);
  }
  _applyGateByName(name, params, qargs, cargs) {
    if (this.gateDefinitions.has(name)) {
      const def = this.gateDefinitions.get(name);
      const argMap = /* @__PURE__ */ new Map();
      for (let i = 0; i < def.qargs.length && i < qargs.length; i++) {
        argMap.set(def.qargs[i], qargs[i]);
      }
      for (const stmt of def.body) {
        const mappedQargs = stmt.qargs.map((q) => argMap.get(q.name) || qargs[0]);
        this._applyGateByName(stmt.gate, stmt.params, mappedQargs, []);
      }
      return;
    }
    const n = name.toLowerCase();
    switch (n) {
      case "h":
        this.circuit.h(qargs[0]);
        break;
      case "x":
        this.circuit.x(qargs[0]);
        break;
      case "y":
        this.circuit.y(qargs[0]);
        break;
      case "z":
        this.circuit.z(qargs[0]);
        break;
      case "s":
        this.circuit.s(qargs[0]);
        break;
      case "sdg":
        this.circuit.sdg(qargs[0]);
        break;
      case "t":
        this.circuit.t(qargs[0]);
        break;
      case "tdg":
        this.circuit.tdg(qargs[0]);
        break;
      case "sx":
        this.circuit.sx(qargs[0]);
        break;
      case "id":
        this.circuit.id(qargs[0]);
        break;
      case "u":
      case "u3":
        this.circuit.u3(params[0], params[1], params[2], qargs[0]);
        break;
      case "u2":
        this.circuit.u2(params[0], params[1], qargs[0]);
        break;
      case "u1":
      case "p":
        this.circuit.u1(params[0], qargs[0]);
        break;
      case "rx":
        this.circuit.rx(params[0], qargs[0]);
        break;
      case "ry":
        this.circuit.ry(params[0], qargs[0]);
        break;
      case "rz":
        this.circuit.rz(params[0], qargs[0]);
        break;
      case "cx":
        this.circuit.cx(qargs[0], qargs[1]);
        break;
      case "cy":
        this.circuit.cy(qargs[0], qargs[1]);
        break;
      case "cz":
        this.circuit.cz(qargs[0], qargs[1]);
        break;
      case "ch":
        this.circuit.ch(qargs[0], qargs[1]);
        break;
      case "swap":
        this.circuit.swap(qargs[0], qargs[1]);
        break;
      case "iswap":
        this.circuit.iswap(qargs[0], qargs[1]);
        break;
      case "ccx":
        this.circuit.ccx(qargs[0], qargs[1], qargs[2]);
        break;
      case "cswap":
        this.circuit.cswap(qargs[0], qargs[1], qargs[2]);
        break;
      case "crx":
        this.circuit.crx(params[0], qargs[0], qargs[1]);
        break;
      case "cry":
        this.circuit.cry(params[0], qargs[0], qargs[1]);
        break;
      case "crz":
        this.circuit.crz(params[0], qargs[0], qargs[1]);
        break;
      case "cp":
      case "cu1":
        this.circuit.cp(params[0], qargs[0], qargs[1]);
        break;
      case "rxx":
        this.circuit.rxx(params[0], qargs[0], qargs[1]);
        break;
      case "ryy":
        this.circuit.ryy(params[0], qargs[0], qargs[1]);
        break;
      case "rzz":
        this.circuit.rzz(params[0], qargs[0], qargs[1]);
        break;
      default: {
        const instr = new Instruction(n, qargs.length, cargs.length, params);
        this.circuit.append(instr, qargs, cargs);
      }
    }
  }
  parseBitRef() {
    const name = this.expect("IDENT").value;
    if (this.accept("LBRACKET")) {
      let idx;
      const t = this.peek();
      if (t.type === "NUMBER") {
        idx = parseInt(this.next().value, 10);
      } else if (t.type === "IDENT" && this.classicalVars.has(t.value)) {
        idx = this.classicalVars.get(t.value);
        this.next();
        if (!Number.isInteger(idx)) {
          throw new Error(`QASM3: register index ${t.value} = ${idx} is not an integer`);
        }
      } else {
        throw new Error(`QASM3: expected NUMBER or classical var in [], got ${t.type} (${t.value})`);
      }
      this.expect("RBRACKET");
      if (this.qregs.has(name)) return this.qregs.get(name)._bits[idx];
      if (this.cregs.has(name)) return this.cregs.get(name)._bits[idx];
      throw new Error(`Unknown register: ${name}`);
    }
    if (this.qregs.has(name)) return this.qregs.get(name)._bits[0];
    if (this.cregs.has(name)) return this.cregs.get(name)._bits[0];
    throw new Error(`Unknown register: ${name}`);
  }
  parseExpression() {
    const t = this.peek();
    if (t.type === "NUMBER") {
      this.next();
      let value = parseFloat(t.value);
      while (this.peek().type === "PLUS" || this.peek().type === "STAR" || this.peek().type === "SLASH") {
        const op = this.next().type;
        const right = this.expect("NUMBER").value;
        if (op === "PLUS") value += parseFloat(right);
        else if (op === "STAR") value *= parseFloat(right);
        else if (op === "SLASH") value /= parseFloat(right);
      }
      return value;
    }
    if (t.type === "IDENT") {
      this.next();
      if (this.classicalVars.has(t.value)) {
        return this.classicalVars.get(t.value);
      }
      return t.value;
    }
    throw new Error(`QASM3 parse error: expected expression, got ${t.type}`);
  }
};
function qasm3_parse(source) {
  const parser = new QASM3Parser(source);
  return parser.parse();
}

// src/primitives/primitives.js
var BaseEstimator = class {
  constructor(options = {}) {
    this.options = options;
    this._circuits = /* @__PURE__ */ new Map();
  }
  run(circuits2, observables, parameterValues = null, options = {}) {
    throw new Error("BaseEstimator.run not implemented");
  }
};
var BaseSampler = class {
  constructor(options = {}) {
    this.options = options;
    this._circuits = /* @__PURE__ */ new Map();
  }
};
var EstimatorResult = class {
  constructor(values, metadata = []) {
    this.values = values;
    this.metadata = metadata;
  }
};
var SamplerResult = class {
  constructor(quasi_dists, metadata = []) {
    this.quasi_dists = quasi_dists;
    this.metadata = metadata;
  }
};
var Estimator = class extends BaseEstimator {
  constructor(options = {}) {
    super(options);
    this.shots = options.shots || 1024;
    this.backend = options.backend || "statevector";
  }
  // Compute expectation values for (circuit, observable) pairs.
  // circuits: QuantumCircuit or array of QuantumCircuit
  // observables: SparsePauliOp or Pauli or array of these (same length as circuits)
  // parameter_values: array of {paramName: number} dicts, or null
  run(circuits2, observables, parameterValues = null, options = {}) {
    const circuitList = Array.isArray(circuits2) ? circuits2 : [circuits2];
    const obsList = Array.isArray(observables) ? observables : [observables];
    if (circuitList.length !== obsList.length) {
      throw new Error("Estimator: number of circuits must match number of observables");
    }
    const values = [];
    const metadata = [];
    for (let i = 0; i < circuitList.length; i++) {
      let circuit = circuitList[i];
      let observable = obsList[i];
      if (parameterValues && parameterValues[i]) {
        circuit = circuit.bind_parameters(parameterValues[i]);
      }
      const sv = Statevector.fromCircuit(circuit);
      let spo;
      if (observable instanceof SparsePauliOp) {
        spo = observable;
      } else if (observable instanceof Pauli) {
        spo = SparsePauliOp.from_list([[observable.label, 1]]);
      } else if (typeof observable === "string") {
        spo = SparsePauliOp.from_list([[observable, 1]]);
      } else {
        throw new TypeError("Observable must be Pauli, SparsePauliOp, or string");
      }
      const expVal = spo.expectation_value(sv);
      values.push(expVal);
      metadata.push({ shots: 0 });
    }
    return new EstimatorResult(values, metadata);
  }
};
var Sampler = class extends BaseSampler {
  constructor(options = {}) {
    super(options);
    this.shots = options.shots || 1024;
    this.backend = options.backend || "statevector";
  }
  // Sample measurement outcomes from circuits.
  run(circuits2, parameterValues = null, options = {}) {
    const circuitList = Array.isArray(circuits2) ? circuits2 : [circuits2];
    const shots = options.shots || this.shots;
    const quasiDists = [];
    const metadata = [];
    for (let i = 0; i < circuitList.length; i++) {
      let circuit = circuitList[i];
      if (parameterValues && parameterValues[i]) {
        circuit = circuit.bind_parameters(parameterValues[i]);
      }
      const hasMeasure = circuit.data.some((ci) => ci.operation.name === "measure");
      if (hasMeasure) {
        const result = simulate(circuit, shots);
        const counts = result.get_counts();
        const quasi = {};
        const total = counts.shots;
        for (const [key, val] of counts.items()) {
          quasi[key] = val / total;
        }
        quasiDists.push(quasi);
        metadata.push({ shots });
      } else {
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
};

// src/opflow/opflow.js
var OperatorBase = class {
  constructor() {
    this.num_qubits = 0;
    this.coeff = Complex.ONE;
  }
  add(other) {
    return new ListOp([this, other], "add");
  }
  sub(other) {
    return new ListOp([this, other.scale(-1)], "add");
  }
  compose(other) {
    return new ListOp([this, other], "compose");
  }
  tensor(other) {
    return new ListOp([this, other], "tensor");
  }
  expand(other) {
    return new ListOp([other, this], "tensor");
  }
  pow(n) {
    return new ListOp([this], "pow", { exponent: n });
  }
  scale(c) {
    if (typeof c === "number") c = new Complex(c, 0);
    return new ListOp([this], "scale", { coeff: c });
  }
  adjoint() {
    return new ListOp([this], "adjoint");
  }
  to_matrix() {
    throw new Error("OperatorBase.to_matrix not implemented");
  }
  to_matrix_op() {
    return new MatrixOp(this.to_matrix());
  }
  to_spmatrix_op() {
    return new PauliSumOp(SparsePauliOp.from_operator(new Operator(this.to_matrix())));
  }
  // Apply to a state
  apply_to(state) {
    const m = this.to_matrix();
    return m.matvec(state.data || state);
  }
  is_hermitian() {
    return this.to_matrix().isHermitian();
  }
};
var PauliOp = class _PauliOp extends OperatorBase {
  constructor(pauli, coeff = Complex.ONE) {
    super();
    this.primitive = pauli instanceof Pauli ? pauli : new Pauli(pauli);
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.num_qubits = this.primitive.num_qubits;
  }
  to_matrix() {
    return this.primitive.to_matrix().scale(this.coeff);
  }
  adjoint() {
    return new _PauliOp(this.primitive, this.coeff.conjugate());
  }
  is_hermitian() {
    return Math.abs(this.coeff.im) < 1e-9;
  }
  toString() {
    return `${this.coeff.toString()} * ${this.primitive.label}`;
  }
};
var PauliSumOp = class _PauliSumOp extends OperatorBase {
  constructor(spmatrixop, coeff = Complex.ONE) {
    super();
    this.primitive = spmatrixop instanceof SparsePauliOp ? spmatrixop : SparsePauliOp.from_list(spmatrixop);
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.num_qubits = this.primitive.num_qubits;
  }
  static from_list(list) {
    return new _PauliSumOp(SparsePauliOp.from_list(list));
  }
  to_matrix() {
    return this.primitive.to_matrix().scale(this.coeff);
  }
  adjoint() {
    return new _PauliSumOp(this.primitive.adjoint(), this.coeff.conjugate());
  }
  is_hermitian() {
    return this.to_matrix().isHermitian();
  }
  // Reduce: combine like terms
  reduce() {
    return new _PauliSumOp(this.primitive.simplify(), this.coeff);
  }
  // Convert to a list of [Pauli, coeff] pairs
  to_list() {
    return this.primitive.to_list();
  }
  toString() {
    return this.primitive.toString();
  }
};
var MatrixOp = class _MatrixOp extends OperatorBase {
  constructor(matrix, coeff = Complex.ONE) {
    super();
    if (matrix instanceof ComplexMatrix2) {
      this.primitive = matrix;
    } else if (matrix instanceof Operator) {
      this.primitive = matrix.data;
    } else {
      this.primitive = matrix;
    }
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    const n = Math.log2(this.primitive.rows);
    this.num_qubits = n;
  }
  to_matrix() {
    return this.primitive.scale(this.coeff);
  }
  adjoint() {
    return new _MatrixOp(this.primitive.dagger(), this.coeff.conjugate());
  }
  is_hermitian() {
    return this.to_matrix().isHermitian();
  }
  toString() {
    return `MatrixOp(num_qubits=${this.num_qubits})`;
  }
};
var CircuitOp = class _CircuitOp extends OperatorBase {
  constructor(circuit, coeff = Complex.ONE) {
    super();
    this.primitive = circuit;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.num_qubits = circuit.num_qubits;
  }
  to_matrix() {
    const op = Operator.fromCircuit(this.primitive);
    return op.data.scale(this.coeff);
  }
  to_circuit() {
    return this.primitive;
  }
  adjoint() {
    return new _CircuitOp(this.primitive.inverse(), this.coeff.conjugate());
  }
  toString() {
    return `CircuitOp(${this.primitive.name})`;
  }
};
var StateFn = class _StateFn extends OperatorBase {
  constructor(primitive, is_measurement = false, coeff = Complex.ONE) {
    super();
    if (primitive instanceof Statevector) {
      this.primitive = primitive;
    } else if (primitive instanceof ComplexVector) {
      this.primitive = new Statevector(primitive);
    } else if (typeof primitive === "string") {
      this.primitive = Statevector.fromLabel(primitive);
    } else {
      this.primitive = primitive;
    }
    this.is_measurement = is_measurement;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.num_qubits = this.primitive.num_qubits;
  }
  static from_label(label, is_measurement = false) {
    return new _StateFn(Statevector.fromLabel(label), is_measurement);
  }
  to_matrix() {
    if (this.is_measurement) {
      const dim2 = this.primitive.dim;
      const m2 = ComplexMatrix2.zeros(1, dim2);
      for (let i = 0; i < dim2; i++) {
        m2.set(0, i, this.primitive.data.get(i).conjugate().mul(this.coeff));
      }
      return m2;
    }
    const dim = this.primitive.dim;
    const m = ComplexMatrix2.zeros(dim, 1);
    for (let i = 0; i < dim; i++) {
      m.set(i, 0, this.primitive.data.get(i).mul(this.coeff));
    }
    return m;
  }
  adjoint() {
    return new _StateFn(this.primitive, !this.is_measurement, this.coeff.conjugate());
  }
  // Apply an operator to the state
  compose(other) {
    if (other instanceof OperatorBase) {
      const newState = this.primitive.evolve(other.to_matrix());
      return new _StateFn(newState, this.is_measurement, this.coeff);
    }
    return super.compose(other);
  }
  toString() {
    const prefix = this.is_measurement ? "~" : "";
    return `${prefix}StateFn(${this.primitive.toString()})`;
  }
};
var CircuitStateFn = class extends StateFn {
  constructor(circuit, is_measurement = false, coeff = Complex.ONE) {
    const sv = Statevector.fromCircuit(circuit);
    super(sv, is_measurement, coeff);
    this.circuit = circuit;
  }
  to_circuit() {
    return this.circuit;
  }
  toString() {
    return `CircuitStateFn(${this.circuit.name})`;
  }
};
var ListOp = class _ListOp extends OperatorBase {
  constructor(oplist, comboFn = "add", auxFields = {}) {
    super();
    this.oplist = oplist;
    this.combo_fn = comboFn;
    this.aux_fields = auxFields;
    if (oplist.length > 0) this.num_qubits = oplist[0].num_qubits;
    this.coeff = auxFields.coeff || Complex.ONE;
  }
  to_matrix() {
    if (this.oplist.length === 0) return ComplexMatrix2.identity(1);
    const matrices = this.oplist.map((op) => op.to_matrix());
    let result;
    switch (this.combo_fn) {
      case "add":
        result = matrices[0];
        for (let i = 1; i < matrices.length; i++) result = result.add(matrices[i]);
        break;
      case "compose":
        result = matrices[0];
        for (let i = 1; i < matrices.length; i++) result = matrices[i].mul(result);
        break;
      case "tensor":
        result = matrices[0];
        for (let i = 1; i < matrices.length; i++) result = result.tensor(matrices[i]);
        break;
      case "adjoint":
        result = matrices[0].dagger();
        break;
      case "scale":
        result = matrices[0].scale(this.aux_fields.coeff || Complex.ONE);
        break;
      case "pow":
        result = matrices[0];
        const n = this.aux_fields.exponent || 1;
        for (let i = 1; i < n; i++) result = result.mul(matrices[0]);
        break;
      default:
        throw new Error(`ListOp: unknown combo function ${this.combo_fn}`);
    }
    return result.scale(this.coeff);
  }
  adjoint() {
    return new _ListOp(
      this.oplist.map((op) => op.adjoint()),
      this.combo_fn,
      Object.assign({}, this.aux_fields, { coeff: this.coeff.conjugate() })
    );
  }
  // Reduce: recursively simplify children
  reduce() {
    const reduced = this.oplist.map((op) => op.reduce ? op.reduce() : op);
    return new _ListOp(reduced, this.combo_fn, this.aux_fields);
  }
  toString() {
    const opStr = this.oplist.map((op) => op.toString()).join(` ${this.combo_fn} `);
    return `(${opStr})`;
  }
};
var EvolvedOp = class _EvolvedOp extends OperatorBase {
  constructor(primitive, coeff = Complex.ONE, time = 1) {
    super();
    this.primitive = primitive;
    this.coeff = coeff instanceof Complex ? coeff : new Complex(coeff, 0);
    this.time = time;
    this.num_qubits = primitive.num_qubits;
  }
  to_matrix() {
    const O = this.primitive.to_matrix();
    const expMatrix = O.scale(new Complex(0, -this.time));
    return expMatrix.expm().scale(this.coeff);
  }
  adjoint() {
    return new _EvolvedOp(this.primitive, this.coeff.conjugate(), -this.time);
  }
  toString() {
    return `EvolvedOp(e^(-i ${this.time} ${this.primitive}))`;
  }
};

// src/opflow/opflow_expec.js
var MatrixExpectation = class {
  constructor() {
  }
  // Compute <ψ|O|ψ> for a CircuitStateFn and operator
  compute(stateFn, operator) {
    const sv = stateFn.primitive || stateFn;
    const statevector = sv instanceof Statevector ? sv : Statevector.fromCircuit(sv);
    const op = operator.primitive || operator;
    const matrix = op instanceof SparsePauliOp ? op.to_matrix() : op.to_matrix();
    const opResult = matrix.matvec(statevector.data);
    const result = statevector.data.inner(opResult);
    return result;
  }
  // Compute expectation for a list of (state, operator) pairs
  compute_list(stateFns, operators) {
    return stateFns.map((sf, i) => this.compute(sf, operators[i]));
  }
};
var PauliExpectation = class {
  constructor(grouping = false) {
    this.grouping = grouping;
  }
  compute(stateFn, operator) {
    const sv = stateFn.primitive || stateFn;
    const statevector = sv instanceof Statevector ? sv : Statevector.fromCircuit(sv);
    const op = operator.primitive || operator;
    const spo = op instanceof SparsePauliOp ? op : SparsePauliOp.from_list(op);
    let result = Complex.ZERO;
    for (let i = 0; i < spo.paulis.size; i++) {
      const pauli = spo.paulis.get(i);
      const coeff = spo.coeffs[i];
      const pauliMatrix = pauli.to_matrix();
      const opResult = pauliMatrix.matvec(statevector.data);
      const expP = statevector.data.inner(opResult);
      result = result.add(expP.mul(coeff));
    }
    return result;
  }
  compute_list(stateFns, operators) {
    return stateFns.map((sf, i) => this.compute(sf, operators[i]));
  }
};
var CircuitSampler = class {
  constructor(backend = null, shots = 1024) {
    this.backend = backend;
    this.shots = shots;
  }
  // Sample a state from a circuit (returns a Statevector or samples)
  sample(circuit, parameters = null) {
    let qc = circuit;
    if (parameters) {
      qc = circuit.bind_parameters(parameters);
    }
    if (this.backend) {
      const result = this.backend.run(qc, this.shots);
      return result;
    }
    return Statevector.fromCircuit(qc);
  }
  // Convert a CircuitStateFn to a DictStateFn (sampled state)
  convert(stateFn, parameters = null) {
    const circuit = stateFn.circuit || stateFn.primitive;
    return this.sample(circuit, parameters);
  }
};
function get_expectation(operator, backend = null) {
  if (backend) {
    return new PauliExpectation();
  }
  return new MatrixExpectation();
}

// src/algorithms/optimizers.js
var OptimizerResult = class {
  constructor(kwargs = {}) {
    this.x = kwargs.x || [];
    this.fun = kwargs.fun || 0;
    this.nfev = kwargs.nfev || 0;
    this.nit = kwargs.nit || 0;
    this.success = kwargs.success !== void 0 ? kwargs.success : true;
    this.message = kwargs.message || "Optimization finished";
    this.history = kwargs.history || [];
  }
};
var GradientDescent = class {
  constructor(options = {}) {
    this.learning_rate = options.learning_rate || 0.01;
    this.tolerance = options.tolerance || 1e-6;
    this.max_iter = options.maxiter || 100;
  }
  minimize(objectiveFn, x0, gradientFn = null) {
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const history = [{ x: x.slice(), fun: fx }];
    let iter = 0;
    while (iter < this.max_iter) {
      let grad;
      if (gradientFn) {
        grad = gradientFn(x);
      } else {
        grad = this._numericalGradient(objectiveFn, x);
        nfev += 2 * x.length;
      }
      const newX = x.map((xi, i) => xi - this.learning_rate * grad[i]);
      const newFx = objectiveFn(newX);
      nfev++;
      history.push({ x: newX.slice(), fun: newFx });
      if (Math.abs(newFx - fx) < this.tolerance) {
        x = newX;
        fx = newFx;
        break;
      }
      x = newX;
      fx = newFx;
      iter++;
    }
    return new OptimizerResult({
      x,
      fun: fx,
      nfev,
      nit: iter,
      success: true,
      message: `GradientDescent finished after ${iter} iterations`,
      history
    });
  }
  _numericalGradient(f, x, eps = 1e-6) {
    const grad = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
      const xp = x.slice();
      xp[i] += eps;
      const xm = x.slice();
      xm[i] -= eps;
      grad[i] = (f(xp) - f(xm)) / (2 * eps);
    }
    return grad;
  }
};
var SPSA = class {
  constructor(options = {}) {
    this.maxiter = options.maxiter || 100;
    this.learning_rate = options.learning_rate || null;
    this.perturbation = options.perturbation || null;
    this.tolerance = options.tolerance || 1e-6;
    this.seed = options.seed != null ? options.seed : null;
    this.trust_region = options.trust_region || false;
    this.max_evals_grouped = options.max_evals_grouped || 1;
  }
  _makeRng() {
    if (this.seed == null) return Math.random;
    let s = this.seed >>> 0;
    return () => {
      s = s + 1831565813 | 0;
      let t = Math.imul(s ^ s >>> 15, 1 | s);
      t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
  }
  minimize(objectiveFn, x0) {
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const rng = this._makeRng();
    const history = [{ x: x.slice(), fun: fx }];
    const a = this.learning_rate || 0.1;
    const c = this.perturbation || 0.1;
    const A = this.maxiter / 10;
    let bestX = x.slice();
    let bestFx = fx;
    for (let k = 0; k < this.maxiter; k++) {
      const delta = x.map(() => rng() < 0.5 ? -1 : 1);
      const ck = c / Math.pow(k + 1 + A, 0.101);
      const xp = x.map((xi, i) => xi + ck * delta[i]);
      const xm = x.map((xi, i) => xi - ck * delta[i]);
      const fp = objectiveFn(xp);
      const fm = objectiveFn(xm);
      nfev += 2;
      const grad = x.map((_, i) => (fp - fm) / (2 * ck * delta[i]));
      const ak = a / Math.pow(k + 1 + A, 0.602);
      const newX = x.map((xi, i) => xi - ak * grad[i]);
      const newFx = objectiveFn(newX);
      nfev++;
      history.push({ x: newX.slice(), fun: newFx });
      if (newFx < bestFx) {
        bestFx = newFx;
        bestX = newX.slice();
      }
      if (Math.abs(newFx - fx) < this.tolerance) {
        x = newX;
        fx = newFx;
        break;
      }
      x = newX;
      fx = newFx;
    }
    return new OptimizerResult({
      x: bestX,
      fun: bestFx,
      nfev,
      nit: history.length - 1,
      success: true,
      message: `SPSA finished after ${history.length - 1} iterations`,
      history
    });
  }
};
var COBYLA = class {
  constructor(options = {}) {
    this.maxiter = options.maxiter || 100;
    this.tolerance = options.tol || 1e-4;
    this.rhobeg = options.rhobeg || 1;
    this.rhoend = options.rhoend || this.tolerance;
    this.constraints = options.constraints || null;
  }
  minimize(objectiveFn, x0) {
    const n = x0.length;
    const rho = this.rhobeg;
    const simplex = [x0.slice()];
    for (let i = 0; i < n; i++) {
      const v = x0.slice();
      v[i] += 0.5 * rho;
      simplex.push(v);
    }
    let fvals = simplex.map((p) => this._eval(objectiveFn, p));
    let nfev = simplex.length;
    const history = [{ x: x0.slice(), fun: Math.min(...fvals) }];
    let curRho = rho;
    let iter = 0;
    while (iter < this.maxiter && curRho > this.rhoend) {
      const order = fvals.map((v, i) => i).sort((a, b) => fvals[a] - fvals[b]);
      const best = order[0];
      const worst = order[n];
      const secondWorst = order[n - 1];
      const centroid = new Array(n).fill(0);
      for (let i = 0; i <= n; i++) {
        if (i === worst) continue;
        for (let j = 0; j < n; j++) centroid[j] += simplex[i][j];
      }
      for (let j = 0; j < n; j++) centroid[j] /= n;
      const reflected = centroid.map((c, j) => c + (c - simplex[worst][j]));
      const fr = this._eval(objectiveFn, reflected);
      nfev++;
      let nextVertex, nextF;
      if (fr < fvals[best]) {
        const expanded = centroid.map((c, j) => c + 2 * (reflected[j] - c));
        const fe = this._eval(objectiveFn, expanded);
        nfev++;
        if (fe < fr) {
          nextVertex = expanded;
          nextF = fe;
        } else {
          nextVertex = reflected;
          nextF = fr;
        }
      } else if (fr < fvals[secondWorst]) {
        nextVertex = reflected;
        nextF = fr;
      } else {
        const contractBest = fr < fvals[worst];
        const contractPoint = contractBest ? reflected : simplex[worst];
        const contracted = centroid.map((c, j) => c + 0.5 * (contractPoint[j] - c));
        const fc = this._eval(objectiveFn, contracted);
        nfev++;
        if (fc < Math.min(fr, fvals[worst])) {
          nextVertex = contracted;
          nextF = fc;
        } else {
          for (let i = 0; i <= n; i++) {
            if (i === best) continue;
            for (let j = 0; j < n; j++) {
              simplex[i][j] = simplex[best][j] + 0.5 * (simplex[i][j] - simplex[best][j]);
            }
            fvals[i] = this._eval(objectiveFn, simplex[i]);
            nfev++;
          }
          nextVertex = null;
          nextF = null;
        }
      }
      if (nextVertex !== null) {
        simplex[worst] = nextVertex;
        fvals[worst] = nextF;
      }
      const spread = Math.max(...fvals) - Math.min(...fvals);
      if (spread < curRho * curRho) {
        curRho = Math.max(this.rhoend, curRho * 0.5);
      }
      const curBest = fvals.indexOf(Math.min(...fvals));
      history.push({ x: simplex[curBest].slice(), fun: fvals[curBest] });
      if (spread < this.tolerance) break;
      iter++;
    }
    const finalBest = fvals.indexOf(Math.min(...fvals));
    return new OptimizerResult({
      x: simplex[finalBest],
      fun: fvals[finalBest],
      nfev,
      nit: iter,
      success: true,
      message: `COBYLA finished after ${iter} iterations (rho=${curRho.toExponential(3)})`,
      history
    });
  }
  // Evaluate the objective, returning +Infinity when constraints are violated.
  _eval(fn, x) {
    const f = fn(x);
    if (this.constraints) {
      for (const c of this.constraints) {
        const cv = c(x);
        if (cv < 0) return f + 1e6 * Math.abs(cv);
      }
    }
    return f;
  }
};

// src/algorithms/optimizers_extra.js
var Adam = class {
  constructor(options = {}) {
    this.learning_rate = options.learning_rate || 1e-3;
    this.beta1 = options.beta1 || 0.9;
    this.beta2 = options.beta2 || 0.999;
    this.epsilon = options.epsilon || 1e-8;
    this.max_iter = options.maxiter || 100;
    this.tolerance = options.tolerance || 1e-6;
  }
  minimize(objectiveFn, x0, gradientFn = null) {
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const m = new Array(x.length).fill(0);
    const v = new Array(x.length).fill(0);
    const history = [{ x: x.slice(), fun: fx }];
    let bestX = x.slice(), bestFx = fx;
    for (let t = 1; t <= this.max_iter; t++) {
      let grad;
      if (gradientFn) {
        grad = gradientFn(x);
      } else {
        grad = this._numericalGradient(objectiveFn, x);
        nfev += 2 * x.length;
      }
      for (let i = 0; i < x.length; i++) {
        m[i] = this.beta1 * m[i] + (1 - this.beta1) * grad[i];
        v[i] = this.beta2 * v[i] + (1 - this.beta2) * grad[i] * grad[i];
      }
      const mHat = m.map((mi) => mi / (1 - Math.pow(this.beta1, t)));
      const vHat = v.map((vi) => vi / (1 - Math.pow(this.beta2, t)));
      const newX = x.map((xi, i) => xi - this.learning_rate * mHat[i] / (Math.sqrt(vHat[i]) + this.epsilon));
      const newFx = objectiveFn(newX);
      nfev++;
      history.push({ x: newX.slice(), fun: newFx });
      if (newFx < bestFx) {
        bestFx = newFx;
        bestX = newX.slice();
      }
      if (Math.abs(newFx - fx) < this.tolerance) {
        x = newX;
        fx = newFx;
        break;
      }
      x = newX;
      fx = newFx;
    }
    return new OptimizerResult({
      x: bestX,
      fun: bestFx,
      nfev,
      nit: history.length - 1,
      success: true,
      message: `Adam finished after ${history.length - 1} iterations`,
      history
    });
  }
  _numericalGradient(f, x, eps = 1e-6) {
    const grad = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
      const xp = x.slice();
      xp[i] += eps;
      const xm = x.slice();
      xm[i] -= eps;
      grad[i] = (f(xp) - f(xm)) / (2 * eps);
    }
    return grad;
  }
};
var LBFGSB = class {
  constructor(options = {}) {
    this.max_iter = options.maxiter || 100;
    this.tolerance = options.tolerance || 1e-6;
    this.memory_size = options.memory_size || 10;
    this.bounds = options.bounds || null;
  }
  minimize(objectiveFn, x0, gradientFn = null) {
    let x = this._clipToBounds(x0.slice());
    let fx = objectiveFn(x);
    let nfev = 1;
    const history = [{ x: x.slice(), fun: fx }];
    const sList = [];
    const yList = [];
    const rhoList = [];
    for (let iter = 0; iter < this.max_iter; iter++) {
      const grad = gradientFn ? gradientFn(x) : this._numericalGradient(objectiveFn, x);
      if (!gradientFn) nfev += 2 * x.length;
      const q = grad.slice();
      const alpha = new Array(sList.length);
      for (let i = sList.length - 1; i >= 0; i--) {
        alpha[i] = rhoList[i] * this._dot(sList[i], q);
        for (let j = 0; j < q.length; j++) q[j] -= alpha[i] * yList[i][j];
      }
      for (let i = 0; i < sList.length; i++) {
        const beta = rhoList[i] * this._dot(yList[i], q);
        for (let j = 0; j < q.length; j++) q[j] += (alpha[i] - beta) * sList[i][j];
      }
      const direction = q.map((qi) => -qi);
      const stepSize = this._lineSearch(objectiveFn, x, direction, fx, grad);
      nfev += 5;
      const newX = this._clipToBounds(x.map((xi, i) => xi + stepSize * direction[i]));
      const newFx = objectiveFn(newX);
      nfev++;
      history.push({ x: newX.slice(), fun: newFx });
      const s = newX.map((ni, i) => ni - x[i]);
      const newGrad = gradientFn ? gradientFn(newX) : this._numericalGradient(objectiveFn, newX);
      if (!gradientFn) nfev += 2 * x.length;
      const y = newGrad.map((gi, i) => gi - grad[i]);
      const sy = this._dot(s, y);
      if (Math.abs(sy) > 1e-15) {
        sList.push(s);
        yList.push(y);
        rhoList.push(1 / sy);
        if (sList.length > this.memory_size) {
          sList.shift();
          yList.shift();
          rhoList.shift();
        }
      }
      if (Math.abs(newFx - fx) < this.tolerance) {
        x = newX;
        fx = newFx;
        break;
      }
      x = newX;
      fx = newFx;
    }
    return new OptimizerResult({
      x,
      fun: fx,
      nfev,
      nit: history.length - 1,
      success: true,
      message: `L-BFGS-B finished after ${history.length - 1} iterations`,
      history
    });
  }
  _clipToBounds(x) {
    if (!this.bounds) return x;
    return x.map((xi, i) => {
      if (this.bounds[i] && this.bounds[i][0] !== null && xi < this.bounds[i][0]) return this.bounds[i][0];
      if (this.bounds[i] && this.bounds[i][1] !== null && xi > this.bounds[i][1]) return this.bounds[i][1];
      return xi;
    });
  }
  _lineSearch(f, x, direction, fx0, grad0) {
    let alpha = 1;
    const c1 = 1e-4;
    for (let i = 0; i < 20; i++) {
      const newX = x.map((xi, j) => xi + alpha * direction[j]);
      const fx = f(newX);
      if (fx <= fx0 + c1 * alpha * this._dot(grad0, direction)) return alpha;
      alpha *= 0.5;
    }
    return alpha;
  }
  _dot(a, b) {
    return a.reduce((s, ai, i) => s + ai * b[i], 0);
  }
  _numericalGradient(f, x, eps = 1e-6) {
    const grad = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
      const xp = x.slice();
      xp[i] += eps;
      const xm = x.slice();
      xm[i] -= eps;
      grad[i] = (f(xp) - f(xm)) / (2 * eps);
    }
    return grad;
  }
};
var SLSQP = class {
  constructor(options = {}) {
    this.max_iter = options.maxiter || 100;
    this.tolerance = options.tolerance || 1e-6;
  }
  minimize(objectiveFn, x0, gradientFn = null) {
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const n = x.length;
    const history = [{ x: x.slice(), fun: fx }];
    let H = new Array(n);
    for (let i = 0; i < n; i++) {
      H[i] = new Array(n).fill(0);
      H[i][i] = 1;
    }
    for (let iter = 0; iter < this.max_iter; iter++) {
      const grad = gradientFn ? gradientFn(x) : this._numericalGradient(objectiveFn, x);
      if (!gradientFn) nfev += 2 * n;
      const direction = new Array(n).fill(0);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) direction[i] -= H[i][j] * grad[j];
      }
      let alpha = 1;
      for (let ls = 0; ls < 20; ls++) {
        const newX = x.map((xi, i) => xi + alpha * direction[i]);
        const newFx = objectiveFn(newX);
        nfev++;
        if (newFx < fx) {
          const newGrad = gradientFn ? gradientFn(newX) : this._numericalGradient(objectiveFn, newX);
          if (!gradientFn) nfev += 2 * n;
          const s = newX.map((ni, i) => ni - x[i]);
          const y = newGrad.map((gi, i) => gi - grad[i]);
          const sy = s.reduce((sum, si, i) => sum + si * y[i], 0);
          if (Math.abs(sy) > 1e-15) {
            const Hy = new Array(n).fill(0);
            for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Hy[i] += H[i][j] * y[j];
            const yHy = y.reduce((sum, yi, i) => sum + yi * Hy[i], 0);
            for (let i = 0; i < n; i++) {
              for (let j = 0; j < n; j++) {
                H[i][j] += (sy + yHy) / (sy * sy) * s[i] * s[j] - (Hy[i] * s[j] + s[i] * Hy[j]) / sy;
              }
            }
          }
          x = newX;
          fx = newFx;
          history.push({ x: x.slice(), fun: fx });
          break;
        }
        alpha *= 0.5;
      }
      if (Math.abs(history[history.length - 1].fun - (history[history.length - 2] || { fun: fx }).fun) < this.tolerance) break;
    }
    return new OptimizerResult({
      x,
      fun: fx,
      nfev,
      nit: history.length - 1,
      success: true,
      message: `SLSQP finished after ${history.length - 1} iterations`,
      history
    });
  }
  _numericalGradient(f, x, eps = 1e-6) {
    const grad = new Array(x.length);
    for (let i = 0; i < x.length; i++) {
      const xp = x.slice();
      xp[i] += eps;
      const xm = x.slice();
      xm[i] -= eps;
      grad[i] = (f(xp) - f(xm)) / (2 * eps);
    }
    return grad;
  }
};
var NFT = class {
  constructor(options = {}) {
    this.max_iter = options.maxiter || 100;
    this.tolerance = options.tolerance || 1e-6;
  }
  minimize(objectiveFn, x0) {
    const n = x0.length;
    let x = x0.slice();
    let fx = objectiveFn(x);
    let nfev = 1;
    const history = [{ x: x.slice(), fun: fx }];
    for (let iter = 0; iter < this.max_iter; iter++) {
      let improved = false;
      for (let i = 0; i < n; i++) {
        const f0 = fx;
        const xp = x.slice();
        xp[i] += Math.PI / 2;
        const fp = objectiveFn(xp);
        nfev++;
        const xp2 = x.slice();
        xp2[i] += Math.PI;
        const fp2 = objectiveFn(xp2);
        nfev++;
        const C = (f0 + fp2) / 2;
        const a = (f0 - fp2) / 2;
        const b = fp - C;
        const R = Math.sqrt(a * a + b * b);
        if (R < 1e-12) continue;
        const deltaTheta = Math.atan2(-b, -a);
        const newX = x.slice();
        newX[i] = x[i] + deltaTheta;
        const newFx = objectiveFn(newX);
        nfev++;
        if (newFx < fx) {
          x = newX;
          fx = newFx;
          improved = true;
          history.push({ x: x.slice(), fun: fx });
        }
      }
      if (!improved || Math.abs(history[history.length - 1].fun - (history[history.length - 2] || { fun: fx }).fun) < this.tolerance) {
        break;
      }
    }
    return new OptimizerResult({
      x,
      fun: fx,
      nfev,
      nit: history.length - 1,
      success: true,
      message: `NFT finished after ${history.length - 1} iterations`,
      history
    });
  }
};
var NelderMead = class {
  constructor(options = {}) {
    this.max_iter = options.maxiter || 200;
    this.tolerance = options.tolerance || 1e-6;
    this.initial_step = options.initial_step || 1;
    this.alpha = 1;
    this.gamma = 2;
    this.rho = 0.5;
    this.sigma = 0.5;
  }
  minimize(objectiveFn, x0) {
    const n = x0.length;
    let simplex = [x0.slice()];
    for (let i = 0; i < n; i++) {
      const point = x0.slice();
      point[i] += this.initial_step;
      simplex.push(point);
    }
    let fvals = simplex.map((p) => objectiveFn(p));
    let nfev = simplex.length;
    const history = [{ x: x0.slice(), fun: Math.min(...fvals) }];
    for (let iter = 0; iter < this.max_iter; iter++) {
      const order = fvals.map((v, i) => i).sort((a, b) => fvals[a] - fvals[b]);
      simplex = order.map((i) => simplex[i]);
      fvals = order.map((i) => fvals[i]);
      if (Math.abs(fvals[fvals.length - 1] - fvals[0]) < this.tolerance) break;
      const centroid = new Array(n).fill(0);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) centroid[j] += simplex[i][j];
      for (let j = 0; j < n; j++) centroid[j] /= n;
      const reflected = centroid.map((c, j) => c + this.alpha * (c - simplex[n][j]));
      const fr = objectiveFn(reflected);
      nfev++;
      if (fr >= fvals[0] && fr < fvals[n - 1]) {
        simplex[n] = reflected;
        fvals[n] = fr;
      } else if (fr < fvals[0]) {
        const expanded = centroid.map((c, j) => c + this.gamma * (c - simplex[n][j]));
        const fe = objectiveFn(expanded);
        nfev++;
        if (fe < fr) {
          simplex[n] = expanded;
          fvals[n] = fe;
        } else {
          simplex[n] = reflected;
          fvals[n] = fr;
        }
      } else {
        const contract = centroid.map((c, j) => c + this.rho * (c - simplex[n][j]));
        const fc = objectiveFn(contract);
        nfev++;
        if (fc < fvals[n]) {
          simplex[n] = contract;
          fvals[n] = fc;
        } else {
          for (let i = 1; i <= n; i++) {
            simplex[i] = simplex[0].map((v, j) => v + this.sigma * (simplex[i][j] - v));
            fvals[i] = objectiveFn(simplex[i]);
            nfev++;
          }
        }
      }
      history.push({ x: simplex[0].slice(), fun: fvals[0] });
    }
    return new OptimizerResult({
      x: simplex[0],
      fun: fvals[0],
      nfev,
      nit: history.length - 1,
      success: true,
      message: `Nelder-Mead finished after ${history.length - 1} iterations`,
      history
    });
  }
};

// src/algorithms/vqe.js
var VQEResult = class {
  constructor(kwargs = {}) {
    this.optimal_parameters = kwargs.optimal_parameters || {};
    this.optimal_value = kwargs.optimal_value || 0;
    this.optimal_circuit = kwargs.optimal_circuit || null;
    this.optimal_state = kwargs.optimal_state || null;
    this.cost_function_evals = kwargs.cost_function_evals || 0;
    this.optimizer_result = kwargs.optimizer_result || null;
    this.eigenvalue = this.optimal_value;
    this.eigenstate = this.optimal_state;
  }
};
var VQE = class {
  constructor(options = {}) {
    this.estimator = options.estimator || new Estimator();
    this.optimizer = options.optimizer || new SPSA({ maxiter: 100 });
    this.ansatz = options.ansatz || null;
    this.initial_point = options.initial_point || null;
    this.callback = options.callback || null;
    this.quantum_instance = options.quantum_instance || null;
  }
  // Compute the ground state of the given Hamiltonian.
  // operator: SparsePauliOp representing the Hamiltonian
  // ansatz: QuantumCircuit with parameters (overrides constructor's ansatz)
  // initial_point: array of numbers, the initial parameter values
  compute_minimum_eigenvalue(operator, options = {}) {
    const ansatz = options.ansatz || this.ansatz;
    if (!ansatz) throw new Error("VQE: ansatz is required");
    if (!operator) throw new Error("VQE: operator (Hamiltonian) is required");
    const paramSet = Array.from(ansatz.parameters);
    const paramNames = paramSet.map((p) => p.name);
    let x0;
    if (options.initial_point) {
      x0 = options.initial_point.slice();
    } else if (this.initial_point) {
      x0 = this.initial_point.slice();
    } else {
      x0 = paramNames.map(() => (Math.random() * 2 - 1) * Math.PI);
    }
    let evalCount = 0;
    const objectiveFn = (x) => {
      const params = {};
      for (let i = 0; i < paramNames.length; i++) {
        params[paramNames[i]] = x[i];
      }
      const boundAnsatz = ansatz.bind_parameters(params);
      const result = this.estimator.run(boundAnsatz, operator);
      const value = result.values[0];
      evalCount++;
      if (this.callback) {
        this.callback(evalCount, x, value);
      }
      return typeof value === "number" ? value : value.re;
    };
    const optResult = this.optimizer.minimize(objectiveFn, x0);
    const optimalParams = {};
    for (let i = 0; i < paramNames.length; i++) {
      optimalParams[paramNames[i]] = optResult.x[i];
    }
    const optimalCircuit = ansatz.bind_parameters(optimalParams);
    const optimalState = Statevector.fromCircuit(optimalCircuit);
    return new VQEResult({
      optimal_parameters: optimalParams,
      optimal_value: optResult.fun,
      optimal_circuit: optimalCircuit,
      optimal_state: optimalState,
      cost_function_evals: evalCount,
      optimizer_result: optResult
    });
  }
};

// src/algorithms/qaoa.js
var QAOAResult = class {
  constructor(kwargs = {}) {
    this.optimal_parameters = kwargs.optimal_parameters || {};
    this.optimal_value = kwargs.optimal_value || 0;
    this.optimal_circuit = kwargs.optimal_circuit || null;
    this.optimal_state = kwargs.optimal_state || null;
    this.cost_function_evals = kwargs.cost_function_evals || 0;
    this.eigenvalue = this.optimal_value;
    this.eigenstate = this.optimal_state;
  }
};
var QAOA = class {
  constructor(options = {}) {
    this.estimator = options.estimator || new Estimator();
    this.optimizer = options.optimizer || new SPSA({ maxiter: 100 });
    this.reps = options.reps || 1;
    this.initial_state = options.initial_state || null;
    this.mixer = options.mixer || null;
    this.callback = options.callback || null;
  }
  // Compute the minimum eigenvalue of `operator` (cost Hamiltonian).
  // operator: SparsePauliOp representing the cost Hamiltonian H_C
  compute_minimum_eigenvalue(operator, options = {}) {
    if (!operator) throw new Error("QAOA: cost operator is required");
    const numQubits = operator.num_qubits;
    let mixerOp;
    if (this.mixer) {
      mixerOp = this.mixer;
    } else {
      const mixerTerms = [];
      for (let i = 0; i < numQubits; i++) {
        const label = "I".repeat(i) + "X" + "I".repeat(numQubits - i - 1);
        mixerTerms.push([label, 1]);
      }
      mixerOp = SparsePauliOp.from_list(mixerTerms);
    }
    const ansatz = new QuantumCircuit(numQubits);
    if (this.initial_state) {
      ansatz.compose(this.initial_state, null, null, false, true);
    } else {
      for (let q = 0; q < numQubits; q++) ansatz.h(q);
    }
    ansatz.barrier();
    const gammaParams = [];
    const betaParams = [];
    for (let p = 0; p < this.reps; p++) {
      gammaParams.push(new Parameter(`gamma_${p}`));
      betaParams.push(new Parameter(`beta_${p}`));
    }
    for (let p = 0; p < this.reps; p++) {
      for (const [label, coeff] of operator.to_list()) {
        if (label === "I".repeat(numQubits)) continue;
        const cr = coeff instanceof Complex ? coeff.re : Number(coeff);
        const ci = coeff instanceof Complex ? coeff.im : 0;
        if (Math.abs(cr) < 1e-15 && Math.abs(ci) < 1e-15) continue;
        const timeExpr = gammaParams[p].mul(cr);
        const evolved = pauliEvolution(label, timeExpr);
        ansatz.compose(evolved, null, null, false, true);
      }
      for (let q = 0; q < numQubits; q++) {
        ansatz.rx(betaParams[p].mul(2), q);
      }
      ansatz.barrier();
    }
    const paramNames = gammaParams.concat(betaParams).map((p) => p.name);
    const x0 = new Array(paramNames.length).fill(0.5);
    let evalCount = 0;
    const objectiveFn = (x) => {
      const params = {};
      for (let i = 0; i < paramNames.length; i++) params[paramNames[i]] = x[i];
      const boundAnsatz = ansatz.bind_parameters(params);
      const result = this.estimator.run(boundAnsatz, operator);
      const value = result.values[0];
      evalCount++;
      if (this.callback) this.callback(evalCount, x, value);
      return typeof value === "number" ? value : value.re;
    };
    const optResult = this.optimizer.minimize(objectiveFn, x0);
    const optimalParams = {};
    for (let i = 0; i < paramNames.length; i++) {
      optimalParams[paramNames[i]] = optResult.x[i];
    }
    const optimalCircuit = ansatz.bind_parameters(optimalParams);
    const optimalState = Statevector.fromCircuit(optimalCircuit);
    return new QAOAResult({
      optimal_parameters: optimalParams,
      optimal_value: optResult.fun,
      optimal_circuit: optimalCircuit,
      optimal_state: optimalState,
      cost_function_evals: evalCount
    });
  }
};

// src/algorithms/grover.js
var GroverResult = class {
  constructor(kwargs = {}) {
    this.top_measurement = kwargs.top_measurement || null;
    this.iterations = kwargs.iterations || 0;
    this.measurement = kwargs.measurement || null;
    this.circuit = kwargs.circuit || null;
  }
};
var Grover = class {
  constructor(options = {}) {
    this.iterations = options.iterations || null;
    this.growth_rate = options.growth_rate || 2;
    this.sample_from_iterations = options.sample_from_iterations || false;
    this.sampler = options.sampler || null;
  }
  // Run Grover's algorithm.
  // oracle: QuantumCircuit that marks the solution(s) by flipping phase
  // numSolutions: number of marked items (for determining iterations)
  amplify(oracle, numSolutions = 1) {
    const numQubits = oracle.num_qubits;
    const N = 1 << numQubits;
    const M = numSolutions;
    let numIters;
    if (this.iterations !== null) {
      numIters = this.iterations;
    } else {
      numIters = Math.floor(Math.PI / 4 * Math.sqrt(N / M));
      if (numIters < 1) numIters = 1;
    }
    const circuit = new QuantumCircuit(numQubits, numQubits);
    for (let q = 0; q < numQubits; q++) circuit.h(q);
    circuit.barrier();
    for (let i = 0; i < numIters; i++) {
      circuit.compose(oracle, null, null, false, true);
      this._addDiffusion(circuit, numQubits);
      circuit.barrier();
    }
    for (let q = 0; q < numQubits; q++) circuit.measure(q, q);
    const result = simulate(circuit, 1024);
    const counts = result.get_counts().to_dict();
    let bestKey = null, bestCount = -1;
    for (const [key, count] of Object.entries(counts)) {
      if (count > bestCount) {
        bestCount = count;
        bestKey = key;
      }
    }
    return new GroverResult({
      top_measurement: bestKey,
      iterations: numIters,
      measurement: counts,
      circuit
    });
  }
  // Diffusion operator: 2|s><s| - I, where |s> = H^n |0...0>.
  // Decomposition: H^n · (2|0><0| - I) · H^n = H^n · X^n · (2|1...1><1...1| - I) · X^n · H^n
  // The middle operator is a multi-controlled Z (MCZ) up to a global phase.
  // For n=1: MCZ = Z. For n=2: MCZ = CZ. For n>=3: MCZ = H(target) · MCX(controls, target) · H(target).
  _addDiffusion(circuit, numQubits) {
    for (let q = 0; q < numQubits; q++) circuit.h(q);
    for (let q = 0; q < numQubits; q++) circuit.x(q);
    if (numQubits === 1) {
      circuit.z(0);
    } else if (numQubits === 2) {
      circuit.cz(0, 1);
    } else {
      circuit.h(numQubits - 1);
      const controls = [];
      for (let q = 0; q < numQubits - 1; q++) controls.push(q);
      circuit.mcx(controls, numQubits - 1);
      circuit.h(numQubits - 1);
    }
    for (let q = 0; q < numQubits; q++) circuit.x(q);
    for (let q = 0; q < numQubits; q++) circuit.h(q);
  }
};

// src/algorithms/phase_estimation.js
var PhaseEstimationResult = class {
  constructor(kwargs = {}) {
    this.phase = kwargs.phase || 0;
    this.measurement = kwargs.measurement || null;
    this.circuit = kwargs.circuit || null;
  }
};
var PhaseEstimation = class {
  constructor(options = {}) {
    this.num_evaluation_qubits = options.num_evaluation_qubits || 3;
    this.sampler = options.sampler || null;
  }
  // Estimate the phase of a unitary.
  // unitary: QuantumCircuit implementing U (any num_qubits >= 1).
  // state_preparation: optional circuit to prepare the eigenstate on U's qubits.
  estimate(unitary, statePreparation = null) {
    const evalQubits = this.num_evaluation_qubits;
    const unitaryQubits = unitary.num_qubits;
    const totalQubits = evalQubits + unitaryQubits;
    const circuit = new QuantumCircuit(totalQubits, evalQubits);
    if (statePreparation) {
      const targetQubits = [];
      for (let q = 0; q < unitaryQubits; q++) targetQubits.push(evalQubits + q);
      circuit.compose(statePreparation, targetQubits, null, false, true);
    }
    for (let q = 0; q < evalQubits; q++) {
      circuit.h(q);
    }
    circuit.barrier();
    for (let k = 0; k < evalQubits; k++) {
      const controlQubit = k;
      const unitaryQubitIndices = [];
      for (let q = 0; q < unitaryQubits; q++) {
        unitaryQubitIndices.push(evalQubits + q);
      }
      for (let p = 0; p < 1 << k; p++) {
        this._applyControlledUnitary(circuit, unitary, controlQubit, unitaryQubitIndices);
      }
    }
    circuit.barrier();
    this._addInverseQFT(circuit, evalQubits);
    for (let q = 0; q < evalQubits; q++) {
      circuit.measure(q, q);
    }
    const result = simulate(circuit, 1024);
    const counts = result.get_counts().to_dict();
    let bestKey = null, bestCount = -1;
    for (const [key, count] of Object.entries(counts)) {
      if (count > bestCount) {
        bestCount = count;
        bestKey = key;
      }
    }
    const measuredValue = parseInt(bestKey, 2);
    const phase = measuredValue / (1 << evalQubits);
    return new PhaseEstimationResult({
      phase,
      measurement: counts,
      circuit
    });
  }
  // Apply the unitary circuit U controlled by `controlQubit`, targeting
  // `unitaryQubits` in order. Each gate in U is wrapped in a ControlledGate
  // and applied to the corresponding target qubit.
  _applyControlledUnitary(circuit, unitary, controlQubit, unitaryQubits) {
    for (const ci of unitary.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (typeof op.to_matrix !== "function") continue;
      const mappedQubits = ci.qubits.map((q) => {
        const idxInUnitary = unitary._qubit_index.get(q);
        return unitaryQubits[idxInUnitary];
      });
      try {
        const cg = new ControlledGate(op, 1);
        cg.name = "c" + op.name;
        circuit.append(cg, [controlQubit].concat(mappedQubits));
      } catch (e) {
        circuit.append(op.copy(), mappedQubits);
      }
    }
  }
  // Add inverse QFT (QFT-dagger) on the first n qubits.
  _addInverseQFT(circuit, n) {
    for (let i = 0; i < Math.floor(n / 2); i++) {
      circuit.swap(i, n - 1 - i);
    }
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < i; j++) {
        const lam = -Math.PI / Math.pow(2, i - j);
        circuit.cp(lam, i, j);
      }
      circuit.h(i);
    }
  }
};
var AmplitudeEstimationResult = class {
  constructor(kwargs = {}) {
    this.estimation = kwargs.estimation || 0;
    this.num_oracle_queries = kwargs.num_oracle_queries || 0;
    this.measurement = kwargs.measurement || null;
    this.circuit = kwargs.circuit || null;
  }
};
var AmplitudeEstimation = class {
  constructor(options = {}) {
    this.num_eval_qubits = options.num_eval_qubits || 3;
    this.sampler = options.sampler || null;
  }
  // Estimate the amplitude of a good state.
  // statePreparation: circuit A that prepares the state with amplitude a.
  // groverOperator: circuit Q = A (2|0><0| - I) A^dagger (2|psi_g><psi_g| - I).
  // The amplitude a is estimated as sin^2(pi * y / 2^m) where y is the
  // measured value and m is the number of evaluation qubits.
  estimate(statePreparation, groverOperator = null) {
    const evalQubits = this.num_eval_qubits;
    const stateQubits = statePreparation.num_qubits;
    const totalQubits = evalQubits + stateQubits;
    const circuit = new QuantumCircuit(totalQubits, evalQubits);
    for (let q = 0; q < evalQubits; q++) circuit.h(q);
    const stateQubitIndices = [];
    for (let q = 0; q < stateQubits; q++) stateQubitIndices.push(evalQubits + q);
    this._applyCircuitToQubits(circuit, statePreparation, stateQubitIndices);
    if (groverOperator) {
      for (let k = 0; k < evalQubits; k++) {
        const controlQubit = k;
        for (let p = 0; p < 1 << k; p++) {
          this._applyControlledCircuit(circuit, groverOperator, controlQubit, stateQubitIndices);
        }
      }
    }
    const pe = new PhaseEstimation({ num_evaluation_qubits: evalQubits });
    pe._addInverseQFT(circuit, evalQubits);
    for (let q = 0; q < evalQubits; q++) circuit.measure(q, q);
    const result = simulate(circuit, 1024);
    const counts = result.get_counts().to_dict();
    let bestKey = null, bestCount = -1;
    for (const [key, count] of Object.entries(counts)) {
      if (count > bestCount) {
        bestCount = count;
        bestKey = key;
      }
    }
    const y = parseInt(bestKey, 2);
    const a = Math.sin(Math.PI * y / (1 << evalQubits)) ** 2;
    return new AmplitudeEstimationResult({
      estimation: a,
      num_oracle_queries: 2 * evalQubits,
      measurement: counts,
      circuit
    });
  }
  // Apply each gate of `sub` onto the specified target qubit indices.
  _applyCircuitToQubits(circuit, sub, targetQubits) {
    for (const ci of sub.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (typeof op.to_matrix !== "function") continue;
      const mappedQubits = ci.qubits.map((q) => {
        const idx = sub._qubit_index.get(q);
        return targetQubits[idx];
      });
      circuit.append(op.copy(), mappedQubits);
    }
  }
  // Apply each gate of `sub` controlled by `controlQubit`.
  _applyControlledCircuit(circuit, sub, controlQubit, targetQubits) {
    for (const ci of sub.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (typeof op.to_matrix !== "function") continue;
      const mappedQubits = ci.qubits.map((q) => {
        const idx = sub._qubit_index.get(q);
        return targetQubits[idx];
      });
      try {
        const cg = new ControlledGate(op, 1);
        cg.name = "c" + op.name;
        circuit.append(cg, [controlQubit].concat(mappedQubits));
      } catch (e) {
        circuit.append(op.copy(), mappedQubits);
      }
    }
  }
};

// src/algorithms/numpy_eigensolver.js
var NumPyMinimumEigensolverResult = class {
  constructor(kwargs = {}) {
    this.eigenvalue = kwargs.eigenvalue || 0;
    this.eigenstate = kwargs.eigenstate || null;
    this.aux_operators_evaluated = kwargs.aux_operators_evaluated || [];
  }
};
var NumPyMinimumEigensolver = class {
  constructor(options = {}) {
    this.aux_operators = options.aux_operators || null;
  }
  // Compute the minimum eigenvalue of `operator` (a SparsePauliOp, Operator,
  // or any object with a to_matrix() method).
  compute_minimum_eigenvalue(operator, auxOperators = null) {
    let mat;
    if (operator instanceof SparsePauliOp) {
      mat = operator.to_matrix();
    } else if (operator instanceof Operator) {
      mat = operator._data;
    } else if (typeof operator.to_matrix === "function") {
      mat = operator.to_matrix();
    } else if (operator instanceof ComplexMatrix) {
      mat = operator;
    } else {
      throw new TypeError("NumPyMinimumEigensolver: operator must have to_matrix()");
    }
    const { eigenvalues, eigenvectors } = mat.eigh();
    let minIdx = 0;
    for (let i = 1; i < eigenvalues.length; i++) {
      if (eigenvalues[i] < eigenvalues[minIdx]) minIdx = i;
    }
    const eigenvalue = eigenvalues[minIdx];
    const dim = mat.rows;
    const nq = Math.log2(dim);
    const eigenvecData = new Array(dim);
    for (let i = 0; i < dim; i++) {
      eigenvecData[i] = eigenvectors.get(i, minIdx);
    }
    const eigenstate = new Statevector(
      // Lazy import to avoid circular deps.
      (function() {
        const ComplexVector2 = globalThis.__ketraComplexVector;
        if (ComplexVector2) return new ComplexVector2(eigenvecData);
        throw new Error("ComplexVector not registered");
      })(),
      nq
    );
    const aux = auxOperators || this.aux_operators;
    const auxResults = [];
    if (aux) {
      for (const op of aux) {
        if (op instanceof SparsePauliOp) {
          auxResults.push([op, op.expectation_value(eigenstate)]);
        } else {
          auxResults.push([op, null]);
        }
      }
    }
    return new NumPyMinimumEigensolverResult({
      eigenvalue,
      eigenstate,
      aux_operators_evaluated: auxResults
    });
  }
};
var NumPyMaximumEigensolver = class extends NumPyMinimumEigensolver {
  compute_maximum_eigenvalue(operator, auxOperators = null) {
    const negOp = operator instanceof SparsePauliOp ? new SparsePauliOp(operator.paulis, operator.coeffs.map((c) => c.scale(-1))) : operator;
    const result = super.compute_minimum_eigenvalue(negOp, auxOperators);
    return new NumPyMinimumEigensolverResult({
      eigenvalue: -result.eigenvalue,
      eigenstate: result.eigenstate,
      aux_operators_evaluated: result.aux_operators_evaluated
    });
  }
};
globalThis.__ketraComplexVector = ComplexVector;

// src/algorithms/gradients.js
var GradientBase = class {
  constructor() {
  }
  // Compute gradient of <ψ(θ)|O|ψ(θ)> with respect to each parameter.
  // circuit: QuantumCircuit with parameters
  // observable: SparsePauliOp
  // parameterValues: { name: number }
  // Returns: array of numbers (gradient w.r.t. each parameter)
  compute(circuit, observable, parameterValues) {
    throw new Error("compute not implemented");
  }
  // Get the list of parameters from a circuit
  _getParamNames(circuit) {
    return Array.from(circuit.parameters).map((p) => p.name);
  }
};
var ParamShift = class extends GradientBase {
  constructor() {
    super();
  }
  compute(circuit, observable, parameterValues) {
    const paramNames = this._getParamNames(circuit);
    const gradient = new Array(paramNames.length);
    for (let i = 0; i < paramNames.length; i++) {
      const name = paramNames[i];
      const plus = Object.assign({}, parameterValues);
      plus[name] += Math.PI / 2;
      const fPlus = this._computeExpectation(circuit, observable, plus);
      const minus = Object.assign({}, parameterValues);
      minus[name] -= Math.PI / 2;
      const fMinus = this._computeExpectation(circuit, observable, minus);
      gradient[i] = (fPlus - fMinus) / 2;
    }
    return gradient;
  }
  _computeExpectation(circuit, observable, parameterValues) {
    const bound = circuit.bind_parameters(parameterValues);
    const sv = Statevector.fromCircuit(bound);
    const expVal = observable.expectation_value(sv);
    return typeof expVal === "number" ? expVal : expVal.re;
  }
};
var FiniteDiff = class extends GradientBase {
  constructor(epsilon = 1e-4) {
    super();
    this.epsilon = epsilon;
  }
  compute(circuit, observable, parameterValues) {
    const paramNames = this._getParamNames(circuit);
    const gradient = new Array(paramNames.length);
    for (let i = 0; i < paramNames.length; i++) {
      const name = paramNames[i];
      const plus = Object.assign({}, parameterValues);
      plus[name] += this.epsilon;
      const fPlus = this._computeExpectation(circuit, observable, plus);
      const minus = Object.assign({}, parameterValues);
      minus[name] -= this.epsilon;
      const fMinus = this._computeExpectation(circuit, observable, minus);
      gradient[i] = (fPlus - fMinus) / (2 * this.epsilon);
    }
    return gradient;
  }
  _computeExpectation(circuit, observable, parameterValues) {
    const bound = circuit.bind_parameters(parameterValues);
    const sv = Statevector.fromCircuit(bound);
    const expVal = observable.expectation_value(sv);
    return typeof expVal === "number" ? expVal : expVal.re;
  }
};
var LinearCombination = class extends GradientBase {
  constructor() {
    super();
  }
  compute(circuit, observable, parameterValues) {
    const paramShift = new ParamShift();
    return paramShift.compute(circuit, observable, parameterValues);
  }
};
var NaturalGradient = class extends GradientBase {
  constructor(regularizer = 1e-3) {
    super();
    this.regularizer = regularizer;
  }
  compute(circuit, observable, parameterValues) {
    const paramNames = this._getParamNames(circuit);
    const N = paramNames.length;
    const grad = new ParamShift().compute(circuit, observable, parameterValues);
    if (N === 0) return grad;
    const psiPlus = new Array(N);
    const psiMinus = new Array(N);
    for (let i = 0; i < N; i++) {
      const name = paramNames[i];
      const plus = Object.assign({}, parameterValues);
      plus[name] += Math.PI / 2;
      const minus = Object.assign({}, parameterValues);
      minus[name] -= Math.PI / 2;
      psiPlus[i] = Statevector.fromCircuit(circuit.bind_parameters(plus));
      psiMinus[i] = Statevector.fromCircuit(circuit.bind_parameters(minus));
    }
    const psi = Statevector.fromCircuit(circuit.bind_parameters(parameterValues));
    const dpsi = new Array(N);
    for (let i = 0; i < N; i++) {
      const dim = psi._data.size;
      const data = new Array(dim);
      for (let k = 0; k < dim; k++) {
        const a = psiPlus[i]._data.get(k);
        const b = psiMinus[i]._data.get(k);
        data[k] = new Complex((a.re - b.re) / 2, (a.im - b.im) / 2);
      }
      dpsi[i] = new Statevector(new ComplexVector(data), psi._numQubits);
    }
    const g = new Array(N);
    for (let i = 0; i < N; i++) {
      g[i] = new Array(N).fill(0);
      for (let j = 0; j < N; j++) {
        const inner_ij = dpsi[i]._data.inner(dpsi[j]._data);
        const inner_i_psi = dpsi[i]._data.inner(psi._data);
        const inner_psi_j = psi._data.inner(dpsi[j]._data);
        g[i][j] = inner_ij.re - inner_i_psi.mul(inner_psi_j).re;
      }
      g[i][i] += this.regularizer;
    }
    const x = _solveLinearSystem(g, grad);
    return x;
  }
};
function _solveLinearSystem(A, b) {
  const n = A.length;
  const M = A.map((row, i) => row.slice().concat([b[i]]));
  for (let col = 0; col < n; col++) {
    let maxRow = col;
    for (let r4 = col + 1; r4 < n; r4++) {
      if (Math.abs(M[r4][col]) > Math.abs(M[maxRow][col])) maxRow = r4;
    }
    [M[col], M[maxRow]] = [M[maxRow], M[col]];
    if (Math.abs(M[col][col]) < 1e-12) {
      return b.slice();
    }
    for (let r4 = 0; r4 < n; r4++) {
      if (r4 === col) continue;
      const factor = M[r4][col] / M[col][col];
      for (let c = col; c <= n; c++) {
        M[r4][c] -= factor * M[col][c];
      }
    }
  }
  const x = new Array(n);
  for (let i = 0; i < n; i++) x[i] = M[i][n] / M[i][i];
  return x;
}

// src/visualization/visualization.js
var CONTROLLED_2Q = ["cx", "cy", "cz", "ch", "csx", "crx", "cry", "crz", "cp", "cu1", "cu3", "cu"];
function draw_circuit(circuit, output = "text", kwargs) {
  if (output === "text" || output === "ascii") return _drawText(circuit);
  if (output === "latex") return _drawLatexSource(circuit);
  if (output === "html") return _drawHtml(circuit);
  return _drawText(circuit);
}
var _touchedCells = /* @__PURE__ */ new Set();
function _setCell(grid, row, col, char) {
  if (!grid[row]) return;
  while (grid[row].length <= col) grid[row].push(" ");
  grid[row][col] = char;
  _touchedCells.add(`${row},${col}`);
}
function _drawText(circuit) {
  const nq = circuit.num_qubits;
  const nc = circuit.num_clbits;
  _touchedCells.clear();
  const allInstructions = [];
  for (let idx = 0; idx < circuit.data.length; idx++) {
    const ci = circuit.data[idx];
    const qubits = ci.qubits.map((q) => circuit._qubit_index.get(q));
    const clbits = ci.clbits.map((c) => circuit._clbit_index.get(c));
    if (qubits.some((q) => q === void 0) || clbits.some((c) => c === void 0)) continue;
    allInstructions.push({ op: ci.operation, qubits, clbits, idx });
  }
  const dagLayers = [];
  const wireBusyUntil = /* @__PURE__ */ new Map();
  for (const instr of allInstructions) {
    const allWires = instr.qubits.concat(instr.clbits);
    let layer = 0;
    for (const w of allWires) {
      if (wireBusyUntil.has(w)) layer = Math.max(layer, wireBusyUntil.get(w) + 1);
    }
    while (dagLayers.length <= layer) dagLayers.push([]);
    dagLayers[layer].push(instr);
    for (const w of allWires) wireBusyUntil.set(w, layer);
  }
  function _gateSpan(instr) {
    if (instr.qubits.length === 0) return [];
    const minQ = Math.min(...instr.qubits);
    const maxQ = Math.max(...instr.qubits);
    const span = [];
    for (let i = minQ; i <= maxQ; i++) span.push(i);
    if (instr.clbits.length > 0 && instr.op.name === "measure") {
      for (let i = maxQ + 1; i < nq; i++) span.push(i);
    }
    return span;
  }
  function _anyCrossover(instr, layerNodes) {
    const span = new Set(_gateSpan(instr));
    for (const existing of layerNodes) {
      const existingSpan = _gateSpan(existing);
      for (const q of existingSpan) {
        if (span.has(q)) return true;
      }
    }
    return false;
  }
  const layers = [];
  for (const dagLayer of dagLayers) {
    const currentIdx = layers.length - 1;
    for (const node of dagLayer) {
      let inserted = false;
      let lastInsertable = -1;
      for (let i = currentIdx; i >= 0; i--) {
        const nodeWires = new Set(node.qubits.concat(node.clbits));
        let found = false;
        for (const existing of layers[i]) {
          for (const w of existing.qubits.concat(existing.clbits)) {
            if (nodeWires.has(w)) {
              found = true;
              break;
            }
          }
          if (found) break;
        }
        if (found) break;
        if (!_anyCrossover(node, layers[i])) {
          lastInsertable = i;
        }
      }
      if (lastInsertable >= 0) {
        layers[lastInsertable].push(node);
        inserted = true;
      } else {
        for (let i = Math.max(0, currentIdx); i < layers.length; i++) {
          if (!_anyCrossover(node, layers[i])) {
            const nodeWires = new Set(node.qubits.concat(node.clbits));
            let conflict = false;
            for (const existing of layers[i]) {
              for (const w of existing.qubits.concat(existing.clbits)) {
                if (nodeWires.has(w)) {
                  conflict = true;
                  break;
                }
              }
              if (conflict) break;
            }
            if (!conflict) {
              layers[i].push(node);
              inserted = true;
              break;
            }
          }
        }
      }
      if (!inserted) {
        layers.push([node]);
      }
    }
  }
  function _instrWidth(instr) {
    if (instr.op.name === "barrier") return 1;
    if (instr.op.name === "measure") return 3;
    if (instr.op.name === "reset") return 5;
    if (instr.op.name === "swap" || instr.op.name === "cswap" || instr.op.name === "cz") return 3;
    return Math.max(3, _gateLabel(instr.op).length) + 2;
  }
  const layerWidths = layers.map((layer) => {
    let maxW = 0;
    for (const instr of layer) {
      const w = _instrWidth(instr);
      if (w > maxW) maxW = w;
    }
    return maxW === 0 ? 5 : maxW;
  });
  const nRows = nc > 0 ? 2 * nq + 3 : 2 * nq + 1;
  const labelWidth = _maxLabelWidth(nq, nc);
  const grid = [];
  for (let r4 = 0; r4 < nRows; r4++) grid.push([]);
  _addLabels(grid, nq, nc, labelWidth);
  let colOffset = labelWidth;
  for (let li = 0; li < layers.length; li++) {
    const layer = layers[li];
    const width = layerWidths[li];
    for (let r4 = 0; r4 < nRows; r4++) {
      while (grid[r4].length < colOffset) grid[r4].push(" ");
    }
    _drawColumn(grid, layer, width, nq, nc, colOffset);
    const layerQubits = /* @__PURE__ */ new Set();
    for (const instr of layer) {
      for (const q of instr.qubits) layerQubits.add(q);
    }
    for (let qi = 0; qi < nq; qi++) {
      if (!layerQubits.has(qi)) {
        const wireRow = 2 * qi + 1;
        for (let c = 0; c < width; c++) {
          if (!_touchedCells.has(`${wireRow},${colOffset + c}`)) {
            _setCell(grid, wireRow, colOffset + c, "\u2500");
          }
        }
      }
    }
    const hasMeasure = layer.some((i) => i.op.name === "measure");
    if (nc > 0 && !hasMeasure) {
      const cWireRow = 2 * nq + 1;
      for (let c = 0; c < width; c++) {
        if (!_touchedCells.has(`${cWireRow},${colOffset + c}`)) {
          _setCell(grid, cWireRow, colOffset + c, "\u2550");
        }
      }
    }
    colOffset += width;
  }
  let maxLen = 0;
  for (let r4 = 0; r4 < nRows; r4++) {
    if (grid[r4].length > maxLen) maxLen = grid[r4].length;
  }
  const targetLen = maxLen;
  for (let r4 = 0; r4 < nRows; r4++) {
    while (grid[r4].length < targetLen) grid[r4].push(" ");
  }
  _fillUntouchedWires(grid, nq, nc, labelWidth, targetLen);
  for (let r4 = 0; r4 < nRows; r4++) {
    while (grid[r4].length < targetLen + 1) grid[r4].push(" ");
  }
  const lines = grid.map((row) => row.join(""));
  const merged = [];
  for (let r4 = 0; r4 < lines.length; r4++) {
    if (merged.length === 0) {
      merged.push(lines[r4]);
      continue;
    }
    const top = merged[merged.length - 1];
    const bot = lines[r4];
    let shouldCompress = true;
    for (let c = 0; c < Math.min(top.length, bot.length); c++) {
      const tc = top[c], bc = bot[c];
      if (tc === "\u2534" && (bc === "\u252C" || bc === "\u2565")) {
        shouldCompress = false;
        break;
      }
      if (tc === "\u2568" && (bc === "\u252C" || bc === "\u2565")) {
        shouldCompress = false;
        break;
      }
      if ((tc >= "a" && tc <= "z" || tc >= "A" && tc <= "Z" || tc >= "0" && tc <= "9") && bc !== " ") {
        shouldCompress = false;
        break;
      }
      if ((bc >= "a" && bc <= "z" || bc >= "A" && bc <= "Z" || bc >= "0" && bc <= "9") && tc !== " ") {
        shouldCompress = false;
        break;
      }
      if (tc === "\u2500" && bc === "\u2550") {
        shouldCompress = false;
        break;
      }
      if ((tc === "\u2518" || tc === "\u2514" || tc === "\u2524" || tc === "\u251C") && bc === "\u2550") {
        shouldCompress = false;
        break;
      }
    }
    if (shouldCompress) {
      merged[merged.length - 1] = _mergeLines(top, bot, "top");
    } else {
      merged.push(bot);
    }
  }
  return merged.join("\n");
}
function _mergeLines(top, bot, icod) {
  icod = icod || "top";
  let ret = "";
  const len = Math.max(top.length, bot.length);
  for (let c = 0; c < len; c++) {
    const tc = c < top.length ? top[c] : " ";
    const bc = c < bot.length ? bot[c] : " ";
    if (tc === bc) {
      ret += tc;
    } else if ((tc === "\u253C" || tc === "\u256A") && bc === " ") {
      ret += "\u2502";
    } else if (tc === " ") {
      ret += bc;
    } else if ((tc === "\u252C" || tc === "\u2565") && (bc === " " || bc === "\u2551" || bc === "\u2502") && icod === "top") {
      ret += tc;
    } else if (tc === "\u252C" && bc === " " && icod === "bot") {
      ret += "\u2502";
    } else if (tc === "\u2565" && bc === " " && icod === "bot") {
      ret += "\u2551";
    } else if ((tc === "\u252C" || tc === "\u2502") && bc === "\u2550") {
      ret += "\u256A";
    } else if ((tc === "\u252C" || tc === "\u2502") && bc === "\u2500") {
      ret += "\u253C";
    } else if ((tc === "\u2514" || tc === "\u2518" || tc === "\u2551" || tc === "\u2502" || tc === "\u2591") && bc === " " && icod === "top") {
      ret += tc;
    } else if ((tc === "\u2500" || tc === "\u2550") && bc === " " && icod === "top") {
      ret += tc;
    } else if ((tc === "\u2500" || tc === "\u2550") && bc === " " && icod === "bot") {
      ret += bc;
    } else if ((tc === "\u2551" || tc === "\u2565") && bc === "\u2550") {
      ret += "\u256C";
    } else if ((tc === "\u2551" || tc === "\u2565") && bc === "\u2500") {
      ret += "\u256B";
    } else if ((tc === "\u2551" || tc === "\u256B" || tc === "\u256C") && bc === " ") {
      ret += "\u2551";
    } else if ((tc === "\u2502" || tc === "\u253C" || tc === "\u256A") && bc === " ") {
      ret += "\u2502";
    } else if (tc === "\u2514" && bc === "\u250C" && icod === "top") {
      ret += "\u251C";
    } else if (tc === "\u2518" && bc === "\u2510" && icod === "top") {
      ret += "\u2524";
    } else if ((bc === "\u250C" || bc === "\u2510") && icod === "top") {
      ret += "\u252C";
    } else if ((tc === "\u2518" || tc === "\u2514") && bc === "\u2500" && icod === "top") {
      ret += "\u2534";
    } else if (bc === " " && icod === "top") {
      ret += tc;
    } else {
      ret += bc;
    }
  }
  return ret;
}
function _gateLabel(op) {
  const name = op.name;
  if (name === "measure") return "M";
  if (name === "reset") return "|0>";
  if (name === "barrier") return "";
  if (name === "delay") return `D${op.params[0] || 0}`;
  const mixedCase = {
    "rx": "Rx",
    "ry": "Ry",
    "rz": "Rz",
    "rxx": "Rxx",
    "ryy": "Ryy",
    "rzz": "Rzz",
    "rzx": "Rzx",
    "crx": "CRx",
    "cry": "CRy",
    "crz": "CRz",
    "p": "P",
    "cp": "P",
    "u1": "U1",
    "u2": "U2",
    "u3": "U3",
    "u": "U",
    "cu1": "CU1",
    "cu3": "CU3",
    "cu": "CU",
    "iswap": "Iswap",
    "dcx": "Dcx",
    "sdg": "Sdg",
    "tdg": "Tdg",
    "sxdg": "Sxdg"
  };
  const params = op.params.length ? `(${op.params.map(_fmtParam).join(",")})` : "";
  const baseName = mixedCase[name] || name.toUpperCase();
  return `${baseName}${params}`;
}
function _controlledTargetLabel(name, op) {
  const paramMap = {
    "cx": "X",
    "cy": "Y",
    "cz": null,
    "ch": "H",
    "csx": "SX",
    "crx": "Rx",
    "cry": "Ry",
    "crz": "Rz",
    "cp": "P",
    "cu1": "U1",
    "cu3": "U3",
    "cu": "U"
  };
  const base = paramMap[name];
  if (base === null) return null;
  if (!base) return name.toUpperCase();
  if (op.params.length > 0 && ["crx", "cry", "crz", "cp", "cu1", "cu3", "cu"].includes(name)) {
    return `${base}(${op.params.map(_fmtParam).join(",")})`;
  }
  return base;
}
function _fmtParam(p) {
  if (typeof p === "number") {
    if (Number.isInteger(p)) return p.toString();
    const rounded = Math.round(p * 1e4) / 1e4;
    let str = rounded.toFixed(4);
    str = str.replace(/\.?0+$/, "");
    return str;
  }
  if (p && typeof p.toString === "function") return p.toString();
  return String(p);
}
function _maxLabelWidth(nq, nc) {
  let maxW = 0;
  for (let i = 0; i < nq; i++) {
    maxW = Math.max(maxW, `q_${i}: `.length);
  }
  if (nc > 0) {
    maxW = Math.max(maxW, `c: ${nc}/`.length);
  }
  return maxW;
}
function _addLabels(grid, nq, nc, labelWidth) {
  for (let i = 0; i < nq; i++) {
    const wireRow = 2 * i + 1;
    const label = `q_${i}: `;
    _setRow(grid, wireRow, label, labelWidth);
    if (wireRow > 0) _setRow(grid, wireRow - 1, "", labelWidth);
    if (wireRow + 1 < grid.length) _setRow(grid, wireRow + 1, "", labelWidth);
  }
  if (nc > 0) {
    const cWireRow = 2 * nq + 1;
    const label = `c: ${nc}/`;
    _setRow(grid, cWireRow, label, labelWidth);
    _setRow(grid, cWireRow - 1, "", labelWidth);
    if (cWireRow + 1 < grid.length) _setRow(grid, cWireRow + 1, "", labelWidth);
  }
}
function _setRow(grid, row, text, width) {
  while (grid[row].length < width) {
    grid[row].push(text[grid[row].length] || " ");
  }
}
function _drawColumn(grid, col, width, nq, nc, colOffset) {
  const measurements = col.filter((i) => i.op.name === "measure");
  const otherInstrs = col.filter((i) => i.op.name !== "measure" && i.op.name !== "barrier");
  const boxQubits = /* @__PURE__ */ new Set();
  const boxLabels = /* @__PURE__ */ new Map();
  for (const instr of otherInstrs) {
    if (instr.op.num_qubits === 1 && instr.op.name !== "reset") {
      boxQubits.add(instr.qubits[0]);
      boxLabels.set(instr.qubits[0], _gateLabel(instr.op));
    } else if (instr.op.name === "reset") {
      boxQubits.add(instr.qubits[0]);
      boxLabels.set(instr.qubits[0], "|0>");
    }
  }
  const controlLines = [];
  const targetBoxes = /* @__PURE__ */ new Set();
  for (const instr of otherInstrs) {
    if (CONTROLLED_2Q.includes(instr.op.name) || instr.op.name === "ccx") {
      const centerCol = colOffset + Math.floor(width / 2);
      if (instr.op.name === "ccx") {
        controlLines.push({ from: Math.min(...instr.qubits), to: Math.max(...instr.qubits), centerCol });
        targetBoxes.add(instr.qubits[2]);
      } else {
        controlLines.push({ from: Math.min(instr.qubits[0], instr.qubits[1]), to: Math.max(instr.qubits[0], instr.qubits[1]), centerCol });
        if (instr.op.name !== "cz") targetBoxes.add(instr.qubits[1]);
      }
    }
  }
  _drawBoxesWithSharedBorders(grid, boxQubits, boxLabels, width, colOffset, nq);
  for (const instr of otherInstrs) {
    if (instr.op.num_qubits >= 2) {
      _drawInstruction(grid, instr, width, nq, nc, colOffset);
    }
  }
  if (measurements.length > 0) {
    _drawMeasurementLayer(grid, measurements, width, nq, nc, colOffset);
  }
  for (const line of controlLines) {
    for (let q = line.from + 1; q < line.to; q++) {
      const wireRow = 2 * q + 1;
      if (boxQubits.has(q) && !targetBoxes.has(q)) {
        _setCell(grid, wireRow, line.centerCol, "\u253C");
      }
    }
  }
}
function _drawBoxesWithSharedBorders(grid, boxQubits, boxLabels, width, colOffset, nq) {
  if (boxQubits.size === 0) return;
  const innerW = Math.max(3, width - 2);
  const sortedQubits = Array.from(boxQubits).sort((a, b) => a - b);
  for (const q of sortedQubits) {
    const wireRow = 2 * q + 1;
    const topRow = wireRow - 1;
    const botRow = wireRow + 1;
    const label = boxLabels.get(q);
    const padL = Math.floor((innerW - label.length) / 2);
    const padR = innerW - label.length - padL;
    const hasBoxAbove = boxQubits.has(q - 1);
    const hasBoxBelow = boxQubits.has(q + 1);
    if (hasBoxAbove) {
      _setCell(grid, topRow, colOffset, "\u251C");
    } else {
      _setCell(grid, topRow, colOffset, "\u250C");
    }
    for (let c = 0; c < innerW; c++) _setCell(grid, topRow, colOffset + 1 + c, "\u2500");
    if (hasBoxAbove) {
      _setCell(grid, topRow, colOffset + 1 + innerW, "\u2524");
    } else {
      _setCell(grid, topRow, colOffset + 1 + innerW, "\u2510");
    }
    _setCell(grid, wireRow, colOffset, "\u2524");
    for (let c = 0; c < padL; c++) _setCell(grid, wireRow, colOffset + 1 + c, " ");
    for (let c = 0; c < label.length; c++) _setCell(grid, wireRow, colOffset + 1 + padL + c, label[c]);
    for (let c = 0; c < padR; c++) _setCell(grid, wireRow, colOffset + 1 + padL + label.length + c, " ");
    _setCell(grid, wireRow, colOffset + 1 + innerW, "\u251C");
    if (hasBoxBelow) {
      _setCell(grid, botRow, colOffset, "\u251C");
    } else {
      _setCell(grid, botRow, colOffset, "\u2514");
    }
    for (let c = 0; c < innerW; c++) _setCell(grid, botRow, colOffset + 1 + c, "\u2500");
    if (hasBoxBelow) {
      _setCell(grid, botRow, colOffset + 1 + innerW, "\u2524");
    } else {
      _setCell(grid, botRow, colOffset + 1 + innerW, "\u2518");
    }
  }
}
function _drawMeasurementLayer(grid, measurements, width, nq, nc, colOffset) {
  const cWireRow = 2 * nq + 1;
  const cLabelRow = 2 * nq + 2;
  const mWidth = 3;
  const totalMeasureWidth = measurements.length * mWidth;
  let mStartOffset = colOffset + Math.floor((width - totalMeasureWidth) / 2);
  if (mStartOffset < colOffset) mStartOffset = colOffset;
  const measurePositions = [];
  for (let mi = 0; mi < measurements.length; mi++) {
    const m = measurements[mi];
    const mOffset = mStartOffset + mi * mWidth;
    const centerCol = mOffset + 1;
    measurePositions.push({
      qubit: m.qubits[0],
      clbit: m.clbits[0],
      mOffset,
      centerCol,
      qWireRow: 2 * m.qubits[0] + 1
    });
  }
  const sorted = measurePositions.slice().sort((a, b) => a.qubit - b.qubit);
  for (const mp of sorted) {
    for (let r4 = mp.qWireRow + 2; r4 < cWireRow; r4++) {
      const isQuantumWire = r4 % 2 === 1 && r4 < 2 * nq + 1;
      if (isQuantumWire) {
        _setCell(grid, r4, mp.centerCol, "\u256B");
      } else {
        _setCell(grid, r4, mp.centerCol, "\u2551");
      }
    }
  }
  for (const mp of sorted) {
    const { qubit, clbit, mOffset, centerCol, qWireRow } = mp;
    let lineFromAbove = false;
    for (const other of sorted) {
      if (other.qubit < qubit && other.centerCol === centerCol) {
        lineFromAbove = true;
        break;
      }
    }
    const topRow = qWireRow - 1;
    const existingLeft = grid[topRow] && grid[topRow][mOffset] || " ";
    if (existingLeft === "\u2500") {
      _setCell(grid, topRow, mOffset, "\u252C");
    } else if (existingLeft === "\u2514") {
      _setCell(grid, topRow, mOffset, "\u251C");
    } else if (existingLeft === "\u2518") {
      _setCell(grid, topRow, mOffset, "\u2524");
    } else {
      _setCell(grid, topRow, mOffset, "\u250C");
    }
    if (lineFromAbove) {
      _setCell(grid, topRow, centerCol, "\u252C");
    } else {
      _setCell(grid, topRow, centerCol, "\u2500");
    }
    const existingRight = grid[topRow] && grid[topRow][mOffset + 2] || " ";
    if (existingRight === "\u2500") {
      _setCell(grid, topRow, mOffset + 2, "\u252C");
    } else if (existingRight === "\u2518") {
      _setCell(grid, topRow, mOffset + 2, "\u2524");
    } else if (existingRight === "\u2514") {
      _setCell(grid, topRow, mOffset + 2, "\u251C");
    } else {
      _setCell(grid, topRow, mOffset + 2, "\u2510");
    }
    _setCell(grid, qWireRow, mOffset, "\u2524");
    _setCell(grid, qWireRow, centerCol, "M");
    _setCell(grid, qWireRow, mOffset + 2, "\u251C");
    _setCell(grid, qWireRow + 1, mOffset, "\u2514");
    _setCell(grid, qWireRow + 1, centerCol, "\u2565");
    _setCell(grid, qWireRow + 1, mOffset + 2, "\u2518");
    _setCell(grid, cWireRow, centerCol, "\u2569");
    _setCell(grid, cLabelRow, centerCol, clbit.toString());
  }
  for (let c = 0; c < width; c++) {
    if (!_touchedCells.has(`${cWireRow},${colOffset + c}`)) {
      _setCell(grid, cWireRow, colOffset + c, "\u2550");
    }
  }
}
function _drawInstruction(grid, instr, width, nq, nc, colOffset) {
  const op = instr.op;
  const name = op.name;
  if (name === "barrier") return;
  if (name === "measure") return;
  if (name === "reset") {
    _drawBox(grid, instr.qubits[0], "|0>", width, colOffset);
    return;
  }
  if (op.num_qubits === 1) {
    _drawBox(grid, instr.qubits[0], _gateLabel(op), width, colOffset);
    return;
  }
  if (op.num_qubits === 2) {
    if (CONTROLLED_2Q.includes(name)) {
      const targetLabel = _controlledTargetLabel(name, op);
      _drawControlled(grid, instr.qubits[0], instr.qubits[1], targetLabel, width, nq, colOffset, name);
    } else {
      _drawSymmetric(grid, instr.qubits[0], instr.qubits[1], name, width, nq, colOffset);
    }
    return;
  }
  if (op.num_qubits === 3) {
    if (name === "ccx") {
      _drawControlled(grid, instr.qubits[0], instr.qubits[2], "X", width, nq, colOffset, "x", instr.qubits[1]);
    } else if (name === "cswap") {
      _drawCSwap(grid, instr.qubits[0], instr.qubits[1], instr.qubits[2], width, nq, colOffset);
    } else {
      _drawControlled(grid, instr.qubits[0], instr.qubits[2], _gateLabel(op), width, nq, colOffset, name, instr.qubits[1]);
    }
    return;
  }
}
function _drawBox(grid, qubit, label, width, colOffset) {
  const wireRow = 2 * qubit + 1;
  const topRow = wireRow - 1;
  const botRow = wireRow + 1;
  const innerW = Math.max(3, width - 2);
  const padL = Math.floor((innerW - label.length) / 2);
  const padR = innerW - label.length - padL;
  _setCell(grid, topRow, colOffset, "\u250C");
  for (let c = 0; c < innerW; c++) _setCell(grid, topRow, colOffset + 1 + c, "\u2500");
  _setCell(grid, topRow, colOffset + 1 + innerW, "\u2510");
  _setCell(grid, wireRow, colOffset, "\u2524");
  for (let c = 0; c < padL; c++) _setCell(grid, wireRow, colOffset + 1 + c, " ");
  for (let c = 0; c < label.length; c++) _setCell(grid, wireRow, colOffset + 1 + padL + c, label[c]);
  for (let c = 0; c < padR; c++) _setCell(grid, wireRow, colOffset + 1 + padL + label.length + c, " ");
  _setCell(grid, wireRow, colOffset + 1 + innerW, "\u251C");
  _setCell(grid, botRow, colOffset, "\u2514");
  for (let c = 0; c < innerW; c++) _setCell(grid, botRow, colOffset + 1 + c, "\u2500");
  _setCell(grid, botRow, colOffset + 1 + innerW, "\u2518");
}
function _drawTargetBox(grid, qubit, label, width, colOffset, lineFromAbove, centerCol) {
  const wireRow = 2 * qubit + 1;
  const topRow = wireRow - 1;
  const botRow = wireRow + 1;
  const innerW = Math.max(3, width - 2);
  const padL = Math.floor((innerW - label.length) / 2);
  const padR = innerW - label.length - padL;
  _setCell(grid, topRow, colOffset, "\u250C");
  for (let c = 0; c < innerW; c++) {
    if (lineFromAbove && colOffset + 1 + c === centerCol) {
      _setCell(grid, topRow, colOffset + 1 + c, "\u2534");
    } else {
      _setCell(grid, topRow, colOffset + 1 + c, "\u2500");
    }
  }
  _setCell(grid, topRow, colOffset + 1 + innerW, "\u2510");
  _setCell(grid, wireRow, colOffset, "\u2524");
  for (let c = 0; c < padL; c++) _setCell(grid, wireRow, colOffset + 1 + c, " ");
  for (let c = 0; c < label.length; c++) _setCell(grid, wireRow, colOffset + 1 + padL + c, label[c]);
  for (let c = 0; c < padR; c++) _setCell(grid, wireRow, colOffset + 1 + padL + label.length + c, " ");
  _setCell(grid, wireRow, colOffset + 1 + innerW, "\u251C");
  _setCell(grid, botRow, colOffset, "\u2514");
  for (let c = 0; c < innerW; c++) {
    if (!lineFromAbove && colOffset + 1 + c === centerCol) {
      _setCell(grid, botRow, colOffset + 1 + c, "\u252C");
    } else {
      _setCell(grid, botRow, colOffset + 1 + c, "\u2500");
    }
  }
  _setCell(grid, botRow, colOffset + 1 + innerW, "\u2518");
}
function _drawControlled(grid, control, target, label, width, nq, colOffset, gateName, extraControl = null) {
  const ctrlRow = 2 * control + 1;
  const centerCol = colOffset + Math.floor(width / 2);
  _setCell(grid, ctrlRow, centerCol, "\u25A0");
  if (extraControl !== null) {
    _setCell(grid, 2 * extraControl + 1, centerCol, "\u25A0");
  }
  const allNodes = extraControl !== null ? [control, extraControl, target] : [control, target];
  allNodes.sort((a, b) => a - b);
  for (let i = 0; i < allNodes.length - 1; i++) {
    for (let r4 = 2 * allNodes[i] + 2; r4 <= 2 * allNodes[i + 1]; r4++) {
      if (r4 % 2 === 1) {
        const qIdx = (r4 - 1) / 2;
        if (qIdx !== target) {
          _setCell(grid, r4, centerCol, "\u253C");
        }
      } else {
        _setCell(grid, r4, centerCol, "\u2502");
      }
    }
  }
  if (gateName === "cz") {
    _setCell(grid, 2 * target + 1, centerCol, "\u25A0");
  } else {
    _drawTargetBox(grid, target, label, width, colOffset, control < target, centerCol);
  }
}
function _drawSymmetric(grid, q0, q1, gateName, width, nq, colOffset) {
  const centerCol = colOffset + Math.floor(width / 2);
  if (gateName === "swap") {
    _setCell(grid, 2 * q0 + 1, centerCol, "X");
    _setCell(grid, 2 * q1 + 1, centerCol, "X");
    const lo = Math.min(q0, q1);
    const hi = Math.max(q0, q1);
    for (let r4 = 2 * lo + 2; r4 <= 2 * hi; r4++) {
      if (r4 % 2 === 1) {
        _setCell(grid, r4, centerCol, "\u253C");
      } else {
        _setCell(grid, r4, centerCol, "\u2502");
      }
    }
  } else {
    const label = _symmetricLabel(gateName);
    _drawBox(grid, q0, label, width, colOffset);
    _drawBox(grid, q1, label, width, colOffset);
  }
}
function _symmetricLabel(name) {
  const map = {
    "iswap": "Iswap",
    "dcx": "Dcx",
    "rxx": "Rxx",
    "ryy": "Ryy",
    "rzz": "ZZ",
    // qiskit shows RZZ as ZZ
    "rzx": "Rzx"
  };
  return map[name] || name.toUpperCase();
}
function _drawCSwap(grid, control, t1, t2, width, nq, colOffset) {
  const centerCol = colOffset + Math.floor(width / 2);
  _setCell(grid, 2 * control + 1, centerCol, "\u25A0");
  _setCell(grid, 2 * t1 + 1, centerCol, "X");
  _setCell(grid, 2 * t2 + 1, centerCol, "X");
  const allQs = [control, t1, t2].sort((a, b) => a - b);
  for (let i = 0; i < allQs.length - 1; i++) {
    for (let r4 = 2 * allQs[i] + 2; r4 <= 2 * allQs[i + 1]; r4++) {
      if (r4 % 2 === 1) {
        _setCell(grid, r4, centerCol, "\u253C");
      } else {
        _setCell(grid, r4, centerCol, "\u2502");
      }
    }
  }
}
function _fillUntouchedWires(grid, nq, nc, startCol, endCol) {
  for (let i = 0; i < nq; i++) {
    const wireRow = 2 * i + 1;
    for (let c = startCol; c < endCol; c++) {
      if (!_touchedCells.has(`${wireRow},${c}`)) {
        grid[wireRow][c] = "\u2500";
      }
    }
  }
  if (nc > 0) {
    const cWireRow = 2 * nq + 1;
    for (let c = startCol; c < endCol; c++) {
      if (!_touchedCells.has(`${cWireRow},${c}`)) {
        grid[cWireRow][c] = "\u2550";
      }
    }
  }
}
function _drawLatexSource(circuit) {
  const nq = circuit.num_qubits;
  const nc = circuit.num_clbits;
  const lines = [];
  lines.push("\\documentclass[border=2pt]{standalone}");
  lines.push("\\usepackage[braket, qm]{qcircuit}");
  lines.push("\\begin{document}");
  lines.push("\\Qcircuit @C=1.0em @R=0.7em {");
  const header = [];
  for (let i = 0; i < nq; i++) header.push(`\\lstick{\\ket{q_${i}}}`);
  for (let i = 0; i < nc; i++) header.push(`\\lstick{c_{i}}`);
  header.push("\\qw");
  lines.push("  " + header.join(" & ") + " \\\\");
  for (const ci of circuit.data) {
    const row = new Array(nq + nc).fill("\\qw");
    const name = ci.operation.name;
    if (name === "barrier") {
      for (const q of ci.qubits) row[circuit._qubit_index.get(q)] = "\\barrier{0}";
    } else if (name === "measure") {
      const qi = circuit._qubit_index.get(ci.qubits[0]);
      const ci2 = circuit._clbit_index.get(ci.clbits[0]);
      row[qi] = "\\meter";
      row[nq + ci2] = "\\cw";
    } else if (ci.operation.num_qubits === 1) {
      const qi = circuit._qubit_index.get(ci.qubits[0]);
      const p = ci.operation.params.length ? `(${ci.operation.params.map(_fmtParam).join(",")})` : "";
      row[qi] = `\\gate{${name.toUpperCase()}${p}}`;
    } else if (ci.operation.num_qubits === 2) {
      const q0 = circuit._qubit_index.get(ci.qubits[0]);
      const q1 = circuit._qubit_index.get(ci.qubits[1]);
      if (CONTROLLED_2Q.includes(name)) {
        row[q0] = "\\ctrl";
        row[q1] = `\\gate{${name.toUpperCase()}}`;
      } else {
        row[q0] = `\\gate{${name.toUpperCase()}}`;
        row[q1] = `\\gate{${name.toUpperCase()}}`;
      }
    } else if (ci.operation.num_qubits === 3) {
      const q0 = circuit._qubit_index.get(ci.qubits[0]);
      const q1 = circuit._qubit_index.get(ci.qubits[1]);
      const q2 = circuit._qubit_index.get(ci.qubits[2]);
      row[q0] = "\\ctrl";
      row[q1] = "\\ctrl";
      row[q2] = `\\gate{${name.toUpperCase()}}`;
    }
    lines.push("  " + row.join(" & ") + " \\\\");
  }
  lines.push("}");
  lines.push("\\end{document}");
  return lines.join("\n");
}
function _drawHtml(circuit) {
  const rows = [];
  const nq = circuit.num_qubits;
  const nc = circuit.num_clbits;
  const cols = [];
  for (const ci of circuit.data) {
    const col = new Array(nq + nc).fill("");
    const name = ci.operation.name;
    if (name === "barrier") {
      for (const q of ci.qubits) col[circuit._qubit_index.get(q)] = "\u2502";
    } else if (name === "measure") {
      col[circuit._qubit_index.get(ci.qubits[0])] = "M";
      col[nq + circuit._clbit_index.get(ci.clbits[0])] = "\u2190";
    } else if (ci.operation.num_qubits === 1) {
      col[circuit._qubit_index.get(ci.qubits[0])] = name.toUpperCase();
    } else if (ci.operation.num_qubits === 2) {
      if (CONTROLLED_2Q.includes(name)) {
        col[circuit._qubit_index.get(ci.qubits[0])] = "\u25A0";
        col[circuit._qubit_index.get(ci.qubits[1])] = name.toUpperCase();
      } else {
        col[circuit._qubit_index.get(ci.qubits[0])] = name.toUpperCase();
        col[circuit._qubit_index.get(ci.qubits[1])] = name.toUpperCase();
      }
    } else if (ci.operation.num_qubits === 3) {
      col[circuit._qubit_index.get(ci.qubits[0])] = "\u25A0";
      col[circuit._qubit_index.get(ci.qubits[1])] = "\u25A0";
      col[circuit._qubit_index.get(ci.qubits[2])] = name.toUpperCase();
    }
    cols.push(col);
  }
  rows.push('<table class="qiskit-circuit" style="border-collapse:collapse;font-family:monospace">');
  for (let r4 = 0; r4 < nq + nc; r4++) {
    rows.push("<tr>");
    const label = r4 < nq ? `q_${r4}` : `c_${r4 - nq}`;
    rows.push(`<td style="padding:4px 8px;border:1px solid #ccc;font-weight:bold">${label}</td>`);
    for (const col of cols) {
      const cell = col[r4] || "";
      rows.push(`<td style="padding:4px 8px;border:1px solid #eee;text-align:center">${cell}</td>`);
    }
    rows.push("</tr>");
  }
  rows.push("</table>");
  return rows.join("");
}
function plot_histogram(data, kwargs = {}) {
  const counts = data && data.to_dict ? data.to_dict() : data;
  const total = Object.values(counts).reduce((a, b) => a + b, 0);
  const normalized = {};
  for (const k in counts) normalized[k] = counts[k] / total;
  return _makePlaceholder("histogram", { counts, normalized, kwargs });
}
function plot_state_city(statevector, kwargs = {}) {
  const dim = statevector.num_qubits ? 1 << statevector.num_qubits : statevector.size;
  const data = statevector.data || statevector;
  const realPart = new Array(dim);
  const imagPart = new Array(dim);
  for (let i = 0; i < dim; i++) {
    realPart[i] = new Array(dim);
    imagPart[i] = new Array(dim);
    for (let j = 0; j < dim; j++) {
      const a = data.get ? data.get(i) : data[i];
      const b = data.get ? data.get(j) : data[j];
      const conjB = b.conjugate ? b.conjugate() : { re: b.re, im: -b.im };
      const prod = a.re !== void 0 ? { re: a.re * conjB.re - a.im * conjB.im, im: a.re * conjB.im + a.im * conjB.re } : a * b;
      realPart[i][j] = prod.re;
      imagPart[i][j] = prod.im;
    }
  }
  return _makePlaceholder("state_city", { realPart, imagPart, numQubits: statevector.num_qubits, kwargs });
}
function plot_bloch_vector(blochVector, kwargs = {}) {
  return _makePlaceholder("bloch", { vector: blochVector, kwargs });
}
function plot_state_hinton(statevector, kwargs = {}) {
  return plot_state_city(statevector, kwargs);
}
function plot_state_qsphere(statevector, kwargs = {}) {
  return _makePlaceholder("qsphere", { statevector, kwargs });
}
function _makePlaceholder(kind, data) {
  if (typeof document !== "undefined") {
    const div = document.createElement("div");
    div.textContent = "Placeholder";
    return div;
  }
  return { kind, data };
}
_setDrawHook(draw_circuit);

// src/visualization/ascii_viz.js
function render_histogram_ascii(counts, options = {}) {
  const maxBarWidth = options.maxBarWidth || 40;
  const data = counts && counts.to_dict ? counts.to_dict() : counts;
  const entries = Object.entries(data).sort((a, b) => b[1] - a[1]);
  const maxCount = Math.max(...entries.map((e) => e[1]));
  const total = entries.reduce((s, e) => s + e[1], 0);
  const barChar = options.barChar || "\u2588";
  const emptyChar = options.emptyChar || " ";
  const labelWidth = Math.max(...entries.map((e) => e[0].length));
  const countWidth = String(maxCount).length;
  const lines = [];
  for (const [key, count] of entries) {
    const barLen = Math.round(count / maxCount * maxBarWidth);
    const bar = barChar.repeat(barLen);
    const pct = (count / total * 100).toFixed(1);
    const paddedKey = key.padEnd(labelWidth);
    const paddedCount = String(count).padStart(countWidth);
    lines.push(`${paddedKey} \u2502${bar} ${paddedCount} (${pct}%)`);
  }
  const axisLine = " ".repeat(labelWidth) + "\u2514" + "\u2500".repeat(maxBarWidth + countWidth + 8);
  lines.push(axisLine);
  return lines.join("\n");
}
function render_bloch_ascii(vector, options = {}) {
  const [x, y, z] = vector;
  const width = options.width || 24;
  const height = options.height || 12;
  const radius = Math.min(width, height * 2) / 2 - 2;
  const px = Math.round(radius * x + width / 2);
  const py = Math.round(-radius * z + height / 2);
  const grid = [];
  for (let r4 = 0; r4 < height; r4++) {
    grid.push(new Array(width).fill(" "));
  }
  for (let theta = 0; theta < 2 * Math.PI; theta += 0.05) {
    const cx = Math.round(radius * Math.cos(theta) + width / 2);
    const cy = Math.round(radius * Math.sin(theta) * 0.5 + height / 2);
    if (cy >= 0 && cy < height && cx >= 0 && cx < width) {
      grid[cy][cx] = "\xB7";
    }
  }
  for (let i = -radius; i <= radius; i++) {
    const cx = Math.round(i + width / 2);
    const cy = Math.round(height / 2);
    if (cy >= 0 && cy < height && cx >= 0 && cx < width && grid[cy][cx] === " ") {
      grid[cy][cx] = "\u2500";
    }
  }
  for (let i = -radius / 2; i <= radius / 2; i++) {
    const cx = Math.round(width / 2);
    const cy = Math.round(i + height / 2);
    if (cy >= 0 && cy < height && cx >= 0 && cx < width && grid[cy][cx] === " ") {
      grid[cy][cx] = "\u2502";
    }
  }
  const cx0 = Math.round(width / 2);
  const cy0 = Math.round(height / 2);
  const steps = Math.max(Math.abs(px - cx0), Math.abs(py - cy0), 1);
  for (let s = 0; s <= steps; s++) {
    const ix = Math.round(cx0 + (px - cx0) * s / steps);
    const iy = Math.round(cy0 + (py - cy0) * s / steps);
    if (iy >= 0 && iy < height && ix >= 0 && ix < width) {
      grid[iy][ix] = "\u2022";
    }
  }
  if (py >= 0 && py < height && px >= 0 && px < width) {
    grid[py][px] = "\u25C9";
  }
  _setGrid(grid, 0, width - 1, "x");
  _setGrid(grid, 0, 0, "-x");
  _setGrid(grid, 0, Math.floor(width / 2), "z");
  return grid.map((row) => row.join("")).join("\n");
}
function render_state_city_ascii(statevector, options = {}) {
  const dim = statevector.num_qubits ? 1 << statevector.num_qubits : statevector.size;
  const data = statevector.data || statevector;
  const realPart = new Array(dim);
  const imagPart = new Array(dim);
  for (let i = 0; i < dim; i++) {
    realPart[i] = new Array(dim);
    imagPart[i] = new Array(dim);
    for (let j = 0; j < dim; j++) {
      const a = data.get ? data.get(i) : data[i];
      const b = data.get ? data.get(j) : data[j];
      const conjB = b.conjugate ? b.conjugate() : { re: b.re, im: -b.im };
      const prod = a.re !== void 0 ? { re: a.re * conjB.re - a.im * conjB.im, im: a.re * conjB.im + a.im * conjB.re } : a * b;
      realPart[i][j] = prod.re;
      imagPart[i][j] = prod.im;
    }
  }
  const lines = [];
  lines.push("=== State City (Real part of \u03C1) ===");
  lines.push(_renderMatrix3D(realPart, options));
  lines.push("");
  lines.push("=== State City (Imaginary part of \u03C1) ===");
  lines.push(_renderMatrix3D(imagPart, options));
  return lines.join("\n");
}
function _renderMatrix3D(matrix, options = {}) {
  const maxBarHeight = options.maxBarHeight || 8;
  const dim = matrix.length;
  const maxVal = Math.max(...matrix.flat().map(Math.abs));
  if (maxVal < 1e-15) return "(all zeros)";
  const lines = [];
  const barChar = "\u2588";
  for (let row = 0; row < dim; row++) {
    let line = "";
    for (let col = 0; col < dim; col++) {
      const val = matrix[row][col];
      const height = Math.round(Math.abs(val) / maxVal * maxBarHeight);
      const sign = val >= 0 ? "+" : "-";
      if (height === 0) {
        line += "  \xB7  ";
      } else {
        line += sign + barChar.repeat(Math.min(height, 4)) + " ";
      }
    }
    lines.push(line);
  }
  lines.push("");
  lines.push("Values:");
  for (let row = 0; row < dim; row++) {
    let line = "";
    for (let col = 0; col < dim; col++) {
      line += matrix[row][col].toFixed(3).padStart(8) + " ";
    }
    lines.push(line);
  }
  return lines.join("\n");
}
function render_qsphere_ascii(statevector, options = {}) {
  const n = statevector.num_qubits;
  const dim = 1 << n;
  const probs = statevector.probabilities();
  const lines = [];
  lines.push("=== QSphere ===");
  const groups = new Array(n + 1);
  for (let i = 0; i <= n; i++) groups[i] = [];
  for (let i = 0; i < dim; i++) {
    const weight = _hammingWeight(i);
    groups[weight].push({ idx: i, prob: probs[i] });
  }
  for (let w = 0; w <= n; w++) {
    const states = groups[w];
    if (states.length === 0) continue;
    const totalProb = states.reduce((s, st) => s + st.prob, 0);
    const ringRadius = w;
    lines.push(`  Ring ${w} (|${"1".repeat(w)}${"0".repeat(n - w)}\u27E9): total prob = ${totalProb.toFixed(4)}`);
    for (const st of states) {
      const bits = st.idx.toString(2).padStart(n, "0");
      const barLen = Math.round(st.prob * 30);
      lines.push(`    |${bits}\u27E9 ${"\u2588".repeat(barLen)} ${st.prob.toFixed(4)}`);
    }
  }
  return lines.join("\n");
}
function _hammingWeight(n) {
  let count = 0;
  while (n) {
    count += n & 1;
    n >>= 1;
  }
  return count;
}
function _setGrid(grid, row, col, char) {
  if (row >= 0 && row < grid.length && col >= 0 && col < grid[0].length) {
    if (char.length === 1) {
      grid[row][col] = char;
    } else {
      for (let i = 0; i < char.length && col + i < grid[0].length; i++) {
        grid[row][col + i] = char[i];
      }
    }
  }
}

// src/index.js
var circuits = circuits_exports;
function execute(circuit, options = {}) {
  const shots = options.shots || 1024;
  return simulate(circuit, shots, options);
}
export {
  Adam,
  AmplitudeEstimation,
  AmplitudeEstimationResult,
  AnalysisPass,
  ApplyLayout,
  BarrierBeforeFinalMeasurements,
  BaseEstimator,
  BaseSampler,
  BasicSwap,
  Bit,
  COBYLA,
  CONSTANTS,
  CXCancellation,
  CheckMap,
  Chi,
  CircuitInstruction,
  CircuitOp,
  CircuitSampler,
  CircuitStateFn,
  ClassicalRegister,
  Clbit,
  Clifford,
  Collect2qBlocks,
  CommutativeCancellation,
  Complex,
  ComplexMatrix2 as ComplexMatrix,
  ComplexVector,
  ConsolidateBlocks,
  ControlledGate,
  CountOps,
  Counts,
  CouplingMap,
  DAGCircuit,
  DAGFixedPointPass,
  DAGInNode,
  DAGOpNode,
  DAGOutNode,
  DAGRegister,
  DECOMP_RULES,
  DEFAULT_BASIS,
  Decompose,
  DenseLayout,
  DensityMatrix,
  Depth,
  DiagonalGate,
  Estimator,
  EvolvedOp,
  FiniteDiff,
  Gate,
  GateDirection,
  GradientBase,
  GradientDescent,
  Grover,
  GroverResult,
  HamiltonianGate,
  Initialize,
  Instruction,
  Kraus,
  LBFGSB,
  Layout,
  LinearCombination,
  ListOp,
  LookaheadSwap,
  MCMTGate,
  MCPhaseGate,
  MCRXGate,
  MCRYGate,
  MCRZGate,
  MatrixExpectation,
  MatrixOp,
  MergeAdjacentBarriers,
  NFT,
  NaturalGradient,
  NelderMead,
  NoiseModel,
  NumPyMaximumEigensolver,
  NumPyMinimumEigensolver,
  NumPyMinimumEigensolverResult,
  Operator,
  OperatorBase,
  Optimize1qGates,
  Optimize1qGatesDecomposition,
  OptimizeSwapBeforeMeasure,
  OptimizerResult,
  PAULI,
  PTM,
  ParamShift,
  Parameter,
  ParameterExpression,
  ParameterVector,
  PassManager,
  PassManagerConfig,
  Pauli,
  PauliExpectation,
  PauliList,
  PauliOp,
  PauliSumOp,
  PermutationGate,
  PhaseEstimation,
  PhaseEstimationResult,
  QAOA,
  QAOAResult,
  QASM3Parser,
  QASMExporter,
  QASMParser,
  QasmSimulator,
  QuantumChannel,
  QuantumCircuit,
  QuantumError,
  QuantumRegister,
  Qubit,
  ReadoutError,
  Register,
  RemoveBarriers,
  RemoveFinalMeasurements,
  RemoveResetInZeroState,
  Result,
  SLSQP,
  SPSA,
  SabreLayout2 as SabreLayout,
  SabreLayout as SabreLayoutReal,
  SabreSwap2 as SabreSwap,
  SabreSwap as SabreSwapReal,
  Sampler,
  ScalarOp,
  SchmidtDecomposition,
  Size,
  SparsePauliOp,
  StabilizerState,
  StateFn,
  Statevector,
  StatevectorSimulator,
  StochasticSwap,
  SuperOp,
  TransformationPass,
  TrivialLayout,
  UnitaryGate,
  Unroll3qOrMore,
  VQE,
  VQEResult,
  Width,
  _getStdGate,
  _registerParamBuilder,
  _registerStd,
  _setDrawHook,
  amplitude_damping_error,
  average_gate_fidelity,
  bellState,
  bit_flip_error,
  circuit_to_dag,
  circuits,
  collect_1q_runs,
  collect_2q_runs,
  combine_errors,
  commutative_cancellation,
  concurrence,
  consolidate_blocks,
  collect_1q_runs2 as dag_collect_1q_runs,
  collect_2q_runs2 as dag_collect_2q_runs,
  consolidate_blocks2 as dag_consolidate_blocks,
  optimize_cliffords2 as dag_optimize_cliffords,
  remove_diagonal_gates_before_measure2 as dag_remove_diag_before_measure,
  dag_to_circuit,
  decomposeGate,
  depolarizing_error,
  diamond_norm,
  draw_circuit,
  efficientSU2,
  elide_permutations,
  entanglement_of_formation,
  execute,
  gate_fidelity,
  generalized_gates_exports as generalizedGates,
  get_expectation,
  ghzState,
  graphState,
  hiddenLinearFunction,
  iqft,
  kraus_error,
  kronMatrices,
  kronVectors,
  mixed_unitary_error,
  mutual_information,
  optimize_cliffords,
  pauliEvolution,
  pauliTwoDesign,
  pauli_x_error,
  pauli_y_error,
  pauli_z_error,
  phase_damping_error,
  phase_flip_error,
  plot_bloch_vector,
  plot_histogram,
  plot_state_city,
  plot_state_hinton,
  plot_state_qsphere,
  presetPassManager,
  process_fidelity,
  purity,
  qasm3_parse,
  qft,
  qftInverse,
  quantumVolume,
  randomUniform,
  random_clifford,
  random_density_matrix,
  random_pauli,
  random_statevector,
  random_unitary,
  realAmplitudes,
  remove_diagonal_gates_before_measure,
  remove_redundant_gates,
  render_bloch_ascii,
  render_histogram_ascii,
  render_qsphere_ascii,
  render_state_city_ascii,
  reset_error,
  sampleDistribution,
  simulate,
  simulate_noisy,
  standard_gates_exports as standardGates,
  state_fidelity,
  template_optimization,
  transpile,
  twoLocal,
  unitarity
};
//# sourceMappingURL=ketra.js.map
