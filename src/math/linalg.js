// Complex number
export class Complex {
  constructor(real, imag) {
    if (typeof real !== "number") {
      throw new TypeError("Complex real part must be a number");
    }
    this.re = real;
    this.im = (typeof imag === "number") ? imag : 0;
  }

  static from(v) {
    if (v instanceof Complex) return v;
    if (typeof v === "number") return new Complex(v, 0);
    throw new TypeError("Cannot coerce value to Complex");
  }

  static get ZERO() { return new Complex(0, 0); }
  static get ONE()  { return new Complex(1, 0); }
  static get I()    { return new Complex(0, 1); }

  add(other) {
    const o = Complex.from(other);
    return new Complex(this.re + o.re, this.im + o.im);
  }

  sub(other) {
    const o = Complex.from(other);
    return new Complex(this.re - o.re, this.im - o.im);
  }

  mul(other) {
    const o = Complex.from(other);
    return new Complex(
      this.re * o.re - this.im * o.im,
      this.re * o.im + this.im * o.re
    );
  }

  div(other) {
    const o = Complex.from(other);
    const denom = o.re * o.re + o.im * o.im;
    if (denom === 0) throw new Error("Complex division by zero");
    return new Complex(
      (this.re * o.re + this.im * o.im) / denom,
      (this.im * o.re - this.re * o.im) / denom
    );
  }

  conjugate() { return new Complex(this.re, -this.im); }
  negate()    { return new Complex(-this.re, -this.im); }

  // Multiply by a real scalar (convenience method)
  scale(s) { return new Complex(this.re * s, this.im * s); }

  abs()   { return Math.hypot(this.re, this.im); }
  abs2()  { return this.re * this.re + this.im * this.im; }
  arg()   { return Math.atan2(this.im, this.re); }

  exp() {
    const r = Math.exp(this.re);
    return new Complex(r * Math.cos(this.im), r * Math.sin(this.im));
  }

  pow(n) {
    const r = Math.pow(this.abs(), n);
    const theta = this.arg() * n;
    return new Complex(r * Math.cos(theta), r * Math.sin(theta));
  }

  sqrt() {
    const r = this.abs();
    const theta = this.arg() / 2;
    const sr = Math.sqrt(r);
    return new Complex(sr * Math.cos(theta), sr * Math.sin(theta));
  }

  log() {
    return new Complex(Math.log(this.abs()), this.arg());
  }

  equals(other, tol) {
    const o = Complex.from(other);
    const eps = (typeof tol === "number") ? tol : 1e-12;
    return Math.abs(this.re - o.re) < eps && Math.abs(this.im - o.im) < eps;
  }

  toString() {
    const sign = this.im < 0 ? "-" : "+";
    return `${this.re}${sign}${Math.abs(this.im)}i`;
  }

  toJSON() { return { re: this.re, im: this.im }; }
}

// Complex vector
export class ComplexVector {
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
    return new ComplexVector(data);
  }

  static basis(n, index) {
    const v = ComplexVector.zeros(n);
    v.data[index] = Complex.ONE;
    return v;
  }

  static fromArray(arr) {
    return new ComplexVector(arr);
  }

  get(i) { return this.data[i]; }
  set(i, v) { this.data[i] = Complex.from(v); }

  add(other) {
    this._checkSameSize(other);
    const out = new Array(this.size);
    for (let i = 0; i < this.size; i++) out[i] = this.data[i].add(other.data[i]);
    return new ComplexVector(out);
  }

  sub(other) {
    this._checkSameSize(other);
    const out = new Array(this.size);
    for (let i = 0; i < this.size; i++) out[i] = this.data[i].sub(other.data[i]);
    return new ComplexVector(out);
  }

  scale(c) {
    const f = Complex.from(c);
    const out = new Array(this.size);
    for (let i = 0; i < this.size; i++) out[i] = this.data[i].mul(f);
    return new ComplexVector(out);
  }

  inner(other) {
    this._checkSameSize(other);
    let acc = Complex.ZERO;
    for (let i = 0; i < this.size; i++) {
      acc = acc.add(this.data[i].conjugate().mul(other.data[i]));
    }
    return acc;
  }

  dot(other) { return this.inner(other); }

  norm() { return Math.sqrt(this.inner(this).re); }

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
    return new ComplexVector(out);
  }

  conjugate() {
    return new ComplexVector(this.data.map(c => c.conjugate()));
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
    if (!(other instanceof ComplexVector) || other.size !== this.size) return false;
    for (let i = 0; i < this.size; i++) {
      if (!this.data[i].equals(other.data[i], tol)) return false;
    }
    return true;
  }

  toArray() { return this.data.slice(); }
  toRealArray() { return this.data.map(c => c.re); }

  _checkSameSize(other) {
    if (!(other instanceof ComplexVector) || other.size !== this.size) {
      throw new Error(`Size mismatch: ${this.size} vs ${other && other.size}`);
    }
  }
}

// Complex matrix (dense, row-major)
export class ComplexMatrix {
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
    const r = rows.length;
    const c = rows[0].length;
    const data = new Array(r * c);
    for (let i = 0; i < r; i++) {
      for (let j = 0; j < c; j++) data[i * c + j] = Complex.from(rows[i][j]);
    }
    return new ComplexMatrix(r, c, data);
  }

  static identity(n) {
    const m = new ComplexMatrix(n, n);
    for (let i = 0; i < n; i++) m.set(i, i, Complex.ONE);
    return m;
  }

  static zeros(n, m) {
    return new ComplexMatrix(n, m || n);
  }

  static fromDiagonal(diag) {
    const n = diag.length;
    const m = new ComplexMatrix(n, n);
    for (let i = 0; i < n; i++) m.set(i, i, diag[i]);
    return m;
  }

  get(i, j) { return this.data[i * this.cols + j]; }
  set(i, j, v) { this.data[i * this.cols + j] = Complex.from(v); }

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
    return new ComplexMatrix(this.rows, this.cols, out);
  }

  // In-place addition (avoids allocating a new matrix for hot paths like
  // Kraus composition and superoperator assembly).
  addInplace(other) {
    this._checkSameShape(other);
    for (let k = 0; k < this.data.length; k++) {
      this.data[k] = this.data[k].add(other.data[k]);
    }
    return this;
  }

  _checkSameShape(other) {
    if (!(other instanceof ComplexMatrix) ||
        other.rows !== this.rows ||
        other.cols !== this.cols) {
      throw new Error(`Matrix shape mismatch: ${this.rows}x${this.cols} vs ${other && other.rows}x${other && other.cols}`);
    }
  }

  sub(other) {
    this._checkSameShape(other);
    const out = new Array(this.rows * this.cols);
    for (let k = 0; k < out.length; k++) out[k] = this.data[k].sub(other.data[k]);
    return new ComplexMatrix(this.rows, this.cols, out);
  }

  scale(c) {
    const f = Complex.from(c);
    const out = new Array(this.rows * this.cols);
    for (let k = 0; k < out.length; k++) out[k] = this.data[k].mul(f);
    return new ComplexMatrix(this.rows, this.cols, out);
  }

  mul(other) {
    if (other instanceof ComplexMatrix) {
      if (this.cols !== other.rows) {
        throw new Error(`Matrix multiply shape mismatch: ${this.rows}x${this.cols} * ${other.rows}x${other.cols}`);
      }
      const out = new ComplexMatrix(this.rows, other.cols);
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
    // Scalar
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
    const out = new ComplexMatrix(this.rows * other.rows, this.cols * other.cols);
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
    return new ComplexMatrix(this.rows, this.cols, out);
  }

  transpose() {
    const out = new ComplexMatrix(this.cols, this.rows);
    for (let i = 0; i < this.rows; i++) {
      for (let j = 0; j < this.cols; j++) {
        out.set(j, i, this.get(i, j));
      }
    }
    return out;
  }

  dagger() {
    const out = new ComplexMatrix(this.cols, this.rows);
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
    // LU decomposition with partial pivoting
    const n = this.rows;
    const A = this._clone();
    let det = Complex.ONE;
    for (let i = 0; i < n; i++) {
      // Find pivot
      let maxRow = i;
      let maxVal = A.get(i, i).abs();
      for (let k = i + 1; k < n; k++) {
        const v = A.get(k, i).abs();
        if (v > maxVal) { maxVal = v; maxRow = k; }
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
    const eps = (typeof tol === "number") ? tol : 1e-9;
    if (this.rows !== this.cols) return false;
    const prod = this.dagger().mul(this);
    const id = ComplexMatrix.identity(this.rows);
    return prod.equals(id, eps);
  }

  isHermitian(tol) {
    const eps = (typeof tol === "number") ? tol : 1e-9;
    return this.equals(this.dagger(), eps);
  }

  equals(other, tol) {
    const eps = (typeof tol === "number") ? tol : 1e-12;
    if (!(other instanceof ComplexMatrix)) return false;
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
    if (this.rows !== this.cols || this.rows !== (1 << qubits)) {
      throw new Error("partialTrace: matrix must be 2^n x 2^n for n qubits");
    }
    const kept = qubits - 1;
    const dim = 1 << kept;
    const out = ComplexMatrix.zeros(dim, dim);
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
    const V = ComplexMatrix.identity(n);
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
      // Sort by descending eigenvalue
      const indices = eigenvalues.map((v, i) => i).sort((a, b) => eigenvalues[b] - eigenvalues[a]);
      const S = indices.map(i => Math.sqrt(Math.max(0, eigenvalues[i])));
      const V = eigenvectors;
      // U = A V S^-1
      const U = ComplexMatrix.zeros(m, Math.min(m, n));
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
      // Vh = V^dagger (rows are eigenvectors)
      const Vh = ComplexMatrix.zeros(k, n);
      for (let i = 0; i < k; i++) {
        for (let j = 0; j < n; j++) {
          Vh.set(i, j, V.get(j, indices[i]).conjugate());
        }
      }
      return { U, S, Vh };
    } else {
      // Transpose trick
      const t = this.transpose();
      const r = t.svd();
      return { U: r.Vh.transpose(), S: r.S, Vh: r.U.transpose() };
    }
  }

  // Matrix inverse via Gauss-Jordan elimination
  inverse() {
    if (this.rows !== this.cols) throw new Error("inverse requires square matrix");
    const n = this.rows;
    const A = this._clone();
    const I = ComplexMatrix.identity(n);
    for (let i = 0; i < n; i++) {
      let maxRow = i;
      let maxVal = A.get(i, i).abs();
      for (let k = i + 1; k < n; k++) {
        const v = A.get(k, i).abs();
        if (v > maxVal) { maxVal = v; maxRow = k; }
      }
      if (maxVal < 1e-15) throw new Error("Matrix is singular");
      if (maxRow !== i) {
        for (let j = 0; j < n; j++) {
          const tmp = A.get(i, j); A.set(i, j, A.get(maxRow, j)); A.set(maxRow, j, tmp);
          const tmp2 = I.get(i, j); I.set(i, j, I.get(maxRow, j)); I.set(maxRow, j, tmp2);
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
    // Scaling: choose s such that ||A/2^s|| < 1/2
    let norm = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) norm = Math.max(norm, this.get(i, j).abs());
    let s = 0;
    while (norm > 0.5) { norm /= 2; s++; }
    const scaled = this.scale(1 / Math.pow(2, s));
    // Taylor series to order 20
    let result = ComplexMatrix.identity(n);
    let term = ComplexMatrix.identity(n);
    for (let k = 1; k <= 20; k++) {
      term = term.mul(scaled).scale(1 / k);
      result = result.add(term);
      if (term.data.every(v => v.abs() < 1e-18)) break;
    }
    // Squaring s times
    for (let i = 0; i < s; i++) result = result.mul(result);
    return result;
  }

  // Kronecker product alias
  kron(other) { return this.tensor(other); }

  _clone() {
    return new ComplexMatrix(this.rows, this.cols, this.data.slice());
  }

  toRows() {
    const out = new Array(this.rows);
    for (let i = 0; i < this.rows; i++) out[i] = this.getRow(i);
    return out;
  }

  toString() {
    const rows = this.toRows().map(r => r.map(c => c.toString()).join("  "));
    return rows.join("\n");
  }
}

// Helpers

function _insertBit(bits, pos, value) {
  const mask = (1 << pos) - 1;
  const low = bits & mask;
  const high = (bits & ~mask) << 1;
  return high | low | (value << pos);
}

export function kronMatrices(matrices) {
  if (matrices.length === 0) return ComplexMatrix.identity(1);
  let acc = matrices[0];
  for (let i = 1; i < matrices.length; i++) acc = acc.tensor(matrices[i]);
  return acc;
}

export function kronVectors(vectors) {
  if (vectors.length === 0) return new ComplexVector([Complex.ONE]);
  let acc = vectors[0];
  for (let i = 1; i < vectors.length; i++) acc = acc.tensor(vectors[i]);
  return acc;
}

export function randomUniform() {
  const c = (typeof globalThis !== "undefined" && globalThis.crypto && globalThis.crypto.getRandomValues);
  if (c) {
    const buf = new Uint32Array(1);
    globalThis.crypto.getRandomValues(buf);
    return buf[0] / 0x100000000;
  }
  return Math.random();
}

export function sampleDistribution(probs, rng) {
  const r = (typeof rng === "function") ? rng() : randomUniform();
  let cum = 0;
  for (let i = 0; i < probs.length; i++) {
    cum += probs[i];
    if (r < cum) return i;
  }
  return probs.length - 1;
}

// Common single-qubit gate matrices
export const PAULI = {
  I: ComplexMatrix.identity(2),
  X: ComplexMatrix.fromRows([[new Complex(0), new Complex(1)], [new Complex(1), new Complex(0)]]),
  Y: ComplexMatrix.fromRows([[new Complex(0), new Complex(0, -1)], [new Complex(0, 1), new Complex(0)]]),
  Z: ComplexMatrix.fromRows([[new Complex(1), new Complex(0)], [new Complex(0), new Complex(-1)]]),
};

export const CONSTANTS = {
  ZERO: Complex.ZERO,
  ONE: Complex.ONE,
  I: Complex.I,
  H: ComplexMatrix.fromRows([
    [new Complex(1 / Math.SQRT2), new Complex(1 / Math.SQRT2)],
    [new Complex(1 / Math.SQRT2), new Complex(-1 / Math.SQRT2)],
  ]),
  S: ComplexMatrix.fromRows([
    [new Complex(1), new Complex(0)],
    [new Complex(0), new Complex(0, 1)],
  ]),
  T: ComplexMatrix.fromRows([
    [new Complex(1), new Complex(0)],
    [new Complex(0), new Complex(Math.cos(Math.PI / 4), Math.sin(Math.PI / 4))],
  ]),
};
