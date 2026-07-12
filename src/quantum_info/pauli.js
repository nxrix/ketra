/**
 * pauli.js - Pauli operators, PauliList, and SparsePauliOp.
 *
 * *
 * A Pauli is represented internally as two bit arrays: z (Z-part) and x (X-part).
 * The full operator is (-i)^(x·z) * X^x · Z^z (symplectic representation).
 */

import { Complex, ComplexMatrix, PAULI } from "./../math/linalg.js";

const PAULI_MATRICES = {
  "I": ComplexMatrix.identity(2),
  "X": PAULI.X,
  "Y": PAULI.Y,
  "Z": PAULI.Z,
};

export class Pauli {
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
      if (ch === "X") { this.x[i] = 1; }
      else if (ch === "Z") { this.z[i] = 1; }
      else if (ch === "Y") { this.x[i] = 1; this.z[i] = 1; }
      else if (ch === "I") { /* both zero */ }
      else throw new Error(`Invalid Pauli label char: ${ch}`);
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
    return new Pauli(label);
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
      // (x1^x2)(z1^z2) is 1 only when both bits differ; subtract that.
      const yzProd = (x1 ^ x2) & (z1 ^ z2);
      phase += x1 * z1 + x2 * z2 - yzProd + 2 * (z1 & x2);
    }
    phase = ((phase % 4) + 4) % 4;
    const result = Pauli.fromSymplectic(newX, newZ);
    const phaseFactor = [new Complex(1, 0), new Complex(0, 1), new Complex(-1, 0), new Complex(0, -1)][phase];
    return { pauli: result, phase: phaseFactor };
  }

  multiply(other) {
    const r = this.compose(other);
    return r.pauli.to_matrix().scale(r.phase);
  }

  to_label() { return this.label; }

  equals(other) {
    if (!(other instanceof Pauli)) return false;
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

  anticommutes(other) { return !this.commutes(other); }

  // Apply X/Y/Z to specific qubits (single-qubit Pauli evolution)
  evolve(other) {
    // P -> P * other (Pauli group multiplication, ignoring phase)
    if (this.num_qubits !== other.num_qubits) {
      throw new Error("evolve: qubit count mismatch");
    }
    const newX = new Array(this.num_qubits).fill(0);
    const newZ = new Array(this.num_qubits).fill(0);
    for (let i = 0; i < this.num_qubits; i++) {
      newX[i] = this.x[i] ^ other.x[i];
      newZ[i] = this.z[i] ^ other.z[i];
    }
    return Pauli.fromSymplectic(newX, newZ);
  }

  toString() { return `Pauli(${this.label})`; }
}

export class PauliList {
  constructor(labels) {
    if (Array.isArray(labels)) {
      this.paulis = labels.map(l => l instanceof Pauli ? l : new Pauli(l));
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

  get size() { return this.paulis.length; }
  get length() { return this.paulis.length; }
  get(i) { return this.paulis[i]; }
  [Symbol.iterator]() { return this.paulis[Symbol.iterator](); }
  to_labels() { return this.paulis.map(p => p.label); }
  to_matrices() { return this.paulis.map(p => p.to_matrix()); }
  weights() { return this.paulis.map(p => p.weight()); }

  append(other) {
    if (other instanceof Pauli) this.paulis.push(other);
    else if (other instanceof PauliList) this.paulis.push(...other.paulis);
    else throw new TypeError("append: expects Pauli or PauliList");
  }
}

export class SparsePauliOp {
  constructor(paulis, coeffs = null) {
    if (paulis instanceof PauliList) this.paulis = paulis;
    else this.paulis = new PauliList(paulis);
    const n = this.paulis.size;
    if (coeffs == null) {
      this.coeffs = new Array(n).fill(0).map(() => new Complex(1, 0));
    } else {
      this.coeffs = coeffs.map(c => c instanceof Complex ? c : new Complex(c, 0));
    }
    if (this.coeffs.length !== n) {
      throw new Error("SparsePauliOp: paulis and coeffs length mismatch");
    }
    this.num_qubits = this.paulis.num_qubits;
  }

  static from_list(list) {
    const paulis = list.map(([p]) => p instanceof Pauli ? p : new Pauli(p));
    const coeffs = list.map(([_, c]) => c instanceof Complex ? c : new Complex(c, 0));
    return new SparsePauliOp(new PauliList(paulis), coeffs);
  }

  static from_sparse_list(zlist, numQubits) {
    const paulis = zlist.map(([label]) => {
      let padded = label;
      while (padded.length < numQubits) padded = "I" + padded;
      return new Pauli(padded);
    });
    const coeffs = zlist.map(([_, c]) => c instanceof Complex ? c : new Complex(c, 0));
    return new SparsePauliOp(new PauliList(paulis), coeffs);
  }

  static from_operator(operator) {
    // Decompose an operator into a sum of Paulis (exponential cost in qubit count)
    const n = operator.num_qubits;
    const dim = 1 << n;
    const pauliLabels = [];
    function genLabels(prefix, depth) {
      if (depth === 0) { pauliLabels.push(prefix); return; }
      for (const c of ["I", "X", "Y", "Z"]) genLabels(prefix + c, depth - 1);
    }
    genLabels("", n);
    const results = [];
    for (const label of pauliLabels) {
      const p = new Pauli(label);
      const m = p.to_matrix();
      // Tr(P * O) / 2^n gives the coefficient
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
    return SparsePauliOp.from_list(results);
  }

  to_matrix() {
    const dim = 1 << this.num_qubits;
    const result = ComplexMatrix.zeros(dim, dim);
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
    if (other instanceof SparsePauliOp) {
      return new SparsePauliOp(
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
    if (other instanceof SparsePauliOp) {
      const newPaulis = [];
      const newCoeffs = [];
      for (let i = 0; i < this.paulis.size; i++) {
        const pi = this.paulis.get(i);
        const ci = this.coeffs[i];
        for (let j = 0; j < other.paulis.size; j++) {
          const pj = other.paulis.get(j);
          const cj = other.coeffs[j];
          const r = pi.compose(pj);
          newPaulis.push(r.pauli);
          newCoeffs.push(ci.mul(r.phase).mul(cj));
        }
      }
      return new SparsePauliOp(new PauliList(newPaulis), newCoeffs).simplify();
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
      // For each Pauli, the non-zero entries are at (i, i XOR mask) with
      // phase iPow * (-1)^popcount((i XOR mask) & z) * coeff.
      let mask = 0, yMask = 0, zOnlyMask = 0, yCount = 0;
      for (let q = 0; q < this.num_qubits; q++) {
        if (x[q]) {
          mask |= (1 << q);
          if (z[q]) { yMask |= (1 << q); yCount++; }
        } else if (z[q]) {
          zOnlyMask |= (1 << q);
        }
      }
      const iPow = [new Complex(1, 0), new Complex(0, 1), new Complex(-1, 0), new Complex(0, -1)][yCount % 4];
      for (let i = 0; i < dim; i++) {
        const j = i ^ mask;
        // Phase from Z-only qubits and Y qubits (depends on j = i XOR mask).
        let parity = 0;
        for (let q = 0; q < this.num_qubits; q++) {
          if (zOnlyMask & (1 << q) & i) parity ^= 1;
          if (yMask & (1 << q) & j) parity ^= 1;
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
        // Check if pi commutes with every Pauli in g.
        let commutes = true;
        for (const pj of g.paulis) {
          if (!pi.commutes(pj)) { commutes = false; break; }
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
    return groups.map(g => new SparsePauliOp(new PauliList(g.paulis), g.coeffs));
  }

  // Sort terms by weight (descending by default) for measurement efficiency.
  sort(reverse = true) {
    const indexed = this.paulis.paulis.map((p, i) => ({ p, c: this.coeffs[i], w: p.weight() }));
    indexed.sort((a, b) => reverse ? b.w - a.w : a.w - b.w);
    return new SparsePauliOp(
      new PauliList(indexed.map(x => x.p)),
      indexed.map(x => x.c),
    );
  }

  // Split into chunks of mutually-commuting terms (alias for noncommutation_groups
  // with a size cap).
  chunk(maxSize = Infinity) {
    const groups = this.noncommutation_groups();
    if (maxSize === Infinity) return groups;
    // Split groups larger than maxSize.
    const result = [];
    for (const g of groups) {
      if (g.paulis.size <= maxSize) {
        result.push(g);
      } else {
        for (let i = 0; i < g.paulis.size; i += maxSize) {
          const slice = g.paulis.paulis.slice(i, i + maxSize);
          const coeffs = g.coeffs.slice(i, i + maxSize);
          result.push(new SparsePauliOp(new PauliList(slice), coeffs));
        }
      }
    }
    return result;
  }

  conjugate() {
    const newCoeffs = this.coeffs.map(c => c.conjugate());
    return new SparsePauliOp(new PauliList(this.paulis.paulis), newCoeffs);
  }

  transpose() {
    const newCoeffs = this.coeffs.map((c, i) => {
      let yCount = 0;
      const p = this.paulis.get(i);
      for (let j = 0; j < p.num_qubits; j++) if (p.x[j] && p.z[j]) yCount++;
      const sign = (yCount % 2 === 0) ? 1 : -1;
      return new Complex(c.re * sign, c.im * sign);
    });
    return new SparsePauliOp(new PauliList(this.paulis.paulis), newCoeffs);
  }

  adjoint() {
    const newCoeffs = this.coeffs.map(c => c.conjugate());
    const adjustedPaulis = this.paulis.paulis.map(p => {
      let yCount = 0;
      for (let i = 0; i < p.num_qubits; i++) if (p.x[i] && p.z[i]) yCount++;
      return { pauli: p, sign: (yCount % 2 === 0) ? 1 : -1 };
    });
    const finalCoeffs = newCoeffs.map((c, i) => new Complex(c.re * adjustedPaulis[i].sign, c.im * adjustedPaulis[i].sign));
    return new SparsePauliOp(new PauliList(this.paulis.paulis), finalCoeffs);
  }

  copy() {
    return new SparsePauliOp(
      new PauliList(this.paulis.paulis.map(p => new Pauli(p.label))),
      this.coeffs.map(c => new Complex(c.re, c.im))
    );
  }

  simplify(tol = 1e-12) {
    const map = new Map();
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
    return new SparsePauliOp(new PauliList(newPaulis), newCoeffs);
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
      // Classify each qubit's contribution:
      //   I: no effect
      //   Z: phase (-1)^(bit q of i), no bit flip
      //   X: bit flip on qubit q, no phase
      //   Y: bit flip on qubit q AND phase (-1)^(bit q of j) * i
      //     (since Y|b> = i * (-1)^b * |1-b>)
      // For the full Pauli P = prod_q P_q applied to |i>:
      //   P|i> = (i^yCount) * (-1)^(popcount(i & zOnlyMask)) * (-1)^(popcount((i^mask) & yMask)) * |i ^ mask>
      // where:
      //   mask     = bits where P_q is X or Y
      //   yMask    = bits where P_q is Y
      //   zOnlyMask = bits where P_q is Z (no X)
      let mask = 0, yMask = 0, zOnlyMask = 0, yCount = 0;
      for (let q = 0; q < n; q++) {
        if (x[q]) {
          mask |= (1 << q);
          if (z[q]) { yMask |= (1 << q); yCount++; }
        } else if (z[q]) {
          zOnlyMask |= (1 << q);
        }
      }
      const iPow = [new Complex(1, 0), new Complex(0, 1), new Complex(-1, 0), new Complex(0, -1)][yCount % 4];
      const iPowRe = iPow.re, iPowIm = iPow.im;
      const cr = coeff.re, ci = coeff.im;
      for (let i = 0; i < dim; i++) {
        const j = i ^ mask;
        // Phase from Z-only qubits (depends on original i): (-1)^popcount(i & zOnlyMask)
        const zPhase = _popcountParity(i & zOnlyMask);
        // Phase from Y qubits (depends on j = i^mask): (-1)^popcount(j & yMask)
        const yPhase = _popcountParity(j & yMask);
        const phase = zPhase * yPhase; // +/- 1
        const psiI = sv.data[i];
        const psiJ = sv.data[j];
        // (P psi)_i = iPow * phase * psiJ
        // contribution to <psi|P|psi> = conj(psiI) * iPow * phase * psiJ
        // Then multiply by coeff.
        const pRe = iPowRe * phase;
        const pIm = iPowIm * phase;
        // conj(psiI) = (psiI.re, -psiI.im)
        const cRe = psiI.re * pRe - (-psiI.im) * pIm;
        const cIm = psiI.re * pIm + (-psiI.im) * pRe;
        // * psiJ
        const tRe = cRe * psiJ.re - cIm * psiJ.im;
        const tIm = cRe * psiJ.im + cIm * psiJ.re;
        // * coeff, accumulate
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
}

// Returns +1 or -1 depending on whether the number of set bits in `x` is
// even or odd. Uses a bit-twiddling trick for speed.
function _popcountParity(x) {
  x ^= x >> 16;
  x ^= x >> 8;
  x ^= x >> 4;
  x &= 0xf;
  // Lookup: parity of each 4-bit value (0..15)
  return ((0x6996 >> x) & 1) ? -1 : 1;
}
