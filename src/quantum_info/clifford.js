/**
 * clifford.js - Clifford and StabilizerState classes.
 *
 * *
 * A Clifford is represented by the symplectic stabilizer table:
 *   - x: n×n binary matrix (X parts of stabilizer generators)
 *   - z: n×n binary matrix (Z parts of stabilizer generators)
 *   - signs: n-bit vector (the +/- sign of each generator)
 *
 * Rows 0..n-1 are the stabilizers, rows n..2n-1 are the destabilizers.
 */

import { Complex, ComplexMatrix, ComplexVector, PAULI } from "./../math/linalg.js";
import { Pauli } from "./pauli.js";
import { Statevector } from "./statevector.js";
import { QuantumCircuit } from "./../core/circuit.js";
import { UnitaryGate } from "./../library/extra_gates.js";

export class Clifford {
  constructor(data) {
    // Accept either:
    //  - { x: number[][], z: number[][], signs: number[] } (symplectic form)
    //  - A label string of basis states like '01' (computational)
    //  - An array of Pauli labels (stabilizer generators)
    if (data && data.x && data.z) {
      this.x = data.x.map(row => row.slice());
      this.z = data.z.map(row => row.slice());
      this.signs = (data.signs || new Array(data.x.length).fill(0)).slice();
      // The symplectic tableau has 2n rows (n destabilizers + n stabilizers),
      // so num_qubits = x.length / 2. (Older code incorrectly set
      // num_qubits = x.length, which broke compose / fromStabilizers /
      // to_matrix whenever the tableau form was used.)
      this.num_qubits = data.x.length / 2;
      if (!Number.isInteger(this.num_qubits)) {
        throw new Error(
          `Clifford: tableau must have 2n rows (got ${data.x.length}); ` +
          `num_qubits would be non-integer.`
        );
      }
    } else if (typeof data === "string") {
      // Build from a computational basis label.
      // The stabilizer tableau has 2n rows: the first n are destabilizers
      // (X_0, X_1, ..., X_{n-1}) and the last n are stabilizers
      // (Z_i if |q_i>=0, -X_i if |q_i>=1).
      this.num_qubits = data.length;
      this.x = [];
      this.z = [];
      this.signs = [];
      // Destabilizers: X_0, X_1, ..., X_{n-1}
      for (let i = 0; i < this.num_qubits; i++) {
        const xRow = new Array(this.num_qubits).fill(0);
        const zRow = new Array(this.num_qubits).fill(0);
        xRow[i] = 1;
        this.x.push(xRow);
        this.z.push(zRow);
        this.signs.push(0);
      }
      // Stabilizers: Z_i for |0>, -Z_i for |1> (since Z|0>=+|0>, Z|1>=-|1>,
      // so |1> is the +1 eigenstate of -Z).
      for (let i = 0; i < this.num_qubits; i++) {
        const xRow = new Array(this.num_qubits).fill(0);
        const zRow = new Array(this.num_qubits).fill(0);
        const bit = parseInt(data[data.length - 1 - i], 10);
        if (bit === 1) {
          // |1>: stabilizer is -Z_i (sign=1 means negative)
          zRow[i] = 1;
          this.signs.push(1);
        } else {
          // |0>: stabilizer is +Z_i
          zRow[i] = 1;
          this.signs.push(0);
        }
        this.x.push(xRow);
        this.z.push(zRow);
      }
    } else if (Array.isArray(data) && data.length > 0 && typeof data[0] === "string") {
      // Array of Pauli labels as stabilizers. Caller is responsible for
      // providing a full set of n mutually commuting stabilizers; we add
      // canonical X_i destabilizers.
      const paulis = data.map(l => new Pauli(l));
      this.num_qubits = paulis[0].num_qubits;
      this.x = [];
      this.z = [];
      this.signs = [];
      // Destabilizers: X_i
      for (let i = 0; i < this.num_qubits; i++) {
        const xRow = new Array(this.num_qubits).fill(0);
        const zRow = new Array(this.num_qubits).fill(0);
        xRow[i] = 1;
        this.x.push(xRow);
        this.z.push(zRow);
        this.signs.push(0);
      }
      // Stabilizers from the input
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
    // Start with |0...0> stabilizers: Z_0, Z_1, ..., Z_{n-1}
    const n = circuit.num_qubits;
    const cliff = Clifford.fromLabel("0".repeat(n));
    for (const ci of circuit.data) {
      const op = ci.operation;
      if (op.name === "barrier" || op.name === "measure" || op.name === "reset") continue;
      if (op.num_qubits === 0) continue;
      const targets = ci.qubits.map(q => circuit._qubit_index.get(q));
      cliff._applyGate(op.name, targets, op.params || []);
    }
    return cliff;
  }

  static fromLabel(label) {
    return new Clifford(label);
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
    // Compute U P U^† for each single-qubit Pauli P on each qubit.
    // The result must be a tensor product of Paulis (up to a phase ±1)
    // for the matrix to be a Clifford.
    const matrixDag = matrix.dagger();
    const paulis = [
      { name: "X", mat: PAULI.X },
      { name: "Z", mat: PAULI.Z },
    ];
    const destabs = []; // for each qubit i: U X_i U^†
    const stabs = [];   // for each qubit i: U Z_i U^†
    for (let q = 0; q < n; q++) {
      for (const { name, mat } of paulis) {
        // Embed the single-qubit Pauli into the full Hilbert space.
        const fullP = _embedPauli(mat, n, q);
        // U P U^†
        const conj = matrix.mul(fullP).mul(matrixDag);
        // Identify the resulting Pauli string.
        const result = _identifyPauli(conj, n);
        if (result === null) {
          throw new Error(
            `Clifford.fromMatrix: U ${name}_${q} U^† is not a Pauli (matrix is not a Clifford)`
          );
        }
        if (name === "X") destabs.push(result);
        else stabs.push(result);
      }
    }
    // Build the tableau. Rows 0..n-1 are destabilizers (U X_i U^†),
    // rows n..2n-1 are stabilizers (U Z_i U^†).
    const x = destabs.map(d => d.x).concat(stabs.map(s => s.x));
    const z = destabs.map(d => d.z).concat(stabs.map(s => s.z));
    const signs = destabs.map(d => d.sign).concat(stabs.map(s => s.sign));
    return new Clifford({ x, z, signs });
  }

  static random(numQubits, seed = null) {
    let rng = Math.random;
    if (seed != null) {
      let s = seed >>> 0;
      rng = () => {
        s = (s + 0x6D2B79F5) | 0;
        let t = Math.imul(s ^ (s >>> 15), 1 | s);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
      };
    }
    // Generate random stabilizer generators via random Clifford actions
    const c = Clifford.fromLabel("0".repeat(numQubits));
    for (let q = 0; q < numQubits; q++) {
      // Apply random single-qubit Clifford
      const r = rng();
      if (r < 0.25) {} // I
      else if (r < 0.5) c._h(q);
      else if (r < 0.75) { c._h(q); c._s(q); }
      else c._s(q);
    }
    // Apply random entangling gates
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
      case "h": case "H": this._h(targets[0]); break;
      case "s": case "S": this._s(targets[0]); break;
      case "sdg": case "SDG": this._s(targets[0]); this._s(targets[0]); this._s(targets[0]); break;
      // X = H · Z · H  (H swaps X and Z, so conjugating Z by H gives X).
      case "x": case "X": this._h(targets[0]); this._s(targets[0]); this._s(targets[0]); this._h(targets[0]); break;
      // Z = Z (identity on the tableau for Z itself is just applying Z,
      // which is S^2 — but to keep the convention consistent we apply
      // S twice, which is exactly the Z action in the H/S basis).
      case "z": case "Z": this._s(targets[0]); this._s(targets[0]); break;
      // Y = i·X·Z = X·Z·(phase). On the tableau, Y acts as (X·Z) up to a
      // phase, which we can realize as Z then X (or equivalently S^2 then
      // HSH). Using X·Z means: apply Z first (S^2), then X (H S^2 H).
      case "y": case "Y": this._s(targets[0]); this._s(targets[0]); this._h(targets[0]); this._s(targets[0]); this._s(targets[0]); this._h(targets[0]); break;
      case "cx": case "CX": this._cx(targets[0], targets[1]); break;
      case "cz": case "CZ": this._h(targets[1]); this._cx(targets[0], targets[1]); this._h(targets[1]); break;
      // CY = (I⊗Sdg) · CX · (I⊗S). Decompose Sdg as S^3.
      case "cy": case "CY":
        this._s(targets[1]); this._s(targets[1]); this._s(targets[1]);
        this._cx(targets[0], targets[1]);
        this._s(targets[1]);
        break;
      case "swap": case "SWAP": this._cx(targets[0], targets[1]); this._cx(targets[1], targets[0]); this._cx(targets[0], targets[1]); break;
      case "sx": case "SX":
        // SX = sqrt(X) = H · sqrt(Z) · H = H · S · H. (Verified: HSH
        // squared is HSH·HSH = HS(HH)SH = HS^2 H = HZH = X, and (SX)^2 = X.)
        this._h(targets[0]); this._s(targets[0]); this._h(targets[0]);
        break;
      case "sxdg": case "SXDG":
        // SX† = H · S† · H = H · S^3 · H.
        this._h(targets[0]); this._s(targets[0]); this._s(targets[0]); this._s(targets[0]); this._h(targets[0]);
        break;
      case "id": case "I": break;
      default:
        throw new Error(`Clifford: cannot apply gate ${name} (not Clifford)`);
    }
  }

  // Apply H to qubit q
  _h(q) {
    for (let i = 0; i < this.num_qubits * 2; i++) {
      const xBit = this.x[i][q];
      const zBit = this.z[i][q];
      // Swap x and z
      this.x[i][q] = zBit;
      this.z[i][q] = xBit;
      // Update sign: -1 if x=1 and z=1 (Y -> -Y under H)
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
      // S: X -> Y, Y -> -X, Z -> Z
      // (x, z) -> (x, x XOR z), sign flip if x=1 and z=1
      this.z[i][q] = xBit ^ zBit;
      if (xBit === 1 && zBit === 1) {
        this.signs[i] ^= 1;
      }
    }
  }

  // Apply CNOT(control, target)
  _cx(control, target) {
    for (let i = 0; i < this.num_qubits * 2; i++) {
      // Update sign
      const xc = this.x[i][control];
      const zc = this.z[i][control];
      const xt = this.x[i][target];
      const zt = this.z[i][target];
      // Phase change: -1 if (xc=1, zc=0, xt=0, zt=1) -> product Y_c * (-X_t * Z_t)
      if (xc === 1 && zc === 0 && xt === 0 && zt === 1) {
        this.signs[i] ^= 1;
      }
      // X_c -> X_c X_t
      this.x[i][target] = xc ^ xt;
      // Z_t -> Z_c Z_t
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

    // Compute U|0> as a stabilizer state.
    const uZero = new StabilizerState(this).to_statevector();
    const uZeroVec = uZero._data.data; // array of Complex

    // Precompute the destabilizer Paulis (top n rows) as {x, z, sign}.
    const destabs = [];
    for (let k = 0; k < n; k++) {
      destabs.push({
        x: this.x[k].slice(),
        z: this.z[k].slice(),
        sign: this.signs[k],
      });
    }

    const m = ComplexMatrix.zeros(dim, dim);

    for (let i = 0; i < dim; i++) {
      // Compute P_i = prod_k destab[k]^{i_k} as (x, z, fullPhase).
      // Start with the identity.
      let pX = new Array(n).fill(0);
      let pZ = new Array(n).fill(0);
      let pPhase = new Complex(1, 0);
      for (let k = 0; k < n; k++) {
        if ((i >> k) & 1) {
          const d = destabs[k];
          // Multiply current (pX, pZ, pPhase) by destab[k] = (-1)^d.sign · i^(d.x·d.z) · X^d.x · Z^d.z.
          let dFullPhase = d.sign ? new Complex(-1, 0) : new Complex(1, 0);
          let dIPhase = 0;
          for (let q = 0; q < n; q++) {
            if (d.x[q] && d.z[q]) dIPhase = (dIPhase + 1) % 4;
          }
          const iPhaseFactor = [
            new Complex(1, 0), new Complex(0, 1),
            new Complex(-1, 0), new Complex(0, -1),
          ][dIPhase];
          dFullPhase = dFullPhase.mul(iPhaseFactor);

          // Anticommutation sign: Z^pZ X^d.x = (-1)^(pZ · d.x) X^d.x Z^pZ.
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

      // Apply P_i = pPhase · X^pX · Z^pZ to U|0>:
      //   (Z^pZ psi)(j)        = (-1)^{popcount(j AND pZ)} psi(j)
      //   (X^pX phi)(j)        = phi(j XOR pX)
      //   (X^pX Z^pZ psi)(j)   = (Z^pZ psi)(j XOR pX)
      //                        = (-1)^{popcount((j XOR pX) AND pZ)} psi(j XOR pX)
      for (let j = 0; j < dim; j++) {
        const srcIdx = j ^ _bitsToNum(pX, n);
        let parity = 0;
        for (let q = 0; q < n; q++) {
          if (pZ[q] && ((srcIdx >> q) & 1)) parity ^= 1;
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

    // Destabilizers of result: U_self * (destabilizer_other[k]) * U_self^†
    for (let k = 0; k < n; k++) {
      const r = this._conjugatePauli(other.x[k], other.z[k], other.signs[k]);
      newX.push(r.x);
      newZ.push(r.z);
      newSigns.push(r.sign);
    }
    // Stabilizers of result: U_self * (stabilizer_other[k]) * U_self^†
    for (let k = 0; k < n; k++) {
      const r = this._conjugatePauli(other.x[n + k], other.z[n + k], other.signs[n + k]);
      newX.push(r.x);
      newZ.push(r.z);
      newSigns.push(r.sign);
    }

    return new Clifford({
      x: newX,
      z: newZ,
      signs: newSigns,
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

    // Track the result as a Pauli (currentPauli) plus an extra phase
    // (currentPhase). Pauli.compose gives us the new Pauli and the phase
    // arising from anticommutation; we fold in the per-row signs from
    // this tableau and the input sign.
    let currentPauli = new Pauli("I".repeat(n)); // P(0, 0) = identity
    let currentPhase = new Complex(1, 0);

    // Input sign: (-1)^inSign
    if (inSign) currentPhase = currentPhase.mul(new Complex(-1, 0));

    // Input i^(inX·inZ) factor (from the P(x,z) convention).
    let inIPhase = 0;
    for (let q = 0; q < n; q++) {
      if (inX[q] && inZ[q]) inIPhase = (inIPhase + 1) % 4;
    }
    const iPhaseFactor = [
      new Complex(1, 0), new Complex(0, 1),
      new Complex(-1, 0), new Complex(0, -1),
    ][inIPhase];
    currentPhase = currentPhase.mul(iPhaseFactor);

    // U X^inX U^† = prod_k (destab_self[k])^{inX_k}
    for (let k = 0; k < n; k++) {
      if (inX[k]) {
        const destab = Pauli.fromSymplectic(this.x[k], this.z[k]);
        const r = currentPauli.compose(destab);
        currentPauli = r.pauli;
        currentPhase = currentPhase.mul(r.phase);
        if (this.signs[k]) currentPhase = currentPhase.mul(new Complex(-1, 0));
      }
    }

    // U Z^inZ U^† = prod_k (stab_self[k])^{inZ_k}
    for (let k = 0; k < n; k++) {
      if (inZ[k]) {
        const stab = Pauli.fromSymplectic(this.x[n + k], this.z[n + k]);
        const r = currentPauli.compose(stab);
        currentPauli = r.pauli;
        currentPhase = currentPhase.mul(r.phase);
        if (this.signs[n + k]) currentPhase = currentPhase.mul(new Complex(-1, 0));
      }
    }

    // The result must be ±P(x, z); the phase must be real.
    if (Math.abs(currentPhase.im) > 1e-9) {
      throw new Error(
        `Clifford._conjugatePauli: result has imaginary phase ${currentPhase}; ` +
        `expected ±1 (input was a Pauli from another Clifford tableau).`
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
        sign: this.signs[i],
      })),
      destabilizer: this.x.slice(this.num_qubits).map((xRow, i) => ({
        x: xRow.slice(),
        z: this.z[this.num_qubits + i].slice(),
        sign: this.signs[this.num_qubits + i],
      })),
    };
  }

  // Adjoint (inverse) of this Clifford: returns the Clifford U^† such that
  // U · U^† = I. Computed by inverting the tableau's symplectic part and
  // recomputing the signs.
  adjoint() {
    const n = this.num_qubits;
    // The Clifford tableau is a 2n × 2n symplectic matrix over GF(2).
    // Its inverse can be computed via Gaussian elimination, but the sign
    // bookkeeping is intricate. The simplest correct approach: build the
    // unitary via to_matrix(), then dagger it, then reconstruct the
    // Clifford from a circuit that applies the daggered unitary.
    // Since to_matrix is exponential, this is only practical for small n.
    // For larger n, we'd implement the symplectic inverse directly.
    // For now, we use the matrix approach.
    const m = this.to_matrix();
    const mDag = m.dagger();
    // Build a circuit that applies mDag and convert to Clifford.
    // Easier: apply the inverse of each generator step.
    // Since we don't have the original gate sequence, we use the matrix.
    return Clifford.fromMatrix(mDag);
  }

  // Alias for adjoint.
  inverse() { return this.adjoint(); }

  // Compose this Clifford with itself n times. For n=0 returns identity.
  power(n) {
    if (n === 0) return Clifford.fromLabel("0".repeat(this.num_qubits));
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
    if (!(other instanceof Clifford)) return false;
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
}

// Build a Clifford from stabilizer generators (Pauli labels)
Clifford.fromStabilizers = function(labels) {
  const paulis = labels.map(l => new Pauli(l));
  const n = paulis[0].num_qubits;
  // Destabilizers are the canonical X_i (or some other mutually unbiased set)
  const x = [];
  const z = [];
  const signs = [];
  // Destabilizers: X_0, X_1, ..., X_{n-1}
  for (let i = 0; i < n; i++) {
    const xRow = new Array(n).fill(0);
    const zRow = new Array(n).fill(0);
    xRow[i] = 1;
    x.push(xRow);
    z.push(zRow);
    signs.push(0);
  }
  // Stabilizers: from labels
  for (const p of paulis) {
    x.push(p.x.slice());
    z.push(p.z.slice());
    signs.push(0);
  }
  return new Clifford({ x, z, signs });
};

export class StabilizerState {
  constructor(clifford) {
    if (clifford instanceof Clifford) {
      this.clifford = clifford;
    } else {
      this.clifford = new Clifford(clifford);
    }
    this.num_qubits = this.clifford.num_qubits;
  }

  static fromLabel(label) {
    return new StabilizerState(Clifford.fromLabel(label));
  }

  static zero(numQubits) {
    return StabilizerState.fromLabel("0".repeat(numQubits));
  }

  // Probabilities of measuring each basis state
  probabilities(qargs = null) {
    // For a stabilizer state, the probabilities can be computed in polynomial
    // time via Gaussian elimination on the stabilizer table.
    // Simple approach: convert to statevector (works for small n only).
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

    // Compute all 2ⁿ stabilizer group elements as {x, z, fullPhase} where
    // the actual operator is fullPhase · X^x · Z^z.
    const groupElements = [];
    for (let subset = 0; subset < dim; subset++) {
      let x = new Array(n).fill(0);
      let z = new Array(n).fill(0);
      let fullPhase = new Complex(1, 0); // identity operator
      for (let k = 0; k < n; k++) {
        if ((subset >> k) & 1) {
          const gx = this.clifford.x[n + k];
          const gz = this.clifford.z[n + k];
          const gSign = this.clifford.signs[n + k];
          // Generator operator: (-1)^gSign · i^(gx·gz) · X^gx · Z^gz.
          let genFullPhase = gSign ? new Complex(-1, 0) : new Complex(1, 0);
          let genIPhase = 0;
          for (let q = 0; q < n; q++) {
            if (gx[q] && gz[q]) genIPhase = (genIPhase + 1) % 4;
          }
          const iPhaseFactor = [
            new Complex(1, 0), new Complex(0, 1),
            new Complex(-1, 0), new Complex(0, -1),
          ][genIPhase];
          genFullPhase = genFullPhase.mul(iPhaseFactor);

          // Multiply current · generator. Z and X anticommute, so
          // X^x Z^z X^gx Z^gz = (-1)^(z·gx) X^(x^gx) Z^(z^gz).
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

    // Find a reference basis state |i*> with ⟨ψ|i*⟩ ≠ 0.
    // ⟨i*|Π|i*⟩ = (1/2ⁿ) Σ_{P with x_P=0} fullPhase_P · (-1)^{popcount(i* AND z_P)}.
    // If this is positive, |i*> has non-zero overlap with |ψ>.
    let refState = -1;
    let refNormSq = 0;
    for (let iStar = 0; iStar < dim; iStar++) {
      let sumRe = 0;
      for (const g of groupElements) {
        // Only include P with x_P = 0 (so that X^x_P doesn't flip the bit).
        let xZero = true;
        for (let q = 0; q < n; q++) {
          if (g.x[q]) { xZero = false; break; }
        }
        if (!xZero) continue;
        // (-1)^{popcount(i* AND z_P)}
        let parity = 0;
        for (let q = 0; q < n; q++) {
          if (g.z[q] && ((iStar >> q) & 1)) parity ^= 1;
        }
        sumRe += parity ? -g.fullPhase.re : g.fullPhase.re;
      }
      // ⟨i*|Π|i*⟩ = sumRe / dim. (Should be real and non-negative.)
      if (sumRe > 1e-9) {
        refState = iStar;
        refNormSq = sumRe / dim;
        break;
      }
    }
    if (refState < 0) {
      throw new Error(
        "StabilizerState.to_statevector: could not find a reference state " +
        "with non-zero overlap (invalid stabilizer state)"
      );
    }
    const refNorm = Math.sqrt(refNormSq);

    // For each basis state |j>, compute the amplitude:
    //   amp_j = (1/2ⁿ) · Σ_{P with x_P = j XOR i*} fullPhase_P · (-1)^{popcount(i* AND z_P)} / refNorm
    const amps = new Array(dim);
    for (let j = 0; j < dim; j++) {
      const targetX = j ^ refState;
      let sumRe = 0, sumIm = 0;
      for (const g of groupElements) {
        // Include only group elements with x_P = targetX.
        let xMatches = true;
        for (let q = 0; q < n; q++) {
          if (((targetX >> q) & 1) !== g.x[q]) { xMatches = false; break; }
        }
        if (!xMatches) continue;
        // (-1)^{popcount(refState AND z_P)}
        let parity = 0;
        for (let q = 0; q < n; q++) {
          if (g.z[q] && ((refState >> q) & 1)) parity ^= 1;
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
    // Target Pauli's x and z parts.
    const targetX = pauli.x;
    const targetZ = pauli.z;

    // Enumerate all 2^n products of stabilizer generators.
    for (let subset = 0; subset < dim; subset++) {
      let x = new Array(n).fill(0);
      let z = new Array(n).fill(0);
      let fullPhase = new Complex(1, 0); // identity operator
      for (let k = 0; k < n; k++) {
        if ((subset >> k) & 1) {
          const gx = this.clifford.x[n + k];
          const gz = this.clifford.z[n + k];
          const gSign = this.clifford.signs[n + k];
          // Generator: (-1)^gSign · i^(gx·gz) · X^gx · Z^gz.
          let genFullPhase = gSign ? new Complex(-1, 0) : new Complex(1, 0);
          let genIPhase = 0;
          for (let q = 0; q < n; q++) {
            if (gx[q] && gz[q]) genIPhase = (genIPhase + 1) % 4;
          }
          const iPhaseFactor = [
            new Complex(1, 0), new Complex(0, 1),
            new Complex(-1, 0), new Complex(0, -1),
          ][genIPhase];
          genFullPhase = genFullPhase.mul(iPhaseFactor);
          // Anticommutation sign.
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
      // Check if (x, z) matches the target Pauli (ignoring phase, since
      // the phase IS what we're trying to recover).
      let matches = true;
      for (let q = 0; q < n; q++) {
        if (x[q] !== targetX[q] || z[q] !== targetZ[q]) { matches = false; break; }
      }
      if (matches) {
        // The product equals fullPhase · X^x · Z^z. But the target Pauli
        // is P(targetX, targetZ) = i^(targetX·targetZ) · X^targetX · Z^targetZ.
        // So the eigenvalue is fullPhase / i^(targetX·targetZ) =
        // fullPhase · i^(-targetX·targetZ).
        let targetIPhase = 0;
        for (let q = 0; q < n; q++) {
          if (targetX[q] && targetZ[q]) targetIPhase = (targetIPhase + 1) % 4;
        }
        // i^(-k) = i^(4-k) mod 4.
        const invIPhase = (4 - targetIPhase) % 4;
        const invIPhaseFactor = [
          new Complex(1, 0), new Complex(0, -1),
          new Complex(-1, 0), new Complex(0, 1),
        ][invIPhase];
        return fullPhase.mul(invIPhaseFactor);
      }
    }
    // No matching product: P is not in the stabilizer group.
    return Complex.ZERO;
  }

  equals(other, tol) {
    if (!(other instanceof StabilizerState)) return false;
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
}

// Convert a bit array (LSB first) to a number.
function _bitsToNum(bits, n) {
  let num = 0;
  for (let q = 0; q < n; q++) {
    if (bits[q]) num |= (1 << q);
  }
  return num;
}

// Embed a single-qubit Pauli matrix into the full n-qubit Hilbert space,
// acting on qubit `q`.
function _embedPauli(pauli2x2, n, q) {
  const dim = 1 << n;
  const out = ComplexMatrix.zeros(dim, dim);
  for (let i = 0; i < dim; i++) {
    for (let j = 0; j < dim; j++) {
      // The bit at position q in i and j selects the row/col of the 2x2 Pauli.
      const ib = (i >> q) & 1;
      const jb = (j >> q) & 1;
      // All other bits must match for the entry to be non-zero.
      let ok = true;
      for (let k = 0; k < n; k++) {
        if (k === q) continue;
        if (((i >> k) & 1) !== ((j >> k) & 1)) { ok = false; break; }
      }
      if (ok) out.set(i, j, pauli2x2.get(ib, jb));
    }
  }
  return out;
}

// Identify a tensor-product Pauli (with phase ±1) from a 2^n x 2^n matrix.
// Returns { x, z, sign } or null if the matrix isn't a Pauli.
function _identifyPauli(m, n) {
  const dim = 1 << n;
  // Find a non-zero entry to determine the global phase.
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
  if (!firstNonZero) return null; // zero matrix
  // The global phase is the value at (firstNonZero.i, firstNonZero.j).
  // For a Pauli, all non-zero entries have the same magnitude (= |phase|).
  const phase = firstNonZero.val;
  const mag = phase.abs();
  if (mag < 1e-9) return null;
  // Determine the X part: for each qubit, check if the matrix has non-zero
  // entries at (i, j) where i and j differ at bit q (X flips the bit).
  const x = new Array(n).fill(0);
  const z = new Array(n).fill(0);
  for (let q = 0; q < n; q++) {
    // Check if there's a non-zero entry where bit q differs between i and j.
    let hasFlip = false;
    let hasSame = false;
    for (let i = 0; i < dim; i++) {
      const j = i ^ (1 << q);
      const v = m.get(i, j);
      if (v.abs() > 1e-9) hasFlip = true;
      const v2 = m.get(i, i);
      if (v2.abs() > 1e-9) hasSame = true;
    }
    if (hasFlip) x[q] = 1;
    // Z is determined by the phase pattern: if entries where bit q = 1 have
    // a relative phase of -1 compared to bit q = 0, then z[q] = 1.
    // We pick a reference row and check the diagonal pattern.
    // For a Z-only qubit (x[q]=0), the diagonal entries at bit q = 1 are -phase.
    // For a Y qubit (x[q]=1, z[q]=1), the off-diagonal entries at (i, i XOR mask_q)
    // have phase ±i * phase depending on bit q of i.
  }
  // Now determine Z part more carefully by examining specific entries.
  // Pick i = 0 (all zeros). The column j = sum of mask bits where x[q] = 1.
  let j0 = 0;
  for (let q = 0; q < n; q++) if (x[q]) j0 |= (1 << q);
  const v00 = m.get(0, j0);
  if (v00.abs() < 1e-9) return null;
  // v00 should equal ±phase (or ±i*phase if there are Y qubits).
  // The relative phase = v00 / phase tells us the Y count parity.
  const rel = v00.div(phase);
  let yCount = 0;
  if (Math.abs(rel.re - 1) < 1e-9 && Math.abs(rel.im) < 1e-9) yCount = 0;
  else if (Math.abs(rel.re) < 1e-9 && Math.abs(rel.im - 1) < 1e-9) yCount = 1;
  else if (Math.abs(rel.re + 1) < 1e-9 && Math.abs(rel.im) < 1e-9) yCount = 2;
  else if (Math.abs(rel.re) < 1e-9 && Math.abs(rel.im + 1) < 1e-9) yCount = 3;
  else return null; // not a Pauli

  // Now figure out which x[q]=1 qubits are actually Y (z[q]=1).
  // For a Y qubit, the entry at (i, i XOR mask) has phase (-1)^{bit q of i} * i^yCount * phase.
  // We can determine z[q] for each qubit with x[q]=1 by checking the sign
  // pattern across different rows.
  // For X-only qubits (z[q]=0), all entries have the same phase (up to the
  // global i^yCount factor).
  // For Y qubits (z[q]=1), entries at rows with bit q = 1 have an extra -1
  // factor (combined with the i factor).
  // We'll determine z[q] by comparing the entry at (i, i XOR mask) for two
  // values of i that differ only at bit q.
  if (yCount === 0) {
    // No Y qubits; z[q] for x[q]=1 qubits is 0.
    // For x[q]=0 qubits, z[q] = 1 iff diagonal entries with bit q = 1 are -phase.
    for (let q = 0; q < n; q++) {
      if (x[q]) continue;
      // Check diagonal: entries (i, i) where bit q = 1.
      let hasNeg = false, hasPos = false;
      for (let i = 0; i < dim; i++) {
        if (((i >> q) & 1) !== 1) continue;
        const v = m.get(i, i);
        if (v.abs() < 1e-9) continue;
        const r = v.div(phase);
        if (r.re < -0.5) hasNeg = true;
        else if (r.re > 0.5) hasPos = true;
      }
      if (hasNeg && !hasPos) z[q] = 1;
      else if (hasNeg && hasPos) return null; // not a clean Pauli
    }
  } else {
    // yCount > 0: we have Y qubits. Identify them by checking the sign
    // pattern across rows.
    // For each qubit with x[q]=1, check if flipping bit q in the row index
    // changes the phase by -1 (Y) or leaves it unchanged (X).
    for (let q = 0; q < n; q++) {
      if (!x[q]) {
        // Check Z part for non-X qubits.
        let hasNeg = false, hasPos = false;
        for (let i = 0; i < dim; i++) {
          if (((i >> q) & 1) !== 1) continue;
          const v = m.get(i, i);
          if (v.abs() < 1e-9) continue;
          const r = v.div(phase);
          if (r.re < -0.5 || r.im < -0.5) hasNeg = true;
          else if (r.re > 0.5 || r.im > 0.5) hasPos = true;
        }
        if (hasNeg && !hasPos) z[q] = 1;
        continue;
      }
      // x[q] = 1: compare entry at (i, i XOR mask) for two rows differing at bit q.
      // We need to find two such rows where the entry is non-zero.
      let found = false;
      for (let i = 0; i < dim && !found; i++) {
        const j = i ^ j0; // column with all X bits flipped
        const i2 = i ^ (1 << q);
        const j2 = i2 ^ j0;
        const v1 = m.get(i, j);
        const v2 = m.get(i2, j2);
        if (v1.abs() < 1e-9 || v2.abs() < 1e-9) continue;
        const r = v2.div(v1);
        // If r ≈ -1, this is a Y qubit (z[q] = 1). If r ≈ 1, it's an X qubit.
        if (Math.abs(r.re + 1) < 1e-9 && Math.abs(r.im) < 1e-9) {
          z[q] = 1;
          found = true;
        } else if (Math.abs(r.re - 1) < 1e-9 && Math.abs(r.im) < 1e-9) {
          found = true; // X qubit, z[q] stays 0
        } else {
          return null; // not a Pauli
        }
      }
    }
  }
  // Determine the overall sign: the phase at (0, j0) divided by i^yCount
  // should be ±1.
  const expectedPhase = [new Complex(1, 0), new Complex(0, 1), new Complex(-1, 0), new Complex(0, -1)][yCount];
  const signPhase = phase.div(expectedPhase);
  if (Math.abs(signPhase.im) > 1e-9) return null;
  const sign = signPhase.re < 0 ? 1 : 0;
  // Verify: re-count Y qubits and check it matches yCount.
  let actualYCount = 0;
  for (let q = 0; q < n; q++) if (x[q] && z[q]) actualYCount++;
  if (actualYCount !== yCount && actualYCount % 4 !== yCount % 4) {
    // The yCount modulo 4 must match (since i^4 = 1).
    return null;
  }
  return { x, z, sign };
}
