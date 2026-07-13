import { Complex, ComplexMatrix } from "./../math/linalg.js";

export class Instruction {
  constructor(name, numQubits, numClbits, params) {
    this.name = name;
    this.numQubits = numQubits;
    this.numClbits = numClbits;
    this.params = params || [];
    this.label = null;
    this.condition = null;
    this._definition = null;
    this._matrixBuilder = null;
  }

  get numParams() { return this.params.length; }

  get definition() {
    if (this._definition) return this._definition(this);
    return null;
  }

  set definition(fn) {
    this._definition = fn;
  }

  toMatrix() {
    if (this._matrixBuilder) {
      return this._matrixBuilder(this.params);
    }
    throw new Error(`Instruction ${this.name} does not define a matrix`);
  }

  hasMatrix() {
    return this._matrixBuilder !== null;
  }

  copy() {
    const c = new Instruction(this.name, this.numQubits, this.numClbits, this.params.slice());
    c.label = this.label;
    c.condition = this.condition;
    c._definition = this._definition;
    c._matrixBuilder = this._matrixBuilder;
    return c;
  }

  control(numCtrlQubits = 1, label = null, ctrlState = null) {
    return new ControlledGate(this, numCtrlQubits, label, ctrlState);
  }

  inverse() {
    // Self-inverse gates: keep the same name and matrix
    const selfInverse = ["h", "x", "y", "z", "cx", "cy", "cz", "swap", "ccx",
                         "cswap", "id", "ch"];
    // Gates with named inverses (swap name, use dagger matrix)
    const inverseMap = {
      "s": "sdg", "sdg": "s",
      "t": "tdg", "tdg": "t",
      "sx": "sxdg", "sxdg": "sx",
    };
    // Parameterized gates whose inverse negates the parameter
    const negateParam = ["rx", "ry", "rz", "p", "u1", "rxx", "ryy", "rzz", "rzx",
                         "cp", "crx", "cry", "crz", "cu1"];
    // Gates that need dagger (not self-inverse, not negatable, not named-inverse)
    const daggerGates = ["iswap", "csx", "dcx", "u2", "u3", "u", "cu3", "cu"];

    const name = this.name.toLowerCase();

    if (selfInverse.includes(name)) {
      return this.copy();
    }

    if (inverseMap[name]) {
      const inv = this.copy();
      inv.name = inverseMap[name];
      if (this._matrixBuilder) {
        const originalBuilder = this._matrixBuilder;
        inv._matrixBuilder = (params) => originalBuilder(params).dagger();
      }
      return inv;
    }

    if (negateParam.includes(name)) {
      const inv = this.copy();
      inv.params = this.params.map(p => {
        if (typeof p === "number") return -p;
        if (p && typeof p.negate === "function") return p.negate();
        if (p && typeof p.mul === "function") return p.mul(-1);
        return p;
      });
      return inv;
    }

    if (daggerGates.includes(name)) {
      const inv = this.copy();
      if (this._matrixBuilder) {
        const originalBuilder = this._matrixBuilder;
        inv._matrixBuilder = (params) => originalBuilder(params).dagger();
      }
      return inv;
    }

    // ControlledGate: invert the base gate and wrap in a new ControlledGate
    if (this instanceof ControlledGate) {
      const invBase = this.baseGate.inverse();
      const inv = new ControlledGate(invBase, this.numCtrlQubits, this.label, this.ctrlState);
      inv.name = this.name;
      return inv;
    }

    // Default: use _dg suffix with dagger matrix
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
}

function _p2s(p) {
  if (p && typeof p.toString === "function") return p.toString();
  return String(p);
}

export class Gate extends Instruction {
  constructor(name, numQubits, params) {
    super(name, numQubits, 0, params);
  }
}

export class ControlledGate extends Gate {
  constructor(baseGate, numCtrlQubits, label, ctrlState) {
    super(
      `c${numCtrlQubits}_${baseGate.name}`,
      baseGate.numQubits + numCtrlQubits,
      baseGate.params.slice()
    );
    this.baseGate = baseGate;
    this.numCtrlQubits = numCtrlQubits;
    this.label = label;
    // Default ctrlState: all-ones (i.e., gate fires when every control is |1>).
    this.ctrlState = (ctrlState !== null && ctrlState !== undefined)
      ? ctrlState
      : ((1 << numCtrlQubits) - 1);
    this._matrixBuilder = this._controlledMatrixBuilder.bind(this);
  }

  _controlledMatrixBuilder(params) {
    // Build the controlled-unitary matrix in the gate's local Hilbert space.
    // Qubit ordering convention (matching _embedGate in statevector.js,
    // density_matrix.js, and the simulators): gate-bit 0 is the LSB.
    // The first `numCtrlQubits` gate-bits are the control qubits, so
    // gate.qubits[0] is the LSB control. The remaining bits are the base
    // gate's qubits, with base qubit 0 at bit position `numCtrlQubits`.
    //
    // Index layout: index = control_bits + (base_bits << numCtrlQubits).
    // So control state c occupies indices {c, c + numCtrl, c + 2*numCtrl, ...}.
    const baseMatrix = this.baseGate.toMatrix();
    const dim = baseMatrix.rows;          // 2^(base.numQubits)
    const numCtrl = 1 << this.numCtrlQubits;
    const total = dim * numCtrl;
    const out = ComplexMatrix.zeros(total, total);
    const ctrlOn = this.ctrlState;
    for (let c = 0; c < numCtrl; c++) {
      if (c === ctrlOn) {
        // Place the base gate at the subspace where control bits == ctrlOn.
        for (let bi = 0; bi < dim; bi++) {
          for (let bj = 0; bj < dim; bj++) {
            const row = c + (bi << this.numCtrlQubits);
            const col = c + (bj << this.numCtrlQubits);
            out.set(row, col, baseMatrix.get(bi, bj));
          }
        }
      } else {
        // Identity: out[i][i] = 1 for all i in this control block.
        for (let bi = 0; bi < dim; bi++) {
          const idx = c + (bi << this.numCtrlQubits);
          out.set(idx, idx, Complex.ONE);
        }
      }
    }
    return out;
  }

  copy() {
    const inv = new ControlledGate(this.baseGate.copy(), this.numCtrlQubits, this.label, this.ctrlState);
    inv.name = this.name;
    inv.params = this.params.slice();
    return inv;
  }

  inverse() {
    const inv = new ControlledGate(this.baseGate.inverse(), this.numCtrlQubits, this.label, this.ctrlState);
    inv.name = this.name + "_dg";
    return inv;
  }
}
