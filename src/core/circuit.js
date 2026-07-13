import { QuantumRegister, ClassicalRegister, Qubit, Clbit } from "./bit.js";
import { Instruction, Gate, ControlledGate } from "./gate.js";
import { Parameter } from "./parameter.js";

// Registry for extra gate classes (UnitaryGate, Initialize, etc.). The
// library/extra_gates.js module calls _registerExtraGateClass() at import
// time to populate this. This avoids the circular-import problem (extra_gates.js
// imports from circuit.js via _registerParamBuilder, so circuit.js can't
// statically import extra_gates.js).
const _extraGateClasses = {};
export function _registerExtraGateClass(name, cls) {
  _extraGateClasses[name] = cls;
}
function _extraGateClass(name) {
  if (!_extraGateClasses[name]) {
    throw new Error(
      `Extra gate class ${name} not registered. Import '@nxrix/ketra/library/extra_gates.js' ` +
      `before using circuit.unitary() / .initialize() / .diagonal() / .permute() / .hamiltonian().`
    );
  }
  return _extraGateClasses[name];
}

// Registry for classes that circuit.js needs but that would create a
// circular import if imported statically. Each is set by its module when
// it loads.
let _OperatorClass = null;
let _StatevectorClass = null;
let _QASMExporterClass = null;
export function _registerOperatorClass(cls) { _OperatorClass = cls; }
export function _registerStatevectorClass(cls) { _StatevectorClass = cls; }
export function _registerQASMExporterClass(cls) { _QASMExporterClass = cls; }

export class CircuitInstruction {
  constructor(operation, qubits, clbits) {
    this.operation = operation;
    this.qubits = qubits || [];
    this.clbits = clbits || [];
  }

  copy() {
    return new CircuitInstruction(this.operation.copy(), this.qubits.slice(), this.clbits.slice());
  }
}

// Standard gate cache (populated by library/standard_gates.js on import)
const _stdCache = {};
const _paramBuilders = {};

export function _registerStd(name, gate) {
  _stdCache[name] = gate;
}

export function _registerParamBuilder(name, builder) {
  _paramBuilders[name] = builder;
}

export function _getStdGate(name) {
  if (!_stdCache[name]) {
    throw new Error(`Standard gate ${name} not registered. Did you import library/standard_gates?`);
  }
  return _stdCache[name].copy();
}

export class QuantumCircuit {
  constructor(...regs) {
    this.qregs = [];
    this.cregs = [];
    this.qubits = [];
    this.clbits = [];
    this.data = [];
    this._qubit_index = new Map();
    this._clbit_index = new Map();
    this._global_phase = 0;
    this.name = "circuit";
    this.metadata = null;

    // Qiskit-style shorthand: QuantumCircuit(n_qubits[, n_clbits])
    if (regs.length > 0 && regs.every(r => typeof r === "number")) {
      const nQ = regs[0];
      const nC = regs.length > 1 ? regs[1] : 0;
      if (nQ > 0) this.addRegister(new QuantumRegister(nQ, "q"));
      if (nC > 0) this.addRegister(new ClassicalRegister(nC, "c"));
      return;
    }

    for (const r of regs) {
      if (r instanceof QuantumRegister) this.addRegister(r);
      else if (r instanceof ClassicalRegister) this.addRegister(r);
      else if (typeof r === "number") this.addRegister(new QuantumRegister(r));
      else if (Array.isArray(r)) {
        for (const sub of r) {
          if (sub instanceof QuantumRegister) this.addRegister(sub);
          else if (sub instanceof ClassicalRegister) this.addRegister(sub);
          else throw new TypeError("Invalid register in list");
        }
      }
      else throw new TypeError("Invalid register argument to QuantumCircuit");
    }
  }

  addRegister(register) {
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
      throw new TypeError("addRegister requires QuantumRegister or ClassicalRegister");
    }
    return register;
  }

  get numQubits() { return this.qubits.length; }
  get numClbits() { return this.clbits.length; }
  get nqubits() { return this.qubits.length; }
  get nclbits() { return this.clbits.length; }

  qubitIndices(qubit) {
    if (qubit instanceof Qubit) {
      const i = this._qubit_index.get(qubit);
      if (i === undefined) throw new Error("Qubit not in circuit");
      return i;
    }
    throw new TypeError("qubitIndices expects a Qubit");
  }

  clbitIndices(clbit) {
    if (clbit instanceof Clbit) {
      const i = this._clbit_index.get(clbit);
      if (i === undefined) throw new Error("Clbit not in circuit");
      return i;
    }
    throw new TypeError("clbitIndices expects a Clbit");
  }

  _resolveQubits(qargs) {
    if (!Array.isArray(qargs)) qargs = [qargs];
    return qargs.map(q => {
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
    return cargs.map(c => {
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
    if (qubits.length !== instruction.numQubits) {
      throw new Error(`Qubit count mismatch: gate ${instruction.name} expects ${instruction.numQubits} but got ${qubits.length}`);
    }
    if (clbits.length !== instruction.numClbits) {
      throw new Error(`Clbit count mismatch: gate ${instruction.name} expects ${instruction.numClbits} but got ${clbits.length}`);
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
    const newData = other.data.map(ci => {
      return new CircuitInstruction(
        ci.operation.copy(),
        ci.qubits.map(q => qMap.get(q)),
        ci.clbits.map(c => cMap.get(c)),
      );
    });
    if (front) target.data = newData.concat(target.data);
    else target.data = target.data.concat(newData);
    return target;
  }

  _composeQubitMap(target, qubits) {
    const map = new Map();
    if (qubits == null) {
      // Default: map this circuit's qubits to target's qubits by index
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
    const map = new Map();
    if (clbits == null) {
      // Default: map this circuit's clbits to target's clbits by index
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
    const c = new QuantumCircuit();
    const qubitMap = new Map();
    const clbitMap = new Map();
    for (const r of this.qregs) {
      const newReg = new QuantumRegister(r.size, r.name);
      c.addRegister(newReg);
      for (let i = 0; i < r.size; i++) qubitMap.set(r._bits[i], newReg._bits[i]);
    }
    for (const r of this.cregs) {
      const newReg = new ClassicalRegister(r.size, r.name);
      c.addRegister(newReg);
      for (let i = 0; i < r.size; i++) clbitMap.set(r._bits[i], newReg._bits[i]);
    }
    c.data = this.data.map(ci => {
      const newOp = ci.operation.copy();
      const newQubits = ci.qubits.map(q => qubitMap.get(q) || q);
      const newClbits = ci.clbits.map(cl => clbitMap.get(cl) || cl);
      return new CircuitInstruction(newOp, newQubits, newClbits);
    });
    c.globalPhase = this.globalPhase;
    c.name = this.name;
    c.metadata = this.metadata;
    return c;
  }

  // Standard gates
  h(q)   { return this.append(_getStdGate("H"),   q); }
  x(q)   { return this.append(_getStdGate("X"),   q); }
  y(q)   { return this.append(_getStdGate("Y"),   q); }
  z(q)   { return this.append(_getStdGate("Z"),   q); }
  s(q)   { return this.append(_getStdGate("S"),   q); }
  sdg(q) { return this.append(_getStdGate("SDG"), q); }
  t(q)   { return this.append(_getStdGate("T"),   q); }
  tdg(q) { return this.append(_getStdGate("TDG"), q); }
  sx(q)  { return this.append(_getStdGate("SX"),  q); }
  sxdg(q){ return this.append(_getStdGate("SXDG"),q); }
  id(q)  { return this.append(_getStdGate("I"),   q); }

  cx(control, target)    { return this.append(_getStdGate("CX"),    [control, target]); }
  cy(control, target)    { return this.append(_getStdGate("CY"),    [control, target]); }
  cz(control, target)    { return this.append(_getStdGate("CZ"),    [control, target]); }
  ch(control, target)    { return this.append(_getStdGate("CH"),    [control, target]); }
  csx(control, target)   { return this.append(_getStdGate("CSX"),   [control, target]); }
  swap(q1, q2)           { return this.append(_getStdGate("SWAP"),  [q1, q2]); }
  iswap(q1, q2)          { return this.append(_getStdGate("ISWAP"), [q1, q2]); }
  dcx(q1, q2)            { return this.append(_getStdGate("DCX"),   [q1, q2]); }

  ccx(c1, c2, target)    { return this.append(_getStdGate("CCX"),    [c1, c2, target]); }
  cswap(c, t1, t2)       { return this.append(_getStdGate("CSWAP"),  [c, t1, t2]); }

  rx(theta, q)           { return this.append(_mkParamGate("RX", 1, [theta], 1), q); }
  ry(theta, q)           { return this.append(_mkParamGate("RY", 1, [theta], 1), q); }
  rz(phi, q)             { return this.append(_mkParamGate("RZ", 1, [phi],   1), q); }
  rxx(theta, q1, q2)     { return this.append(_mkParamGate("RXX", 2, [theta], 1), [q1, q2]); }
  ryy(theta, q1, q2)     { return this.append(_mkParamGate("RYY", 2, [theta], 1), [q1, q2]); }
  rzz(theta, q1, q2)     { return this.append(_mkParamGate("RZZ", 2, [theta], 1), [q1, q2]); }
  rzx(theta, q1, q2)     { return this.append(_mkParamGate("RZX", 2, [theta], 1), [q1, q2]); }
  p(theta, q)            { return this.append(_mkParamGate("P",  1, [theta], 1), q); }
  u(theta, phi, lam, q)  { return this.append(_mkParamGate("U", 1, [theta, phi, lam], 1), q); }
  u1(lam, q)             { return this.append(_mkParamGate("U1", 1, [lam], 1), q); }
  u2(phi, lam, q)        { return this.append(_mkParamGate("U2", 1, [phi, lam], 1), q); }
  u3(theta, phi, lam, q) { return this.append(_mkParamGate("U3", 1, [theta, phi, lam], 1), q); }

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
    // Use the dedicated 4-parameter CU gate builder so that the global
    // implementation built a ControlledGate around a 3-parameter U gate,
    // silently dropping γ.
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
    const UnitaryGate = _extraGateClass("UnitaryGate");
    const g = new UnitaryGate(matrix, label);
    const qs = Array.isArray(qubits) ? qubits : [qubits];
    return this.append(g, qs);
  }

  // Append a state-preparation instruction. Amplitudes can be an array of
  // Complex/numbers or a state label like "01" or "+-".
  initialize(amplitudes, qubits) {
    const Initialize = _extraGateClass("Initialize");
    // qubits is required: it must be an array or a single qubit index/Qubit.
    if (qubits == null) {
      // Default: apply to the first n qubits, where n = log2(len(amplitudes)).
      let n;
      if (typeof amplitudes === "string") n = amplitudes.length;
      else if (Array.isArray(amplitudes)) n = Math.log2(amplitudes.length);
      else throw new TypeError("initialize requires an array or label string");
      qubits = [];
      for (let i = 0; i < n; i++) qubits.push(i);
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
    // The Initialize instruction acts on `numQubits` qubits determined by
    // the amplitudes array length.
    const n = Math.log2(amps.length);
    if (!Number.isInteger(n)) {
      throw new Error(`initialize: amplitudes length ${amps.length} is not 2^n`);
    }
    const instr = new Initialize(amps, n);
    return this.append(instr, qs);
  }

  // Append a diagonal unitary, given as a list of phases (radians).
  diagonal(diag, qubits) {
    const DiagonalGate = _extraGateClass("DiagonalGate");
    const g = new DiagonalGate(diag);
    const qs = Array.isArray(qubits) ? qubits : [qubits];
    return this.append(g, qs);
  }

  // Append a qubit-permutation gate. `pattern[i]` is the qubit index that
  // ends up at position i.
  permute(pattern, qubits) {
    const PermutationGate = _extraGateClass("PermutationGate");
    const g = new PermutationGate(pattern);
    const qs = Array.isArray(qubits) ? qubits : qubits !== undefined ? [qubits] : null;
    if (qs === null) {
      // Default: apply to the first n qubits.
      const qs2 = [];
      for (let i = 0; i < pattern.length; i++) qs2.push(i);
      return this.append(g, qs2);
    }
    return this.append(g, qs);
  }

  // Append a Hamiltonian evolution gate e^{-i*t*H}.
  hamiltonian(operator, time, qubits, label = null) {
    const HamiltonianGate = _extraGateClass("HamiltonianGate");
    const g = new HamiltonianGate(operator, time, label);
    const qs = Array.isArray(qubits) ? qubits : [qubits];
    return this.append(g, qs);
  }

  // Convert this circuit to a reusable Gate (with the same qubit count).
  toGate(label = null) {
    const g = new Gate(this.name || "circuit", this.numQubits, []);
    g.label = label;
    g._matrixBuilder = () => {
      // Build the unitary by composing each instruction's matrix.
      // Use the Operator registry (set when quantum_info/operator.js loads).
      if (!_OperatorClass) {
        throw new Error("Operator class not registered. Import quantum_info/operator.js first.");
      }
      return _OperatorClass.fromCircuit(this)._data;
    };
    // Carry a definition so that decomposition / transpilation can work.
    g._definition = () => {
      return this.data.map(ci => [ci.operation, ci.qubits, ci.clbits]);
    };
    return g;
  }

  // Convert this circuit to an Instruction (can have clbits).
  toInstruction(label = null) {
    const instr = new Instruction(this.name || "circuit", this.numQubits, this.numClbits, []);
    instr.label = label;
    instr._matrixBuilder = () => {
      if (!_OperatorClass) {
        throw new Error("Operator class not registered. Import quantum_info/operator.js first.");
      }
      return _OperatorClass.fromCircuit(this)._data;
    };
    instr._definition = () => {
      return this.data.map(ci => [ci.operation, ci.qubits, ci.clbits]);
    };
    return instr;
  }

  // Repeat the circuit `n` times on the same qubits (in series).
  repeat(n) {
    if (n < 0) throw new Error("repeat: n must be non-negative");
    const c = new QuantumCircuit();
    // Copy registers
    for (const r of this.qregs) c.addRegister(new QuantumRegister(r.size, r.name));
    for (const r of this.cregs) c.addRegister(new ClassicalRegister(r.size, r.name));
    c.globalPhase = this.globalPhase;
    c.name = this.name;
    // Copy qubit/clbit mappings: this.qubits[i] -> c.qubits[i]
    const qMap = new Map();
    const cMap = new Map();
    for (let i = 0; i < this.qubits.length; i++) qMap.set(this.qubits[i], c.qubits[i]);
    for (let i = 0; i < this.clbits.length; i++) cMap.set(this.clbits[i], c.clbits[i]);
    for (let rep = 0; rep < n; rep++) {
      for (const ci of this.data) {
        c.append(
          ci.operation.copy(),
          ci.qubits.map(q => qMap.get(q)),
          ci.clbits.map(cl => cMap.get(cl)),
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
      // Negative power = repeat the inverse.
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
  set globalPhase(value) { this._global_phase = value; }
  get globalPhase() { return this._global_phase || 0; }

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

  measureAll(inplace = true) {
    const target = inplace ? this : this.copy();
    if (target.clbits.length === 0) {
      target.addRegister(new ClassicalRegister(target.qubits.length));
    }
    for (let i = 0; i < target.qubits.length; i++) {
      target.measure(target.qubits[i], target.clbits[i]);
    }
    return target;
  }

  measureActive(inplace = true) {
    // Add measurements for every qubit that participates in at least one
    // gate (excluding barriers, which don't count as "active" for this
    // purpose). All qubits of multi-qubit gates are inspected.
    const target = inplace ? this : this.copy();
    const active = new Set();
    for (const ci of target.data) {
      if (ci.operation.name === "barrier") continue;
      for (const q of ci.qubits) active.add(q);
    }
    const activeBits = target.qubits.filter(q => active.has(q));
    if (activeBits.length === 0) return target;
    const startClbit = target.clbits.length;
    target.addRegister(new ClassicalRegister(activeBits.length));
    activeBits.forEach((q, i) => target.measure(q, target.clbits[startClbit + i]));
    return target;
  }

  reset(qubit) {
    const qubits = this._resolveQubits(Array.isArray(qubit) ? qubit : [qubit]);
    const out = [];
    for (const q of qubits) {
      const r = new Instruction("reset", 1, 0, []);
      out.push(this.append(r, q));
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

  ifTest(condition, true_predicate, false_predicate = null) {
    const instr = new Instruction("if_else", 0, 0, [condition, true_predicate, false_predicate]);
    this.data.push(new CircuitInstruction(instr, [], []));
    return instr;
  }

  whileLoop(condition, body, qubits, clbits) {
    const instr = new Instruction("whileLoop", 0, 0, [condition, body]);
    this.data.push(new CircuitInstruction(instr, [], []));
    return instr;
  }

  // Parameter binding
  bindParameters(values) {
    const c = this.copy();
    c.data = c.data.map(ci => {
      const newOp = ci.operation.copy();
      if (newOp.params && newOp.params.length) {
        newOp.params = newOp.params.map(p => {
          if (p && typeof p.bind === "function") return p.bind(values);
          return p;
        });
      }
      return new CircuitInstruction(newOp, ci.qubits.slice(), ci.clbits.slice());
    });
    return c;
  }

  get parameters() {
    // Collect every (leaf) Parameter referenced by any instruction in the
    // circuit. Both Parameter and ParameterExpression expose `parameters`
    // as a getter (a Set), NOT as a method. We accept either a getter
    // (Set) or a method (returns Set), and we only add leaf Parameters
    // (not ParameterExpression intermediates) to match qiskit's convention:
    // `circuit.parameters` is the set of unbound Parameters that need
    // values supplied to `bindParameters`.
    const set = new Set();
    const seen = [];
    for (const ci of this.data) {
      for (const p of ci.operation.params || []) {
        if (p == null) continue;
        // Plain numeric parameter (already bound).
        if (typeof p === "number") continue;
        // Leaf Parameter: add directly.
        if (p instanceof Parameter) {
          _addToSet(set, seen, p);
          continue;
        }
        // ParameterExpression: collect its sub-parameters (which are the
        // actual unbound Parameters).
        let subParams = null;
        if (p && p.parameters !== undefined) {
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

  countOps() {
    const counts = {};
    for (const ci of this.data) {
      counts[ci.operation.name] = (counts[ci.operation.name] || 0) + 1;
    }
    return counts;
  }

  numNonlocalGates() {
    let n = 0;
    for (const ci of this.data) {
      if (ci.operation.numQubits > 1 && ci.operation.name !== "barrier") n++;
    }
    return n;
  }

  numTensorFactors() {
    const parent = new Array(this.qubits.length).fill(0).map((_, i) => i);
    const find = (x) => parent[x] === x ? x : (parent[x] = find(parent[x]));
    const union = (a, b) => { parent[find(a)] = find(b); };
    for (const ci of this.data) {
      if (ci.operation.numQubits < 2) continue;
      const idx = ci.qubits.map(q => this._qubit_index.get(q));
      for (let i = 1; i < idx.length; i++) union(idx[0], idx[i]);
    }
    const roots = new Set();
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
  reverseOps() {
    const c = this.copy();
    c.data.reverse();
    return c;
  }

  inverse() {
    const c = this.copy();
    c.data = c.data
      .filter(ci => !["measure", "reset", "barrier"].includes(ci.operation.name))
      .reverse()
      .map(ci => new CircuitInstruction(ci.operation.inverse(), ci.qubits.slice(), ci.clbits.slice()));
    return c;
  }

  removeFinalMeasurements(inplace = false) {
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
    // The visualization module hooks into QuantumCircuit via a registration
    // function. If not registered, use the inline ASCII renderer.
    if (_drawHook) {
      return _drawHook(this, output, kwargs);
    }
    return _fallbackDraw(this);
  }

  toString() {
    return _fallbackDraw(this);
  }
}

// Drawing hook - set by visualization.js when imported
let _drawHook = null;
export function _setDrawHook(fn) { _drawHook = fn; }

function _mkParamGate(name, numQubits, params, numClbits) {
  if (!_paramBuilders[name]) {
    throw new Error(`Parameterized gate ${name} not registered. Did you import library/generalized_gates?`);
  }
  return _paramBuilders[name].apply(null, params);
}

function _fallbackDraw(circuit) {
  const lines = [];
  lines.push(`QuantumCircuit "${circuit.name}" (qubits=${circuit.numQubits}, clbits=${circuit.numClbits})`);
  for (const ci of circuit.data) {
    const qs = ci.qubits.map(q => q.toString()).join(",");
    const cs = ci.clbits.map(c => c.toString()).join(",");
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

// Decompose one level: replace any instruction that has a definition with
// the instructions in its definition. Returns a new QuantumCircuit.
function _decomposeOnce(circuit) {
  const c = circuit.copy();
  const newData = [];
  for (const ci of c.data) {
    const op = ci.operation;
    // Check if the operation has a definition (a function that returns a
    // list of [instruction, qubits, clbits] tuples).
    const defFn = op._definition;
    if (typeof defFn === "function") {
      try {
        const def = defFn(op);
        if (def && def.length > 0) {
          // Build a map from definition-qubits (which are the SUB-circuit's
          // Qubit objects) to the parent instruction's qargs, using the
          // order of first appearance in the definition.
          const defUniqueQubits = [];
          const seen = new Set();
          for (const [_, sqs, scs] of def) {
            for (const sq of sqs) {
              if (!seen.has(sq)) { seen.add(sq); defUniqueQubits.push(sq); }
            }
          }
          const defUniqueClbits = [];
          const seenC = new Set();
          for (const [_, sqs, scs] of def) {
            for (const sc of scs) {
              if (!seenC.has(sc)) { seenC.add(sc); defUniqueClbits.push(sc); }
            }
          }
          for (const [subOp, subQubits, subClbits] of def) {
            newData.push(new CircuitInstruction(
              subOp.copy(),
              subQubits.map(sq => {
                const idx = defUniqueQubits.indexOf(sq);
                return idx >= 0 && idx < ci.qubits.length ? ci.qubits[idx] : sq;
              }),
              subClbits.map(sc => {
                const idx = defUniqueClbits.indexOf(sc);
                return idx >= 0 && idx < ci.clbits.length ? ci.clbits[idx] : sc;
              }),
            ));
          }
          continue;
        }
      } catch (e) {
        // Definition threw — fall through to keep the original instruction.
      }
    }
    newData.push(ci);
  }
  c.data = newData;
  return c;
}
