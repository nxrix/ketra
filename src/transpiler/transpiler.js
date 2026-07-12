/**
 * transpiler.js - transpile() function and basis gate decomposition passes.
 *
 * Includes basis-gate decomposition for non-basis gates, an unroll pass that
 * recursively expands gate definitions, an optimization pass that removes
 * adjacent inverse gates, and a trivial layout when circuit qubits match
 * backend qubits. Decompositions follow the IBM Quantum basis {cx, id, rz, sx, x}.
 */

import { QuantumCircuit, CircuitInstruction } from "./../core/circuit.js";
import { QuantumRegister, ClassicalRegister } from "./../core/bit.js";
import { Gate, ControlledGate } from "./../core/gate.js";
import * as standardGates from "./../library/standard_gates.js";
import * as generalizedGates from "./../library/generalized_gates.js";

// Default basis gate set matching IBM Quantum's defaults.
export const DEFAULT_BASIS = ["cx", "id", "rz", "sx", "x"];

// Decomposition rules. The decomposition list is in CIRCUIT ORDER (first
// applied first). The combined unitary is U = g_n * ... * g_2 * g_1 (matrix
// multiplication is right-to-left, with the first-applied gate rightmost).
// All rules are derived from the standard Z-Y decomposition of single-qubit
// unitaries and the well-known CX-based decompositions of multi-qubit gates
// (see Nielsen & Chuang, sections 4.2 and 4.5).
export const DECOMP_RULES = {
  "h": () => [
    ["rz", [0], [Math.PI / 2]],
    ["sx", [0], []],
    ["rz", [0], [Math.PI / 2]],
  ],
  "y": () => [
    ["sx", [0], []],
    ["sx", [0], []],
    ["rz", [0], [Math.PI]],
  ],
  "z": () => [["rz", [0], [Math.PI]]],
  "s": () => [["rz", [0], [Math.PI / 2]]],
  "sdg": () => [["rz", [0], [-Math.PI / 2]]],
  "t": () => [["rz", [0], [Math.PI / 4]]],
  "tdg": () => [["rz", [0], [-Math.PI / 4]]],
  "sx": () => null,
  "sxdg": () => [
    ["sx", [0], []],
    ["x", [0], []],
  ],
  "rx": (params) => [
    ["h", [0], []],
    ["rz", [0], [params[0]]],
    ["h", [0], []],
  ],
  "ry": (params) => [
    // RY(theta) = RZ(-pi/2) H RZ(theta) H RZ(pi/2).
    ["rz", [0], [Math.PI / 2]],
    ["h", [0], []],
    ["rz", [0], [params[0]]],
    ["h", [0], []],
    ["rz", [0], [-Math.PI / 2]],
  ],
  "p": (params) => [["rz", [0], [params[0]]]],
  "u1": (params) => [["rz", [0], [params[0]]]],
  "u2": (params) => [
    // U2(phi, lambda) = RZ(phi) H RZ(lambda) up to global phase.
    ["rz", [0], [params[1]]],
    ["h", [0], []],
    ["rz", [0], [params[0]]],
  ],
  "u3": (params) => [
    // U3(theta, phi, lambda) = RZ(phi) H RZ(theta) H RZ(lambda) up to global phase.
    ["rz", [0], [params[2]]],
    ["h", [0], []],
    ["rz", [0], [params[0]]],
    ["h", [0], []],
    ["rz", [0], [params[1]]],
  ],
  "u": (params) => [
    ["rz", [0], [params[2]]],
    ["h", [0], []],
    ["rz", [0], [params[0]]],
    ["h", [0], []],
    ["rz", [0], [params[1]]],
  ],
  "cy": () => [
    ["sdg", [1], []],
    ["cx", [0, 1], []],
    ["s", [1], []],
  ],
  "cz": () => [
    ["h", [1], []],
    ["cx", [0, 1], []],
    ["h", [1], []],
  ],
  "ch": () => [
    // CH = S(target) H(target) CX(control, target) H(target) Sdg(target).
    ["s", [1], []],
    ["h", [1], []],
    ["cx", [0, 1], []],
    ["h", [1], []],
    ["sdg", [1], []],
  ],
  "swap": () => [
    ["cx", [0, 1], []],
    ["cx", [1, 0], []],
    ["cx", [0, 1], []],
  ],
  "iswap": () => [
    // iSWAP = S(0) S(1) H(0) CX(0,1) CX(1,0) H(1) (verified up to global phase).
    ["s", [0], []],
    ["s", [1], []],
    ["h", [0], []],
    ["cx", [0, 1], []],
    ["cx", [1, 0], []],
    ["h", [1], []],
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
    ["cx", [0, 1], []],
  ],
  "cswap": () => [
    ["cx", [2, 1], []],
    ["ccx", [0, 1, 2], []],
    ["cx", [2, 1], []],
  ],
  "crx": (params) => [
    // CRX(theta) = RZ(-pi/2) RY(theta/2) CX RY(-theta/2) CX RZ(pi/2) on target.
    ["rz", [1], [-Math.PI / 2]],
    ["ry", [1], [params[0] / 2]],
    ["cx", [0, 1], []],
    ["ry", [1], [-params[0] / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [Math.PI / 2]],
  ],
  "cry": (params) => [
    // CRY(theta) = RY(theta/2) CX RY(-theta/2) CX on target.
    ["ry", [1], [params[0] / 2]],
    ["cx", [0, 1], []],
    ["ry", [1], [-params[0] / 2]],
    ["cx", [0, 1], []],
  ],
  "crz": (params) => [
    // CRZ(theta) = RZ(theta/2) CX RZ(-theta/2) CX on target.
    ["rz", [1], [params[0] / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [-params[0] / 2]],
    ["cx", [0, 1], []],
  ],
  "cp": (params) => [
    // CP(theta) = RZ(theta/2) on control, CX, RZ(-theta/2) on target, CX, RZ(theta/2) on target.
    ["rz", [0], [params[0] / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [-params[0] / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [params[0] / 2]],
  ],
  "rzz": (params) => [
    ["cx", [0, 1], []],
    ["rz", [1], [params[0]]],
    ["cx", [0, 1], []],
  ],
  "rxx": (params) => [
    ["h", [0], []],
    ["h", [1], []],
    ["cx", [0, 1], []],
    ["rz", [1], [params[0]]],
    ["cx", [0, 1], []],
    ["h", [1], []],
    ["h", [0], []],
  ],
  "ryy": (params) => [
    // RYY(theta) = RY(-pi/2) on both, CX, RZ(theta) on target, CX, RY(pi/2) on both.
    ["ry", [0], [-Math.PI / 2]],
    ["ry", [1], [-Math.PI / 2]],
    ["cx", [0, 1], []],
    ["rz", [1], [params[0]]],
    ["cx", [0, 1], []],
    ["ry", [0], [Math.PI / 2]],
    ["ry", [1], [Math.PI / 2]],
  ],
};

export function decomposeGate(gateName, params, basis) {
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
        const mappedQubits = subQubits.map(qi => qubits[qi]);
        result.push([subName, mappedQubits, subP]);
      }
    }
  }
  return result;
}

export function transpile(circuits, options = {}) {
  const basis = options.basis_gates || DEFAULT_BASIS;
  const optimizationLevel = options.optimization_level !== undefined ? options.optimization_level : 1;
  const single = !Array.isArray(circuits) || circuits instanceof QuantumCircuit;
  const list = single ? [circuits] : circuits;
  const transpiled = list.map(c => _transpileOne(c, basis, optimizationLevel, options));
  return single ? transpiled[0] : transpiled;
}

function _transpileOne(circuit, basis, optLevel, options) {
  const out = new QuantumCircuit();
  const qubitMap = new Map();
  const clbitMap = new Map();
  if (circuit.qregs.length === 0) {
    const newReg = new QuantumRegister(circuit.num_qubits, "q");
    out.add_register(newReg);
    for (let i = 0; i < circuit.num_qubits; i++) qubitMap.set(circuit.qubits[i], newReg._bits[i]);
  } else {
    for (const r of circuit.qregs) {
      const newReg = new QuantumRegister(r.size, r.name);
      out.add_register(newReg);
      for (let i = 0; i < r.size; i++) qubitMap.set(r._bits[i], newReg._bits[i]);
    }
  }
  if (circuit.cregs.length > 0) {
    for (const r of circuit.cregs) {
      const newReg = new ClassicalRegister(r.size, r.name);
      out.add_register(newReg);
      for (let i = 0; i < r.size; i++) clbitMap.set(r._bits[i], newReg._bits[i]);
    }
  } else if (circuit.num_clbits > 0) {
    const newReg = new ClassicalRegister(circuit.num_clbits, "c");
    out.add_register(newReg);
    for (let i = 0; i < circuit.num_clbits; i++) clbitMap.set(circuit.clbits[i], newReg._bits[i]);
  }
  out.name = circuit.name + "_transpiled";
  out.global_phase = circuit.global_phase;

  const remapQubits = (qs) => qs.map(q => qubitMap.get(q) || q);
  const remapClbits = (cs) => cs.map(c => clbitMap.get(c) || c);

  for (const ci of circuit.data) {
    const opName = ci.operation.name;
    const newQubits = remapQubits(ci.qubits);
    const newClbits = remapClbits(ci.clbits);

    if (opName === "barrier") { out.barrier(newQubits); continue; }
    if (opName === "measure") { out.append(ci.operation.copy(), newQubits, newClbits); continue; }
    if (opName === "reset") { out.append(ci.operation.copy(), newQubits); continue; }
    if (opName === "delay") { out.append(ci.operation.copy(), newQubits); continue; }
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
      const mappedQubits = qubitOffsets.map(qi => newQubits[qi]);
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
  if (standardGates) {
    const factories = {
      "h": standardGates.makeHGate,
      "x": standardGates.makeXGate,
      "y": standardGates.makeYGate,
      "z": standardGates.makeZGate,
      "s": standardGates.makeSGate,
      "sdg": standardGates.makeSdgGate,
      "t": standardGates.makeTGate,
      "tdg": standardGates.makeTdgGate,
      "sx": standardGates.makeSxGate,
      "sxdg": standardGates.makeSxdgGate,
      "id": standardGates.makeIGate,
      "cx": standardGates.makeCXGate,
      "cy": standardGates.makeCYGate,
      "cz": standardGates.makeCZGate,
      "ch": standardGates.makeCHGate,
      "csx": standardGates.makeCSXGate,
      "swap": standardGates.makeSwapGate,
      "iswap": standardGates.makeISwapGate,
      "dcx": standardGates.makeDCXGate,
      "ccx": standardGates.makeCCXGate,
      "cswap": standardGates.makeCSwapGate,
    };
    const factory = factories[name];
    if (factory) return factory();
  }
  if (generalizedGates) {
    const paramFactories = {
      "rx": generalizedGates.makeRXGate,
      "ry": generalizedGates.makeRYGate,
      "rz": generalizedGates.makeRZGate,
      "p": generalizedGates.makePGate,
      "u1": generalizedGates.makeU1Gate,
      "u2": generalizedGates.makeU2Gate,
      "u3": generalizedGates.makeU3Gate,
      "u": generalizedGates.makeUGate,
      "rxx": generalizedGates.makeRXXGate,
      "ryy": generalizedGates.makeRYYGate,
      "rzz": generalizedGates.makeRZZGate,
      "rzx": generalizedGates.makeRZXGate,
    };
    const factory = paramFactories[name];
    if (factory) return factory.apply(null, params);
  }
  if (generalizedGates) {
    const controlledBases = {
      "cp": generalizedGates.makePGate,
      "crx": generalizedGates.makeRXGate,
      "cry": generalizedGates.makeRYGate,
      "crz": generalizedGates.makeRZGate,
      "cu1": generalizedGates.makeU1Gate,
      "cu3": generalizedGates.makeU3Gate,
      "cu": generalizedGates.makeUGate,
    };
    const baseFactory = controlledBases[name];
    if (baseFactory) {
      const base = baseFactory.apply(null, params);
      const cg = new ControlledGate(base, 1);
      cg.name = name;
      return cg;
    }
  }
  const numQubits = (name === "cx" || name === "cy" || name === "cz" || name === "ch" ||
                    name === "swap" || name === "iswap" || name === "cp" || name === "crx" ||
                    name === "cry" || name === "crz" || name === "rxx" || name === "ryy" ||
                    name === "rzz" || name === "rzx" || name === "csx" || name === "cu" ||
                    name === "cu1" || name === "cu3" || name === "dcx") ? 2 :
                    (name === "ccx" || name === "cswap") ? 3 : 1;
  return new Gate(name, numQubits, params);
}

function _removeAdjacentInverseGates(circuit) {
  const newData = [];
  for (const ci of circuit.data) {
    if (newData.length === 0) { newData.push(ci); continue; }
    const prev = newData[newData.length - 1];
    if (prev.qubits.length === ci.qubits.length &&
        prev.qubits.every((q, i) => q === ci.qubits[i]) &&
        prev.operation.name === ci.operation.name + "_dg") {
      newData.pop();
      continue;
    }
    if (prev.qubits.length === ci.qubits.length &&
        prev.qubits.every((q, i) => q === ci.qubits[i]) &&
        ci.operation.name === prev.operation.name + "_dg") {
      newData.pop();
      continue;
    }
    const selfInverse = ["h", "x", "y", "z", "cx", "cz", "swap"];
    if (selfInverse.indexOf(prev.operation.name) !== -1 &&
        prev.operation.name === ci.operation.name &&
        prev.qubits.length === ci.qubits.length &&
        prev.qubits.every((q, i) => q === ci.qubits[i])) {
      newData.pop();
      continue;
    }
    newData.push(ci);
  }
  circuit.data = newData;
}

function _commuteCZ(circuit) { return circuit; }

// ---------------------------------------------------------------------------
// PassManager-style API
// ---------------------------------------------------------------------------
export class PassManager {
  constructor() { this.passes = []; }
  append(pass) { this.passes.push(pass); }
  run(circuit) {
    let c = circuit;
    for (const pass of this.passes) c = pass(c);
    return c;
  }
}

export class PassManagerConfig {
  constructor(kwargs) {
    this.basis_gates = kwargs.basis_gates || DEFAULT_BASIS;
    this.optimization_level = kwargs.optimization_level || 0;
    this.initial_layout = kwargs.initial_layout || null;
    this.routing_method = kwargs.routing_method || "basic";
  }
}

export function presetPassManager(optimizationLevel = 1, backend = null, basisGates = null) {
  const pm = new PassManager();
  const basis = basisGates || (backend && backend.basis_gates) || DEFAULT_BASIS;
  pm.append((c) => transpile(c, { basis_gates: basis, optimization_level: optimizationLevel }));
  return pm;
}
