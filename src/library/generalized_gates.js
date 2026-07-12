/**
 * generalized_gates.js - Parameterized single- and two-qubit gates.
 *
 * generalized_gates: RX, RY, RZ, P (Phase),
 * U1, U2, U3, U, RXX, RYY, RZZ, RZX.
 */

import { Complex, ComplexMatrix } from "./../math/linalg.js";
import { Gate } from "./../core/gate.js";
import { _registerParamBuilder } from "./../core/circuit.js";

const r = (re, im = 0) => new Complex(re, im);

function _val(p) {
  if (typeof p === "number") return p;
  if (p && typeof p.bind === "function") {
    throw new Error("Cannot build matrix with unbound parameter");
  }
  throw new TypeError(`Unsupported parameter type: ${typeof p}`);
}

export function makeRXGate(theta) {
  const g = new Gate("rx", 1, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix.fromRows([
      [r(ct),   r(0, -st)],
      [r(0, -st), r(ct)],
    ]);
  };
  return g;
}

export function makeRYGate(theta) {
  const g = new Gate("ry", 1, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix.fromRows([
      [r(ct),   r(st)],
      [r(-st),  r(ct)],
    ]);
  };
  return g;
}

export function makeRZGate(phi) {
  const g = new Gate("rz", 1, [phi]);
  g._matrixBuilder = (params) => {
    const p = _val(params[0]);
    return ComplexMatrix.fromRows([
      [r(Math.cos(p / 2), -Math.sin(p / 2)), r(0)],
      [r(0), r(Math.cos(p / 2), Math.sin(p / 2))],
    ]);
  };
  return g;
}

export function makePGate(theta) {
  const g = new Gate("p", 1, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    return ComplexMatrix.fromRows([
      [r(1), r(0)],
      [r(0), r(Math.cos(t), Math.sin(t))],
    ]);
  };
  return g;
}

export function makeU1Gate(lambda) {
  const g = new Gate("u1", 1, [lambda]);
  g._matrixBuilder = (params) => {
    const l = _val(params[0]);
    return ComplexMatrix.fromRows([
      [r(1), r(0)],
      [r(0), r(Math.cos(l), Math.sin(l))],
    ]);
  };
  return g;
}

export function makeU2Gate(phi, lambda) {
  const g = new Gate("u2", 1, [phi, lambda]);
  g._matrixBuilder = (params) => {
    const phi = _val(params[0]);
    const lam = _val(params[1]);
    const v = 1 / Math.SQRT2;
    return ComplexMatrix.fromRows([
      [r(v),                   r(-v * Math.cos(lam), -v * Math.sin(lam))],
      [r(v * Math.cos(phi), v * Math.sin(phi)), r(v * Math.cos(lam + phi), v * Math.sin(lam + phi))],
    ]);
  };
  return g;
}

export function makeU3Gate(theta, phi, lambda) {
  const g = new Gate("u3", 1, [theta, phi, lambda]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const p = _val(params[1]);
    const l = _val(params[2]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix.fromRows([
      [r(ct), r(-st * Math.cos(l), -st * Math.sin(l))],
      [r(st * Math.cos(p), st * Math.sin(p)), r(ct * Math.cos(l + p), ct * Math.sin(l + p))],
    ]);
  };
  return g;
}

export function makeUGate(theta, phi, lambda) {
  const g = makeU3Gate(theta, phi, lambda);
  g.name = "u";
  return g;
}

// Controlled-U gate with global phase.
// Qiskit's CU(θ, φ, λ, γ) = |0><0| ⊗ I + |1><1| ⊗ (e^{iγ} U3(θ, φ, λ)).
// The γ parameter is a global phase applied to the U3 block when the
// control is |1> (it does NOT affect the |0> branch). The previous
// implementation dropped γ entirely, which made any circuit using cu(…,
// γ≠0) produce wrong unitaries.
export function makeCUGate(theta, phi, lambda, gamma) {
  const g = new Gate("cu", 2, [theta, phi, lambda, gamma]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const p = _val(params[1]);
    const l = _val(params[2]);
    const gm = _val(params[3]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    // e^{iγ} * U3(θ, φ, λ) block:
    //   [[ e^{iγ} ct,                              -e^{iγ} e^{iλ} st ],
    //    [  e^{iγ} e^{iφ} st,   e^{iγ} e^{i(φ+λ)} ct ]]
    const gRe = Math.cos(gm);
    const gIm = Math.sin(gm);
    // Helper: (gRe + i gIm) * (a + i b) = (gRe*a - gIm*b) + i (gRe*b + gIm*a)
    const cmul = (a, b) => [gRe * a - gIm * b, gRe * b + gIm * a];
    // Top-left: e^{iγ} * ct
    const m00 = cmul(ct, 0);
    // Top-right: -e^{iγ} * (cos(l) + i sin(l)) * st
    const e_il = [Math.cos(l), Math.sin(l)];
    const tr_prod = cmul(e_il[0] * st, e_il[1] * st);
    const m01 = [-tr_prod[0], -tr_prod[1]];
    // Bottom-left: e^{iγ} * (cos(p) + i sin(p)) * st
    const e_ip = [Math.cos(p), Math.sin(p)];
    const bl_prod = cmul(e_ip[0] * st, e_ip[1] * st);
    const m10 = bl_prod;
    // Bottom-right: e^{iγ} * (cos(p+l) + i sin(p+l)) * ct
    const e_ipl = [Math.cos(p + l), Math.sin(p + l)];
    const br_prod = cmul(e_ipl[0] * ct, e_ipl[1] * ct);
    const m11 = br_prod;
    // Full 4x4 controlled-U: identity on |0> branch, e^{iγ} U3 on |1> branch.
    // Qubit ordering (control, target) with control as MSB:
    //   |00> -> |00>, |01> -> |01>,
    //   |10> -> e^{iγ} U3 |0> on target = e^{iγ} (ct |0> + e^{iφ} st |1>) on target
    //   |11> -> e^{iγ} U3 |1> on target = e^{iγ} (-e^{iλ} st |0> + e^{i(φ+λ)} ct |1>) on target
    return ComplexMatrix.fromRows([
      [r(1),              r(0),       r(0),       r(0)],
      [r(0),              r(1),       r(0),       r(0)],
      [r(0),              r(0),       r(m00[0], m00[1]),  r(m01[0], m01[1])],
      [r(0),              r(0),       r(m10[0], m10[1]),  r(m11[0], m11[1])],
    ]);
  };
  return g;
}

// Two-qubit parameterized gates
export function makeRXXGate(theta) {
  const g = new Gate("rxx", 2, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix.fromRows([
      [r(ct),     r(0),       r(0),       r(0, -st)],
      [r(0),      r(ct),      r(0, -st),  r(0)],
      [r(0),      r(0, -st),  r(ct),      r(0)],
      [r(0, -st), r(0),       r(0),       r(ct)],
    ]);
  };
  return g;
}

export function makeRYYGate(theta) {
  const g = new Gate("ryy", 2, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix.fromRows([
      [r(ct),    r(0),      r(0),      r(0, st)],
      [r(0),     r(ct),     r(0, -st), r(0)],
      [r(0),     r(0, -st), r(ct),     r(0)],
      [r(0, st), r(0),      r(0),      r(ct)],
    ]);
  };
  return g;
}

export function makeRZZGate(theta) {
  const g = new Gate("rzz", 2, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix.fromRows([
      [r(ct, -st), r(0),       r(0),       r(0)],
      [r(0),       r(ct, st),  r(0),       r(0)],
      [r(0),       r(0),       r(ct, st),  r(0)],
      [r(0),       r(0),       r(0),       r(ct, -st)],
    ]);
  };
  return g;
}

export function makeRZXGate(theta) {
  const g = new Gate("rzx", 2, [theta]);
  g._matrixBuilder = (params) => {
    const t = _val(params[0]);
    const ct = Math.cos(t / 2);
    const st = Math.sin(t / 2);
    return ComplexMatrix.fromRows([
      [r(ct),     r(0, -st), r(0),      r(0)],
      [r(0, -st), r(ct),     r(0),      r(0)],
      [r(0),      r(0),      r(ct),     r(0, st)],
      [r(0),      r(0),      r(0, st),  r(ct)],
    ]);
  };
  return g;
}

// Register all parameterized builders
_registerParamBuilder("RX", makeRXGate);
_registerParamBuilder("RY", makeRYGate);
_registerParamBuilder("RZ", makeRZGate);
_registerParamBuilder("P",  makePGate);
_registerParamBuilder("U1", makeU1Gate);
_registerParamBuilder("U2", makeU2Gate);
_registerParamBuilder("U3", makeU3Gate);
_registerParamBuilder("U",  makeUGate);
_registerParamBuilder("CU", makeCUGate);
_registerParamBuilder("RXX", makeRXXGate);
_registerParamBuilder("RYY", makeRYYGate);
_registerParamBuilder("RZZ", makeRZZGate);
_registerParamBuilder("RZX", makeRZXGate);
