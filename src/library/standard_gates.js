import { Complex, ComplexMatrix } from "./../math/linalg.js";
import { Gate } from "./../core/gate.js";
import { _registerStd } from "./../core/circuit.js";

const PI = Math.PI;
const SQRT2 = Math.SQRT2;
const r = (re, im = 0) => new Complex(re, im);

// Single-qubit gates
export function makeHGate() {
  const g = new Gate("h", 1, []);
  const v = 1 / SQRT2;
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(v), r(v)],
    [r(v), r(-v)],
  ]);
  return g;
}

export function makeXGate() {
  const g = new Gate("x", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(0), r(1)],
    [r(1), r(0)],
  ]);
  return g;
}

export function makeYGate() {
  const g = new Gate("y", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(0), r(0, -1)],
    [r(0, 1), r(0)],
  ]);
  return g;
}

export function makeZGate() {
  const g = new Gate("z", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0)],
    [r(0), r(-1)],
  ]);
  return g;
}

export function makeSGate() {
  const g = new Gate("s", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0)],
    [r(0), r(0, 1)],
  ]);
  return g;
}

export function makeSdgGate() {
  const g = new Gate("sdg", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0)],
    [r(0), r(0, -1)],
  ]);
  return g;
}

export function makeTGate() {
  const g = new Gate("t", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0)],
    [r(0), r(Math.cos(PI / 4), Math.sin(PI / 4))],
  ]);
  return g;
}

export function makeTdgGate() {
  const g = new Gate("tdg", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0)],
    [r(0), r(Math.cos(PI / 4), -Math.sin(PI / 4))],
  ]);
  return g;
}

export function makeSxGate() {
  const g = new Gate("sx", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(0.5, 0.5), r(0.5, -0.5)],
    [r(0.5, -0.5), r(0.5, 0.5)],
  ]);
  return g;
}

export function makeSxdgGate() {
  const g = new Gate("sxdg", 1, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(0.5, -0.5), r(0.5, 0.5)],
    [r(0.5, 0.5), r(0.5, -0.5)],
  ]);
  return g;
}

export function makeIGate() {
  const g = new Gate("id", 1, []);
  g._matrixBuilder = () => ComplexMatrix.identity(2);
  return g;
}

// Two-qubit gates
export function makeCXGate() {
  const g = new Gate("cx", 2, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0), r(0), r(1)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(1), r(0), r(0)],
  ]);
  return g;
}

export function makeCYGate() {
  const g = new Gate("cy", 2, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0),       r(0), r(0)],
    [r(0), r(0),       r(0), r(0, -1)],
    [r(0), r(0),       r(1), r(0)],
    [r(0), r(0, 1),    r(0), r(0)],
  ]);
  return g;
}

export function makeCZGate() {
  const g = new Gate("cz", 2, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(1), r(0), r(0)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(0), r(0), r(-1)],
  ]);
  return g;
}

export function makeCHGate() {
  const g = new Gate("ch", 2, []);
  const v = 1 / SQRT2;
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0),  r(0),  r(0)],
    [r(0), r(v),  r(0),  r(v)],
    [r(0), r(0),  r(1),  r(0)],
    [r(0), r(v),  r(0),  r(-v)],
  ]);
  return g;
}

export function makeCSXGate() {
  const g = new Gate("csx", 2, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0),         r(0),          r(0)],
    [r(0), r(0.5, 0.5),  r(0),          r(0.5, -0.5)],
    [r(0), r(0),         r(1),          r(0)],
    [r(0), r(0.5, -0.5), r(0),          r(0.5, 0.5)],
  ]);
  return g;
}

export function makeSwapGate() {
  const g = new Gate("swap", 2, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0), r(1), r(0)],
    [r(0), r(1), r(0), r(0)],
    [r(0), r(0), r(0), r(1)],
  ]);
  return g;
}

export function makeISwapGate() {
  const g = new Gate("iswap", 2, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0),      r(0),      r(0)],
    [r(0), r(0),      r(0, 1),   r(0)],
    [r(0), r(0, 1),   r(0),      r(0)],
    [r(0), r(0),      r(0),      r(1)],
  ]);
  return g;
}

export function makeDCXGate() {
  const g = new Gate("dcx", 2, []);
  g._matrixBuilder = () => ComplexMatrix.fromRows([
    [r(1), r(0), r(0), r(0)],
    [r(0), r(0), r(0), r(1)],
    [r(0), r(1), r(0), r(0)],
    [r(0), r(0), r(1), r(0)],
  ]);
  return g;
}

// Three-qubit gates
export function makeCCXGate() {
  const g = new Gate("ccx", 3, []);
  g._matrixBuilder = () => {
    const m = ComplexMatrix.identity(8);
    m.set(3, 3, r(0));
    m.set(7, 7, r(0));
    m.set(3, 7, r(1));
    m.set(7, 3, r(1));
    return m;
  };
  return g;
}

export function makeCSwapGate() {
  const g = new Gate("cswap", 3, []);
  g._matrixBuilder = () => {
    const m = ComplexMatrix.identity(8);
    m.set(3, 3, r(0));
    m.set(5, 5, r(0));
    m.set(3, 5, r(1));
    m.set(5, 3, r(1));
    return m;
  };
  return g;
}

// Registration
export const HGate = makeHGate();
export const XGate = makeXGate();
export const YGate = makeYGate();
export const ZGate = makeZGate();
export const SGate = makeSGate();
export const SdgGate = makeSdgGate();
export const TGate = makeTGate();
export const TdgGate = makeTdgGate();
export const SxGate = makeSxGate();
export const SxdgGate = makeSxdgGate();
export const IGate = makeIGate();

export const CXGate = makeCXGate();
export const CYGate = makeCYGate();
export const CZGate = makeCZGate();
export const CHGate = makeCHGate();
export const CSXGate = makeCSXGate();
export const SwapGate = makeSwapGate();
export const ISwapGate = makeISwapGate();
export const DCXGate = makeDCXGate();

export const CCXGate = makeCCXGate();
export const CSwapGate = makeCSwapGate();

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
