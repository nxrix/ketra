import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("Quantum channels", () => {
  test("Kraus from unitary", () => {
    const X = K.PAULI.X;
    const k = new K.Kraus([X]);
    expect(k.numQubits).toBe(1);
  });

  test("Kraus toSuperop", () => {
    const X = K.PAULI.X;
    const k = new K.Kraus([X]);
    const s = k.toSuperop();
    expect(s.rows).toBe(4);
  });

  test("Kraus toChoi", () => {
    const X = K.PAULI.X;
    const k = new K.Kraus([X]);
    const choi = k.toChoi();
    expect(choi.rows).toBe(4);
  });

  test("Kraus toPtm", () => {
    const X = K.PAULI.X;
    const k = new K.Kraus([X]);
    const ptm = k.toPtm();
    expect(ptm.rows).toBe(4);
  });

  test("SuperOp fromMatrix", () => {
    const X = K.PAULI.X;
    const k = new K.Kraus([X]);
    const s = SuperOp.fromMatrix(k.toSuperop(), 1);
    expect(s.numQubits).toBe(1);
  });

  test("stateFidelity", () => {
    const sv0 = K.Statevector.fromLabel("0");
    const sv1 = K.Statevector.fromLabel("1");
    expect(approxEq(K.stateFidelity(sv0, sv0), 1, 1e-9)).toBe(true);
    expect(approxEq(K.stateFidelity(sv0, sv1), 0, 1e-9)).toBe(true);
  });

  test("processFidelity", () => {
    const id = new K.Kraus([K.ComplexMatrix.identity(2)]);
    expect(K.processFidelity(id, K.ComplexMatrix.identity(2))).toBeGreaterThan(0);
  });

  test("averageGateFidelity", () => {
    const id = new K.Kraus([K.ComplexMatrix.identity(2)]);
    const agf = K.averageGateFidelity(id, K.ComplexMatrix.identity(2));
    expect(agf).toBeGreaterThan(0.9);
  });

  test("diamondNorm", () => {
    const id = new K.Kraus([K.ComplexMatrix.identity(2)]);
    const dn = K.diamondNorm(id);
    expect(dn).toBeGreaterThan(0);
  });
});

import { SuperOp } from "../src/quantum_info/channels.js";

describe("Noise models", () => {
  test("depolarizingError", () => {
    const e = K.depolarizingError(0.1, 1);
    expect(e).toBeDefined();
  });

  test("bitFlipError", () => {
    const e = K.bitFlipError(0.1);
    expect(e).toBeDefined();
  });

  test("phaseFlipError", () => {
    expect(K.phaseFlipError(0.1)).toBeDefined();
  });

  test("amplitudeDampingError", () => {
    expect(K.amplitudeDampingError(0.1)).toBeDefined();
  });

  test("phaseDampingError", () => {
    expect(K.phaseDampingError(0.1)).toBeDefined();
  });

  test("NoiseModel", () => {
    const nm = new K.NoiseModel();
    nm.addAllQubitQuantumError(K.bitFlipError(0.1), ["x"]);
    expect(nm).toBeDefined();
  });

  test("ReadoutError", () => {
    const re = new K.ReadoutError([[0.9, 0.1], [0.1, 0.9]]);
    expect(re).toBeDefined();
  });
});
