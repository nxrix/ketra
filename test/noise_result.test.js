import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("Boolean logic gates", () => {
  test("ANDGate matrix", () => {
    const g = new K.ANDGate(2);
    const m = g.toMatrix();
    expect(m.rows).toBe(8);
    // Verify the gate is unitary (permutation matrix).
    expect(m.isUnitary(1e-9)).toBe(true);
  });

  test("ORGate matrix", () => {
    const g = new K.ORGate(2);
    const m = g.toMatrix();
    expect(m.rows).toBe(8);
    expect(m.isUnitary(1e-9)).toBe(true);
  });

  test("XORGate matrix", () => {
    const g = new K.XORGate(2);
    const m = g.toMatrix();
    expect(m.rows).toBe(8);
    expect(m.isUnitary(1e-9)).toBe(true);
  });

  test("NANDGate matrix", () => {
    const g = new K.NANDGate(2);
    const m = g.toMatrix();
    expect(m.rows).toBe(8);
    expect(m.isUnitary(1e-9)).toBe(true);
  });

  test("NORGate matrix", () => {
    const g = new K.NORGate(2);
    const m = g.toMatrix();
    expect(m.rows).toBe(8);
    expect(m.isUnitary(1e-9)).toBe(true);
  });

  test("XNORGate matrix", () => {
    const g = new K.XNORGate(2);
    const m = g.toMatrix();
    expect(m.rows).toBe(8);
    expect(m.isUnitary(1e-9)).toBe(true);
  });

  test("ANDGate copy", () => {
    const g = new K.ANDGate(2);
    const c = g.copy();
    expect(c.numQubits).toBe(3);
  });

  test("ANDGate with 3 inputs", () => {
    const g = new K.ANDGate(3);
    expect(g.numQubits).toBe(4);
  });

  test("mcxVChain", () => {
    const qc = new K.QuantumCircuit(5);
    K.mcxVChain(qc, [0, 1, 2], 3, [4]);
    expect(qc.data.length).toBeGreaterThan(0);
  });

  test("mcxVChain with 2 controls", () => {
    const qc = new K.QuantumCircuit(3);
    K.mcxVChain(qc, [0, 1], 2, []);
    expect(qc.data.length).toBe(1);
  });

  test("mcxRecursive with 3 controls", () => {
    const qc = new K.QuantumCircuit(5);
    K.mcxRecursive(qc, [0, 1, 2], 3, 4);
    expect(qc.data.length).toBeGreaterThan(0);
  });

  test("mcxRecursive with 2 controls", () => {
    const qc = new K.QuantumCircuit(3);
    K.mcxRecursive(qc, [0, 1], 2);
    expect(qc.data.length).toBe(1);
  });

  test("mcxRecursive with 1 control", () => {
    const qc = new K.QuantumCircuit(2);
    K.mcxRecursive(qc, [0], 1);
    expect(qc.data.length).toBe(1);
  });

  test("mcxNoAncilla with 1 control", () => {
    const qc = new K.QuantumCircuit(2);
    K.mcxNoAncilla(qc, [0], 1);
    expect(qc.data.length).toBe(1);
  });

  test("mcxNoAncilla with 2 controls", () => {
    const qc = new K.QuantumCircuit(3);
    K.mcxNoAncilla(qc, [0, 1], 2);
    expect(qc.data.length).toBe(1);
  });

  test("mcxNoAncilla with 3 controls", () => {
    const qc = new K.QuantumCircuit(4);
    K.mcxNoAncilla(qc, [0, 1, 2], 3);
    expect(qc.data.length).toBe(1);
  });
});

describe("Noise models", () => {
  test("QuantumError construction", () => {
    const e = new K.QuantumError([K.PAULI.X], 1);
    expect(e).toBeDefined();
  });

  test("depolarizingError", () => {
    const e = K.depolarizingError(0.1, 1);
    expect(e).toBeDefined();
    const e2 = K.depolarizingError(0.1, 2);
    expect(e2).toBeDefined();
  });

  test("bitFlipError", () => {
    const e = K.bitFlipError(0.1);
    expect(e).toBeDefined();
  });

  test("phaseFlipError", () => {
    const e = K.phaseFlipError(0.1);
    expect(e).toBeDefined();
  });

  test("amplitudeDampingError", () => {
    const e = K.amplitudeDampingError(0.1);
    expect(e).toBeDefined();
  });

  test("phaseDampingError", () => {
    const e = K.phaseDampingError(0.1);
    expect(e).toBeDefined();
  });

  test("pauliXError, pauliYError, pauliZError", () => {
    expect(K.pauliXError(0.1)).toBeDefined();
    expect(K.pauliYError(0.1)).toBeDefined();
    expect(K.pauliZError(0.1)).toBeDefined();
  });

  test("resetError", () => {
    expect(K.resetError(0.1)).toBeDefined();
  });

  test("krausError", () => {
    expect(K.krausError([K.PAULI.X])).toBeDefined();
  });

  test("mixedUnitaryError", () => {
    expect(K.mixedUnitaryError([[K.PAULI.X, 0.5], [K.PAULI.Z, 0.5]])).toBeDefined();
  });

  test("combineErrors", () => {
    const e1 = K.bitFlipError(0.1);
    const e2 = K.phaseFlipError(0.1);
    const combined = K.combineErrors(e1, e2);
    expect(combined).toBeDefined();
  });

  test("NoiseModel addQuantumError with qubits", () => {
    const nm = new K.NoiseModel();
    const e = K.bitFlipError(0.1);
    nm.addQuantumError(e, ["x"], [0]);
    expect(nm).toBeDefined();
  });

  test("NoiseModel basisGates", () => {
    const nm = new K.NoiseModel();
    nm.addAllQubitQuantumError(K.bitFlipError(0.1), ["x"]);
    expect(nm.basisGates).toContain("x");
  });

  test("NoiseModel getQuantumError", () => {
    const nm = new K.NoiseModel();
    const e = K.bitFlipError(0.1);
    nm.addAllQubitQuantumError(e, ["x"]);
    const retrieved = nm.getQuantumError("x", [0]);
    expect(retrieved).toBeDefined();
  });

  test("NoiseModel addReadoutError", () => {
    const nm = new K.NoiseModel();
    const re = new K.ReadoutError([[0.9, 0.1], [0.1, 0.9]]);
    nm.addReadoutError(re, [0]);
    expect(nm).toBeDefined();
  });

  test("NoiseModel addAllQubitReadoutError", () => {
    const nm = new K.NoiseModel();
    const re = new K.ReadoutError([[0.9, 0.1], [0.1, 0.9]]);
    nm.addAllQubitReadoutError(re);
    expect(nm).toBeDefined();
  });

  test("ReadoutError", () => {
    const re = new K.ReadoutError([[0.9, 0.1], [0.1, 0.9]]);
    expect(re).toBeDefined();
  });
});

describe("Result and Counts", () => {
  test("Counts basic operations", () => {
    const c = new K.Counts({ "00": 512, "11": 512 });
    expect(c.get("00")).toBe(512);
    expect(c.get("01")).toBe(0);
    expect(c.get("11")).toBe(512);
  });

  test("Counts items", () => {
    const c = new K.Counts({ "00": 512, "11": 512 });
    const items = Array.from(c.items());
    expect(items.length).toBe(2);
  });

  test("Counts from simulation", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.x(0);
    qc.measure(0, 0);
    const result = K.simulate(qc, 100);
    const counts = result.getCounts();
    expect(counts.get("1")).toBe(100);
  });

  test("Result properties", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.x(0);
    qc.measure(0, 0);
    const result = K.simulate(qc, 100);
    expect(result.backendName).toBeDefined();
    expect(result.success).toBe(true);
  });
});
