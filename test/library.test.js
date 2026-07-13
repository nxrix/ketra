import * as K from "../src/index.js";

describe("Library circuits", () => {
  test("bellState", () => {
    const qc = K.circuits.bellState();
    expect(qc.numQubits).toBe(2);
  });

  test("ghzState", () => {
    const qc = K.circuits.ghzState(3);
    expect(qc.numQubits).toBe(3);
  });

  test("qft", () => {
    const qc = K.circuits.qft(3);
    expect(qc.numQubits).toBe(3);
  });

  test("realAmplitudes has parameters", () => {
    const qc = K.circuits.realAmplitudes(3, 2);
    expect(qc.parameters.size).toBeGreaterThan(0);
  });

  test("efficientSU2 has parameters", () => {
    const qc = K.circuits.efficientSU2(3, 2);
    expect(qc.parameters.size).toBeGreaterThan(0);
  });

  test("twoLocal has unique parameter names", () => {
    const qc = K.circuits.twoLocal(3, ["ry"], ["cx"], 2);
    const paramNames = Array.from(qc.parameters).map(p => p.name);
    const unique = new Set(paramNames);
    expect(unique.size).toBe(paramNames.length);
  });

  test("graphState", () => {
    const qc = K.circuits.graphState(3, [[0, 1], [1, 2]]);
    expect(qc.numQubits).toBe(3);
  });

  test("pauliEvolution", () => {
    const qc = K.circuits.pauliEvolution("ZZ", 1.0);
    expect(qc.numQubits).toBe(2);
  });

  test("hiddenLinearFunction", () => {
    const qc = K.circuits.hiddenLinearFunction(2, [[0, 1]]);
    expect(qc.numQubits).toBe(2);
  });
});

describe("Extra gates", () => {
  test("UnitaryGate", () => {
    expect(K.UnitaryGate).toBeDefined();
    const g = new K.UnitaryGate(K.PAULI.X);
    expect(g.numQubits).toBe(1);
  });

  test("DiagonalGate", () => {
    expect(K.DiagonalGate).toBeDefined();
    const g = new K.DiagonalGate([0, Math.PI]);
    expect(g.numQubits).toBe(1);
  });

  test("PermutationGate", () => {
    expect(K.PermutationGate).toBeDefined();
    const g = new K.PermutationGate([1, 0]);
    expect(g.numQubits).toBe(2);
  });

  test("HamiltonianGate", () => {
    expect(K.HamiltonianGate).toBeDefined();
    const g = new K.HamiltonianGate(K.PAULI.X, Math.PI / 2);
    expect(g.numQubits).toBe(1);
  });

  test("Initialize", () => {
    expect(K.Initialize).toBeDefined();
    const i = new K.Initialize([1, 0]);
    expect(i.numQubits).toBe(1);
  });

  test("MCPhaseGate, MCRXGate, MCRYGate, MCRZGate, MCMTGate", () => {
    expect(K.MCPhaseGate).toBeDefined();
    expect(K.MCRXGate).toBeDefined();
    expect(K.MCRYGate).toBeDefined();
    expect(K.MCRZGate).toBeDefined();
    expect(K.MCMTGate).toBeDefined();
  });
});

describe("Boolean logic gates", () => {
  test("ANDGate", () => {
    const g = new K.ANDGate(2);
    expect(g.numQubits).toBe(3);
  });

  test("ORGate", () => {
    expect(new K.ORGate(2).numQubits).toBe(3);
  });

  test("XORGate", () => {
    expect(new K.XORGate(2).numQubits).toBe(3);
  });

  test("NANDGate, NORGate, XNORGate", () => {
    expect(K.NANDGate).toBeDefined();
    expect(K.NORGate).toBeDefined();
    expect(K.XNORGate).toBeDefined();
  });

  test("AND(1,1) -> output=1", () => {
    const andGate = new K.ANDGate(2);
    const qc = new K.QuantumCircuit(3);
    qc.x(0); qc.x(1);
    qc.append(andGate, [0, 1, 2]);
    const sv = K.Statevector.fromCircuit(qc);
    expect(sv.probabilities()[7]).toBeCloseTo(1, 9);
  });
});

describe("Arithmetic circuits", () => {
  test("linearPauliRotations", () => {
    const qc = K.arithmetic.linearPauliRotations(2, [0.5, 0.3], 0.1);
    expect(qc.numQubits).toBe(3);
  });

  test("weightedAdder", () => {
    const qc = K.arithmetic.weightedAdder(2, [1, 2]);
    expect(qc.numQubits).toBeGreaterThanOrEqual(2);
  });

  test("draperQFTAdder is unitary", () => {
    const qc = K.arithmetic.draperQFTAdder(2);
    expect(qc.numQubits).toBe(4);
    const op = K.Operator.fromCircuit(qc);
    expect(op.isUnitary(1e-9)).toBe(true);
  });

  test("quadraticForm, integerComparator, etc.", () => {
    expect(K.arithmetic.quadraticForm(2)).toBeDefined();
    expect(K.arithmetic.integerComparator(2, 1)).toBeDefined();
    expect(K.arithmetic.cdkmRippleCarryAdder(2)).toBeDefined();
    expect(K.arithmetic.hrsCumulativeMultiplier(2)).toBeDefined();
    expect(K.arithmetic.functionalPauliRotations(2, [0], [0.5], [0])).toBeDefined();
  });
});
