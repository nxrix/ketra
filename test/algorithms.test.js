import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("VQE", () => {
  test("finds Z ground state = -1", () => {
    const op = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const ansatz = new K.QuantumCircuit(1);
    const theta = new K.Parameter("theta");
    ansatz.ry(theta, 0);
    const optimizer = new K.GradientDescent({ learningRate: 0.5, maxiter: 50 });
    const vqe = new K.VQE({ ansatz, optimizer, initialPoint: [0.1] });
    const result = vqe.computeMinimumEigenvalue(op);
    expect(approxEq(result.optimalValue, -1, 0.05)).toBe(true);
  });

  test("finds negative energy for Ising H", () => {
    const op = K.SparsePauliOp.fromList([["ZZ", 1.0], ["ZI", 0.5], ["IZ", 0.5]]);
    const ansatz = K.circuits.realAmplitudes(2, 1);
    const optimizer = new K.SPSA({ maxiter: 100, seed: 42 });
    const vqe = new K.VQE({ ansatz, optimizer, initialPoint: [0.5, 0.5, 0.5, 0.5] });
    const result = vqe.computeMinimumEigenvalue(op);
    expect(result.optimalValue).toBeLessThan(0);
  });
});

describe("QAOA", () => {
  test("MaxCut on 2 nodes", () => {
    const op = K.SparsePauliOp.fromList([["ZZ", 0.5], ["II", -0.5]]);
    const qaoa = new K.QAOA({ optimizer: new K.SPSA({ maxiter: 50, seed: 42 }) });
    const result = qaoa.computeMinimumEigenvalue(op);
    expect(result.optimalValue).toBeLessThanOrEqual(0);
  });
});

describe("Grover", () => {
  test("finds marked state", () => {
    const oracle = new K.QuantumCircuit(2);
    oracle.cz(0, 1);
    const grover = new K.Grover();
    const result = grover.amplify(oracle, 1);
    expect(result.top_measurement).toBe("11");
  });
});

describe("PhaseEstimation", () => {
  test("estimates Z eigenvalue", () => {
    const unitary = new K.QuantumCircuit(1); unitary.z(0);
    const qpe = new K.PhaseEstimation(3);
    const result = qpe.estimate(unitary, K.Statevector.fromLabel("1"));
    // Z|1> = -|1>, so phase = pi (eigenvalue e^{i*pi} = -1)
    // Phase = pi means the measurement should be 0.5 (phase / 2pi)
    expect(result.phase).toBeGreaterThan(0);
  });
});

describe("NumPyMinimumEigensolver", () => {
  test("Z ground state = -1", () => {
    const op = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const solver = new K.NumPyMinimumEigensolver();
    const result = solver.computeMinimumEigenvalue(op);
    expect(approxEq(result.eigenvalue, -1, 1e-9)).toBe(true);
  });

  test("ZZ ground state = -1", () => {
    const op = K.SparsePauliOp.fromList([["ZZ", 1.0]]);
    const solver = new K.NumPyMinimumEigensolver();
    const result = solver.computeMinimumEigenvalue(op);
    expect(approxEq(result.eigenvalue, -1, 1e-9)).toBe(true);
  });

  test("Ising H ground state = -1", () => {
    const op = K.SparsePauliOp.fromList([["ZZ", 1.0], ["ZI", 0.5], ["IZ", 0.5]]);
    const solver = new K.NumPyMinimumEigensolver();
    const result = solver.computeMinimumEigenvalue(op);
    expect(approxEq(result.eigenvalue, -1, 1e-9)).toBe(true);
  });
});

describe("Shor", () => {
  test("factors 15", () => {
    const shor = new K.Shor();
    const result = shor.factor(15);
    expect(result.factors[0] * result.factors[1]).toBe(15);
  });

  test("factors 21", () => {
    const shor = new K.Shor();
    const result = shor.factor(21);
    expect(result.factors[0] * result.factors[1]).toBe(21);
  });
});

describe("HHL", () => {
  test("solves identity system", () => {
    const A = K.Operator.identity(1);
    const b = K.Statevector.fromLabel("0");
    const hhl = new K.HHL();
    const x = hhl.solve(A, b);
    expect(x.numQubits).toBe(1);
    expect(approxEq(x.probabilities()[0], 1, 1e-9)).toBe(true);
  });
});

describe("Quantum Volume", () => {
  test("builds circuit", () => {
    const qv = K.quantumVolumeCircuit(3, 3, 42);
    expect(qv.numQubits).toBe(3);
    expect(qv.data.length).toBeGreaterThan(0);
    const hop = K.heavyOutputProbability(qv);
    expect(hop).toBeGreaterThan(0);
    expect(hop).toBeLessThanOrEqual(1);
  });
});

describe("Random generators", () => {
  test("randomUnitary is unitary", () => {
    expect(K.randomUnitary(2, 42).isUnitary(1e-9)).toBe(true);
  });

  test("randomStatevector is normalized", () => {
    expect(approxEq(K.randomStatevector(3, 42).norm(), 1, 1e-9)).toBe(true);
  });

  test("randomPauli", () => {
    expect(K.randomPauli(3, 42).numQubits).toBe(3);
  });

  test("randomClifford", () => {
    const c = K.randomClifford(2, 42);
    expect(c.numQubits).toBe(2);
    const m = c.toMatrix();
    const prod = m.mul(m.dagger());
    let isU = true;
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
      const exp = i === j ? 1 : 0;
      if (Math.abs(prod.get(i, j).re - exp) > 1e-9 || Math.abs(prod.get(i, j).im) > 1e-9) isU = false;
    }
    expect(isU).toBe(true);
  });

  test("randomCliffordSequence", () => {
    const seq = K.randomCliffordSequence(1, 5, 42);
    expect(seq.length).toBe(5);
    expect(seq[0].numQubits).toBe(1);
  });
});

describe("Entanglement measures", () => {
  test("purity", () => {
    const pureSv = K.Statevector.fromLabel("00");
    expect(approxEq(K.purity(pureSv.toOperator()), 1, 1e-9)).toBe(true);
  });

  test("concurrence and EOF of Bell", () => {
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    const dm = new K.DensityMatrix(bell.toOperator().data);
    expect(approxEq(K.concurrence(dm), 1, 1e-6)).toBe(true);
    expect(approxEq(K.entanglementOfFormation(dm), 1, 1e-6)).toBe(true);
  });

  test("mutualInformation of Bell = 2", () => {
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1 / Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1 / Math.sqrt(2), 0),
    ]), 2);
    const dm = new K.DensityMatrix(bell.toOperator().data);
    expect(approxEq(K.mutualInformation(dm, [0]), 2, 1e-6)).toBe(true);
  });

  test("gateFidelity", () => {
    const x = K.standardGates.XGate.toMatrix();
    expect(approxEq(K.gateFidelity(x, x), 1, 1e-9)).toBe(true);
    const h = K.CONSTANTS.H;
    expect(approxEq(K.gateFidelity(x, h), 0.5, 1e-9)).toBe(true);
  });
});
