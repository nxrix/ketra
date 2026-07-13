import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("Optimizers", () => {
  test("GradientDescent", () => {
    const opt = new K.GradientDescent({ learningRate: 0.1, maxiter: 100 });
    const result = opt.minimize(x => (x[0] - 2) ** 2, [0]);
    expect(Math.abs(result.x[0] - 2)).toBeLessThan(0.5);
  });

  test("COBYLA", () => {
    const opt = new K.COBYLA({ maxiter: 200, rhobeg: 1.0, tol: 1e-6 });
    const result = opt.minimize(x => (x[0] - 3) ** 2 + (x[1] + 2) ** 2, [0, 0]);
    expect(Math.abs(result.x[0] - 3)).toBeLessThan(0.2);
    expect(Math.abs(result.x[1] - (-2))).toBeLessThan(0.2);
  });

  test("SPSA", () => {
    const opt = new K.SPSA({ maxiter: 100, seed: 42 });
    const result = opt.minimize(x => (x[0] - 1) ** 2, [0]);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(0.5);
  });

  test("NFT", () => {
    const nft = new K.NFT({ maxiter: 20 });
    const result = nft.minimize(x => Math.cos(x[0]) - 1, [1.0]);
    expect(approxEq(result.fun, -2, 0.01)).toBe(true);
  });

  test("Adam", () => {
    const opt = new K.Adam({ learningRate: 0.05, maxiter: 200 });
    const result = opt.minimize(x => (x[0] - 1) ** 2, [0]);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(0.5);
  });

  test("LBFGSB", () => {
    const opt = new K.LBFGSB({ maxiter: 100 });
    const result = opt.minimize(x => (x[0] - 1) ** 2 + (x[1] - 2) ** 2, [0, 0]);
    expect(Math.abs(result.x[0] - 1)).toBeLessThan(0.5);
  });

  test("NelderMead", () => {
    const opt = new K.NelderMead({ maxiter: 200 });
    const result = opt.minimize(x => (x[0] - 2) ** 2, [0]);
    expect(Math.abs(result.x[0] - 2)).toBeLessThan(0.5);
  });
});

describe("Gradients", () => {
  test("ParamShift", () => {
    const qc = new K.QuantumCircuit(1);
    const p = new K.Parameter("theta");
    qc.ry(p, 0);
    const obs = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const grad = new K.ParamShift().compute(qc, obs, { theta: 0.5 });
    expect(grad.length).toBe(1);
    expect(Math.abs(grad[0])).toBeGreaterThan(0);
  });

  test("FiniteDiff", () => {
    const qc = new K.QuantumCircuit(1);
    const p = new K.Parameter("theta");
    qc.ry(p, 0);
    const obs = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const grad = new K.FiniteDiff().compute(qc, obs, { theta: 0.5 });
    expect(grad.length).toBe(1);
  });

  test("NaturalGradient", () => {
    const qc = new K.QuantumCircuit(1);
    const p = new K.Parameter("theta");
    qc.ry(p, 0);
    const obs = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const regGrad = new K.ParamShift().compute(qc, obs, { theta: 0.5 });
    const natGrad = new K.NaturalGradient(1e-3).compute(qc, obs, { theta: 0.5 });
    expect(Math.abs(natGrad[0])).toBeGreaterThan(Math.abs(regGrad[0]));
  });

  test("LinearCombination", () => {
    const qc = new K.QuantumCircuit(1);
    const p = new K.Parameter("theta");
    qc.ry(p, 0);
    const obs = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const grad = new K.LinearCombination().compute(qc, obs, { theta: 0.5 });
    expect(grad.length).toBe(1);
  });
});
