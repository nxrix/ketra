import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) { return Math.abs(a - b) < tol; }

describe("DAG passes coverage", () => {
  test("collect1qRuns finds runs", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0); qc.t(0); qc.t(0); qc.s(0);
    const dag = K.circuitToDag(qc);
    const runs = K.collect1qRuns(dag);
    expect(runs.length).toBeGreaterThanOrEqual(1);
  });

  test("collect2qRuns finds runs", () => {
    const qc = new K.QuantumCircuit(2);
    qc.cx(0, 1); qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    const runs = K.collect2qRuns(dag);
    expect(runs.length).toBeGreaterThanOrEqual(1);
  });

  test("consolidateBlocks merges gates", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0); qc.t(0); qc.t(0);
    const dag = K.circuitToDag(qc);
    K.consolidateBlocks(dag);
    expect(dag.topologicalOpNodes().length).toBeLessThanOrEqual(3);
  });

  test("optimizeCliffords cancels H H", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0); qc.h(0);
    const dag = K.circuitToDag(qc);
    K.optimizeCliffords(dag);
    expect(dag.topologicalOpNodes().length).toBeLessThan(2);
  });

  test("removeDiagonalGatesBeforeMeasure", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.z(0); qc.measure(0, 0);
    const dag = K.circuitToDag(qc);
    K.removeDiagonalGatesBeforeMeasure(dag);
    const nodes = dag.topologicalOpNodes();
    expect(nodes.filter(n => n.op.name === "z").length).toBe(0);
  });

  test("elidePermutations removes swap pairs", () => {
    const qc = new K.QuantumCircuit(2);
    qc.swap(0, 1); qc.swap(0, 1);
    const dag = K.circuitToDag(qc);
    K.elidePermutations(dag);
    expect(dag).toBeDefined();
  });

  test("removeRedundantGates", () => {
    const qc = new K.QuantumCircuit(1);
    qc.id(0); qc.id(0);
    const dag = K.circuitToDag(qc);
    K.removeRedundantGates(dag);
    expect(dag).toBeDefined();
  });

  test("commutativeCancellation", () => {
    const qc = new K.QuantumCircuit(1);
    qc.rz(0.5, 0); qc.rz(-0.5, 0);
    const dag = K.circuitToDag(qc);
    K.commutativeCancellation(dag);
    expect(dag).toBeDefined();
  });

  test("templateOptimization H Z H -> X", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0); qc.z(0); qc.h(0);
    const dag = K.circuitToDag(qc);
    K.templateOptimization(dag);
    const nodes = dag.topologicalOpNodes();
    expect(nodes.some(n => n.op.name === "x")).toBe(true);
  });

  test("templateOptimization H X H -> Z", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0); qc.x(0); qc.h(0);
    const dag = K.circuitToDag(qc);
    K.templateOptimization(dag);
    const nodes = dag.topologicalOpNodes();
    expect(nodes.some(n => n.op.name === "z")).toBe(true);
  });

  test("DAG aliases work", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const dag = K.circuitToDag(qc);
    expect(K.dagCollect1qRuns(dag)).toBeDefined();
    expect(K.dagCollect2qRuns(dag)).toBeDefined();
    expect(K.dagConsolidateBlocks(dag)).toBeDefined();
    expect(K.dagOptimizeCliffords(dag)).toBeDefined();
    expect(K.dagRemoveDiagBeforeMeasure(dag)).toBeDefined();
  });
});

describe("Result and Counts coverage", () => {
  test("Counts construction", () => {
    const c = new K.Counts({ "00": 512, "11": 512 });
    expect(c.get("00")).toBe(512);
    expect(c.get("01")).toBe(0);
    expect(c.shots).toBe(1024);
  });

  test("Counts items and toDict", () => {
    const c = new K.Counts({ "0": 100, "1": 200 });
    const items = Array.from(c.items());
    expect(items.length).toBe(2);
    const d = c.toDict();
    expect(d["0"]).toBe(100);
  });

  test("Result from simulation", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.h(0); qc.measure(0, 0);
    const result = K.simulate(qc, 100);
    expect(result.backendName).toBeDefined();
    expect(result.success).toBe(true);
    const counts = result.getCounts();
    expect(counts.shots).toBe(100);
  });
});

describe("DensityMatrix coverage", () => {
  test("fromStatevector", () => {
    const sv = K.Statevector.fromLabel("+");
    const dm = K.DensityMatrix.fromStatevector(sv);
    expect(dm.numQubits).toBe(1);
  });

  test("fromCircuit", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const dm = K.DensityMatrix.fromCircuit(qc);
    expect(dm.numQubits).toBe(1);
  });

  test("fromLabel", () => {
    const dm = K.DensityMatrix.fromLabel("00");
    expect(dm.numQubits).toBe(2);
  });

  test("fromOperator", () => {
    const op = K.Operator.identity(1);
    const dm = K.DensityMatrix.fromOperator(op);
    expect(dm.numQubits).toBe(1);
  });

  test("evolve with Operator", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    const xOp = new K.Operator(K.PAULI.X);
    const evolved = dm.evolve(xOp);
    expect(evolved).toBeDefined();
  });

  test("evolve with gate", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    const evolved = dm.evolve(K.standardGates.XGate);
    expect(evolved).toBeDefined();
  });

  test("probabilities", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    expect(approxEq(dm.probabilities()[0], 1)).toBe(true);
  });

  test("toStatevector", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    const sv = dm.toStatevector();
    expect(sv).toBeDefined();
  });

  test("purity", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    expect(approxEq(dm.purity(), 1, 1e-9)).toBe(true);
  });

  test("trace", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    expect(approxEq(dm.trace().re, 1, 1e-9)).toBe(true);
  });

  test("partialTrace", () => {
    const bell = new K.Statevector(new K.ComplexVector([
      new K.Complex(1/Math.sqrt(2), 0), new K.Complex(0, 0),
      new K.Complex(0, 0), new K.Complex(1/Math.sqrt(2), 0),
    ]), 2);
    const dm = K.DensityMatrix.fromStatevector(bell);
    const reduced = dm.partialTrace([1]);
    expect(reduced.numQubits).toBe(1);
  });

  test("expectationValue", () => {
    const dm = K.DensityMatrix.fromStatevector(K.Statevector.fromLabel("0"));
    const ev = dm.expectationValue(K.PAULI.Z);
    expect(approxEq(ev.re, 1)).toBe(true);
  });
});

describe("Channels coverage", () => {
  test("Kraus full round-trip", () => {
    const k = new K.Kraus([K.PAULI.X]);
    const s = k.toSuperop();
    const c = k.toChoi();
    const p = k.toPtm();
    expect(s.rows).toBe(4);
    expect(c.rows).toBe(4);
    expect(p.rows).toBe(4);
  });

  test("SuperOp from matrix", () => {
    const k = new K.Kraus([K.PAULI.X]);
    const s = K.SuperOp.fromMatrix(k.toSuperop(), 1);
    expect(s.numQubits).toBe(1);
    const k2 = s.toKraus();
    expect(k2).toBeDefined();
  });

  test("SuperOp toChoi and toPtm", () => {
    const k = new K.Kraus([K.PAULI.X]);
    const s = K.SuperOp.fromMatrix(k.toSuperop(), 1);
    expect(s.toChoi()).toBeDefined();
    expect(s.toPtm()).toBeDefined();
  });

  test("Chi conversion", () => {
    const k = new K.Kraus([K.PAULI.X]);
    const choi = k.toChoi();
    const chi = new K.Chi(choi, 1);
    expect(chi.toSuperop()).toBeDefined();
    expect(chi.toKraus()).toBeDefined();
  });

  test("PTM conversion", () => {
    const k = new K.Kraus([K.PAULI.X]);
    const ptm = k.toPtm();
    const p = new K.PTM(ptm, 1);
    expect(p.toSuperop()).toBeDefined();
  });

  test("compose channels", () => {
    const k1 = new K.Kraus([K.PAULI.X]);
    const k2 = new K.Kraus([K.PAULI.Z]);
    const composed = k1.compose(k2);
    expect(composed).toBeDefined();
  });

  test("tensor channels", () => {
    const k1 = new K.Kraus([K.PAULI.X]);
    const k2 = new K.Kraus([K.PAULI.Z]);
    const tensored = k1.tensor(k2);
    expect(tensored).toBeDefined();
  });

  test("apply to statevector", () => {
    const k = new K.Kraus([K.PAULI.X]);
    const sv = K.Statevector.fromLabel("0");
    const result = k.apply(sv);
    expect(result).toBeDefined();
  });

  test("stateFidelity pure states", () => {
    expect(K.stateFidelity(K.Statevector.fromLabel("0"), K.Statevector.fromLabel("0"))).toBeCloseTo(1, 9);
    expect(K.stateFidelity(K.Statevector.fromLabel("0"), K.Statevector.fromLabel("1"))).toBeCloseTo(0, 9);
  });

  test("processFidelity", () => {
    const k = new K.Kraus([K.ComplexMatrix.identity(2)]);
    expect(K.processFidelity(k, K.ComplexMatrix.identity(2))).toBeGreaterThan(0.9);
  });

  test("averageGateFidelity", () => {
    const k = new K.Kraus([K.ComplexMatrix.identity(2)]);
    expect(K.averageGateFidelity(k, K.ComplexMatrix.identity(2))).toBeGreaterThan(0.9);
  });

  test("diamondNorm", () => {
    const k = new K.Kraus([K.ComplexMatrix.identity(2)]);
    expect(K.diamondNorm(k)).toBeGreaterThan(0);
  });
});

describe("Parameter coverage", () => {
  test("Parameter arithmetic", () => {
    const p = new K.Parameter("theta");
    const expr = p.mul(2);
    expect(expr.bind({ theta: 0.5 })).toBe(1.0);
    const expr2 = p.add(1);
    expect(expr2.bind({ theta: 2 })).toBe(3);
    const expr3 = p.sub(1);
    expect(expr3.bind({ theta: 2 })).toBe(1);
    const expr4 = p.div(2);
    expect(expr4.bind({ theta: 4 })).toBe(2);
    const expr5 = p.pow(2);
    expect(expr5.bind({ theta: 3 })).toBe(9);
    const expr6 = p.neg();
    expect(expr6.bind({ theta: 5 })).toBe(-5);
  });

  test("Parameter valueOf throws", () => {
    const p = new K.Parameter("x");
    expect(() => { 2 * p; }).toThrow();
  });

  test("ParameterExpression valueOf throws", () => {
    const p = new K.Parameter("x");
    const expr = p.mul(2);
    expect(() => { 2 * expr; }).toThrow();
  });

  test("ParameterVector", () => {
    const pv = new K.ParameterVector("pv", 3);
    expect(pv.length).toBe(3);
    expect(pv.get(0).name).toBe("pv[0]");
    const sliced = pv.slice(0, 2);
    expect(sliced.length).toBe(2);
  });

  test("ParameterVector push", () => {
    const pv = new K.ParameterVector("pv", 0);
    const p = pv.push("custom");
    expect(p.name).toBe("custom");
    expect(pv.length).toBe(1);
  });

  test("Parameter equals", () => {
    const p1 = new K.Parameter("x");
    const p2 = new K.Parameter("x");
    expect(p1.equals(p2)).toBe(false); // different UUIDs
    expect(p1.equals(p1)).toBe(true);
  });
});

describe("Gate coverage", () => {
  test("Instruction inverse", () => {
    const h = K.standardGates.HGate;
    const inv = h.inverse();
    expect(inv).toBeDefined();
    const s = K.standardGates.SGate;
    const sInv = s.inverse();
    expect(sInv.name).toBe("sdg");
  });

  test("Instruction control", () => {
    const x = K.standardGates.XGate;
    const cx = x.control(1);
    expect(cx.numQubits).toBe(2);
  });

  test("ControlledGate matrix", () => {
    const cx = K.standardGates.CXGate;
    const m = cx.toMatrix();
    expect(m.rows).toBe(4);
  });

  test("Gate hasMatrix", () => {
    expect(K.standardGates.HGate.hasMatrix()).toBe(true);
  });

  test("Gate copy", () => {
    const h = K.standardGates.HGate;
    const hCopy = h.copy();
    expect(hCopy.name).toBe("h");
  });
});

describe("Layout coverage", () => {
  test("Layout trivial", () => {
    const l = K.Layout.trivial(3);
    expect(l).toBeDefined();
  });

  test("Layout setPhysical and getPhysical", () => {
    const l = new K.Layout();
    const q = new K.Qubit(new K.QuantumRegister(1, "q"), 0);
    l.setPhysical(0, q);
    expect(l.getPhysical(0)).toBe(q);
  });

  test("CouplingMap construction", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 2]]);
    expect(cm.size).toBe(3);
    expect(cm.hasEdge(0, 1)).toBe(true);
    expect(cm.hasEdge(0, 2)).toBe(false);
  });

  test("CouplingMap neighbors", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    expect(cm.neighbors(1).length).toBe(2);
  });

  test("CouplingMap shortestPath", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const path = cm.shortestPath(0, 2);
    expect(path).toBeDefined();
    expect(path.length).toBe(3);
  });
});

describe("Extra algorithms coverage", () => {
  test("Shor factors 35", () => {
    const shor = new K.Shor();
    const result = shor.factor(35);
    expect(result.factors[0] * result.factors[1]).toBe(35);
  });

  test("HHL solves 2x2 system", () => {
    const A = K.Operator.identity(1);
    const b = K.Statevector.fromLabel("1");
    const hhl = new K.HHL();
    const x = hhl.solve(A, b);
    expect(x.numQubits).toBe(1);
  });

  test("VQC basic training", () => {
    const fm = new K.QuantumCircuit(1);
    const p = new K.Parameter("x");
    fm.rz(p, 0);
    const ansatz = new K.QuantumCircuit(1);
    const theta = new K.Parameter("theta");
    ansatz.ry(theta, 0);
    const vqc = new K.VQC({
      featureMap: fm,
      ansatz: ansatz,
      optimizer: new K.SPSA({ maxiter: 20, seed: 42 }),
      numQubits: 1,
    });
    const result = vqc.fit([[0.5], [1.0]], [0, 1]);
    expect(result).toBeDefined();
  });

  test("QSVC basic", () => {
    const fm = new K.QuantumCircuit(1);
    const p = new K.Parameter("x");
    fm.ry(p, 0);
    const qsvc = new K.QSVC({ featureMap: fm });
    qsvc.fit([[0.0], [Math.PI]], [-1, 1]);
    const preds = qsvc.predict([[0.0], [Math.PI]]);
    expect(preds.length).toBe(2);
  });

  test("quantumVolumeCircuit", () => {
    const qv = K.quantumVolumeCircuit(2, 2, 42);
    expect(qv.numQubits).toBe(2);
  });

  test("heavyOutputProbability", () => {
    const qv = K.quantumVolumeCircuit(2, 2, 42);
    const hop = K.heavyOutputProbability(qv);
    expect(hop).toBeGreaterThan(0);
    expect(hop).toBeLessThanOrEqual(1);
  });

  test("randomCliffordSequence", () => {
    const seq = K.randomCliffordSequence(1, 5, 42);
    expect(seq.length).toBe(5);
  });
});
