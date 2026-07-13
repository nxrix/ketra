import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("Layout", () => {
  test("trivial layout", () => {
    const layout = K.Layout.trivial(3);
    expect(layout).toBeDefined();
  });

  test("Layout constructor and setPhysical", () => {
    const layout = new K.Layout();
    layout.setPhysical(0, 0);
    layout.setPhysical(1, 1);
    expect(layout).toBeDefined();
  });

  test("CouplingMap", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    expect(cm.size).toBe(3);
    expect(cm.hasEdge(0, 1)).toBe(true);
    expect(cm.hasEdge(0, 2)).toBe(false);
    expect(cm.neighbors(1).length).toBe(2);
  });

  test("CouplingMap shortest path", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const path = cm.shortestPath(0, 2);
    expect(path).toBeDefined();
    expect(path.length).toBe(3);
  });
});

describe("Transpiler layout passes", () => {
  test("TrivialLayout", () => {
    const dag = K.circuitToDag(new K.QuantumCircuit(2));
    const layout = new K.TrivialLayout().run(dag);
    expect(dag._layout).toBeDefined();
  });

  test("DenseLayout", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const dag = K.circuitToDag(new K.QuantumCircuit(2));
    new K.DenseLayout(cm).run(dag);
    expect(dag._layout).toBeDefined();
  });

  test("SabreLayout", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const qc = new K.QuantumCircuit(3);
    qc.cx(0, 1);
    qc.cx(1, 2);
    const dag = K.circuitToDag(qc);
    new K.SabreLayout(cm, 42).run(dag);
    expect(dag._layout).toBeDefined();
  });

  test("ApplyLayout", () => {
    const dag = K.circuitToDag(new K.QuantumCircuit(2));
    new K.TrivialLayout().run(dag);
    new K.ApplyLayout().run(dag);
    expect(dag).toBeDefined();
  });
});

describe("Transpiler routing passes", () => {
  test("LookaheadSwap", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const qc = new K.QuantumCircuit(3);
    qc.cx(0, 2);
    const dag = K.circuitToDag(qc);
    const routed = new K.LookaheadSwap(cm).run(dag);
    expect(routed).toBeDefined();
  });

  test("StochasticSwap", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const qc = new K.QuantumCircuit(3);
    qc.cx(0, 2);
    const dag = K.circuitToDag(qc);
    const routed = new K.StochasticSwap(cm, 42, 3).run(dag);
    expect(routed).toBeDefined();
  });

  test("SabreSwap", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const qc = new K.QuantumCircuit(3);
    qc.cx(0, 2);
    const dag = K.circuitToDag(qc);
    const routed = new K.SabreSwap(cm, "lookahead", 42).run(dag);
    expect(routed).toBeDefined();
  });
});

describe("Transpiler optimization passes", () => {
  test("CommutativeCancellation", () => {
    const qc = new K.QuantumCircuit(1);
    qc.rz(0.5, 0);
    qc.rz(-0.5, 0);
    const dag = K.circuitToDag(qc);
    new K.CommutativeCancellation().run(dag);
    expect(dag).toBeDefined();
  });

  test("Optimize1qGates", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    qc.h(0);
    const dag = K.circuitToDag(qc);
    new K.Optimize1qGates().run(dag);
    expect(dag).toBeDefined();
  });

  test("OptimizeSwapBeforeMeasure", () => {
    const qc = new K.QuantumCircuit(2, 2);
    qc.swap(0, 1);
    qc.measure(0, 0);
    qc.measure(1, 1);
    const dag = K.circuitToDag(qc);
    new K.OptimizeSwapBeforeMeasure().run(dag);
    expect(dag).toBeDefined();
  });

  test("RemoveResetInZeroState", () => {
    const qc = new K.QuantumCircuit(1);
    qc.reset(0);
    qc.h(0);
    const dag = K.circuitToDag(qc);
    new K.RemoveResetInZeroState().run(dag);
    expect(dag).toBeDefined();
  });

  test("DAGFixedPointPass", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    const dag = K.circuitToDag(qc);
    new K.DAGFixedPointPass().run(dag);
    expect(dag).toBeDefined();
  });
});

describe("Transpiler extra passes", () => {
  test("BarrierBeforeFinalMeasurements", () => {
    const qc = new K.QuantumCircuit(2, 2);
    qc.h(0);
    qc.measure(0, 0);
    qc.measure(1, 1);
    const dag = K.circuitToDag(qc);
    new K.BarrierBeforeFinalMeasurements().run(dag);
    expect(dag).toBeDefined();
  });

  test("Collect2qBlocks", () => {
    const qc = new K.QuantumCircuit(2);
    qc.cx(0, 1);
    qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    const pass = new K.Collect2qBlocks();
    pass.run(dag);
    expect(pass.propertySet["2q_blocks"]).toBeDefined();
  });

  test("ConsolidateBlocks", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.t(0);
    qc.t(0);
    const dag = K.circuitToDag(qc);
    new K.ConsolidateBlocks().run(dag);
    expect(dag).toBeDefined();
  });

  test("Unroll3qOrMore", () => {
    const qc = new K.QuantumCircuit(3);
    qc.ccx(0, 1, 2);
    const dag = K.circuitToDag(qc);
    new K.Unroll3qOrMore().run(dag);
    expect(dag).toBeDefined();
  });

  test("GateDirection", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0]]);
    const qc = new K.QuantumCircuit(2);
    qc.cx(1, 0);
    const dag = K.circuitToDag(qc);
    new K.GateDirection(cm).run(dag);
    expect(dag).toBeDefined();
  });

  test("Size and Width", () => {
    const qc = new K.QuantumCircuit(2, 2);
    qc.h(0);
    qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    const sizePass = new K.Size();
    sizePass.run(dag);
    expect(sizePass.propertySet.size).toBeGreaterThan(0);
    const widthPass = new K.Width();
    widthPass.run(dag);
    expect(widthPass.propertySet.width).toBe(4);
  });

  test("Optimize1qGatesDecomposition", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    qc.h(0);
    const dag = K.circuitToDag(qc);
    new K.Optimize1qGatesDecomposition().run(dag);
    expect(dag).toBeDefined();
  });
});

describe("Transpiler core", () => {
  test("transpile with basisGates", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.cx(0, 1);
    const transpiled = K.transpile(qc, { basisGates: ["cx", "u3"] });
    expect(transpiled.data.length).toBeGreaterThan(0);
  });

  test("transpile with coupling map", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const qc = new K.QuantumCircuit(3);
    qc.cx(0, 2);
    const transpiled = K.transpile(qc, { couplingMap: cm });
    expect(transpiled.data.length).toBeGreaterThan(0);
  });

  test("decomposeGate", () => {
    const h = K.standardGates.HGate;
    const result = K.decomposeGate(h, ["u3", "cx"]);
    expect(result).toBeDefined();
  });

  test("PassManager", () => {
    const pm = new K.PassManager();
    expect(pm).toBeDefined();
  });

  test("PassManagerConfig", () => {
    const config = new K.PassManagerConfig({ basisGates: ["cx", "u3"] });
    expect(config).toBeDefined();
  });
});
