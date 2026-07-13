import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("Transpiler passes", () => {
  test("BasicSwap inserts SWAPs", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const qc = new K.QuantumCircuit(3);
    qc.cx(0, 2);
    const dag = K.circuitToDag(qc);
    const routed = new K.BasicSwap(cm).run(dag);
    const routedQc = K.dagToCircuit(routed);
    const swapCount = routedQc.data.filter(ci => ci.operation.name === "swap").length;
    expect(swapCount).toBeGreaterThanOrEqual(2);
  });

  test("Decompose expands gates", () => {
    const sub = new K.QuantumCircuit(2); sub.h(0); sub.cx(0, 1);
    const qc = new K.QuantumCircuit(2, 2);
    qc.append(sub.toGate("bell"), [0, 1]);
    qc.measure(0, 0);
    qc.measure(1, 1);
    const dag = K.circuitToDag(qc);
    const decDag = new K.Decompose().run(dag);
    const decQc = K.dagToCircuit(decDag);
    const hCount = decQc.data.filter(ci => ci.operation.name === "h").length;
    expect(hCount).toBe(1);
  });

  test("RemoveFinalMeasurements", () => {
    const qc = new K.QuantumCircuit(2, 2);
    qc.h(0);
    qc.measure(0, 0);
    qc.measure(1, 1);
    const dag = K.circuitToDag(qc);
    const remDag = new K.RemoveFinalMeasurements().run(dag);
    const remQc = K.dagToCircuit(remDag);
    const measCount = remQc.data.filter(ci => ci.operation.name === "measure").length;
    expect(measCount).toBe(0);
  });

  test("CountOps", () => {
    const qc = new K.QuantumCircuit(2, 2);
    qc.append(sub_gate(), [0, 1]);
    qc.measure(0, 0);
    qc.measure(1, 1);
    const dag = K.circuitToDag(qc);
    const countPass = new K.CountOps();
    countPass.run(dag);
    expect(countPass.propertySet.countOps).toBeDefined();
  });

  test("Depth", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    const depthPass = new K.Depth();
    depthPass.run(dag);
    expect(depthPass.propertySet.depth).toBeGreaterThan(0);
  });

  test("CheckMap", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const qc1 = new K.QuantumCircuit(3); qc1.cx(0, 1);
    const dag1 = K.circuitToDag(qc1);
    const check1 = new K.CheckMap(cm);
    check1.run(dag1);
    expect(check1.propertySet.is_mapped).toBe(true);
    const qc2 = new K.QuantumCircuit(3); qc2.cx(0, 2);
    const dag2 = K.circuitToDag(qc2);
    const check2 = new K.CheckMap(cm);
    check2.run(dag2);
    expect(check2.propertySet.is_mapped).toBe(false);
  });

  test("MergeAdjacentBarriers", () => {
    const qc = new K.QuantumCircuit(2);
    qc.barrier();
    qc.barrier();
    qc.h(0);
    const dag = K.circuitToDag(qc);
    const mergedDag = new K.MergeAdjacentBarriers().run(dag);
    const mergedQc = K.dagToCircuit(mergedDag);
    const barrierCount = mergedQc.data.filter(ci => ci.operation.name === "barrier").length;
    expect(barrierCount).toBe(1);
  });

  test("BasisTranslator: H -> U3", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const dag = K.circuitToDag(qc);
    const bt = new K.BasisTranslator(["u3", "cx"]);
    const translatedDag = bt.run(dag);
    const translatedQc = K.dagToCircuit(translatedDag);
    const hCount = translatedQc.data.filter(ci => ci.operation.name === "h").length;
    expect(hCount).toBe(0);
    const u3Count = translatedQc.data.filter(ci => ci.operation.name === "u3").length;
    expect(u3Count).toBeGreaterThanOrEqual(1);
  });

  test("BasisTranslator: SWAP -> 3 CX", () => {
    const qc = new K.QuantumCircuit(2); qc.swap(0, 1);
    const dag = K.circuitToDag(qc);
    const bt = new K.BasisTranslator(["cx", "u3"]);
    const translatedDag = bt.run(dag);
    const translatedQc = K.dagToCircuit(translatedDag);
    const swapCount = translatedQc.data.filter(ci => ci.operation.name === "swap").length;
    expect(swapCount).toBe(0);
    const cxCount = translatedQc.data.filter(ci => ci.operation.name === "cx").length;
    expect(cxCount).toBe(3);
  });

  test("CXCancellation", () => {
    const qc = new K.QuantumCircuit(2);
    qc.cx(0, 1);
    qc.cx(0, 1);
    const dag = K.circuitToDag(qc);
    const cancelledDag = new K.CXCancellation().run(dag);
    const cancelledQc = K.dagToCircuit(cancelledDag);
    const cxCount = cancelledQc.data.filter(ci => ci.operation.name === "cx").length;
    expect(cxCount).toBe(0);
  });

  test("RemoveBarriers", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.barrier();
    qc.x(1);
    const dag = K.circuitToDag(qc);
    const noBarriersDag = new K.RemoveBarriers().run(dag);
    const noBarriersQc = K.dagToCircuit(noBarriersDag);
    const barrierCount = noBarriersQc.data.filter(ci => ci.operation.name === "barrier").length;
    expect(barrierCount).toBe(0);
  });

  test("transpile", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0);
    qc.cx(0, 1);
    const transpiled = K.transpile(qc, { basisGates: ["cx", "u3"] });
    expect(transpiled.data.length).toBeGreaterThan(0);
  });
});

function sub_gate() {
  const sub = new K.QuantumCircuit(2);
  sub.h(0);
  sub.cx(0, 1);
  return sub.toGate("bell");
}
