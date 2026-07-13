import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) { return Math.abs(a - b) < tol; }

describe("Transpiler core coverage", () => {
  test("transpile with optimization level 0", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0); qc.cx(0, 1);
    const t = K.transpile(qc, { optimizationLevel: 0 });
    expect(t.data.length).toBeGreaterThan(0);
  });

  test("transpile with optimization level 1", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0); qc.cx(0, 1);
    const t = K.transpile(qc, { optimizationLevel: 1, basisGates: ["cx", "u3"] });
    expect(t.data.length).toBeGreaterThan(0);
  });

  test("transpile with optimization level 2", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0); qc.h(0); qc.cx(0, 1);
    const t = K.transpile(qc, { optimizationLevel: 2, basisGates: ["cx", "u3"] });
    expect(t).toBeDefined();
  });

  test("transpile with coupling map", () => {
    const cm = new K.CouplingMap([[0, 1], [1, 0], [1, 2], [2, 1]]);
    const qc = new K.QuantumCircuit(3);
    qc.cx(0, 2);
    const t = K.transpile(qc, { couplingMap: cm, basisGates: ["cx", "u3"] });
    expect(t.data.length).toBeGreaterThan(0);
  });

  test("decomposeGate", () => {
    const h = K.standardGates.HGate;
    const result = K.decomposeGate(h, ["u3", "cx"]);
    expect(result).toBeDefined();
  });

  test("PassManager", () => {
    const pm = new K.PassManager();
    expect(pm.passes.length).toBe(0);
    const qc = new K.QuantumCircuit(2);
    qc.cx(0, 1); qc.cx(0, 1);
    pm.append(new K.CXCancellation());
    const dag = K.circuitToDag(qc);
    pm.run(dag);
    expect(dag).toBeDefined();
  });

  test("PassManagerConfig", () => {
    const config = new K.PassManagerConfig({ basisGates: ["cx", "u3"] });
    expect(config.basisGates).toContain("cx");
  });

  test("presetPassManager", () => {
    const qc = new K.QuantumCircuit(2);
    qc.h(0); qc.cx(0, 1);
    const pm = K.presetPassManager(qc, { basisGates: ["cx", "u3"] });
    expect(pm).toBeDefined();
  });

  test("DECOMP_RULES and DEFAULT_BASIS", () => {
    expect(K.DECOMP_RULES).toBeDefined();
    expect(K.DEFAULT_BASIS).toBeDefined();
  });
});

describe("Extra gates coverage", () => {
  test("UnitaryGate from Operator", () => {
    const op = new K.Operator(K.PAULI.X);
    const g = new K.UnitaryGate(op);
    expect(g.numQubits).toBe(1);
  });

  test("UnitaryGate from 2D array", () => {
    const g = new K.UnitaryGate([[1, 0], [0, 1]]);
    expect(g.numQubits).toBe(1);
  });

  test("UnitaryGate copy", () => {
    const g = new K.UnitaryGate(K.PAULI.X);
    const c = g.copy();
    expect(c.numQubits).toBe(1);
  });

  test("UnitaryGate control", () => {
    const g = new K.UnitaryGate(K.PAULI.X);
    const cg = g.control(1);
    expect(cg.numQubits).toBe(2);
  });

  test("DiagonalGate copy", () => {
    const g = new K.DiagonalGate([0, Math.PI]);
    const c = g.copy();
    expect(c.numQubits).toBe(1);
  });

  test("PermutationGate copy", () => {
    const g = new K.PermutationGate([1, 0]);
    const c = g.copy();
    expect(c.numQubits).toBe(2);
  });

  test("HamiltonianGate from SparsePauliOp", () => {
    const h = K.SparsePauliOp.fromList([["Z", 1.0]]);
    const g = new K.HamiltonianGate(h, Math.PI / 2);
    expect(g.numQubits).toBe(1);
  });

  test("HamiltonianGate copy", () => {
    const g = new K.HamiltonianGate(K.PAULI.X, 1.0);
    const c = g.copy();
    expect(c.numQubits).toBe(1);
  });

  test("Initialize with label", () => {
    const i = new K.Initialize("01");
    expect(i.numQubits).toBe(2);
  });

  test("Initialize copy", () => {
    const i = new K.Initialize([1, 0]);
    const c = i.copy();
    expect(c.numQubits).toBe(1);
  });

  test("MCPhaseGate", () => {
    const g = new K.MCPhaseGate(Math.PI / 4, 2);
    const m = g.toMatrix();
    expect(m.rows).toBe(8);
  });

  test("MCRXGate, MCRYGate, MCRZGate", () => {
    expect(new K.MCRXGate(Math.PI / 2, 1).toMatrix().rows).toBe(4);
    expect(new K.MCRYGate(Math.PI / 2, 1).toMatrix().rows).toBe(4);
    expect(new K.MCRZGate(Math.PI / 2, 1).toMatrix().rows).toBe(4);
  });

  test("MCMTGate single target", () => {
    const base = K.standardGates.XGate;
    const g = new K.MCMTGate(base, 1, 1);
    const m = g.toMatrix();
    expect(m.rows).toBe(4);
  });

  test("MCMTGate multi target", () => {
    const base = K.standardGates.XGate;
    const g = new K.MCMTGate(base, 1, 2);
    const m = g.toMatrix();
    expect(m.rows).toBe(8);
  });
});

describe("Statevector extended coverage", () => {
  test("fromInstruction", () => {
    const h = K.standardGates.HGate;
    const sv = K.Statevector.fromInstruction(h);
    expect(sv.numQubits).toBe(1);
  });

  test("fromInt", () => {
    const sv = K.Statevector.fromInt(3, 2);
    expect(sv.probabilities()[3]).toBe(1);
  });

  test("zero and one", () => {
    expect(K.Statevector.zero(2).probabilities()[0]).toBe(1);
    expect(K.Statevector.one(2).probabilities()[3]).toBe(1);
  });

  test("conjugate", () => {
    const sv = K.Statevector.fromLabel("+");
    const conj = sv.conjugate();
    expect(conj).toBeDefined();
  });

  test("inner and dot", () => {
    const sv1 = K.Statevector.fromLabel("0");
    const sv2 = K.Statevector.fromLabel("0");
    expect(sv1.inner(sv2).re).toBe(1);
    expect(sv1.dot(sv2).re).toBe(1);
  });

  test("norm", () => {
    expect(approxEq(K.Statevector.fromLabel("0").norm(), 1)).toBe(true);
  });

  test("tensor and expand", () => {
    const sv1 = K.Statevector.fromLabel("0");
    const sv2 = K.Statevector.fromLabel("1");
    const t = sv1.tensor(sv2);
    expect(t.numQubits).toBe(2);
    const e = sv1.expand(sv2);
    expect(e.numQubits).toBe(2);
  });

  test("toOperator", () => {
    const sv = K.Statevector.fromLabel("0");
    const op = sv.toOperator();
    expect(op.numQubits).toBe(1);
  });

  test("equiv", () => {
    const sv1 = K.Statevector.fromLabel("+");
    const sv2 = new K.Statevector(new K.ComplexVector([
      new K.Complex(-1/Math.sqrt(2), 0), new K.Complex(-1/Math.sqrt(2), 0),
    ]), 1);
    expect(sv1.equiv(sv2)).toBe(true);
  });

  test("toDict", () => {
    const sv = K.Statevector.fromLabel("0");
    const d = sv.toDict();
    expect(d["0"]).toBeDefined();
  });

  test("sampleMemory", () => {
    const sv = K.Statevector.fromLabel("0");
    const mem = sv.sampleMemory(10);
    expect(mem.length).toBe(10);
  });

  test("measure all", () => {
    const sv = K.Statevector.fromLabel("00");
    const meas = sv.measure(null, () => 0);
    expect(meas.bits).toBe("00");
  });
});

describe("Operator extended coverage", () => {
  test("fromMatrix", () => {
    const op = K.Operator.fromMatrix(K.PAULI.X);
    expect(op.numQubits).toBe(1);
  });

  test("fromLabel", () => {
    const op = K.Operator.fromLabel("X");
    expect(op.numQubits).toBe(1);
  });

  test("fromGate", () => {
    const op = K.Operator.fromGate(K.standardGates.HGate);
    expect(op.numQubits).toBe(1);
  });

  test("identity and zero", () => {
    expect(K.Operator.identity(1).numQubits).toBe(1);
    expect(K.Operator.zero(1).numQubits).toBe(1);
  });

  test("fromPhase", () => {
    const op = K.Operator.fromPhase(new K.Complex(0, 1), 1);
    expect(op.numQubits).toBe(1);
  });

  test("compose, tensor, expand", () => {
    const x = new K.Operator(K.PAULI.X);
    const z = new K.Operator(K.PAULI.Z);
    expect(x.compose(z).numQubits).toBe(1);
    expect(x.tensor(z).numQubits).toBe(2);
    expect(x.expand(z).numQubits).toBe(2);
  });

  test("adjoint, conjugate, transpose", () => {
    const x = new K.Operator(K.PAULI.X);
    expect(x.adjoint()).toBeDefined();
    expect(x.conjugate()).toBeDefined();
    expect(x.transpose()).toBeDefined();
  });

  test("trace and det", () => {
    const x = new K.Operator(K.PAULI.X);
    expect(x.trace().re).toBe(0);
    expect(x.det().re).toBe(-1);
  });

  test("pow", () => {
    const x = new K.Operator(K.PAULI.X);
    const x2 = x.pow(2);
    expect(x2).toBeDefined();
  });

  test("exp", () => {
    const x = new K.Operator(K.PAULI.X);
    const exp = x.exp();
    expect(exp).toBeDefined();
  });

  test("equalsUpToPhase", () => {
    const x = new K.Operator(K.PAULI.X);
    const xScaled = new K.Operator(K.PAULI.X.scale(new K.Complex(0, 1)));
    expect(x.equalsUpToPhase(xScaled)).toBe(true);
  });

  test("eigvals and svd", () => {
    const z = new K.Operator(K.PAULI.Z);
    const ev = z.eigvals();
    expect(ev.length).toBe(2);
    const svd = z.svd();
    expect(svd).toBeDefined();
  });

  test("toDict and toString", () => {
    const x = new K.Operator(K.PAULI.X);
    expect(x.toDict()).toBeDefined();
    expect(typeof x.toString()).toBe("string");
  });
});

describe("Clifford extended coverage", () => {
  test("fromStabilizers", () => {
    const c = K.Clifford.fromStabilizers(["ZI", "IZ"]);
    expect(c.numQubits).toBe(2);
  });

  test("random", () => {
    const c = K.Clifford.random(2, 42);
    expect(c.numQubits).toBe(2);
  });

  test("toDict", () => {
    const c = K.Clifford.fromLabel("00");
    const d = c.toDict();
    expect(d).toBeDefined();
  });

  test("fromMatrix", () => {
    const qc = new K.QuantumCircuit(1); qc.h(0);
    const cliff = K.Clifford.fromCircuit(qc);
    const m = cliff.toMatrix();
    const c2 = K.Clifford.fromMatrix(m);
    expect(c2.numQubits).toBe(1);
  });

  test("toString", () => {
    const c = K.Clifford.fromLabel("0");
    expect(typeof c.toString()).toBe("string");
  });
});

describe("Noise extended coverage", () => {
  test("QuantumError", () => {
    const e = new K.QuantumError([K.PAULI.X], 1);
    expect(e).toBeDefined();
  });

  test("NoiseModel is_empty", () => {
    const nm = new K.NoiseModel();
    expect(nm.isEmpty()).toBe(true);
  });

  test("NoiseModel with readout error", () => {
    const nm = new K.NoiseModel();
    const re = new K.ReadoutError([[0.9, 0.1], [0.1, 0.9]]);
    nm.addReadoutError(re, [0]);
    expect(nm.isEmpty()).toBe(false);
  });

  test("ReadoutError apply", () => {
    const re = new K.ReadoutError([[1, 0], [0, 1]]);
    const result = re.apply([0, 1]);
    expect(result).toEqual([0, 1]);
  });

  test("ReadoutError copy", () => {
    const re = new K.ReadoutError([[0.9, 0.1], [0.1, 0.9]]);
    const c = re.copy();
    expect(c.probabilities).toEqual(re.probabilities);
  });
});

describe("QASM extended coverage", () => {
  test("QASMExporter with options", () => {
    const qc = new K.QuantumCircuit(1);
    qc.h(0);
    const exporter = new K.QASMExporter({ includes: [], basisGates: ["cx", "u3"] });
    const qasm = exporter.export(qc);
    expect(qasm.includes("OPENQASM")).toBe(true);
  });

  test("QASM3Exporter with options", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.h(0);
    qc.measure(0, 0);
    const exporter = new K.QASM3Exporter({ useModernSyntax: false });
    const qasm = exporter.export(qc);
    expect(qasm.includes("qreg")).toBe(true);
  });

  test("QASM2 round-trip", () => {
    const qasm = `OPENQASM 2.0; include "qelib1.inc"; qreg q[1]; h q[0];`;
    const parsed = K.qasm2Parse(qasm);
    expect(parsed.numQubits).toBe(1);
    const exported = new K.QASMExporter().export(parsed);
    expect(exported.includes("h")).toBe(true);
  });

  test("QASM3 round-trip", () => {
    const qasm = `OPENQASM 3.0; include "stdgates.inc"; qubit[1] q; h q[0];`;
    const parsed = K.qasm3Parse(qasm);
    expect(parsed.numQubits).toBe(1);
    const exported = K.qasm3Export(parsed);
    expect(exported.includes("h")).toBe(true);
  });
});
