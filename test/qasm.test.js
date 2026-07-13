import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("QASM2 parser", () => {
  test("basic circuit", () => {
    const qasm = `
      OPENQASM 2.0;
      include "qelib1.inc";
      qreg q[2];
      creg c[2];
      h q[0];
      cx q[0], q[1];
      measure q[0] -> c[0];
      measure q[1] -> c[1];
    `;
    const parsed = new K.QASMParser(qasm).parse();
    expect(parsed.numQubits).toBe(2);
    expect(parsed.data.length).toBe(4);
  });

  test("if conditional", () => {
    const qasm = `
      OPENQASM 2.0;
      include "qelib1.inc";
      qreg q[1];
      creg c[1];
      if (c[0] == 1) x q[0];
    `;
    const parsed = new K.QASMParser(qasm).parse();
    const xGate = parsed.data.find(ci => ci.operation.name === "x");
    expect(xGate).toBeDefined();
    expect(xGate.operation.condition).not.toBeNull();
    expect(xGate.operation.condition.register).toBe("c");
    expect(xGate.operation.condition.value).toBe(1);
  });

  test("parameterized gates", () => {
    const qasm = `
      OPENQASM 2.0;
      include "qelib1.inc";
      qreg q[1];
      rx(pi/4) q[0];
    `;
    const parsed = new K.QASMParser(qasm).parse();
    expect(parsed.data.length).toBe(1);
    expect(parsed.data[0].operation.name).toBe("rx");
  });
});

describe("QASM2 exporter", () => {
  test("export circuit", () => {
    const qc = new K.QuantumCircuit(2, 2);
    qc.h(0);
    qc.cx(0, 1);
    qc.measure(0, 0);
    qc.measure(1, 1);
    const qasm = new K.QASMExporter().export(qc);
    expect(qasm.includes("OPENQASM 2.0")).toBe(true);
    expect(qasm.includes("h")).toBe(true);
    expect(qasm.includes("cx")).toBe(true);
  });
});

describe("QASM3 parser", () => {
  test("basic circuit with modern syntax", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[2] q;
      bit[2] c;
      h q[0];
      cx q[0], q[1];
      c[0] = measure q[0];
      c[1] = measure q[1];
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(2);
  });

  test("for loop unrolling", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[3] q;
      for (int i in [0:3]) {
        x q[i];
      }
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    const xGates = parsed.data.filter(ci => ci.operation.name === "x");
    expect(xGates.length).toBe(3);
  });
});

describe("QASM3 exporter", () => {
  test("export circuit", () => {
    const qc = new K.QuantumCircuit(2, 2);
    qc.h(0);
    qc.cx(0, 1);
    qc.measure(0, 0);
    qc.measure(1, 1);
    const qasm = K.qasm3Export(qc);
    expect(qasm.includes("OPENQASM 3.0")).toBe(true);
    expect(qasm.includes("qubit[2]")).toBe(true);
    expect(qasm.includes("bit[2]")).toBe(true);
  });
});
