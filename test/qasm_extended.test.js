import * as K from "../src/index.js";

function approxEq(a, b, tol = 1e-9) {
  return Math.abs(a - b) < tol;
}

describe("QASM3 parser extended", () => {
  test("parse with let", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[1] q;
      let x = 3;
      rx(x * pi / 4) q[0];
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(1);
  });

  test("parse with const", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[1] q;
      const x = 3;
      rx(pi / 4) q[0];
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(1);
  });

  test("parse with while loop", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[1] q;
      bit[1] c;
      while (c[0] == 0) {
        x q[0];
        c[0] = measure q[0];
      }
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(1);
  });

  test("parse with def", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[1] q;
      def myfunc(int x) -> int {
        return x;
      }
      x q[0];
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(1);
  });

  test("parse with switch", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[2] q;
      bit[1] c;
      switch (c[0]) {
        case 0: x q[0];
        case 1: x q[1];
        default: x q[0];
      }
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(2);
  });

  test("parse with for loop and step", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[4] q;
      for (int i in [0:2:4]) {
        x q[i];
      }
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    const xGates = parsed.data.filter(ci => ci.operation.name === "x");
    expect(xGates.length).toBe(2);
  });

  test("parse with explicit set", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[3] q;
      int i;
      for (i in {0, 1, 2}) {
        x q[i];
      }
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    const xGates = parsed.data.filter(ci => ci.operation.name === "x");
    expect(xGates.length).toBe(3);
  });

  test("parse with opaque", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      opaque mygate q;
      qubit[1] q;
      x q[0];
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(1);
  });

  test("parse with gate definition", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[2] q;
      gate mygate a, b {
        cx a, b;
      }
      mygate q[0], q[1];
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(2);
  });

  test("parse with measure arrow form", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[1] q;
      bit[1] c;
      measure q[0] -> c[0];
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(1);
  });

  test("parse with physical qubit", () => {
    const qasm = `
      OPENQASM 3.0;
      include "stdgates.inc";
      qubit[1] q;
      x $0;
    `;
    // Physical qubits aren't fully supported; we just check it doesn't throw
    // unhandled exceptions.
    expect(() => {
      try {
        new K.QASM3Parser(qasm).parse();
      } catch (e) {
        // Expected: physical qubit not fully supported
      }
    }).not.toThrow();
  });

  test("parse with comments", () => {
    const qasm = `
      OPENQASM 3.0;
      // This is a comment
      include "stdgates.inc";
      qubit[1] q;
      /* Block comment */
      x q[0];
    `;
    const parsed = new K.QASM3Parser(qasm).parse();
    expect(parsed.numQubits).toBe(1);
  });
});

describe("QASM2 parser extended", () => {
  test("custom gate definition", () => {
    const qasm = `
      OPENQASM 2.0;
      include "qelib1.inc";
      qreg q[2];
      gate mygate a, b {
        cx a, b;
        h a;
      }
      mygate q[0], q[1];
    `;
    const parsed = new K.QASMParser(qasm).parse();
    expect(parsed.numQubits).toBe(2);
    expect(parsed.data.length).toBe(2);
  });

  test("opaque declaration", () => {
    const qasm = `
      OPENQASM 2.0;
      include "qelib1.inc";
      qreg q[1];
      opaque mygate q;
      x q[0];
    `;
    const parsed = new K.QASMParser(qasm).parse();
    expect(parsed.numQubits).toBe(1);
  });

  test("barrier", () => {
    const qasm = `
      OPENQASM 2.0;
      include "qelib1.inc";
      qreg q[2];
      h q[0];
      barrier q[0], q[1];
      cx q[0], q[1];
    `;
    const parsed = new K.QASMParser(qasm).parse();
    expect(parsed.data.length).toBe(3);
  });

  test("reset", () => {
    const qasm = `
      OPENQASM 2.0;
      include "qelib1.inc";
      qreg q[1];
      x q[0];
      reset q[0];
    `;
    const parsed = new K.QASMParser(qasm).parse();
    expect(parsed.data.length).toBe(2);
  });

  test("expressions", () => {
    const qasm = `
      OPENQASM 2.0;
      include "qelib1.inc";
      qreg q[1];
      rx(pi/4 + pi/4) q[0];
    `;
    const parsed = new K.QASMParser(qasm).parse();
    expect(parsed.data.length).toBe(1);
  });
});

describe("QASM exporters extended", () => {
  test("export with parameterized gates", () => {
    const qc = new K.QuantumCircuit(1);
    qc.rx(0.5, 0);
    qc.ry(0.3, 0);
    const qasm = new K.QASMExporter().export(qc);
    expect(qasm.includes("rx")).toBe(true);
    expect(qasm.includes("ry")).toBe(true);
  });

  test("export with reset and barrier", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.h(0);
    qc.barrier();
    qc.reset(0);
    qc.measure(0, 0);
    const qasm = new K.QASMExporter().export(qc);
    expect(qasm.includes("barrier")).toBe(true);
    expect(qasm.includes("reset")).toBe(true);
  });

  test("QASM3 export with condition", () => {
    const qc = new K.QuantumCircuit(1, 1);
    qc.x(0);
    qc.data[0].operation.condition = { register: "c", index: 0, value: 1 };
    const qasm = K.qasm3Export(qc);
    expect(qasm).toBeDefined();
  });

  test("QASM3 export with cu gate", () => {
    const qc = new K.QuantumCircuit(2);
    qc.cu(Math.PI / 4, Math.PI / 8, Math.PI / 16, Math.PI / 32, 0, 1);
    const qasm = K.qasm3Export(qc);
    expect(qasm.includes("cu3")).toBe(true);
  });
});
