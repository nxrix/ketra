import { _registerQASMExporterClass } from "./../core/circuit.js";

const PI = Math.PI;

function _formatNumber(n) {
  if (Math.abs(n - Math.PI) < 1e-12) return "pi";
  if (Math.abs(n + Math.PI) < 1e-12) return "-pi";
  if (Math.abs(n - Math.PI / 2) < 1e-12) return "pi/2";
  if (Math.abs(n + Math.PI / 2) < 1e-12) return "-pi/2";
  if (Math.abs(n - Math.PI / 4) < 1e-12) return "pi/4";
  if (Math.abs(n + Math.PI / 4) < 1e-12) return "-pi/4";
  if (Math.abs(n - Math.PI / 8) < 1e-12) return "pi/8";
  if (Math.abs(n + Math.PI / 8) < 1e-12) return "-pi/8";
  if (Math.abs(n - 3 * Math.PI / 2) < 1e-12) return "3*pi/2";
  if (Math.abs(n - 2 * Math.PI) < 1e-12) return "2*pi";
  if (Number.isInteger(n)) return n.toString();
  const piMultiple = n / Math.PI;
  if (Math.abs(piMultiple - Math.round(piMultiple)) < 1e-9) {
    return Math.round(piMultiple) === 1 ? "pi" : `${Math.round(piMultiple)}*pi`;
  }
  return n.toString();
}

function _formatParam(p) {
  if (typeof p === "number") return _formatNumber(p);
  if (p && typeof p.toString === "function") return p.toString();
  return String(p);
}

export class QASMExporter {
  constructor(options = {}) {
    this.includes = options.includes || ["qelib1.inc"];
    this.basisGates = options.basisGates || ["cx", "u3", "u1", "u2"];
  }

  export(circuit) {
    const lines = [];
    lines.push("OPENQASM 2.0;");
    for (const inc of this.includes) {
      lines.push(`include "${inc}";`);
    }

    // Register declarations
    const qrName = circuit.qregs.length > 0 ? circuit.qregs[0].name : "q";
    const crName = circuit.cregs.length > 0 ? circuit.cregs[0].name : "c";
    if (circuit.numQubits > 0) {
      // Combine all qregs into a single register for export simplicity
      lines.push(`qreg ${qrName}[${circuit.numQubits}];`);
    }
    if (circuit.numClbits > 0) {
      lines.push(`creg ${crName}[${circuit.numClbits}];`);
    }

    // Helper: get qubit index as `qrName[i]`
    const qRef = (q) => `${qrName}[${circuit.qubitIndices(q)}]`;
    const cRef = (c) => `${crName}[${circuit.clbitIndices(c)}]`;

    // Track gate definitions we need to emit
    const customGates = new Map();

    // Walk circuit instructions
    for (const ci of circuit.data) {
      const op = ci.operation;
      const name = op.name;
      const params = op.params || [];
      const qargs = ci.qubits.map(qRef).join(",");
      const cargs = ci.clbits.map(cRef).join(",");

      if (name === "barrier") {
        // Barrier on the specified qubits
        const qList = ci.qubits.map(qRef).join(",");
        lines.push(`barrier ${qList};`);
        continue;
      }
      if (name === "measure") {
        // measure q[i] -> c[j];
        lines.push(`measure ${qRef(ci.qubits[0])} -> ${cRef(ci.clbits[0])};`);
        continue;
      }
      if (name === "reset") {
        lines.push(`reset ${qRef(ci.qubits[0])};`);
        continue;
      }
      if (name === "delay") {
        // delay isn't standard QASM 2; skip or emit as a comment
        lines.push(`// delay ${params[0]} ${params[1] || 'dt'} on ${qargs}`);
        continue;
      }
      if (name === "if_else" || name === "whileLoop") {
        // Control-flow nodes are skipped; conditions are
        // attached to the individual gated instructions below.
        continue;
      }

      // If this instruction has a classical condition, emit `if (c[i] == v)`.
      const cond = op.condition;
      let condPrefix = "";
      if (cond) {
        condPrefix = `if (${cond.register}[${cond.index}] == ${cond.value}) `;
      }

      // Format params
      const paramStr = params.length > 0
        ? `(${params.map(p => _formatParam(p)).join(",")})`
        : "";

      // Map gate names to QASM canonical
      let qasmName = name;
      if (name === "u") qasmName = "u3"; // qiskit's U is U3 in QASM 2.0
      if (name === "p") qasmName = "u1"; // phase gate is u1 in QASM 2.0
      if (name === "rx" || name === "ry" || name === "rz" || name === "rxx" || name === "ryy" || name === "rzz" || name === "rzx") {
        // These aren't in qelib1; emit as custom gate definitions.
      }

      // Check if gate is parameterized
      if (params.length > 0) {
        lines.push(`${condPrefix}${qasmName}${paramStr} ${qargs};`);
      } else {
        const stdNames = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "sx", "id",
                          "cx", "cy", "cz", "ch", "swap", "ccx", "cswap", "iswap", "crx", "cry", "crz"];
        if (stdNames.indexOf(name) !== -1) {
          lines.push(`${condPrefix}${qasmName} ${qargs};`);
        } else {
          lines.push(`${condPrefix}${qasmName} ${qargs};`);
        }
      }
    }

    return lines.join("\n");
  }
}

export function qasm2Export(circuit, options = {}) {
  const exporter = new QASMExporter(options);
  return exporter.export(circuit);
}

// Register the QASMExporter class so QuantumCircuit.qasm() can find it.
_registerQASMExporterClass(QASMExporter);
