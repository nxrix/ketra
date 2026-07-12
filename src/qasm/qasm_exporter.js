/**
 * qasm_exporter.js - Export a QuantumCircuit to OpenQASM 2.0 source code.
 *
 * */

import { Parameter } from "./../core/parameter.js";
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
  // Try fraction with pi
  const piMultiple = n / Math.PI;
  if (Math.abs(piMultiple - Math.round(piMultiple)) < 1e-9) {
    return Math.round(piMultiple) === 1 ? "pi" : `${Math.round(piMultiple)}*pi`;
  }
  return n.toString();
}

export class QASMExporter {
  constructor(options = {}) {
    this.includes = options.includes || ["qelib1.inc"];
    this.basis_gates = options.basis_gates || ["cx", "u3", "u1", "u2"];
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
    if (circuit.num_qubits > 0) {
      // Combine all qregs into a single register for export simplicity
      lines.push(`qreg ${qrName}[${circuit.num_qubits}];`);
    }
    if (circuit.num_clbits > 0) {
      lines.push(`creg ${crName}[${circuit.num_clbits}];`);
    }

    // Helper: get qubit index as `qrName[i]`
    const qRef = (q) => `${qrName}[${circuit.qubit_indices(q)}]`;
    const cRef = (c) => `${crName}[${circuit.clbit_indices(c)}]`;

    // Track gate definitions we need to emit
    const customGates = new Map(); // name -> definition source

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
      if (name === "if_else" || name === "while_loop") {
        // Emit a comment for control-flow nodes (their condition is
        // already attached to the individual gated instructions, so we
        // don't need to emit anything here). Skipping is correct — the
        // per-instruction `if (...)` is emitted below where the condition
        // is set.
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
        // These aren't in qelib1; emit as custom gate definitions
        // For simplicity, decompose them: rx(θ) = u3(θ, -π/2, π/2) etc.
        // Or just emit them as-is and assume the user has them defined.
      }

      // Check if gate is parameterized
      if (params.length > 0) {
        lines.push(`${condPrefix}${qasmName}${paramStr} ${qargs};`);
      } else {
        // Standard non-parameterized gates
        const stdNames = ["h", "x", "y", "z", "s", "sdg", "t", "tdg", "sx", "id",
                          "cx", "cy", "cz", "ch", "swap", "ccx", "cswap", "iswap", "crx", "cry", "crz"];
        if (stdNames.indexOf(name) !== -1) {
          lines.push(`${condPrefix}${qasmName} ${qargs};`);
        } else {
          // Unknown gate: emit as-is
          lines.push(`${condPrefix}${qasmName} ${qargs};`);
        }
      }
    }

    return lines.join("\n");
  }
}

// Convenience function
export function qasm2_export(circuit, options = {}) {
  const exporter = new QASMExporter(options);
  return exporter.export(circuit);
}

// Register the QASMExporter class so QuantumCircuit.qasm() can find it.
_registerQASMExporterClass(QASMExporter);
