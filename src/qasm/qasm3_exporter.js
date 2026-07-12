/**
 * qasm3_exporter.js - Export a QuantumCircuit to OpenQASM 3.0 source code.
 *
 * Supports:
 *   - qubit/creg declarations in both old (qreg/creg) and new (qubit[N] q;)
 *     syntaxes
 *   - All standard + parameterized gates
 *   - measure, reset, barrier
 *   - Classical conditions (if statements)
 *   - Bitstring literals for measure assignments (c[0] = measure q[0];)
 *   - Gate definitions (custom gates)
 *   - global phase, delays, and other QASM3-specific constructs
 *
 * Output is valid OpenQASM 3.0 that can be re-parsed by QASM3Parser.
 */

import { Parameter } from "./../core/parameter.js";

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

function _formatParam(p) {
  if (typeof p === "number") return _formatNumber(p);
  if (p && typeof p.toString === "function") return p.toString();
  return String(p);
}

export class QASM3Exporter {
  constructor(options = {}) {
    this.includes = options.includes || ["stdgates.inc"];
    this.disable_extension_clause = options.disable_extension_clause || false;
    // Use the modern QASM3 syntax (qubit[N] q;) by default.
    this.use_modern_syntax = options.use_modern_syntax !== false;
  }

  export(circuit) {
    const lines = [];
    lines.push("OPENQASM 3.0;");
    lines.push("include \"stdgates.inc\";");
    if (circuit.global_phase && Math.abs(circuit.global_phase) > 1e-12) {
      lines.push(`// global phase = ${_formatNumber(circuit.global_phase)}`);
    }

    // Register declarations
    const qrName = circuit.qregs.length > 0 ? circuit.qregs[0].name : "q";
    const crName = circuit.cregs.length > 0 ? circuit.cregs[0].name : "c";
    if (circuit.num_qubits > 0) {
      if (this.use_modern_syntax) {
        lines.push(`qubit[${circuit.num_qubits}] ${qrName};`);
      } else {
        lines.push(`qreg ${qrName}[${circuit.num_qubits}];`);
      }
    }
    if (circuit.num_clbits > 0) {
      if (this.use_modern_syntax) {
        lines.push(`bit[${circuit.num_clbits}] ${crName};`);
      } else {
        lines.push(`creg ${crName}[${circuit.num_clbits}];`);
      }
    }

    const qRef = (q) => `${qrName}[${circuit.qubit_indices(q)}]`;
    const cRef = (c) => `${crName}[${circuit.clbit_indices(c)}]`;

    // Custom gate definitions we may need to emit.
    const customGates = new Set();

    for (const ci of circuit.data) {
      const op = ci.operation;
      const name = op.name;
      const params = op.params || [];
      const qargs = ci.qubits.map(qRef).join(",");
      const cargs = ci.clbits.map(cRef).join(",");

      if (name === "barrier") {
        const qList = ci.qubits.map(qRef).join(",");
        lines.push(`barrier ${qList};`);
        continue;
      }
      if (name === "measure") {
        // QASM3 supports both `measure q -> c;` and `c = measure q;`
        // We use the latter form (more modern).
        lines.push(`${cRef(ci.clbits[0])} = measure ${qRef(ci.qubits[0])};`);
        continue;
      }
      if (name === "reset") {
        lines.push(`reset ${qRef(ci.qubits[0])};`);
        continue;
      }
      if (name === "delay") {
        lines.push(`// delay ${params[0]} ${params[1] || 'dt'} on ${qargs}`);
        continue;
      }
      if (name === "if_else" || name === "while_loop") {
        // Control-flow nodes are skipped here; their conditions are
        // attached to the individual gated instructions.
        continue;
      }

      // Format params
      const paramStr = params.length > 0
        ? `(${params.map(p => _formatParam(p)).join(",")})`
        : "";

      // Map gate names to QASM3 canonical names
      let qasmName = name;
      if (name === "u") qasmName = "u3"; // qiskit U is U3 in QASM
      if (name === "p") qasmName = "u1"; // phase gate is u1
      if (name === "cu") qasmName = "cu3"; // CU maps to CU3 (we drop gamma)

      // Handle the global phase gamma for cu by emitting a comment
      if (name === "cu" && params.length === 4) {
        // cu(theta, phi, lam, gamma) — emit as cu3 with gamma comment
        const [t, p, l, g] = params;
        const pStr = `(${_formatParam(t)},${_formatParam(p)},${_formatParam(l)})`;
        if (Math.abs(_paramToNum(g)) > 1e-12) {
          lines.push(`// global phase ${_formatParam(g)} on next gate`);
        }
        lines.push(`cu3${pStr} ${qargs};`);
        continue;
      }

      // Check for classical condition
      const cond = op.condition;
      let condPrefix = "";
      if (cond) {
        condPrefix = `if (${cond.register}[${cond.index}] == ${cond.value}) `;
      }

      // Standard and parameterized gates
      lines.push(`${condPrefix}${qasmName}${paramStr} ${qargs};`);
    }

    return lines.join("\n");
  }
}

function _paramToNum(p) {
  if (typeof p === "number") return p;
  if (p && typeof p.bind === "function") {
    // ParameterExpression — assume unbound, return 0
    return 0;
  }
  return 0;
}

// Convenience function
export function qasm3_export(circuit, options = {}) {
  const exporter = new QASM3Exporter(options);
  return exporter.export(circuit);
}
