
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
    this.disableExtensionClause = options.disableExtensionClause || false;
    // Use the modern QASM3 syntax (qubit[N] q;) by default.
    this.useModernSyntax = options.useModernSyntax !== false;
  }

  export(circuit) {
    const lines = [];
    lines.push("OPENQASM 3.0;");
    lines.push("include \"stdgates.inc\";");
    if (circuit.globalPhase && Math.abs(circuit.globalPhase) > 1e-12) {
      lines.push(`// global phase = ${_formatNumber(circuit.globalPhase)}`);
    }

    // Register declarations
    const qrName = circuit.qregs.length > 0 ? circuit.qregs[0].name : "q";
    const crName = circuit.cregs.length > 0 ? circuit.cregs[0].name : "c";
    if (circuit.numQubits > 0) {
      if (this.useModernSyntax) {
        lines.push(`qubit[${circuit.numQubits}] ${qrName};`);
      } else {
        lines.push(`qreg ${qrName}[${circuit.numQubits}];`);
      }
    }
    if (circuit.numClbits > 0) {
      if (this.useModernSyntax) {
        lines.push(`bit[${circuit.numClbits}] ${crName};`);
      } else {
        lines.push(`creg ${crName}[${circuit.numClbits}];`);
      }
    }

    const qRef = (q) => `${qrName}[${circuit.qubitIndices(q)}]`;
    const cRef = (c) => `${crName}[${circuit.clbitIndices(c)}]`;

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
        // QASM3 supports both `measure q -> c;` and `c = measure q;`.
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
      if (name === "if_else" || name === "whileLoop") {
        // Control-flow nodes are skipped; conditions are
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

export function qasm3Export(circuit, options = {}) {
  const exporter = new QASM3Exporter(options);
  return exporter.export(circuit);
}
