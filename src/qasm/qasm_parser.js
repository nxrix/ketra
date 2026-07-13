import { QuantumCircuit } from "./../core/circuit.js";
import { QuantumRegister, ClassicalRegister, Clbit } from "./../core/bit.js";
import { Instruction } from "./../core/gate.js";

class Token {
  constructor(type, value, pos) {
    this.type = type;
    this.value = value;
    this.pos = pos;
  }
  toString() { return `${this.type}(${this.value})`; }
}

const KEYWORDS = new Set([
  "OPENQASM", "include", "qreg", "creg", "gate", "measure",
  "barrier", "reset", "if", "opaque", "format"
]);

function tokenize(source) {
  const tokens = [];
  let i = 0;
  let line = 1, col = 1;
  while (i < source.length) {
    const c = source[i];
    if (c === " " || c === "\t") { i++; col++; continue; }
    if (c === "\n") { i++; line++; col = 1; continue; }
    if (c === "\r") { i++; continue; }
    if (c === "/" && source[i + 1] === "/") {
      while (i < source.length && source[i] !== "\n") { i++; col++; }
      continue;
    }
    if (c === "\"") {
      let str = "";
      i++; col++;
      while (i < source.length && source[i] !== "\"") {
        str += source[i]; i++; col++;
      }
      i++; col++; // closing quote
      tokens.push(new Token("STRING", str, { line, col }));
      continue;
    }
    // Number: integer or decimal, with optional sign (but not if followed by >)
    if ((c === "-" && source[i + 1] !== ">") || (c >= "0" && c <= "9")) {
      let num = "";
      if (c === "-") { num += c; i++; col++; }
      while (i < source.length && ((source[i] >= "0" && source[i] <= "9") || source[i] === ".")) {
        num += source[i]; i++; col++;
      }
      // Check for "*pi" suffix
      if (source[i] === "*") {
        const afterStar = source.slice(i + 1, i + 3);
        if (afterStar === "pi") {
          num = String(parseFloat(num) * Math.PI);
          i += 3; col += 3;
        }
      }
      // Check for trailing pi, pi/2, pi/4 etc.
      if (source.slice(i, i + 2) === "pi") {
        const piVal = parseFloat(num);
        if (source[i + 2] === "/") {
          const denom = parseInt(source.slice(i + 3), 10);
          num = String(piVal * Math.PI / denom);
          i += 3 + String(denom).length;
          col += 3 + String(denom).length;
        } else {
          num = String(piVal * Math.PI);
          i += 2; col += 2;
        }
      }
      tokens.push(new Token("NUMBER", num, { line, col }));
      continue;
    }
    if ((c >= "a" && c <= "z") || (c >= "A" && c <= "Z") || c === "_") {
      let id = "";
      while (i < source.length && ((source[i] >= "a" && source[i] <= "z") ||
                                    (source[i] >= "A" && source[i] <= "Z") ||
                                    (source[i] >= "0" && source[i] <= "9") ||
                                    source[i] === "_")) {
        id += source[i]; i++; col++;
      }
      if (KEYWORDS.has(id)) {
        tokens.push(new Token("KEYWORD", id, { line, col }));
      } else {
        tokens.push(new Token("IDENT", id, { line, col }));
      }
      // Check for "pi" as standalone identifier (e.g., `u1(pi)`)
      if (id === "pi") {
        tokens.pop();
        tokens.push(new Token("NUMBER", String(Math.PI), { line, col }));
      }
      continue;
    }
    if ("()[]{};,->=+*/".indexOf(c) !== -1) {
      if (c === "-" && source[i + 1] === ">") {
        tokens.push(new Token("ARROW", "->", { line, col }));
        i += 2; col += 2; continue;
      }
      if (c === "=" && source[i + 1] === "=") {
        tokens.push(new Token("EQ", "==", { line, col }));
        i += 2; col += 2; continue;
      }
      const typeMap = {
        "(": "LPAREN", ")": "RPAREN", "[": "LBRACKET", "]": "RBRACKET",
        "{": "LBRACE", "}": "RBRACE", ";": "SEMICOLON", ",": "COMMA",
        "=": "ASSIGN", "+": "PLUS", "*": "STAR", "/": "SLASH"
      };
      tokens.push(new Token(typeMap[c] || c, c, { line, col }));
      i++; col++;
      continue;
    }
    i++; col++;
  }
  tokens.push(new Token("EOF", null, { line, col }));
  return tokens;
}

export class QASMParser {
  constructor(source) {
    this.source = source;
    this.tokens = tokenize(source);
    this.pos = 0;
    this.circuit = null;
    this.gateDefinitions = new Map();
    this.qregs = new Map();
    this.cregs = new Map();
    this.allQubits = [];
    this.allClbits = [];
  }

  peek() { return this.tokens[this.pos]; }
  next() { return this.tokens[this.pos++]; }
  expect(type) {
    const t = this.tokens[this.pos];
    if (t.type !== type) {
      throw new Error(`QASM parse error at line ${t.pos.line}: expected ${type} but got ${t.type} (${t.value})`);
    }
    this.pos++;
    return t;
  }
  accept(type) {
    if (this.tokens[this.pos].type === type) {
      return this.tokens[this.pos++];
    }
    return null;
  }

  parse() {
    const versionTok = this.expect("KEYWORD");
    if (versionTok.value !== "OPENQASM") {
      throw new Error(`Expected OPENQASM, got ${versionTok.value}`);
    }
    this.expect("NUMBER");
    this.expect("SEMICOLON");

    this.circuit = new QuantumCircuit();

    while (this.peek().type !== "EOF") {
      this.parseStatement();
    }

    return this.circuit;
  }

  parseStatement() {
    const t = this.peek();
    if (t.type === "KEYWORD") {
      switch (t.value) {
        case "include": this.parseInclude(); return;
        case "qreg": this.parseQregDecl(); return;
        case "creg": this.parseCregDecl(); return;
        case "gate": this.parseGateDef(); return;
        case "measure": this.parseMeasure(); return;
        case "barrier": this.parseBarrier(); return;
        case "reset": this.parseReset(); return;
        case "if": this.parseIf(); return;
        case "opaque": this.parseOpaque(); return;
      }
    }
    if (t.type === "IDENT") {
      this.parseGateCall();
      return;
    }
    throw new Error(`QASM parse error at line ${t.pos.line}: unexpected token ${t.type} (${t.value})`);
  }

  parseInclude() {
    this.expect("KEYWORD");
    this.expect("STRING");
    this.expect("SEMICOLON");
    // Includes are no-ops; we assume qelib1 is built-in.
  }

  parseQregDecl() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const size = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("SEMICOLON");
    const reg = new QuantumRegister(size, name);
    this.circuit.addRegister(reg);
    this.qregs.set(name, reg);
    for (const q of reg._bits) this.allQubits.push({ reg: name, index: q.index, bit: q });
  }

  parseCregDecl() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const size = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("SEMICOLON");
    const reg = new ClassicalRegister(size, name);
    this.circuit.addRegister(reg);
    this.cregs.set(name, reg);
    for (const c of reg._bits) this.allClbits.push({ reg: name, index: c.index, bit: c });
  }

  parseGateDef() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    let params = [];
    if (this.accept("LBRACKET")) {
      while (this.peek().type !== "RBRACKET") {
        params.push(this.expect("IDENT").value);
        this.accept("COMMA");
      }
      this.expect("RBRACKET");
    }
    let qargs = [];
    if (this.peek().type === "IDENT") {
      while (this.peek().type === "IDENT") {
        qargs.push(this.expect("IDENT").value);
        this.accept("COMMA");
      }
    }
    this.expect("LBRACE");
    const body = [];
    while (this.peek().type !== "RBRACE") {
      const gateName = this.expect("IDENT").value;
      let gateParams = [];
      if (this.accept("LPAREN")) {
        while (this.peek().type !== "RPAREN") {
          gateParams.push(this.parseExpression());
          this.accept("COMMA");
        }
        this.expect("RPAREN");
      }
      const gateQargs = [];
      while (this.peek().type === "IDENT") {
        const qName = this.expect("IDENT").value;
        if (this.accept("LBRACKET")) {
          const qIdx = parseInt(this.expect("NUMBER").value, 10);
          this.expect("RBRACKET");
          gateQargs.push({ name: qName, index: qIdx });
        } else {
          gateQargs.push({ name: qName, index: 0 });
        }
        this.accept("COMMA");
      }
      this.expect("SEMICOLON");
      body.push({ gate: gateName, params: gateParams, qargs: gateQargs });
    }
    this.expect("RBRACE");
    this.gateDefinitions.set(name, { params, qargs, body });
  }

  parseOpaque() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    while (this.peek().type !== "SEMICOLON") this.next();
    this.expect("SEMICOLON");
  }

  parseMeasure() {
    this.expect("KEYWORD");
    const qubit = this.parseBitRef();
    this.expect("ARROW");
    const clbit = this.parseBitRef();
    this.expect("SEMICOLON");
    this.circuit.measure(qubit, clbit);
  }

  parseBarrier() {
    this.expect("KEYWORD");
    const qubits = [];
    if (this.peek().type === "IDENT") {
      do {
        qubits.push(this.parseBitRef());
      } while (this.accept("COMMA"));
    }
    this.expect("SEMICOLON");
    if (qubits.length === 0) this.circuit.barrier();
    else this.circuit.barrier(qubits);
  }

  parseReset() {
    this.expect("KEYWORD");
    const qubit = this.parseBitRef();
    this.expect("SEMICOLON");
    this.circuit.reset(qubit);
  }

  parseIf() {
    // if (creg[idx] == value) gate ...;
    // Attach a classical condition to the next instruction. The
    // instruction's `condition` field is `{ register, index, value }`; the
    // simulator / statevector path is responsible for skipping the
    // instruction when the condition is false.
    this.expect("KEYWORD");
    this.expect("LPAREN");
    const cregName = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const cregIdx = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("EQ");
    const value = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RPAREN");

    // Capture the position just before the gated call so we can attach
    // the condition to the resulting CircuitInstruction.
    const dataLenBefore = this.circuit.data.length;
    this.parseGateCall();
    const creg = this.cregs.get(cregName);
    for (let i = dataLenBefore; i < this.circuit.data.length; i++) {
      const ci = this.circuit.data[i];
      ci.operation.condition = {
        register: cregName,
        index: cregIdx,
        value: value,
        bit: creg ? creg._bits[cregIdx] : null,
      };
    }
  }

  parseGateCall() {
    const name = this.expect("IDENT").value;
    let params = [];
    if (this.accept("LPAREN")) {
      while (this.peek().type !== "RPAREN") {
        params.push(this.parseExpression());
        this.accept("COMMA");
      }
      this.expect("RPAREN");
    }
    const qargs = [];
    const cargs = [];
    while (this.peek().type === "IDENT" || this.peek().type === "KEYWORD") {
      // Some keywords might be used as identifiers here
      const ref = this.parseBitRef();
      if (ref instanceof Clbit || (ref && ref.register instanceof ClassicalRegister)) {
        cargs.push(ref);
      } else {
        qargs.push(ref);
      }
      if (!this.accept("COMMA")) break;
    }
    this.expect("SEMICOLON");

    this._applyGateByName(name, params, qargs, cargs);
  }

  _applyGateByName(name, params, qargs, cargs) {
    if (this.gateDefinitions.has(name)) {
      const def = this.gateDefinitions.get(name);
      const argMap = new Map();
      for (let i = 0; i < def.qargs.length && i < qargs.length; i++) {
        argMap.set(def.qargs[i], qargs[i]);
      }
      for (const stmt of def.body) {
        const mappedQargs = stmt.qargs.map(q => argMap.get(q.name) || qargs[0]);
        this._applyGateByName(stmt.gate, stmt.params, mappedQargs, []);
      }
      return;
    }
    const n = name.toLowerCase();
    switch (n) {
      case "h": this.circuit.h(qargs[0]); break;
      case "x": this.circuit.x(qargs[0]); break;
      case "y": this.circuit.y(qargs[0]); break;
      case "z": this.circuit.z(qargs[0]); break;
      case "s": this.circuit.s(qargs[0]); break;
      case "sdg": this.circuit.sdg(qargs[0]); break;
      case "t": this.circuit.t(qargs[0]); break;
      case "tdg": this.circuit.tdg(qargs[0]); break;
      case "sx": this.circuit.sx(qargs[0]); break;
      case "id": this.circuit.id(qargs[0]); break;
      case "u":
      case "u3": this.circuit.u3(params[0], params[1], params[2], qargs[0]); break;
      case "u2": this.circuit.u2(params[0], params[1], qargs[0]); break;
      case "u1":
      case "p": this.circuit.u1(params[0], qargs[0]); break;
      case "rx": this.circuit.rx(params[0], qargs[0]); break;
      case "ry": this.circuit.ry(params[0], qargs[0]); break;
      case "rz": this.circuit.rz(params[0], qargs[0]); break;
      case "cx": this.circuit.cx(qargs[0], qargs[1]); break;
      case "cy": this.circuit.cy(qargs[0], qargs[1]); break;
      case "cz": this.circuit.cz(qargs[0], qargs[1]); break;
      case "ch": this.circuit.ch(qargs[0], qargs[1]); break;
      case "csx": this.circuit.csx(qargs[0], qargs[1]); break;
      case "swap": this.circuit.swap(qargs[0], qargs[1]); break;
      case "iswap": this.circuit.iswap(qargs[0], qargs[1]); break;
      case "ccx": this.circuit.ccx(qargs[0], qargs[1], qargs[2]); break;
      case "cswap": this.circuit.cswap(qargs[0], qargs[1], qargs[2]); break;
      case "crx": this.circuit.crx(params[0], qargs[0], qargs[1]); break;
      case "cry": this.circuit.cry(params[0], qargs[0], qargs[1]); break;
      case "crz": this.circuit.crz(params[0], qargs[0], qargs[1]); break;
      case "cu1":
      case "cp": this.circuit.cp(params[0], qargs[0], qargs[1]); break;
      case "cu3": this.circuit.cu3(params[0], params[1], params[2], qargs[0], qargs[1]); break;
      case "rxx": this.circuit.rxx(params[0], qargs[0], qargs[1]); break;
      case "ryy": this.circuit.ryy(params[0], qargs[0], qargs[1]); break;
      case "rzz": this.circuit.rzz(params[0], qargs[0], qargs[1]); break;
      case "rzx": this.circuit.rzx(params[0], qargs[0], qargs[1]); break;
      default:
        const numQubits = qargs.length;
        const numClbits = cargs.length;
        const instr = new Instruction(n, numQubits, numClbits, params);
        this.circuit.append(instr, qargs, cargs);
    }
  }

  parseBitRef() {
    const name = this.expect("IDENT").value;
    if (this.accept("LBRACKET")) {
      const idx = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
      if (this.qregs.has(name)) {
        return this.qregs.get(name)._bits[idx];
      }
      if (this.cregs.has(name)) {
        return this.cregs.get(name)._bits[idx];
      }
      throw new Error(`Unknown register: ${name}`);
    }
    if (this.qregs.has(name)) return this.qregs.get(name)._bits[0];
    if (this.cregs.has(name)) return this.cregs.get(name)._bits[0];
    throw new Error(`Unknown register: ${name}`);
  }

  parseExpression() {
    const t = this.peek();
    if (t.type === "NUMBER") {
      this.next();
      let value = parseFloat(t.value);
      while (this.peek().type === "PLUS" || this.peek().type === "STAR" || this.peek().type === "SLASH") {
        const op = this.next().type;
        const right = this.expect("NUMBER").value;
        if (op === "PLUS") value += parseFloat(right);
        else if (op === "STAR") value *= parseFloat(right);
        else if (op === "SLASH") value /= parseFloat(right);
      }
      return value;
    }
    if (t.type === "IDENT") {
      this.next();
      if (t.value === "pi") {
        let value = Math.PI;
        while (this.peek().type === "STAR" || this.peek().type === "SLASH") {
          const op = this.next().type;
          const right = this.expect("NUMBER").value;
          if (op === "STAR") value *= parseFloat(right);
          else if (op === "SLASH") value /= parseFloat(right);
        }
        return value;
      }
      return t.value;
    }
    throw new Error(`QASM parse error: expected expression, got ${t.type}`);
  }
}

export function qasm2Parse(source) {
  const parser = new QASMParser(source);
  return parser.parse();
}
