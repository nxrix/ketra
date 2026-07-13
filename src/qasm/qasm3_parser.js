import { QuantumCircuit } from "../core/circuit.js";
import { QuantumRegister, ClassicalRegister, Clbit } from "../core/bit.js";
import { Instruction } from "../core/gate.js";

class Token3 {
  constructor(type, value, pos) {
    this.type = type;
    this.value = value;
    this.pos = pos;
  }
}

const KEYWORDS3 = new Set([
  "OPENQASM", "include", "qreg", "creg", "gate", "measure",
  "barrier", "reset", "if", "opaque", "format",
  // QASM 3.0 additions
  "let", "for", "while", "switch", "case", "default", "break", "continue",
  "def", "return", "cal", "defcal", "stretch", "durationof", "extern",
  "const", "int", "uint", "float", "bool", "bit", "qubit", "qubits",
  "in", "array", "void", "input", "output",
]);

function tokenize3(source) {
  const tokens = [];
  let i = 0, line = 1, col = 1;
  while (i < source.length) {
    const c = source[i];
    if (c === " " || c === "\t") { i++; col++; continue; }
    if (c === "\n") { i++; line++; col = 1; continue; }
    if (c === "\r") { i++; continue; }
    if (c === "/" && source[i + 1] === "/") {
      while (i < source.length && source[i] !== "\n") { i++; col++; }
      continue;
    }
    if (c === "/" && source[i + 1] === "*") {
      while (i < source.length && !(source[i] === "*" && source[i + 1] === "/")) {
        if (source[i] === "\n") { line++; col = 1; } else col++;
        i++;
      }
      i += 2; col += 2;
      continue;
    }
    // String literal (also bitstring)
    if (c === "\"") {
      let str = "";
      i++; col++;
      while (i < source.length && source[i] !== "\"") {
        str += source[i]; i++; col++;
      }
      i++; col++;
      tokens.push(new Token3("STRING", str, { line, col }));
      continue;
    }
    if ((c === "-" && source[i + 1] !== ">") || (c >= "0" && c <= "9")) {
      let num = "";
      if (c === "-") { num += c; i++; col++; }
      while (i < source.length && ((source[i] >= "0" && source[i] <= "9") || source[i] === ".")) {
        num += source[i]; i++; col++;
      }
      if (source.slice(i, i + 2) === "pi") {
        const val = parseFloat(num);
        if (source[i + 2] === "/") {
          const denom = parseInt(source.slice(i + 3), 10);
          num = String(val * Math.PI / denom);
          i += 3 + String(denom).length; col += 3 + String(denom).length;
        } else {
          num = String(val * Math.PI);
          i += 2; col += 2;
        }
      }
      tokens.push(new Token3("NUMBER", num, { line, col }));
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
      if (KEYWORDS3.has(id)) tokens.push(new Token3("KEYWORD", id, { line, col }));
      else if (id === "pi") tokens.push(new Token3("NUMBER", String(Math.PI), { line, col }));
      else tokens.push(new Token3("IDENT", id, { line, col }));
      continue;
    }
    // $-prefixed physical qubit
    if (c === "$") {
      let num = "";
      i++; col++;
      while (i < source.length && source[i] >= "0" && source[i] <= "9") {
        num += source[i]; i++; col++;
      }
      tokens.push(new Token3("PHYSICAL_QUBIT", num, { line, col }));
      continue;
    }
    if (c === "-" && source[i + 1] === ">") {
      tokens.push(new Token3("ARROW", "->", { line, col }));
      i += 2; col += 2; continue;
    }
    if (c === "=" && source[i + 1] === "=") {
      tokens.push(new Token3("EQ", "==", { line, col }));
      i += 2; col += 2; continue;
    }
    if (c === "!" && source[i + 1] === "=") {
      tokens.push(new Token3("NEQ", "!=", { line, col }));
      i += 2; col += 2; continue;
    }
    if (c === "<" && source[i + 1] === "=") {
      tokens.push(new Token3("LE", "<=", { line, col }));
      i += 2; col += 2; continue;
    }
    if (c === ">" && source[i + 1] === "=") {
      tokens.push(new Token3("GE", ">=", { line, col }));
      i += 2; col += 2; continue;
    }
    if (c === "+" && source[i + 1] === "+") {
      tokens.push(new Token3("INC", "++", { line, col }));
      i += 2; col += 2; continue;
    }
    if (c === "-" && source[i + 1] === "-") {
      tokens.push(new Token3("DEC", "--", { line, col }));
      i += 2; col += 2; continue;
    }
    const typeMap = {
      "(": "LPAREN", ")": "RPAREN", "[": "LBRACKET", "]": "RBRACKET",
      "{": "LBRACE", "}": "RBRACE", ";": "SEMICOLON", ",": "COMMA",
      "=": "ASSIGN", "+": "PLUS", "*": "STAR", "/": "SLASH",
      "<": "LT", ">": "GT", "&": "AMP", "|": "PIPE", "^": "CARET",
      "~": "TILDE", "%": "PERCENT", "?": "QUESTION", ":": "COLON",
    };
    if (typeMap[c]) {
      tokens.push(new Token3(typeMap[c], c, { line, col }));
      i++; col++; continue;
    }
    i++; col++;
  }
  tokens.push(new Token3("EOF", null, { line, col }));
  return tokens;
}

export class QASM3Parser {
  constructor(source) {
    this.source = source;
    this.tokens = tokenize3(source);
    this.pos = 0;
    this.circuit = null;
    this.qregs = new Map();
    this.cregs = new Map();
    this.classicalVars = new Map(); // QASM 3.0 classical variables
    this.gateDefinitions = new Map();
  }

  peek() { return this.tokens[this.pos]; }
  next() { return this.tokens[this.pos++]; }
  expect(type) {
    const t = this.tokens[this.pos];
    if (t.type !== type) {
      throw new Error(`QASM3 parse error at line ${t.pos.line}: expected ${type} but got ${t.type} (${t.value})`);
    }
    this.pos++;
    return t;
  }
  accept(type) {
    if (this.tokens[this.pos].type === type) return this.tokens[this.pos++];
    return null;
  }

  parse() {
    const versionTok = this.expect("KEYWORD");
    if (versionTok.value !== "OPENQASM") {
      throw new Error(`Expected OPENQASM, got ${versionTok.value}`);
    }
    // Version can be "3.0" or "3" or a number
    const verTok = this.next();
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
        case "qreg": case "qubit": case "qubits": this.parseQregDecl(); return;
        case "creg": case "bit": this.parseCregDecl(); return;
        case "gate": this.parseGateDef(); return;
        case "def": this.parseDef(); return;
        case "measure": this.parseMeasure(); return;
        case "barrier": this.parseBarrier(); return;
        case "reset": this.parseReset(); return;
        case "if": this.parseIf(); return;
        case "for": this.parseFor(); return;
        case "while": this.parseWhile(); return;
        case "switch": this.parseSwitch(); return;
        case "let": this.parseLet(); return;
        case "const": this.parseConst(); return;
        case "int": case "uint": case "float": case "bool": case "array":
          this.parseClassicalDecl(); return;
        case "opaque": this.parseOpaque(); return;
      }
    }
    if (t.type === "IDENT") {
      // QASM 3.0: c[0] = measure q[0]; (assignment form)
      const savedPos = this.pos;
      const name = this.next().value;
      if (this.accept("LBRACKET")) {
        const idx = parseInt(this.expect("NUMBER").value, 10);
        this.expect("RBRACKET");
        if (this.accept("ASSIGN") && this.peek().type === "KEYWORD" && this.peek().value === "measure") {
          this.expect("KEYWORD");
          const qubit = this.parseBitRef();
          this.expect("SEMICOLON");
          if (this.cregs.has(name)) {
            this.circuit.measure(qubit, this.cregs.get(name)._bits[idx]);
          }
          return;
        }
      }
      // Not an assignment — rewind and parse as gate call.
      this.pos = savedPos;
      this.parseGateCall();
      return;
    }
    throw new Error(`QASM3 parse error at line ${t.pos.line}: unexpected ${t.type} (${t.value})`);
  }

  parseInclude() {
    this.expect("KEYWORD");
    this.expect("STRING");
    this.expect("SEMICOLON");
  }

  parseQregDecl() {
    this.expect("KEYWORD");
    // QASM 3.0: qubit[2] q;  OR  qreg q[2];
    let name, size;
    if (this.peek().type === "LBRACKET") {
      this.expect("LBRACKET");
      size = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
      name = this.expect("IDENT").value;
    } else {
      name = this.expect("IDENT").value;
      this.expect("LBRACKET");
      size = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
    }
    this.expect("SEMICOLON");
    const reg = new QuantumRegister(size, name);
    this.circuit.addRegister(reg);
    this.qregs.set(name, reg);
  }

  parseCregDecl() {
    this.expect("KEYWORD");
    let name, size;
    if (this.peek().type === "LBRACKET") {
      // bit[2] c;
      this.expect("LBRACKET");
      size = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
      name = this.expect("IDENT").value;
    } else {
      // creg c[2];
      name = this.expect("IDENT").value;
      this.expect("LBRACKET");
      size = parseInt(this.expect("NUMBER").value, 10);
      this.expect("RBRACKET");
    }
    this.expect("SEMICOLON");
    const reg = new ClassicalRegister(size, name);
    this.circuit.addRegister(reg);
    this.cregs.set(name, reg);
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
    while (this.peek().type === "IDENT") {
      qargs.push(this.expect("IDENT").value);
      this.accept("COMMA");
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

  parseDef() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("LPAREN");
    const args = [];
    while (this.peek().type !== "RPAREN") {
      if (this.peek().type === "KEYWORD") this.next();
      args.push(this.expect("IDENT").value);
      this.accept("COMMA");
    }
    this.expect("RPAREN");
    if (this.accept("ARROW")) {
      this.next();
    }
    this.expect("LBRACE");
    // Parse the body as a list of tokens for later evaluation; quantum
    // subroutines would have gate calls in the body.
    const bodyTokens = [];
    let depth = 1;
    while (depth > 0) {
      const t = this.peek();
      if (t.type === "EOF") break;
      if (t.type === "LBRACE") depth++;
      else if (t.type === "RBRACE") { depth--; if (depth === 0) { this.next(); break; } }
      bodyTokens.push(this.next());
    }
    // Register the def name and arity so it can be referenced.
    this.gateDefinitions.set(name, {
      params: args,
      qargs: [],
      body: [],
      isDef: true,
      bodyTokens,
    });
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
      do { qubits.push(this.parseBitRef()); } while (this.accept("COMMA"));
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
    // if (condition) { body } [else { body }]
    // We attach a classical condition to every instruction in the body.
    // The condition is parsed as a comparison `creg[idx] == value` (the
    // most common QASM3 form). Other comparisons (!=, <, >, <=, >=) are
    // parsed but only == and != are honored at runtime.
    this.expect("KEYWORD");
    this.expect("LPAREN");
    // Parse the condition: creg[idx] OP value
    const cregName = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const cregIdx = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    const opTok = this.next();
    const op = opTok.type === "EQ" ? "==" : opTok.type === "NEQ" ? "!=" :
               opTok.type === "LT" ? "<" : opTok.type === "GT" ? ">" :
               opTok.type === "LE" ? "<=" : opTok.type === "GE" ? ">=" : "==";
    const value = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RPAREN");

    const creg = this.cregs.get(cregName);
    const condition = {
      register: cregName,
      index: cregIdx,
      op,
      value,
      bit: creg ? creg._bits[cregIdx] : null,
    };

    // Body can be a single statement or a { ... } block.
    const dataLenBefore = this.circuit.data.length;
    if (this.accept("LBRACE")) {
      while (this.peek().type !== "RBRACE") {
        this.parseStatement();
      }
      this.expect("RBRACE");
    } else {
      this.parseStatement();
    }
    // Attach the condition to every instruction in the body.
    for (let i = dataLenBefore; i < this.circuit.data.length; i++) {
      this.circuit.data[i].operation.condition = condition;
    }

    // Optional else clause: parse but mark with inverted condition.
    if (this.peek().type === "IDENT" && this.peek().value === "else") {
      this.next();
      const elseCondition = { ...condition, invert: true };
      const elseBefore = this.circuit.data.length;
      if (this.accept("LBRACE")) {
        while (this.peek().type !== "RBRACE") this.parseStatement();
        this.expect("RBRACE");
      } else {
        this.parseStatement();
      }
      for (let i = elseBefore; i < this.circuit.data.length; i++) {
        this.circuit.data[i].operation.condition = elseCondition;
      }
    }
  }

  parseFor() {
    // for (type var in [start:stop]) { body }       — range with step 1
    // for (type var in [start:step:stop]) { body }  — range with explicit step
    // for (type var in {expr1, expr2, ...}) { body } — explicit set
    //
    // The loop is unrolled at parse time: each iteration re-parses the
    // body with `varName` bound to the current value.
    this.expect("KEYWORD");
    this.expect("LPAREN");
    if (this.peek().type === "KEYWORD") this.next();
    const varName = this.expect("IDENT").value;
    this.expect("KEYWORD");

    const values = [];
    if (this.accept("LBRACKET")) {
      // Range [start[:step]:stop]
      const start = parseFloat(this.expect("NUMBER").value);
      let stop, step = 1;
      if (this.accept("COLON")) {
        const next = parseFloat(this.expect("NUMBER").value);
        if (this.accept("COLON")) {
          step = next;
          stop = parseFloat(this.expect("NUMBER").value);
        } else {
          stop = next;
        }
      }
      this.expect("RBRACKET");
      // Generate values. For integer loops, use integer arithmetic.
      if (Number.isInteger(start) && Number.isInteger(step) && Number.isInteger(stop)) {
        if (step > 0) {
          for (let v = start; v < stop; v += step) values.push(v);
        } else {
          for (let v = start; v > stop; v += step) values.push(v);
        }
      } else {
        const n = Math.ceil(Math.abs((stop - start) / step));
        for (let i = 0; i < n; i++) values.push(start + i * step);
      }
    } else if (this.accept("LBRACE")) {
      // Explicit set { expr, expr, ... }
      while (this.peek().type !== "RBRACE") {
        values.push(this.parseExpression());
        if (!this.accept("COMMA")) break;
      }
      this.expect("RBRACE");
    } else {
      // Single expression (e.g., a previously-declared array).
      values.push(this.parseExpression());
    }
    this.expect("RPAREN");

    // Body can be a single statement or a { ... } block.
    // To unroll, we need to re-parse the body for each value. We do this
    // by saving the token position and restoring it for each iteration.
    const bodyStartPos = this.pos;
    if (this.accept("LBRACE")) {
      const braceStartPos = this.pos;
      for (const v of values) {
        this.classicalVars.set(varName, v);
        if (v !== values[0]) {
          this.pos = braceStartPos;
        }
        while (this.peek().type !== "RBRACE") {
          this.parseStatement();
        }
      }
      this.expect("RBRACE");
    } else {
      const stmtStartPos = this.pos;
      for (const v of values) {
        this.classicalVars.set(varName, v);
        if (v !== values[0]) {
          this.pos = stmtStartPos;
        }
        this.parseStatement();
      }
    }
    // Clean up the loop variable after the loop ends.
    this.classicalVars.delete(varName);
  }

  parseWhile() {
    // while (creg[idx] == value) { body }
    // We unroll a bounded number of iterations (max 1000) to avoid
    // infinite loops in malformed QASM. The condition is checked by
    // re-evaluating it against the current classical state.
    this.expect("KEYWORD");
    this.expect("LPAREN");
    const cregName = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const cregIdx = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("EQ");
    const value = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RPAREN");

    // For QASM3 `while`, we need a runtime condition. Since we don't have
    // a runtime classical state at parse time, we unroll a fixed number
    // of iterations (the user is responsible for ensuring the loop
    // terminates; we cap at 1000 to avoid runaway loops).
    const MAX_ITERATIONS = 1000;
    const bodyStartPos = this.pos;
    if (this.accept("LBRACE")) {
      const braceStartPos = this.pos;
      for (let iter = 0; iter < MAX_ITERATIONS; iter++) {
        // Attach a "while" condition to the body instructions so the
        // simulator can skip them when the condition becomes false.
        const dataLenBefore = this.circuit.data.length;
        if (iter > 0) this.pos = braceStartPos;
        while (this.peek().type !== "RBRACE") {
          this.parseStatement();
        }
        const creg = this.cregs.get(cregName);
        const condition = {
          register: cregName,
          index: cregIdx,
          op: "==",
          value,
          bit: creg ? creg._bits[cregIdx] : null,
          isWhile: true,
        };
        for (let i = dataLenBefore; i < this.circuit.data.length; i++) {
          this.circuit.data[i].operation.condition = condition;
        }
      }
      this.expect("RBRACE");
    } else {
      const stmtStartPos = this.pos;
      for (let iter = 0; iter < MAX_ITERATIONS; iter++) {
        if (iter > 0) this.pos = stmtStartPos;
        this.parseStatement();
      }
    }
  }

  parseSwitch() {
    this.expect("KEYWORD"); // switch
    this.expect("LPAREN");
    const switchVar = this.expect("IDENT").value;
    this.expect("LBRACKET");
    const switchIdx = parseInt(this.expect("NUMBER").value, 10);
    this.expect("RBRACKET");
    this.expect("RPAREN");
    this.expect("LBRACE");
    while (this.peek().type !== "RBRACE") {
      const t = this.peek();
      if (t.type === "KEYWORD" && t.value === "case") {
        this.next();
        const caseValue = parseInt(this.expect("NUMBER").value, 10);
        this.accept("COLON");
        // Parse statements until the next case/default/break/}.
        const dataLenBefore = this.circuit.data.length;
        while (this.peek().type !== "RBRACE" &&
               !(this.peek().type === "KEYWORD" && (this.peek().value === "case" || this.peek().value === "default" || this.peek().value === "break"))) {
          this.parseStatement();
        }
        // Attach a switch condition to each statement we just parsed.
        for (let i = dataLenBefore; i < this.circuit.data.length; i++) {
          this.circuit.data[i].operation.condition = {
            register: switchVar,
            index: switchIdx,
            op: "==",
            value: caseValue,
            isSwitch: true,
          };
        }
        if (this.peek().type === "KEYWORD" && this.peek().value === "break") this.next();
      } else if (t.type === "KEYWORD" && t.value === "default") {
        this.next();
        this.accept("COLON");
        // Default case: parse statements (no condition attached, since
        // default applies when no other case matches).
        while (this.peek().type !== "RBRACE" &&
               !(this.peek().type === "KEYWORD" && (this.peek().value === "case" || this.peek().value === "break"))) {
          this.parseStatement();
        }
        if (this.peek().type === "KEYWORD" && this.peek().value === "break") this.next();
      } else {
        this.next();
      }
    }
    this.expect("RBRACE");
  }

  parseLet() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("ASSIGN");
    const value = this.parseExpression();
    this.expect("SEMICOLON");
    this.classicalVars.set(name, value);
  }

  parseClassicalDecl() {
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    if (this.accept("ASSIGN")) {
      const value = this.parseExpression();
      this.classicalVars.set(name, value);
    } else {
      this.classicalVars.set(name, 0);
    }
    this.expect("SEMICOLON");
  }

  parseConst() {
    // const name = value;
    this.expect("KEYWORD");
    const name = this.expect("IDENT").value;
    this.expect("ASSIGN");
    const value = this.parseExpression();
    this.expect("SEMICOLON");
    this.classicalVars.set(name, value);
  }

  parseOpaque() {
    this.expect("KEYWORD");
    while (this.peek().type !== "SEMICOLON") this.next();
    this.expect("SEMICOLON");
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
      case "u": case "u3": this.circuit.u3(params[0], params[1], params[2], qargs[0]); break;
      case "u2": this.circuit.u2(params[0], params[1], qargs[0]); break;
      case "u1": case "p": this.circuit.u1(params[0], qargs[0]); break;
      case "rx": this.circuit.rx(params[0], qargs[0]); break;
      case "ry": this.circuit.ry(params[0], qargs[0]); break;
      case "rz": this.circuit.rz(params[0], qargs[0]); break;
      case "cx": this.circuit.cx(qargs[0], qargs[1]); break;
      case "cy": this.circuit.cy(qargs[0], qargs[1]); break;
      case "cz": this.circuit.cz(qargs[0], qargs[1]); break;
      case "ch": this.circuit.ch(qargs[0], qargs[1]); break;
      case "swap": this.circuit.swap(qargs[0], qargs[1]); break;
      case "iswap": this.circuit.iswap(qargs[0], qargs[1]); break;
      case "ccx": this.circuit.ccx(qargs[0], qargs[1], qargs[2]); break;
      case "cswap": this.circuit.cswap(qargs[0], qargs[1], qargs[2]); break;
      case "crx": this.circuit.crx(params[0], qargs[0], qargs[1]); break;
      case "cry": this.circuit.cry(params[0], qargs[0], qargs[1]); break;
      case "crz": this.circuit.crz(params[0], qargs[0], qargs[1]); break;
      case "cp": case "cu1": this.circuit.cp(params[0], qargs[0], qargs[1]); break;
      case "rxx": this.circuit.rxx(params[0], qargs[0], qargs[1]); break;
      case "ryy": this.circuit.ryy(params[0], qargs[0], qargs[1]); break;
      case "rzz": this.circuit.rzz(params[0], qargs[0], qargs[1]); break;
      default: {
        const instr = new Instruction(n, qargs.length, cargs.length, params);
        this.circuit.append(instr, qargs, cargs);
      }
    }
  }

  parseBitRef() {
    const name = this.expect("IDENT").value;
    if (this.accept("LBRACKET")) {
      // Index can be a NUMBER or a classical variable (e.g., a loop
      // variable like `i` in `for (int i in [0:3]) { x q[i]; }`).
      let idx;
      const t = this.peek();
      if (t.type === "NUMBER") {
        idx = parseInt(this.next().value, 10);
      } else if (t.type === "IDENT" && this.classicalVars.has(t.value)) {
        idx = this.classicalVars.get(t.value);
        this.next();
        if (!Number.isInteger(idx)) {
          throw new Error(`QASM3: register index ${t.value} = ${idx} is not an integer`);
        }
      } else {
        throw new Error(`QASM3: expected NUMBER or classical var in [], got ${t.type} (${t.value})`);
      }
      this.expect("RBRACKET");
      if (this.qregs.has(name)) return this.qregs.get(name)._bits[idx];
      if (this.cregs.has(name)) return this.cregs.get(name)._bits[idx];
      throw new Error(`Unknown register: ${name}`);
    }
    if (this.qregs.has(name)) return this.qregs.get(name)._bits[0];
    if (this.cregs.has(name)) return this.cregs.get(name)._bits[0];
    throw new Error(`Unknown register: ${name}`);
  }

  parseExpression() {
    const t = this.peek();
    let value;
    if (t.type === "NUMBER") {
      this.next();
      value = parseFloat(t.value);
    } else if (t.type === "IDENT") {
      this.next();
      if (this.classicalVars.has(t.value)) {
        value = this.classicalVars.get(t.value);
      } else {
        // Unknown identifier — return its name as a string.
        return t.value;
      }
    } else {
      throw new Error(`QASM3 parse error: expected expression, got ${t.type}`);
    }
    // Handle arithmetic operators (+, *, /).
    while (this.peek().type === "PLUS" || this.peek().type === "STAR" || this.peek().type === "SLASH") {
      const op = this.next().type;
      const rightTok = this.peek();
      let rightVal;
      if (rightTok.type === "NUMBER") {
        rightVal = parseFloat(this.next().value);
      } else if (rightTok.type === "IDENT" && this.classicalVars.has(rightTok.value)) {
        rightVal = this.classicalVars.get(rightTok.value);
        this.next();
      } else {
        break;
      }
      if (op === "PLUS") value += rightVal;
      else if (op === "STAR") value *= rightVal;
      else if (op === "SLASH") value /= rightVal;
    }
    return value;
  }
}

export function qasm3Parse(source) {
  const parser = new QASM3Parser(source);
  return parser.parse();
}
