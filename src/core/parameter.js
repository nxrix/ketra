/**
 * parameter.js - Parameter and ParameterExpression for parameterized circuits.
 *
 * Symbolic placeholders for circuits that need to be bound to numeric values
 * before simulation (e.g. VQE / QAOA ansaetze).
 */

let _paramCounter = 0;

export class Parameter {
  constructor(name) {
    this.name = name || `param_${_paramCounter++}`;
    this._uuid = ++_paramCounter;
  }

  add(other) { return new ParameterExpression("+", this, _coerce(other)); }
  sub(other) { return new ParameterExpression("-", this, _coerce(other)); }
  mul(other) { return new ParameterExpression("*", this, _coerce(other)); }
  div(other) { return new ParameterExpression("/", this, _coerce(other)); }
  neg()      { return new ParameterExpression("*", this, _coerce(-1)); }
  pow(n)     { return new ParameterExpression("**", this, _coerce(n)); }

  bind(values) {
    if (typeof values === "number") return values;
    if (!(this.name in values)) {
      throw new Error(`Missing value for parameter ${this.name}`);
    }
    return Number(values[this.name]);
  }

  get parameters() {
    return new Set([this]);
  }

  toString() { return this.name; }

  equals(other) {
    return other instanceof Parameter && other._uuid === this._uuid;
  }

  // Reject implicit numeric coercion so that bugs like `2 * param.mul(c)`
  // fail loudly instead of silently producing NaN.
  valueOf() {
    throw new TypeError(
      `Parameter "${this.name}" cannot be implicitly coerced to a number. ` +
      `Call .bind(values) first, or chain ParameterExpression methods: ` +
      `use param.mul(2) instead of 2 * param.`
    );
  }
}

export class ParameterExpression {
  constructor(op, left, right) {
    this.op = op;
    this.left = left;
    this.right = right;
  }

  add(other) { return new ParameterExpression("+", this, _coerce(other)); }
  sub(other) { return new ParameterExpression("-", this, _coerce(other)); }
  mul(other) { return new ParameterExpression("*", this, _coerce(other)); }
  div(other) { return new ParameterExpression("/", this, _coerce(other)); }
  pow(n)     { return new ParameterExpression("**", this, _coerce(n)); }

  bind(values) {
    const l = _bindNode(this.left, values);
    const r = _bindNode(this.right, values);
    switch (this.op) {
      case "+": return l + r;
      case "-": return l - r;
      case "*": return l * r;
      case "/": return l / r;
      case "**": return Math.pow(l, r);
      default: throw new Error(`Unknown op ${this.op}`);
    }
  }

  get parameters() {
    const set = new Set();
    _collectParams(this, set);
    return set;
  }

  toString() {
    return `(${_nodeToString(this.left)} ${this.op} ${_nodeToString(this.right)})`;
  }

  // Same rationale as Parameter.valueOf: catch implicit coercion bugs.
  valueOf() {
    throw new TypeError(
      `ParameterExpression "${this.toString()}" cannot be implicitly coerced to a number. ` +
      `Call .bind(values) first.`
    );
  }
}

function _coerce(v) {
  if (v instanceof Parameter || v instanceof ParameterExpression) return v;
  if (typeof v === "number") return v;
  throw new TypeError("Cannot coerce value to ParameterExpression");
}

function _bindNode(node, values) {
  if (typeof node === "number") return node;
  if (node instanceof Parameter) return node.bind(values);
  if (node instanceof ParameterExpression) return node.bind(values);
  throw new TypeError("Unbindable parameter node");
}

function _collectParams(node, set) {
  if (node instanceof Parameter) set.add(node);
  else if (node instanceof ParameterExpression) {
    _collectParams(node.left, set);
    _collectParams(node.right, set);
  }
}

function _nodeToString(node) {
  if (typeof node === "number") return String(node);
  return node.toString();
}

export class ParameterVector {
  constructor(name, length) {
    this.name = name || `pv_${_paramCounter++}`;
    this._params = new Array(length || 0);
    for (let i = 0; i < this._params.length; i++) {
      this._params[i] = new Parameter(`${this.name}[${i}]`);
    }
  }

  get length() { return this._params.length; }
  get(i) { return this._params[i]; }
  [Symbol.iterator]() { return this._params[Symbol.iterator](); }
  slice(start, end) { return this._params.slice(start, end); }
  push(name) {
    const p = new Parameter(name || `${this.name}[${this._params.length}]`);
    this._params.push(p);
    return p;
  }
}
