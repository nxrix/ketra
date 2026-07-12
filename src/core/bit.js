/**
 * bit.js - Bit / Register primitives.
 *
 * Bit / Register primitives: Bit, Qubit, Clbit, Register, QuantumRegister,
 * ClassicalRegister. Used to construct QuantumCircuit instances.
 */

let _bitCounter = 0;

export class Bit {
  constructor(register, index) {
    this.register = register;
    this.index = index;
    this._uuid = ++_bitCounter;
  }

  equals(other) {
    return other instanceof Bit &&
      this.register === other.register &&
      this.index === other.index;
  }

  toString() {
    const name = (this.register && this.register.name) || "?";
    return `${name}[${this.index}]`;
  }
}

export class Qubit extends Bit {}
export class Clbit extends Bit {}

export class Register {
  constructor(size, name, bitClass) {
    if (typeof size !== "number" || size <= 0) {
      throw new TypeError("Register size must be a positive integer");
    }
    this._size = size;
    this.name = name;
    this._bits = new Array(size);
    for (let i = 0; i < size; i++) {
      this._bits[i] = new bitClass(this, i);
    }
  }

  get size() { return this._size; }
  set size(v) { this._size = v; }

  get nbits() { return this._size; }

  [Symbol.iterator]() { return this._bits[Symbol.iterator](); }

  get(i) {
    if (i < 0 || i >= this._size) throw new RangeError(`Register index ${i} out of range`);
    return this._bits[i];
  }

  slice(start, end) {
    return this._bits.slice(start, end);
  }

  len() { return this._size; }
}

let _qrCounter = 0;
let _crCounter = 0;

export class QuantumRegister extends Register {
  constructor(size, name) {
    super(size, name || null, Qubit);
    if (!name) this.name = `q${_qrCounter++}`;
  }

  get nqubits() { return this._size; }
}

export class ClassicalRegister extends Register {
  constructor(size, name) {
    super(size, name || null, Clbit);
    if (!name) this.name = `c${_crCounter++}`;
  }

  get nclbits() { return this._size; }
}
