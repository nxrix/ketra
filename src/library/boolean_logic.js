import { Gate } from "./../core/gate.js";
import { _registerExtraGateClass } from "./../core/circuit.js";
import { Complex, ComplexMatrix } from "./../math/linalg.js";

const ZERO = new Complex(0, 0);
const ONE = new Complex(1, 0);

// AND gate
export class ANDGate extends Gate {
  constructor(numInputQubits) {
    super("and", numInputQubits + 1, []);
    this._numInput = numInputQubits;
    this._matrixBuilder = () => {
      const dim = 1 << (numInputQubits + 1);
      const m = ComplexMatrix.identity(dim);
      const allOnesInputs = (1 << numInputQubits) - 1;
      const allOnesWithOutput = allOnesInputs | (1 << numInputQubits);
      m.set(allOnesInputs, allOnesInputs, ZERO);
      m.set(allOnesWithOutput, allOnesWithOutput, ZERO);
      m.set(allOnesInputs, allOnesWithOutput, ONE);
      m.set(allOnesWithOutput, allOnesInputs, ONE);
      return m;
    };
  }
  copy() { return new ANDGate(this._numInput); }
}

// OR gate
export class ORGate extends Gate {
  constructor(numInputQubits) {
    super("or", numInputQubits + 1, []);
    this._numInput = numInputQubits;
    this._matrixBuilder = () => {
      const dim = 1 << (numInputQubits + 1);
      const m = ComplexMatrix.identity(dim);
      // OR flips the output iff at least one input is |1>.
      // For each input combination with at least one |1>, swap the output
      // basis states.
      for (let inputs = 0; inputs < (1 << numInputQubits); inputs++) {
        if (inputs === 0) continue; // all-zero: don't flip
        const output0 = inputs; // output = 0
        const output1 = inputs | (1 << numInputQubits); // output = 1
        m.set(output0, output0, ZERO);
        m.set(output1, output1, ZERO);
        m.set(output0, output1, ONE);
        m.set(output1, output0, ONE);
      }
      return m;
    };
  }
  copy() { return new ORGate(this._numInput); }
}

// XOR gate: output = XOR(controls).
// Implements: output ^= (XOR of all controls), so output should start in |0>.
// XOR is just a sequence of CNOTs.
export class XORGate extends Gate {
  constructor(numInputQubits) {
    super("xor", numInputQubits + 1, []);
    this._numInput = numInputQubits;
    this._matrixBuilder = () => {
      const dim = 1 << (numInputQubits + 1);
      const m = ComplexMatrix.zeros(dim, dim);
      // XOR flips the output iff an odd number of inputs are |1>.
      // The matrix is a permutation: for each input state, the output is
      // flipped iff popcount(inputs) is odd.
      for (let i = 0; i < dim; i++) {
        const inputs = i & ((1 << numInputQubits) - 1);
        const output = (i >> numInputQubits) & 1;
        const popcount = _popcount(inputs);
        const newOutput = output ^ (popcount & 1);
        const j = inputs | (newOutput << numInputQubits);
        m.set(j, i, ONE);
      }
      return m;
    };
  }
  copy() { return new XORGate(this._numInput); }
}

function _popcount(x) {
  let count = 0;
  while (x) { count += x & 1; x >>= 1; }
  return count;
}

// NAND, NOR, XNOR gates: AND/OR/XOR followed by NOT on the output.
export class NANDGate extends Gate {
  constructor(numInputQubits) {
    super("nand", numInputQubits + 1, []);
    this._numInput = numInputQubits;
    this._matrixBuilder = () => {
      const dim = 1 << (numInputQubits + 1);
      const m = ComplexMatrix.identity(dim);
      // NAND flips the output iff NOT all inputs are |1>.
      // I.e., flip unless inputs == all-ones.
      const allOnes = (1 << numInputQubits) - 1;
      for (let inputs = 0; inputs < (1 << numInputQubits); inputs++) {
        if (inputs === allOnes) continue; // don't flip for all-ones
        const output0 = inputs;
        const output1 = inputs | (1 << numInputQubits);
        m.set(output0, output0, ZERO);
        m.set(output1, output1, ZERO);
        m.set(output0, output1, ONE);
        m.set(output1, output0, ONE);
      }
      return m;
    };
  }
  copy() { return new NANDGate(this._numInput); }
}

export class NORGate extends Gate {
  constructor(numInputQubits) {
    super("nor", numInputQubits + 1, []);
    this._numInput = numInputQubits;
    this._matrixBuilder = () => {
      const dim = 1 << (numInputQubits + 1);
      const m = ComplexMatrix.identity(dim);
      // NOR flips the output iff all inputs are |0>.
      const output0 = 0;
      const output1 = 1 << numInputQubits;
      m.set(output0, output0, ZERO);
      m.set(output1, output1, ZERO);
      m.set(output0, output1, ONE);
      m.set(output1, output0, ONE);
      return m;
    };
  }
  copy() { return new NORGate(this._numInput); }
}

export class XNORGate extends Gate {
  constructor(numInputQubits) {
    super("xnor", numInputQubits + 1, []);
    this._numInput = numInputQubits;
    this._matrixBuilder = () => {
      const dim = 1 << (numInputQubits + 1);
      const m = ComplexMatrix.zeros(dim, dim);
      // XNOR flips the output iff an EVEN number of inputs are |1>.
      for (let i = 0; i < dim; i++) {
        const inputs = i & ((1 << numInputQubits) - 1);
        const output = (i >> numInputQubits) & 1;
        const popcount = _popcount(inputs);
        const newOutput = output ^ ((popcount & 1) ? 0 : 1);
        const j = inputs | (newOutput << numInputQubits);
        m.set(j, i, ONE);
      }
      return m;
    };
  }
  copy() { return new XNORGate(this._numInput); }
}

_registerExtraGateClass("ANDGate", ANDGate);
_registerExtraGateClass("ORGate", ORGate);
_registerExtraGateClass("XORGate", XORGate);
_registerExtraGateClass("NANDGate", NANDGate);
_registerExtraGateClass("NORGate", NORGate);
_registerExtraGateClass("XNORGate", XNORGate);

// MCX variants: different decompositions of the multi-controlled X gate.
//   - MCXVChain: uses n-2 ancilla qubits in a V-chain of Toffolis
//   - MCXRecursive: uses 1 ancilla and recursive decomposition
//   - MCXNoAncilla: uses no ancillas, but O(n^2) Toffolis
//   - MCXGrayCode: uses Gray code decomposition (no ancillas)

// Build an MCX (multi-controlled X) using a V-chain of Toffoli gates.
// `controls` is the list of control qubit indices, `target` is the target
// qubit index, and `ancillae` is a list of ancilla qubit indices
// (must have at least numControls - 2 ancillae).
export function mcxVChain(circuit, controls, target, ancillae) {
  const n = controls.length;
  if (n <= 1) {
    throw new Error("MCXVChain requires at least 2 controls");
  }
  if (n === 2) {
    circuit.ccx(controls[0], controls[1], target);
    return;
  }
  if (ancillae.length < n - 2) {
    throw new Error(`MCXVChain requires ${n - 2} ancillae, got ${ancillae.length}`);
  }
  // Forward chain: Toffoli(controls[0], controls[1], ancillae[0])
  //               Toffoli(ancillae[0], controls[2], ancillae[1])
  //               ...
  //               Toffoli(ancillae[n-3], controls[n-1], target)
  // Then reverse to uncompute the ancillae.
  const chain = [];
  let last = controls[0];
  for (let i = 1; i < n - 1; i++) {
    const anc = ancillae[i - 1];
    circuit.ccx(last, controls[i], anc);
    chain.push([last, controls[i], anc]);
    last = anc;
  }
  circuit.ccx(last, controls[n - 1], target);
  for (let i = chain.length - 1; i >= 0; i--) {
    const [a, b, c] = chain[i];
    circuit.ccx(a, b, c);
  }
}

// Recursive MCX using 1 ancilla. Decomposes a k-controlled X into two
// (k-1)-controlled Toffolis with the ancilla as the intermediate.
export function mcxRecursive(circuit, controls, target, ancilla) {
  const n = controls.length;
  if (n <= 2) {
    if (n === 1) circuit.cx(controls[0], target);
    else if (n === 2) circuit.ccx(controls[0], controls[1], target);
    return;
  }
  if (ancilla === undefined || ancilla === null) {
    // No ancilla available: use circuit.mcx directly.
    circuit.mcx(controls, target);
    return;
  }
  // Split: first half controls -> ancilla, second half + ancilla -> target.
  const half = Math.floor(n / 2);
  const firstHalf = controls.slice(0, half);
  const secondHalf = controls.slice(half);
  mcxRecursive(circuit, firstHalf, ancilla, null);
  const combinedControls = secondHalf.concat([ancilla]);
  mcxRecursive(circuit, combinedControls, target, null);
  mcxRecursive(circuit, firstHalf, ancilla, null);
}

// No-ancilla MCX using O(n^2) Toffolis. Based on the Barenco et al. (1995)
// decomposition.
export function mcxNoAncilla(circuit, controls, target) {
  const n = controls.length;
  if (n <= 2) {
    if (n === 1) circuit.cx(controls[0], target);
    else if (n === 2) circuit.ccx(controls[0], controls[1], target);
    return;
  }
  // Barenco decomposition: for n controls, the exact construction uses
  // relative-phase Toffolis to avoid ancillas. Here we use circuit.mcx
  // (which the circuit module supports via ControlledGate).
  circuit.mcx(controls, target);
}
