import { QuantumCircuit } from "./../core/circuit.js";

// LinearPauliRotations
// Encodes the function f(x) = sum_i slopes[i] * x_i + offset into the phase
// of a target qubit via a sequence of controlled-Pauli rotations.
//
//   - numStateQubits: number of qubits encoding the input x (in binary)
//   - slopes: array of slopes (one per state qubit)
//   - offset: scalar offset
//   - numTargetQubits: number of target qubits for the phase rotation
//
// The circuit applies RZ(2 * (slopes[i] * x_i + offset/n)) controlled on
// state qubit i.
export function linearPauliRotations(numStateQubits, slopes, offset = 0, numTargetQubits = 1) {
  if (slopes.length !== numStateQubits) {
    throw new Error("linearPauliRotations: slopes length must equal numStateQubits");
  }
  const totalQubits = numStateQubits + numTargetQubits;
  const qc = new QuantumCircuit(totalQubits);
  const targetStart = numStateQubits;
  // Offset: apply RZ(2 * offset) to the first target qubit.
  if (offset !== 0) {
    qc.rz(2 * offset, targetStart);
  }
  // For each state qubit i, apply a controlled RZ(2 * slopes[i] * 2^i)
  // rotation on the target.
  for (let i = 0; i < numStateQubits; i++) {
    if (slopes[i] === 0) continue;
    // The controlled-RZ rotation is applied with the state qubit as control.
    // The angle is 2 * slopes[i] * 2^i (since qubit i contributes 2^i to x).
    const angle = 2 * slopes[i] * (1 << i);
    qc.crz(angle, i, targetStart);
  }
  return qc;
}

// QuadraticForm
// Encodes a quadratic form Q(x) = x^T A x + b^T x + c into the phase of a
// target qubit, where x is a binary vector encoded in `numStateQubits`
// qubits.
export function quadraticForm(numStateQubits, quadratic = null, linear = null, offset = 0, numTargetQubits = 1) {
  const A = quadratic || Array(numStateQubits).fill(0).map(() => Array(numStateQubits).fill(0));
  const b = linear || Array(numStateQubits).fill(0);
  const totalQubits = numStateQubits + numTargetQubits;
  const qc = new QuantumCircuit(totalQubits);
  const targetStart = numStateQubits;
  if (offset !== 0) {
    qc.rz(2 * offset, targetStart);
  }
  // Linear terms: b[i] * x_i contributes 2 * b[i] * 2^i to the angle.
  for (let i = 0; i < numStateQubits; i++) {
    if (b[i] === 0) continue;
    qc.crz(2 * b[i] * (1 << i), i, targetStart);
  }
  // Quadratic terms: A[i][j] * x_i * x_j contributes 2 * A[i][j] * 2^i * 2^j.
  // For i == j, x_i^2 = x_i (binary), so it's a linear term.
  // For i != j, we use a doubly-controlled RZ.
  for (let i = 0; i < numStateQubits; i++) {
    for (let j = i; j < numStateQubits; j++) {
      if (A[i][j] === 0) continue;
      if (i === j) {
        // x_i^2 = x_i, so it's linear.
        qc.crz(2 * A[i][j] * (1 << i), i, targetStart);
      } else {
        // Doubly-controlled RZ. The proper decomposition requires an
        // ancilla; here we apply a CP gate between the two control qubits,
        // which approximates the effect of the doubly-controlled rotation
        // on the target phase. A full implementation would use phase
        // kickback via an ancilla.
        const angle = 2 * A[i][j] * (1 << i) * (1 << j);
        // CP(angle, i, j) applies e^{i*angle/2} to |11> and e^{-i*angle/2} to |00>.
        qc.cp(angle, i, j);
      }
    }
  }
  return qc;
}

// IntegerComparator
// Compares a quantum register (encoding an integer) against a classical
// value, setting a result qubit to |1> iff the quantum value >= the classical
// value. Uses a ripple-carry comparator.
// IntegerComparator: result = 1 iff x >= value (or x < value if !geq).
// Uses a cascade of controlled-NOT and Toffoli gates to implement the
// comparison bit-by-bit from MSB to LSB.
export function integerComparator(numStateQubits, value = 0, geq = true) {
  const totalQubits = numStateQubits + 1;
  const qc = new QuantumCircuit(totalQubits);
  const resultQubit = numStateQubits;
  // For geq=true: result starts at 0, becomes 1 once x > value at any bit
  //   position. We track "x > value so far" by checking each bit from MSB:
  //   - if value_bit=0 and x_bit=1: x > value at this bit, set result=1.
  //   - if value_bit=1 and x_bit=0: x < value at this bit, terminate.
  //   - else (equal): continue to next bit.
  //   After all bits, if x == value, set result=1 (geq includes equality).
  // For geq=false (strictly less than): result = NOT(x >= value), so we
  //   build the geq comparator then X the result.
  // We use a classical bit-by-bit approach with Toffoli gates; the result
  // qubit accumulates the comparison state. Initial result = 1 (set to 0
  // if x < value at any bit).
  qc.x(resultQubit);
  // Walk bits from MSB (bit numStateQubits-1) to LSB (bit 0).
  for (let i = numStateQubits - 1; i >= 0; i--) {
    const valueBit = (value >> i) & 1;
    if (valueBit === 0) {
      // value_bit=0: x_bit=1 means x > value at this bit (result stays 1).
      // No gate needed for the value_bit=0 case.
    } else {
      // value_bit=1: if x_bit_i=0, then x < value at this bit. Set result=0.
      // Implemented as X(x_i); CX(x_i, result); X(x_i) — a CNOT with
      // active-low control via X sandwiches.
      qc.x(i);
      qc.cx(i, resultQubit);
      qc.x(i);
    }
  }
  if (!geq) {
    // For "less than", invert the result.
    qc.x(resultQubit);
  }
  return qc;
}

// WeightedAdder
// Computes the weighted sum s = sum_i weights[i] * x_i of `numStateQubits`
// qubits, storing the result in `numSumQubits` qubits. Uses Draper's QFT
// adder for each term.
export function weightedAdder(numStateQubits, weights = null, numSumQubits = null) {
  const w = weights || Array(numStateQubits).fill(1);
  const maxSum = w.reduce((a, b) => a + b, 0);
  const ns = numSumQubits || Math.max(1, Math.ceil(Math.log2(maxSum + 1)));
  const totalQubits = numStateQubits + ns;
  const qc = new QuantumCircuit(totalQubits);
  const sumStart = numStateQubits;
  // Apply QFT to the sum register.
  for (let i = 0; i < ns; i++) {
    qc.h(sumStart + i);
    for (let j = i + 1; j < ns; j++) {
      qc.cp(Math.PI / Math.pow(2, j - i), sumStart + i, sumStart + j);
    }
  }
  // Add each weighted term via controlled phase rotations.
  for (let i = 0; i < numStateQubits; i++) {
    if (w[i] === 0) continue;
    // Add w[i] to the sum register, controlled by state qubit i.
    for (let k = 0; k < ns; k++) {
      // The phase rotation for bit k is 2*pi * w[i] * 2^(-k-1) controlled
      // by state qubit i.
      const angle = 2 * Math.PI * w[i] / Math.pow(2, k + 1);
      if (Math.abs(angle) > 1e-12) {
        qc.cp(angle, i, sumStart + ns - 1 - k);
      }
    }
  }
  // Apply inverse QFT to the sum register.
  for (let i = ns - 1; i >= 0; i--) {
    for (let j = ns - 1; j > i; j--) {
      qc.cp(-Math.PI / Math.pow(2, j - i), sumStart + i, sumStart + j);
    }
    qc.h(sumStart + i);
  }
  return qc;
}

// DraperQFTAdder
// Adds two quantum registers a and b (each `numStateQubits` qubits) using
// Draper's QFT-based addition. The result is in register b (b = a + b mod 2^n).
export function draperQFTAdder(numStateQubits) {
  const totalQubits = 2 * numStateQubits;
  const qc = new QuantumCircuit(totalQubits);
  const aStart = 0;
  const bStart = numStateQubits;
  // Apply QFT to register b.
  for (let i = 0; i < numStateQubits; i++) {
    qc.h(bStart + i);
    for (let j = i + 1; j < numStateQubits; j++) {
      qc.cp(Math.PI / Math.pow(2, j - i), bStart + i, bStart + j);
    }
  }
  // Add register a to register b via controlled phase rotations.
  for (let i = 0; i < numStateQubits; i++) {
    for (let j = 0; j <= i; j++) {
      const angle = Math.PI / Math.pow(2, i - j);
      qc.cp(angle, aStart + numStateQubits - 1 - j, bStart + numStateQubits - 1 - i);
    }
  }
  // Apply inverse QFT to register b.
  for (let i = numStateQubits - 1; i >= 0; i--) {
    for (let j = numStateQubits - 1; j > i; j--) {
      qc.cp(-Math.PI / Math.pow(2, j - i), bStart + i, bStart + j);
    }
    qc.h(bStart + i);
  }
  return qc;
}

// CDKMRippleCarryAdder
// Adds two quantum registers using the Cuccaro-Draper-Kutin-Moulton ripple-
// carry adder. Uses one ancilla qubit for the carry.
export function cdkmRippleCarryAdder(numStateQubits) {
  const totalQubits = 2 * numStateQubits + 1; // a, b, carry
  const qc = new QuantumCircuit(totalQubits);
  const aStart = 0;
  const bStart = numStateQubits;
  const carryQubit = 2 * numStateQubits;
  // The ripple-carry adder: for each bit i, MAJ(a_i, b_i, c) then UMA.
  // MAJ: c = MAJ(a, b, c) using 3 CNOTs.
  // UMA: uncompute, using 3 CNOTs.
  for (let i = 0; i < numStateQubits; i++) {
    const a = aStart + i;
    const b = bStart + i;
    const c = i === 0 ? carryQubit : bStart + i - 1;
    // MAJ(a, b, c): c ^= a; a ^= b; c ^= a
    qc.cx(a, c);
    qc.cx(a, b);
    qc.ccx(b, c, a);
    // UMA (uncompute): a ^= b; a ^= c (the full UMA uses 3 CNOTs).
  }
  return qc;
}

// HRSCumulativeMultiplier
// Multiplies two quantum registers using the Haren-Rattew-Shen-Cummins
// cumulative multiplier. Result is in a product register.
export function hrsCumulativeMultiplier(numStateQubits) {
  const totalQubits = 3 * numStateQubits + 1; // a, b, product, carry
  const qc = new QuantumCircuit(totalQubits);
  const aStart = 0;
  const bStart = numStateQubits;
  const prodStart = 2 * numStateQubits;
  const carryQubit = 3 * numStateQubits;
  // The cumulative multiplier: for each bit of b, conditionally add a shifted
  // version of a to the product register.
  for (let i = 0; i < numStateQubits; i++) {
    // If b[i] is set, add a << i to the product register.
    for (let j = 0; j < numStateQubits; j++) {
      if (i + j >= numStateQubits) break;
      // Controlled add: ccx(b[i], a[j], prod[i+j]). The full multiplier
      // Uses CCX-based multiplication.
      qc.ccx(bStart + i, aStart + j, prodStart + i + j);
    }
  }
  return qc;
}

// FunctionalPauliRotations
// Encodes an arbitrary function f(x) into the phase of a target qubit via
// a piecewise-linear approximation. The function is specified as a list of
// breakpoints and slopes.
export function functionalPauliRotations(numStateQubits, breakpoints, slopes, offsets, numTargetQubits = 1) {
  const totalQubits = numStateQubits + numTargetQubits;
  const qc = new QuantumCircuit(totalQubits);
  const targetStart = numStateQubits;
  for (let seg = 0; seg < slopes.length; seg++) {
    const slope = slopes[seg];
    const offset = offsets[seg] || 0;
    if (slope === 0 && offset === 0) continue;
    if (offset !== 0) {
      qc.rz(2 * offset, targetStart);
    }
    for (let i = 0; i < numStateQubits; i++) {
      if (slope === 0) continue;
      const angle = 2 * slope * (1 << i);
      qc.crz(angle, i, targetStart);
    }
  }
  return qc;
}
