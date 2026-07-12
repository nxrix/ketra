/**
 * converters.js - Convert between QuantumCircuit and DAGCircuit.
 *
 * Converters between QuantumCircuit and DAGCircuit.
 */

import { QuantumCircuit } from "./../core/circuit.js";
import { QuantumRegister, ClassicalRegister } from "./../core/bit.js";
import { DAGCircuit } from "./../dagcircuit/dagcircuit.js";

export function circuit_to_dag(circuit) {
  const dag = new DAGCircuit();
  dag.name = circuit.name;
  dag.global_phase = circuit.global_phase;
  dag.metadata = circuit.metadata;

  // Add registers
  for (const r of circuit.qregs) {
    dag.add_qreg(new QuantumRegister(r.size, r.name));
  }
  for (const r of circuit.cregs) {
    dag.add_creg(new ClassicalRegister(r.size, r.name));
  }

  // If the circuit was created with the (n_qubits, n_clbits) shorthand,
  // it has no registers but does have qubits/clbits. Add them with default names.
  if (dag.qubits.length < circuit.qubits.length) {
    const qr = new QuantumRegister(circuit.qubits.length - dag.qubits.length, "q");
    dag.add_qreg(qr);
  }
  if (dag.clbits.length < circuit.clbits.length) {
    const cr = new ClassicalRegister(circuit.clbits.length - dag.clbits.length, "c");
    dag.add_creg(cr);
  }

  // Map circuit qubits/clbits to DAG qubits/clbits (by index)
  const qubitMap = new Map();
  for (let i = 0; i < circuit.qubits.length; i++) {
    qubitMap.set(circuit.qubits[i], dag.qubits[i]);
  }
  const clbitMap = new Map();
  for (let i = 0; i < circuit.clbits.length; i++) {
    clbitMap.set(circuit.clbits[i], dag.clbits[i]);
  }

  // Apply each operation in order
  for (const ci of circuit.data) {
    const op = ci.operation.copy();
    const qargs = ci.qubits.map(q => qubitMap.get(q));
    const cargs = ci.clbits.map(c => clbitMap.get(c));
    dag.apply_operation(op, qargs, cargs);
  }

  return dag;
}

export function dag_to_circuit(dag) {
  // Build a circuit with matching registers
  const regs = [];
  for (const [name, reg] of dag.qregs.entries()) {
    regs.push(new QuantumRegister(reg.bits.length, name));
  }
  for (const [name, reg] of dag.cregs.entries()) {
    regs.push(new ClassicalRegister(reg.bits.length, name));
  }

  const circuit = new QuantumCircuit(...regs);
  circuit.name = dag.name;
  circuit.global_phase = dag.global_phase;
  circuit.metadata = dag.metadata;

  // Map DAG qubits/clbits to circuit qubits/clbits (by index)
  const qubitMap = new Map();
  for (let i = 0; i < dag.qubits.length; i++) {
    qubitMap.set(dag.qubits[i], circuit.qubits[i]);
  }
  const clbitMap = new Map();
  for (let i = 0; i < dag.clbits.length; i++) {
    clbitMap.set(dag.clbits[i], circuit.clbits[i]);
  }

  // Append each op in topological order
  for (const node of dag.topological_op_nodes()) {
    const op = node.op.copy();
    const qargs = node.qargs.map(q => qubitMap.get(q));
    const cargs = node.cargs.map(c => clbitMap.get(c));
    circuit.append(op, qargs, cargs);
  }

  return circuit;
}
