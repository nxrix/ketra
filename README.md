# Ketra

[![npm](https://img.shields.io/npm/v/@nxrix/ketra)](https://www.npmjs.com/package/@nxrix/ketra)
[![Minzipped Size](https://badgen.net/bundlephobia/minzip/@nxrix/ketra)](https://bundlephobia.com/package/@nxrix/ketra)
[![Downloads](https://img.shields.io/npm/dm/@nxrix/ketra)](https://www.npmjs.com/package/@nxrix/ketra)
[![GitHub Stars](https://img.shields.io/github/stars/nxrix/ketra?style=social)](https://github.com/nxrix/ketra)
[![GitHub Issues](https://img.shields.io/github/issues/nxrix/ketra)](https://github.com/nxrix/ketra/issues)
<!--[![jsDelivr](https://data.jsdelivr.com/v1/package/npm/@nxrix/ketra/badge?style=rounded)](https://www.jsdelivr.com/package/npm/@nxrix/ketra)-->
<!--[![License](https://img.shields.io/npm/l/@nxrix/ketra)](LICENSE)-->

A JavaScript quantum computing framework inspired by Qiskit

## Features

- Circuits & Simulation
- Algorithms (Grover, VQE, QAOA, HHL, Shor...)
- Statevector, Density Matrix, Sparse Pauli Operators...
- Noise models
- OpenQASM 2 & 3
- Transpiler & DAG circuits
- Visualization (experimental), optimizers & primitives

## Quick start

```
npm install "@nxrix/ketra"
```
or
```html
<script type="module">
  import * as ketra from "https://cdn.jsdelivr.net/npm/@nxrix/ketra@1.2.1/+esm";
</script>
```

### Bell State

```js
import { QuantumCircuit, simulate } from "@nxrix/ketra";

const qc = new QuantumCircuit(2, 2);
qc.h(0);                               // Hadamard on qubit 0
qc.cx(0, 1);                           // CNOT: qubit 0 controls qubit 1
qc.measure(0, 0);                      // measure qubit 0 -> classical bit 0
qc.measure(1, 1);                      // measure qubit 1 -> classical bit 1

const result = simulate(qc, 1024);     // 1024 shots
console.log(result.get_counts().to_dict());
console.log(qc.draw());
```
Output:
```
{ "00": ~512, "11": ~512 }
     ┌───┐     ┌─┐    
q_0: ┤ H ├──■──┤M├─── 
     └───┘┌─┴─┐└╥┘┌─┐ 
q_1: ─────┤ X ├─╫─┤M├ 
          └───┘ ║ └╥┘ 
c: 2/═══════════╩══╩═ 
                0  1  
```

### Grover's Search

```js
import { QuantumCircuit, Grover } from "@nxrix/ketra";

// Oracle: mark the |11> state with a phase flip
const oracle = new QuantumCircuit(2);
oracle.cz(0, 1);

const grover = new Grover();
const result = grover.amplify(oracle, 1);   // 1 marked item
console.log(result.measurement);            // => { "11": 1024 }
```

### Variational Quantum Eigensolver (VQE)

```js
import {
  QuantumCircuit, Parameter, VQE, GradientDescent,
  SparsePauliOp,
} from "@nxrix/ketra";

// Hamiltonian H = Z (ground state energy = -1)
const hamiltonian = SparsePauliOp.from_list([["Z", 1.0]]);

// Ansatz: RY(theta) on qubit 0
const ansatz = new QuantumCircuit(1);
const theta = new Parameter("theta");
ansatz.ry(theta, 0);

const vqe = new VQE({
  ansatz,
  optimizer: new GradientDescent({ learning_rate: 0.5, maxiter: 50 }),
  initial_point: [0.1],
});

const result = vqe.compute_minimum_eigenvalue(hamiltonian);
console.log(result.optimal_value);   // => -1.0 (within optimizer tolerance)
console.log(result.optimal_parameters); // => { theta: ~3.14159 }
```
