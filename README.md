# Ketra

A JavaScript quantum computing framework inspired by Qiskit

## Features

- Quantum circuits & simulation
- Quantum algorithms (Grover, VQE, QAOA, HHL, Shor...)
- Statevector, Density Matrix & quantum information
- Noise models
- OpenQASM 2 & 3 support
- Transpiler & DAG circuits
- Visualization, optimizers & primitives

## Quick start

```
npm install "@nxrix/ketra"
```
or
```js
import * as ketra from "https://cdn.jsdelivr.net/npm/@nxrix/ketra@1.2.0/+esm"
```

### Bell state

```js
import { QuantumCircuit, simulate } from "@nxrix/ketra";

const qc = new QuantumCircuit(2, 2);   // 2 qubits, 2 classical bits
qc.h(0);                               // Hadamard on qubit 0
qc.cx(0, 1);                           // CNOT: qubit 0 controls qubit 1
qc.measure(0, 0);                      // measure qubit 0 -> classical bit 0
qc.measure(1, 1);                      // measure qubit 1 -> classical bit 1

const result = simulate(qc, 1024);     // 1024 shots
console.log(result.get_counts().to_dict());
// => { "00": 512, "11": 512 }   (50/50 Bell state, varies by shot)
```

### Grover's search

```js
import { QuantumCircuit, Grover } from "@nxrix/ketra";

// Oracle: mark the |11> state with a phase flip
const oracle = new QuantumCircuit(2);
oracle.cz(0, 1);

const grover = new Grover();
const result = grover.amplify(oracle, 1);   // 1 marked item
console.log(result.top_measurement);        // => "11"
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
