# Ketra

[![npm](https://img.shields.io/npm/v/@nxrix/ketra?logo=npm&color=cb0000)](https://www.npmjs.com/package/@nxrix/ketra)
[![Minzipped Size](https://img.shields.io/bundlejs/size/@nxrix/ketra?color=86f)](https://bundlephobia.com/package/@nxrix/ketra)
[![License](https://img.shields.io/npm/l/@nxrix/ketra?color=86f)](LICENSE)
[![Downloads](https://img.shields.io/npm/dm/@nxrix/ketra?color=cb0000)](https://www.npmjs.com/package/@nxrix/ketra)
[![jsDelivr](https://data.jsdelivr.com/v1/package/npm/@nxrix/ketra/badge?style=rounded)](https://www.jsdelivr.com/package/npm/@nxrix/ketra)
[![GitHub Issues](https://img.shields.io/github/issues/nxrix/ketra)](https://github.com/nxrix/ketra/issues)
[![GitHub Stars](https://img.shields.io/github/stars/nxrix/ketra?style=social)](https://github.com/nxrix/ketra)
<!--
[![Minzipped Size](https://badgen.net/bundlephobia/minzip/@nxrix/ketra)](https://bundlephobia.com/package/@nxrix/ketra)
-->

A JavaScript quantum computing framework inspired by Qiskit.

## Features

- Circuits & Simulation (statevector + noisy)
- Algorithms: Grover, VQE, QAOA, HHL, Shor, VQC, QSVC, NumPy eigensolver
- Quantum info: Statevector, DensityMatrix, Operator, Pauli, SparsePauliOp, Clifford, StabilizerState, channels
- Noise models: depolarizing, bit flip, phase flip, amplitude/phase damping, Kraus, readout
- OpenQASM 2 & 3 (parser + exporter)
- Transpiler: layout (Trivial, Dense, Sabre), routing (Basic, Lookahead, Stochastic, Sabre), optimization passes
- DAG circuit, converters, analysis passes
- Primitives: Estimator/Sampler (V1 + V2)
- Optimizers: COBYLA, SPSA, Adam, L-BFGS-B, SLSQP, Nelder-Mead, NFT, GradientDescent
- Gradients: ParamShift, FiniteDiff, LinearCombination, NaturalGradient
- Arithmetic circuits: adders, multipliers, comparators, QFT adder
- Boolean logic gates: AND, OR, XOR, NAND, NOR, XNOR
- Extra gates: UnitaryGate, DiagonalGate, PermutationGate, HamiltonianGate, Initialize, MCMT
- Random generators: randomUnitary, randomStatevector, randomPauli, randomClifford
- Entanglement measures: purity, concurrence, entanglementOfFormation, mutualInformation, gateFidelity
- Circuit library: Bell, GHZ, QFT, QuantumVolume, RealAmplitudes, EfficientSU2, TwoLocal, GraphState
- Visualization: text/ASCII circuit drawer

## Quick start

```bash
npm install "@nxrix/ketra"
```
or
```html
<script type="module">
  import * as ketra from "https://cdn.jsdelivr.net/npm/@nxrix/ketra@1.3.1/+esm";
</script>
```

### Bell State

```js
const qc = new ketra.QuantumCircuit(2, 2);
qc.h(0);
qc.cx(0, 1);
qc.measure(0, 0);
qc.measure(1, 1);

const result = ketra.simulate(qc, 1024);
console.log(result.getCounts().toDict());
console.log(qc.draw());
```
Output:
```js
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
const oracle = new ketra.QuantumCircuit(2);
oracle.cz(0, 1);

const grover = new ketra.Grover();
const result = grover.amplify(oracle, 1);
console.log(result.measurement);
``````
Output:
```js
{ "11": 1024 }
```

### Variational Quantum Eigensolver (VQE)

```js
const hamiltonian = ketra.SparsePauliOp.fromList([["Z", 1.0]]);

const ansatz = new ketra.QuantumCircuit(1);
const theta = new ketra.Parameter("theta");
ansatz.ry(theta, 0);

const vqe = new ketra.VQE({
  ansatz,
  optimizer: new ketra.GradientDescent({ learningRate: 0.5, maxIter: 50 }),
  initialPoint: [0.1],
});

const result = vqe.computeMinimumEigenvalue(hamiltonian);
console.log(result.optimalValue);
console.log(result.optimalParameters);
```
Output:
```js
-1
{ theta: ~3.14159 }
```
