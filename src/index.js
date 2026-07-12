/**
 * index.js
 */

// Math
export {
  Complex, ComplexVector, ComplexMatrix,
  kronMatrices, kronVectors, sampleDistribution, randomUniform,
  PAULI, CONSTANTS,
} from "./math/linalg.js";

// Core
export { Bit, Qubit, Clbit, Register, QuantumRegister, ClassicalRegister } from "./core/bit.js";
export { Parameter, ParameterExpression, ParameterVector } from "./core/parameter.js";
export { Instruction, Gate, ControlledGate } from "./core/gate.js";
export {
  QuantumCircuit, CircuitInstruction,
  _registerStd, _registerParamBuilder, _getStdGate, _setDrawHook,
} from "./core/circuit.js";

// Library
export * as standardGates from "./library/standard_gates.js";
export * as generalizedGates from "./library/generalized_gates.js";
export {
  UnitaryGate, DiagonalGate, PermutationGate, HamiltonianGate,
  Initialize, MCPhaseGate, MCRXGate, MCRYGate, MCRZGate, MCMTGate,
} from "./library/extra_gates.js";
import "./library/extra_gates.js"; // side-effect: registers extra gate classes
export {
  ANDGate, ORGate, XORGate, NANDGate, NORGate, XNORGate,
  mcxVChain, mcxRecursive, mcxNoAncilla,
} from "./library/boolean_logic.js";
import "./library/boolean_logic.js"; // side-effect: registers Boolean gates
export {
  linearPauliRotations, quadraticForm, integerComparator,
  weightedAdder, draperQFTAdder, cdkmRippleCarryAdder,
  hrsCumulativeMultiplier, functionalPauliRotations,
} from "./library/arithmetic.js";
import * as _arithmetic from "./library/arithmetic.js";
export const arithmetic = _arithmetic;
export {
  bellState, ghzState, qft, qftInverse, iqft,
  quantumVolume, realAmplitudes, efficientSU2,
  pauliTwoDesign, twoLocal, graphState, pauliEvolution,
  hiddenLinearFunction,
} from "./library/circuits.js";
import * as _circuits from "./library/circuits.js";
export const circuits = _circuits;

// Quantum info
export { Statevector } from "./quantum_info/statevector.js";
export { Operator } from "./quantum_info/operator.js";
export { Pauli, PauliList, SparsePauliOp } from "./quantum_info/pauli.js";
export { DensityMatrix } from "./quantum_info/density_matrix.js";
export { Clifford, StabilizerState } from "./quantum_info/clifford.js";
export { ScalarOp } from "./quantum_info/scalar_op.js";
export { SchmidtDecomposition } from "./quantum_info/schmidt.js";

// Quantum channels
export {
  QuantumChannel, Kraus, SuperOp, Chi, PTM,
  state_fidelity, process_fidelity, average_gate_fidelity, diamond_norm,
} from "./quantum_info/channels.js";

// Additional quantum_info functions (random generators, entanglement measures)
export {
  random_unitary, random_statevector, random_pauli, random_clifford,
  random_density_matrix,
  purity, concurrence, entanglement_of_formation,
  mutual_information, gate_fidelity, unitarity,
} from "./quantum_info/quantum_info_extra.js";

// Simulator
export { StatevectorSimulator, simulate } from "./simulator/statevector_simulator.js";

// Noisy simulator
export { QasmSimulator, simulate_noisy } from "./simulator/qasm_simulator.js";

// Noise models
export {
  QuantumError, NoiseModel, ReadoutError,
  depolarizing_error, bit_flip_error, phase_flip_error,
  amplitude_damping_error, phase_damping_error,
  pauli_x_error, pauli_y_error, pauli_z_error,
  reset_error, kraus_error, mixed_unitary_error,
  combine_errors,
} from "./noise/noise_models.js";

// Transpiler
export {
  transpile, decomposeGate, PassManager, PassManagerConfig,
  presetPassManager, DECOMP_RULES, DEFAULT_BASIS,
} from "./transpiler/transpiler.js";
export { Layout, CouplingMap, AnalysisPass, TransformationPass } from "./transpiler/layout.js";
export {
  SabreLayout, SabreSwap, TrivialLayout, DenseLayout,
  ApplyLayout, BasicSwap, LookaheadSwap, StochasticSwap,
  CommutativeCancellation, CXCancellation, Optimize1qGates, OptimizeSwapBeforeMeasure,
  RemoveBarriers, RemoveResetInZeroState, DAGFixedPointPass,
} from "./transpiler/passes.js";
export {
  Decompose, RemoveFinalMeasurements, BarrierBeforeFinalMeasurements,
  Collect2qBlocks, ConsolidateBlocks, Unroll3qOrMore, MergeAdjacentBarriers,
  CheckMap, GateDirection, CountOps, Depth, Size, Width,
  Optimize1qGatesDecomposition, BasisTranslator,
} from "./transpiler/extra_passes.js";

// Real Sabre routing
export {
  SabreLayout as SabreLayoutReal, SabreSwap as SabreSwapReal,
  collect_1q_runs, collect_2q_runs, consolidate_blocks,
  optimize_cliffords, remove_diagonal_gates_before_measure,
} from "./transpiler/sabre.js";

// Result
export { Result, Counts } from "./result/result.js";

// DAGCircuit + converters + DAG passes
export { DAGCircuit, DAGOpNode, DAGInNode, DAGOutNode, DAGRegister } from "./dagcircuit/dagcircuit.js";
export {
  collect_1q_runs as dag_collect_1q_runs, collect_2q_runs as dag_collect_2q_runs,
  consolidate_blocks as dag_consolidate_blocks, optimize_cliffords as dag_optimize_cliffords,
  remove_diagonal_gates_before_measure as dag_remove_diag_before_measure,
  elide_permutations, remove_redundant_gates, commutative_cancellation, template_optimization,
} from "./dagcircuit/dag_passes.js";
export { circuit_to_dag, dag_to_circuit } from "./converters/converters.js";

// QASM
export { QASMParser } from "./qasm/qasm_parser.js";
export { QASMExporter } from "./qasm/qasm_exporter.js";
export { QASM3Parser, qasm3_parse } from "./qasm/qasm3_parser.js";
export { QASM3Exporter, qasm3_export } from "./qasm/qasm3_exporter.js";

// Primitives
export { Estimator, Sampler, BaseEstimator, BaseSampler } from "./primitives/primitives.js";
export {
  EstimatorV2, SamplerV2,
  PrimitiveResultV2, EstimatorResultV2, SamplerResultV2, PrimitivePubResult,
} from "./primitives/primitives.js";

// Opflow
export {
  PauliSumOp, CircuitStateFn, StateFn, ListOp, OperatorBase,
  PauliOp, MatrixOp, CircuitOp, EvolvedOp,
} from "./opflow/opflow.js";

// Opflow expectations
export {
  MatrixExpectation, PauliExpectation, CircuitSampler, get_expectation,
} from "./opflow/opflow_expec.js";

// Algorithms
export { SPSA, COBYLA, GradientDescent, OptimizerResult } from "./algorithms/optimizers.js";
export { Adam, LBFGSB, SLSQP, NelderMead, NFT } from "./algorithms/optimizers_extra.js";
export { VQE, VQEResult } from "./algorithms/vqe.js";
export { QAOA, QAOAResult } from "./algorithms/qaoa.js";
export { Grover, GroverResult } from "./algorithms/grover.js";
export { PhaseEstimation, PhaseEstimationResult, AmplitudeEstimation, AmplitudeEstimationResult } from "./algorithms/phase_estimation.js";
export { NumPyMinimumEigensolver, NumPyMaximumEigensolver, NumPyMinimumEigensolverResult } from "./algorithms/numpy_eigensolver.js";
export { Shor, HHL, VQC, QSVC, quantumVolumeCircuit, heavyOutputProbability, randomCliffordSequence } from "./algorithms/extra_algorithms.js";

// Gradients
export { GradientBase, ParamShift, FiniteDiff, LinearCombination, NaturalGradient } from "./algorithms/gradients.js";

// Visualization
export {
  draw_circuit, plot_histogram, plot_state_city,
  plot_bloch_vector, plot_state_hinton, plot_state_qsphere,
} from "./visualization/visualization.js";

// ASCII renderers
export {
  render_histogram_ascii, render_bloch_ascii,
  render_state_city_ascii, render_qsphere_ascii,
} from "./visualization/ascii_viz.js";

// Convenience: top-level execute (qiskit.execute equivalent)
import { simulate as _simulate } from "./simulator/statevector_simulator.js";
export function execute(circuit, options = {}) {
  const shots = options.shots || 1024;
  return _simulate(circuit, shots, options);
}
