import { QuantumRegister, ClassicalRegister, Qubit, Clbit } from "./../core/bit.js";
import { CircuitInstruction } from "./../core/circuit.js";

let _nodeCounter = 0;

export class DAGNode {
  constructor() {
    this._id = ++_nodeCounter;
    this._node_id = this._id;
    this.name = "";
    this.wires = []; // qubits/clbits this node touches
  }

  is_op_node() { return false; }
  is_in_node() { return false; }
  is_out_node() { return false; }
}

export class DAGOpNode extends DAGNode {
  constructor(op, qargs = [], cargs = []) {
    super();
    this.op = op;
    this.qargs = qargs;
    this.cargs = cargs;
    this.name = op.name;
    this.wires = qargs.concat(cargs);
  }

  is_op_node() { return true; }

  semantic_eq(other) {
    return other instanceof DAGOpNode &&
      this.op.name === other.op.name &&
      this.qargs.length === other.qargs.length &&
      this.qargs.every((q, i) => q.equals(other.qargs[i]));
  }
}

export class DAGInNode extends DAGNode {
  constructor(wire) {
    super();
    this.wire = wire;
    this.name = wire.toString();
    this.wires = [wire];
  }

  is_in_node() { return true; }
}

export class DAGOutNode extends DAGNode {
  constructor(wire) {
    super();
    this.wire = wire;
    this.name = wire.toString();
    this.wires = [wire];
  }

  is_out_node() { return true; }
}

export class DAGRegister {
  constructor(name, bits, bitClass) {
    this.name = name;
    this.bits = bits;
    this.bitClass = bitClass;
  }
}

export class DAGCircuit {
  constructor() {
    this.qregs = new Map(); // name -> DAGRegister
    this.cregs = new Map();
    this.qubits = [];
    this.clbits = [];
    this._input_nodes = new Map(); // wire -> DAGInNode
    this._output_nodes = new Map(); // wire -> DAGOutNode
    this._multi_graph = new Map(); // node -> { successors: Set, predecessors: Set }
    this._wire_to_in_edges = new Map();
    this._wire_to_out_edges = new Map();
    this.globalPhase = 0;
    this.name = "dag";
    this.metadata = null;
    this.calibrations = [];
    this.duration = null;
    this.unit = "dt";
  }

  addQreg(register) {
    if (!(register instanceof QuantumRegister)) {
      throw new TypeError("addQreg requires a QuantumRegister");
    }
    if (this.qregs.has(register.name)) {
      throw new Error(`Quantum register ${register.name} already exists`);
    }
    const dagReg = new DAGRegister(register.name, register._bits.slice(), Qubit);
    this.qregs.set(register.name, dagReg);
    for (const q of register._bits) {
      this.qubits.push(q);
      const inNode = new DAGInNode(q);
      const outNode = new DAGOutNode(q);
      this._input_nodes.set(q, inNode);
      this._output_nodes.set(q, outNode);
      this._multi_graph.set(inNode, { successors: new Set([outNode]), predecessors: new Set() });
      this._multi_graph.set(outNode, { successors: new Set(), predecessors: new Set([inNode]) });
      this._wire_to_in_edges.set(q, [{ from: inNode, to: outNode }]);
      this._wire_to_out_edges.set(q, []);
    }
    return dagReg;
  }

  addCreg(register) {
    if (!(register instanceof ClassicalRegister)) {
      throw new TypeError("addCreg requires a ClassicalRegister");
    }
    if (this.cregs.has(register.name)) {
      throw new Error(`Classical register ${register.name} already exists`);
    }
    const dagReg = new DAGRegister(register.name, register._bits.slice(), Clbit);
    this.cregs.set(register.name, dagReg);
    for (const c of register._bits) {
      this.clbits.push(c);
      const inNode = new DAGInNode(c);
      const outNode = new DAGOutNode(c);
      this._input_nodes.set(c, inNode);
      this._output_nodes.set(c, outNode);
      this._multi_graph.set(inNode, { successors: new Set([outNode]), predecessors: new Set() });
      this._multi_graph.set(outNode, { successors: new Set(), predecessors: new Set([inNode]) });
      this._wire_to_in_edges.set(c, [{ from: inNode, to: outNode }]);
      this._wire_to_out_edges.set(c, []);
    }
    return dagReg;
  }

  get numQubits() { return this.qubits.length; }
  get numClbits() { return this.clbits.length; }

  // Apply an operation to qubits and clbits
  applyOperation(op, qargs = [], cargs = []) {
    if (!Array.isArray(qargs)) qargs = [qargs];
    if (!Array.isArray(cargs)) cargs = [cargs];
    const node = new DAGOpNode(op, qargs, cargs);
    this._multi_graph.set(node, { successors: new Set(), predecessors: new Set() });

    // For each wire, connect the previous output (the current "end" of the wire's
    // chain) to this new node, and this node to the wire's final output node.
    for (const wire of qargs.concat(cargs)) {
      // Find the current "last producer" on this wire (the most recent op or
      // the input node if no ops yet)
      const inEdges = this._wire_to_in_edges.get(wire);
      if (!inEdges) {
        throw new Error(`Wire ${wire} not in DAG`);
      }
      // Current last edge is the last in the inEdges list
      const lastEdge = inEdges[inEdges.length - 1];
      // Disconnect lastEdge: from -> to, replace with from -> node -> to
      const fromNode = lastEdge.from;
      const toNode = lastEdge.to;
      this._removeEdge(fromNode, toNode);
      this._addEdge(fromNode, node);
      this._addEdge(node, toNode);
      // Update the in-edge list: the new "current end" is from node -> to
      inEdges[inEdges.length - 1] = { from: node, to: toNode };
      // Track this op as a producer on the wire
      this._wire_to_out_edges.get(wire).push({ from: fromNode, to: node });
    }
    return node;
  }

  _addEdge(from, to) {
    this._multi_graph.get(from).successors.add(to);
    this._multi_graph.get(to).predecessors.add(from);
  }

  _removeEdge(from, to) {
    this._multi_graph.get(from).successors.delete(to);
    this._multi_graph.get(to).predecessors.delete(from);
  }

  // Topologically iterate over operation nodes
  topologicalOpNodes() {
    // Kahn's algorithm
    const inDegree = new Map();
    const queue = [];
    for (const [node, links] of this._multi_graph.entries()) {
      let deg = 0;
      // Count op-node predecessors
      for (const pred of links.predecessors) {
        if (pred.is_op_node()) deg++;
      }
      inDegree.set(node, deg);
      if (deg === 0 && node.is_op_node()) queue.push(node);
    }
    const result = [];
    while (queue.length > 0) {
      const node = queue.shift();
      if (node.is_op_node()) result.push(node);
      const links = this._multi_graph.get(node);
      for (const succ of links.successors) {
        const newDeg = inDegree.get(succ) - 1;
        inDegree.set(succ, newDeg);
        if (newDeg === 0 && succ.is_op_node()) queue.push(succ);
      }
    }
    return result;
  }

  // Get all operation nodes (unsorted)
  op_nodes() {
    return Array.from(this._multi_graph.keys()).filter(n => n.is_op_node());
  }

  // Get all nodes (input, output, ops)
  allNodes() {
    return Array.from(this._multi_graph.keys());
  }

  // Count operations by name
  countOps() {
    const counts = {};
    for (const node of this.op_nodes()) {
      counts[node.op.name] = (counts[node.op.name] || 0) + 1;
    }
    return counts;
  }

  // Depth (longest path through op nodes)
  depth() {
    const opNodes = this.topologicalOpNodes();
    const longestPath = new Map();
    for (const node of opNodes) {
      let maxPred = 0;
      const links = this._multi_graph.get(node);
      for (const pred of links.predecessors) {
        if (pred.is_op_node()) {
          const predLen = longestPath.get(pred) || 0;
          if (predLen > maxPred) maxPred = predLen;
        }
      }
      longestPath.set(node, maxPred + 1);
    }
    let max = 0;
    for (const v of longestPath.values()) if (v > max) max = v;
    return max;
  }

  // Get edges on a wire
  edgesOnWire(wire) {
    const result = [];
    const inEdges = this._wire_to_in_edges.get(wire);
    if (!inEdges) return result;
    // Walk the chain from input to output
    let current = this._input_nodes.get(wire);
    const outNode = this._output_nodes.get(wire);
    while (current !== outNode) {
      const links = this._multi_graph.get(current);
      // Find the successor that leads to outNode on this wire
      let next = null;
      for (const succ of links.successors) {
        // Check if this successor eventually reaches outNode
                if (succ.wires.includes(wire) || succ === outNode) {
          next = succ;
          break;
        }
      }
      if (!next) break;
      result.push({ from: current, to: next });
      current = next;
    }
    return result;
  }

  // Remove an operation node
  removeOpNode(node) {
    if (!node.is_op_node()) throw new Error("removeOpNode requires an op node");
    // Reconnect predecessors to successors
    const links = this._multi_graph.get(node);
    const preds = Array.from(links.predecessors);
    const succs = Array.from(links.successors);
    // For each wire, find the matching pred/succ pair
    for (const wire of node.wires) {
      let pred = null, succ = null;
      for (const p of preds) if (p.wires.includes(wire)) { pred = p; break; }
      for (const s of succs) if (s.wires.includes(wire)) { succ = s; break; }
      if (pred && succ) {
        this._removeEdge(pred, node);
        this._removeEdge(node, succ);
        this._addEdge(pred, succ);
        // Update in-edge list
        const inEdges = this._wire_to_in_edges.get(wire);
        for (let i = 0; i < inEdges.length; i++) {
          if (inEdges[i].from === node && inEdges[i].to === succ) {
            inEdges[i] = { from: pred, to: succ };
            break;
          }
        }
      }
    }
    this._multi_graph.delete(node);
  }

  // Substitute one node with another (preserves edges)
  substituteNode(old_node, new_node) {
    if (!old_node.is_op_node()) throw new Error("substituteNode: old must be op node");
    const links = this._multi_graph.get(old_node);
    this._multi_graph.set(new_node, {
      successors: new Set(links.successors),
      predecessors: new Set(links.predecessors),
    });
    // Update predecessors' successor lists
    for (const pred of links.predecessors) {
      const predLinks = this._multi_graph.get(pred);
      predLinks.successors.delete(old_node);
      predLinks.successors.add(new_node);
    }
    // Update successors' predecessor lists
    for (const succ of links.successors) {
      const succLinks = this._multi_graph.get(succ);
      succLinks.predecessors.delete(old_node);
      succLinks.predecessors.add(new_node);
    }
    this._multi_graph.delete(old_node);
  }

  // Replace a node with a subcircuit (DAGCircuit)
  substituteNodeWithDag(node, sub_dag, wires_map = null) {
    // wires_map: { old_wire: new_wire } for the sub-DAG's wires
    if (!node.is_op_node()) throw new Error("substituteNodeWithDag: old must be op node");
    // Map sub-DAG wires to this DAG's wires
    const subQubits = sub_dag.qubits;
    const subClbits = sub_dag.clbits;
    const defaultMap = new Map();
    if (!wires_map) {
      // Assume 1-to-1 mapping in order
      for (let i = 0; i < subQubits.length; i++) defaultMap.set(subQubits[i], node.qargs[i]);
      for (let i = 0; i < subClbits.length; i++) defaultMap.set(subClbits[i], node.cargs[i]);
    } else {
      for (const [k, v] of wires_map.entries()) defaultMap.set(k, v);
    }
    // Apply each op from sub_dag in topological order
    const predLinks = this._multi_graph.get(node).predecessors;
    const succLinks = this._multi_graph.get(node).successors;
    // Get the wire-to-pred and wire-to-succ mapping
    const wirePred = new Map();
    const wireSucc = new Map();
    for (const pred of predLinks) for (const w of pred.wires) wirePred.set(w, pred);
    for (const succ of succLinks) for (const w of succ.wires) wireSucc.set(w, succ);
    // Disconnect node
    for (const w of node.wires) {
      const p = wirePred.get(w);
      const s = wireSucc.get(w);
      if (p && s) {
        this._removeEdge(p, node);
        this._removeEdge(node, s);
      }
    }
    this._multi_graph.delete(node);
    // Apply sub-DAG ops
    let lastNodeOnWire = new Map(wirePred);
    for (const opNode of sub_dag.topologicalOpNodes()) {
      const mappedQargs = opNode.qargs.map(q => defaultMap.get(q) || q);
      const mappedCargs = opNode.cargs.map(c => defaultMap.get(c) || c);
      const newNode = new DAGOpNode(opNode.op, mappedQargs, mappedCargs);
      this._multi_graph.set(newNode, { successors: new Set(), predecessors: new Set() });
      for (const w of mappedQargs.concat(mappedCargs)) {
        const pred = lastNodeOnWire.get(w) || this._input_nodes.get(w);
        this._addEdge(pred, newNode);
        lastNodeOnWire.set(w, newNode);
      }
    }
    // Connect last nodes to successors
    for (const w of node.wires) {
      const last = lastNodeOnWire.get(w) || wirePred.get(w);
      const succ = wireSucc.get(w);
      if (last && succ) this._addEdge(last, succ);
    }
  }

  // Get node predecessors
  predecessors(node) {
    const links = this._multi_graph.get(node);
    return links ? Array.from(links.predecessors) : [];
  }

  // Get node successors
  successors(node) {
    const links = this._multi_graph.get(node);
    return links ? Array.from(links.successors) : [];
  }

  // Get quantum register by name
  qreg(name) {
    return this.qregs.get(name);
  }

  creg(name) {
    return this.cregs.get(name);
  }

  // Copy the DAG
  copy() {
    const newDag = new DAGCircuit();
    for (const [name, reg] of this.qregs.entries()) {
      newDag.addQreg(new QuantumRegister(reg.bits.length, name));
    }
    for (const [name, reg] of this.cregs.entries()) {
      newDag.addCreg(new ClassicalRegister(reg.bits.length, name));
    }
    newDag.name = this.name;
    newDag.globalPhase = this.globalPhase;
    newDag.metadata = this.metadata;
    // Copy op nodes in topological order
    for (const node of this.topologicalOpNodes()) {
      newDag.applyOperation(node.op.copy(), node.qargs.slice(), node.cargs.slice());
    }
    return newDag;
  }

  // Convert DAG to a list of CircuitInstructions (in topological order)
  toInstructions() {
    return this.topologicalOpNodes().map(node => new CircuitInstruction(
      node.op.copy(),
      node.qargs.slice(),
      node.cargs.slice()
    ));
  }
}
