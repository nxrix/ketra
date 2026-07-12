/**
 * layout.js - Layout, CouplingMap, and transpiler pass base classes.
 *
 * Layout, CouplingMap, AnalysisPass,
 * TransformationPass.
 *
 * A Layout maps virtual (circuit) qubits to physical (backend) qubits.
 * A CouplingMap describes which physical qubit pairs can have 2-qubit gates.
 */

export class Layout {
  constructor(mapping = null) {
    // mapping: { physical: virtual } or { virtual: physical }
    // Internally store as two maps for O(1) both directions.
    this._p2v = new Map(); // physical -> virtual
    this._v2p = new Map(); // virtual -> physical
    if (mapping) {
      for (const [k, v] of Object.entries(mapping)) {
        // Treat keys as physical (numbers), values as virtual (any)
        const phys = parseInt(k, 10);
        if (!isNaN(phys)) {
          this.setPhysical(phys, v);
        } else {
          // Virtual -> physical mapping
          this.setVirtual(k, parseInt(v, 10));
        }
      }
    }
  }

  static trivial(numQubits) {
    const layout = new Layout();
    for (let i = 0; i < numQubits; i++) {
      layout.setPhysical(i, i);
    }
    return layout;
  }

  static from_dict(d) {
    const layout = new Layout();
    for (const [k, v] of Object.entries(d)) {
      layout.setPhysical(parseInt(k, 10), v);
    }
    return layout;
  }

  static generate_trivial(numQubits, qubits) {
    return Layout.trivial(numQubits);
  }

  static generate_from_integers(intList) {
    const layout = new Layout();
    intList.forEach((v, p) => layout.setPhysical(p, v));
    return layout;
  }

  setPhysical(phys, virtual) {
    // Remove old mapping
    if (this._p2v.has(phys)) {
      const oldVirt = this._p2v.get(phys);
      this._v2p.delete(oldVirt);
    }
    if (this._v2p.has(virtual)) {
      const oldPhys = this._v2p.get(virtual);
      this._p2v.delete(oldPhys);
    }
    this._p2v.set(phys, virtual);
    this._v2p.set(virtual, phys);
  }

  setVirtual(virtual, phys) { this.setPhysical(phys, virtual); }

  getPhysical(phys) { return this._p2v.get(phys); }
  getVirtual(virtual) { return this._v2p.get(virtual); }

  // physical -> virtual mapping (object form)
  get_physical_bits() {
    const out = {};
    for (const [k, v] of this._p2v.entries()) out[k] = v;
    return out;
  }

  get_virtual_bits() {
    const out = {};
    for (const [k, v] of this._v2p.entries()) out[k] = v;
    return out;
  }

  // Get all physical qubits in use
  get_physical_bits_list() {
    return Array.from(this._p2v.keys()).sort((a, b) => a - b);
  }

  // Add an ancilla qubit
  add(phys, virtual) {
    this.setPhysical(phys, virtual);
  }

  // Swap two physical qubits' mappings
  swap(phys1, phys2) {
    const v1 = this._p2v.get(phys1);
    const v2 = this._p2v.get(phys2);
    this.setPhysical(phys1, v2);
    this.setPhysical(phys2, v1);
  }

  // Combine two layouts
  combine(other) {
    const result = new Layout();
    for (const [p, v] of this._p2v.entries()) result.setPhysical(p, v);
    for (const [p, v] of other._p2v.entries()) result.setPhysical(p, v);
    return result;
  }

  copy() {
    const result = new Layout();
    for (const [p, v] of this._p2v.entries()) result.setPhysical(p, v);
    return result;
  }

  size() { return this._p2v.size; }

  to_dict() { return this.get_physical_bits(); }

  toString() {
    const pairs = [];
    for (const [p, v] of this._p2v.entries()) pairs.push(`${p}: ${v}`);
    return `Layout({${pairs.join(", ")}})`;
  }
}

export class CouplingMap {
  constructor(connections = null, numQubits = null) {
    // connections: list of [a, b] directed edges (a -> b)
    this._edges = [];
    this._graph = new Map(); // node -> Set of neighbors
    this._numQubits = numQubits;

    if (connections) {
      if (typeof connections === "string") {
        // Parse coupling map string like "0-1,1-2,2-3"
        const parts = connections.split(",");
        for (const part of parts) {
          const [a, b] = part.trim().split("-").map(s => parseInt(s.trim(), 10));
          this.addEdge(a, b);
        }
      } else if (Array.isArray(connections)) {
        for (const [a, b] of connections) this.addEdge(a, b);
      } else if (connections instanceof Map) {
        for (const [k, neighbors] of connections.entries()) {
          for (const n of neighbors) this.addEdge(k, n);
        }
      } else {
        throw new TypeError("CouplingMap: invalid connections argument");
      }
    }
  }

  static fromLine(numQubits) {
    // Linear chain: 0-1-2-...-n-1
    const c = new CouplingMap(null, numQubits);
    for (let i = 0; i + 1 < numQubits; i++) {
      c.addEdge(i, i + 1);
      c.addEdge(i + 1, i);
    }
    return c;
  }

  static fromRing(numQubits) {
    const c = new CouplingMap(null, numQubits);
    for (let i = 0; i < numQubits; i++) {
      c.addEdge(i, (i + 1) % numQubits);
      c.addEdge((i + 1) % numQubits, i);
    }
    return c;
  }

  static fromGrid(rows, cols) {
    const n = rows * cols;
    const c = new CouplingMap(null, n);
    for (let r = 0; r < rows; r++) {
      for (let col = 0; col < cols; col++) {
        const i = r * cols + col;
        if (col + 1 < cols) {
          c.addEdge(i, i + 1);
          c.addEdge(i + 1, i);
        }
        if (r + 1 < rows) {
          c.addEdge(i, i + cols);
          c.addEdge(i + cols, i);
        }
      }
    }
    return c;
  }

  static fullyConnected(numQubits) {
    const c = new CouplingMap(null, numQubits);
    for (let i = 0; i < numQubits; i++) {
      for (let j = 0; j < numQubits; j++) {
        if (i !== j) c.addEdge(i, j);
      }
    }
    return c;
  }

  addEdge(a, b) {
    this._edges.push([a, b]);
    if (!this._graph.has(a)) this._graph.set(a, new Set());
    if (!this._graph.has(b)) this._graph.set(b, new Set());
    this._graph.get(a).add(b);
    if (this._numQubits == null) {
      this._numQubits = Math.max(a, b) + 1;
    } else {
      this._numQubits = Math.max(this._numQubits, a + 1, b + 1);
    }
  }

  get edges() { return this._edges.slice(); }
  get size() { return this._numQubits; }
  get num_qubits() { return this._numQubits; }

  neighbors(node) {
    return Array.from(this._graph.get(node) || []);
  }

  hasEdge(a, b) {
    return this._graph.has(a) && this._graph.get(a).has(b);
  }

  is_connected() {
    if (this._numQubits <= 1) return true;
    // BFS
    const visited = new Set([0]);
    const queue = [0];
    while (queue.length > 0) {
      const n = queue.shift();
      for (const m of this.neighbors(n)) {
        if (!visited.has(m)) {
          visited.add(m);
          queue.push(m);
        }
      }
    }
    return visited.size === this._numQubits;
  }

  // Shortest path between two physical qubits (BFS)
  shortestPath(a, b) {
    if (a === b) return [a];
    const visited = new Set([a]);
    const queue = [[a]];
    while (queue.length > 0) {
      const path = queue.shift();
      const node = path[path.length - 1];
      for (const next of this.neighbors(node)) {
        if (visited.has(next)) continue;
        if (next === b) return path.concat([next]);
        visited.add(next);
        queue.push(path.concat([next]));
      }
    }
    return null;
  }

  // Distance between two physical qubits
  distance(a, b) {
    const path = this.shortestPath(a, b);
    return path ? path.length - 1 : Infinity;
  }

  // Get the undirected graph (each edge once)
  get_undirected_edges() {
    const seen = new Set();
    const result = [];
    for (const [a, b] of this._edges) {
      const key = a < b ? `${a}-${b}` : `${b}-${a}`;
      if (!seen.has(key)) {
        seen.add(key);
        result.push([Math.min(a, b), Math.max(a, b)]);
      }
    }
    return result;
  }

  copy() {
    const c = new CouplingMap(null, this._numQubits);
    for (const [a, b] of this._edges) c.addEdge(a, b);
    return c;
  }

  to_dict() {
    const out = {};
    for (const [k, neighbors] of this._graph.entries()) {
      out[k] = Array.from(neighbors);
    }
    return out;
  }
}

// ---------------------------------------------------------------------------
// Pass base classes
// ---------------------------------------------------------------------------
export class AnalysisPass {
  constructor() {
    this.requires = [];
    this.preserves = [];
    this.analysis_name = this.constructor.name;
    this.property_set = {};
  }

  run(dag) {
    throw new Error("AnalysisPass.run not implemented");
  }
}

export class TransformationPass {
  constructor() {
    this.requires = [];
    this.preserves = [];
    this.transformation_name = this.constructor.name;
  }

  run(dag) {
    throw new Error("TransformationPass.run not implemented");
  }
}
