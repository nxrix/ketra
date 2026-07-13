import { _setDrawHook } from "./../core/circuit.js";

const CONTROLLED_2Q = ["cx", "cy", "cz", "ch", "csx", "crx", "cry", "crz", "cp", "cu1", "cu3", "cu"];
const SYMMETRIC_2Q = ["swap", "iswap", "dcx", "rxx", "ryy", "rzz", "rzx"];

export function drawCircuit(circuit, output = "text", kwargs) {
  if (output === "text" || output === "ascii") return _drawText(circuit);
  if (output === "latex") return _drawLatexSource(circuit);
  if (output === "html") return _drawHtml(circuit);
  return _drawText(circuit);
}

// Track which cells have been explicitly set by drawing functions
const _touchedCells = new Set();

function _setCell(grid, row, col, char) {
  if (!grid[row]) return;
  while (grid[row].length <= col) grid[row].push(" ");
  grid[row][col] = char;
  _touchedCells.add(`${row},${col}`);
}

// Text circuit drawer — matches qiskit's text drawer format exactly

function _drawText(circuit) {
  const nq = circuit.numQubits;
  const nc = circuit.numClbits;
  _touchedCells.clear();

  // Build instruction list with resolved qubit/clbit indices
  const allInstructions = [];
  for (let idx = 0; idx < circuit.data.length; idx++) {
    const ci = circuit.data[idx];
    const qubits = ci.qubits.map(q => circuit._qubit_index.get(q));
    const clbits = ci.clbits.map(c => circuit._clbit_index.get(c));
    if (qubits.some(q => q === undefined) || clbits.some(c => c === undefined)) continue;
    allInstructions.push({ op: ci.operation, qubits, clbits, idx });
  }

  // DAG LAYER ASSIGNMENT (matches qiskit's _LayerSpooler)
  // Step 1: Compute DAG layers (topological layers via Kahn's algorithm)
  const dagLayers = [];
  const wireBusyUntil = new Map(); // wire -> last layer index
  for (const instr of allInstructions) {
    const allWires = instr.qubits.concat(instr.clbits);
    let layer = 0;
    for (const w of allWires) {
      if (wireBusyUntil.has(w)) layer = Math.max(layer, wireBusyUntil.get(w) + 1);
    }
    while (dagLayers.length <= layer) dagLayers.push([]);
    dagLayers[layer].push(instr);
    for (const w of allWires) wireBusyUntil.set(w, layer);
  }

  // Step 2: Slide each node left (qiskit's slide_from_left algorithm)
  // Each node tries to slide into the earliest existing layer where it doesn't
  // crossover with any existing node. Crossover = overlapping qubit spans.
  function _gateSpan(instr) {
    if (instr.qubits.length === 0) return [];
    const minQ = Math.min(...instr.qubits);
    const maxQ = Math.max(...instr.qubits);
    const span = [];
    for (let i = minQ; i <= maxQ; i++) span.push(i);
    // For measurements, also include classical wire span
    if (instr.clbits.length > 0 && instr.op.name === "measure") {
      // Measurement span includes all qubits from the measured qubit down to
      // the classical wire (so it crosses all qubits below it)
      for (let i = maxQ + 1; i < nq; i++) span.push(i);
    }
    return span;
  }

  function _anyCrossover(instr, layerNodes) {
    const span = new Set(_gateSpan(instr));
    for (const existing of layerNodes) {
      const existingSpan = _gateSpan(existing);
      for (const q of existingSpan) {
        if (span.has(q)) return true;
      }
    }
    return false;
  }

  // Slide nodes left
  const layers = []; // final layers after sliding
  for (const dagLayer of dagLayers) {
    const currentIdx = layers.length - 1;
    for (const node of dagLayer) {
      let inserted = false;
      // Try to slide left from currentIdx down to 0
      let lastInsertable = -1;
      for (let i = currentIdx; i >= 0; i--) {
        // Check if any qubit is already used in this layer
        const nodeWires = new Set(node.qubits.concat(node.clbits));
        let found = false;
        for (const existing of layers[i]) {
          for (const w of existing.qubits.concat(existing.clbits)) {
            if (nodeWires.has(w)) { found = true; break; }
          }
          if (found) break;
        }
        if (found) break;
        if (!_anyCrossover(node, layers[i])) {
          lastInsertable = i;
        }
      }
      if (lastInsertable >= 0) {
        layers[lastInsertable].push(node);
        inserted = true;
      } else {
        // Try from currentIdx onwards
        for (let i = Math.max(0, currentIdx); i < layers.length; i++) {
          if (!_anyCrossover(node, layers[i])) {
            const nodeWires = new Set(node.qubits.concat(node.clbits));
            let conflict = false;
            for (const existing of layers[i]) {
              for (const w of existing.qubits.concat(existing.clbits)) {
                if (nodeWires.has(w)) { conflict = true; break; }
              }
              if (conflict) break;
            }
            if (!conflict) {
              layers[i].push(node);
              inserted = true;
              break;
            }
          }
        }
      }
      if (!inserted) {
        layers.push([node]);
      }
    }
  }

  // LAYER WIDTH (matches qiskit's normalize_width: MAX of all widths)
  function _instrWidth(instr) {
    if (instr.op.name === "barrier") return 1;
    if (instr.op.name === "measure") return 3;
    if (instr.op.name === "reset") return 5;
    if (instr.op.name === "swap" || instr.op.name === "cswap" || instr.op.name === "cz") return 3;
    return Math.max(3, _gateLabel(instr.op).length) + 2;
  }

  const layerWidths = layers.map(layer => {
    let maxW = 0;
    for (const instr of layer) {
      const w = _instrWidth(instr);
      if (w > maxW) maxW = w;
    }
    return maxW === 0 ? 5 : maxW;
  });

  // GRID DRAWING
  const nRows = nc > 0 ? 2 * nq + 3 : 2 * nq + 1;
  const labelWidth = _maxLabelWidth(nq, nc);
  const grid = [];
  for (let r = 0; r < nRows; r++) grid.push([]);
  _addLabels(grid, nq, nc, labelWidth);

  let colOffset = labelWidth;
  for (let li = 0; li < layers.length; li++) {
    const layer = layers[li];
    const width = layerWidths[li];
    for (let r = 0; r < nRows; r++) {
      while (grid[r].length < colOffset) grid[r].push(" ");
    }
    _drawColumn(grid, layer, width, nq, nc, colOffset);

    // Fill empty wires for qubits not in this layer
    const layerQubits = new Set();
    for (const instr of layer) { for (const q of instr.qubits) layerQubits.add(q); }
    for (let qi = 0; qi < nq; qi++) {
      if (!layerQubits.has(qi)) {
        const wireRow = 2 * qi + 1;
        for (let c = 0; c < width; c++) {
          if (!_touchedCells.has(`${wireRow},${colOffset + c}`)) {
            _setCell(grid, wireRow, colOffset + c, "\u2500");
          }
        }
      }
    }
    // Fill classical wire if no measurement
    const hasMeasure = layer.some(i => i.op.name === "measure");
    if (nc > 0 && !hasMeasure) {
      const cWireRow = 2 * nq + 1;
      for (let c = 0; c < width; c++) {
        if (!_touchedCells.has(`${cWireRow},${colOffset + c}`)) {
          _setCell(grid, cWireRow, colOffset + c, "\u2550");
        }
      }
    }
    colOffset += width;
  }

  // Extend all rows to the same length
  let maxLen = 0;
  for (let r = 0; r < nRows; r++) { if (grid[r].length > maxLen) maxLen = grid[r].length; }
  // Qiskit: wire rows get 1 trailing wire char, non-wire rows get spaces
  const targetLen = maxLen;
  for (let r = 0; r < nRows; r++) {
    while (grid[r].length < targetLen) grid[r].push(" ");
  }
  _fillUntouchedWires(grid, nq, nc, labelWidth, targetLen);
  // No trailing wire chars — qiskit pads all rows to the same length with spaces.
  // Wire rows that already end with wire chars keep them; rows that end with box
  // borders get spaces.
  for (let r = 0; r < nRows; r++) {
    while (grid[r].length < targetLen + 1) grid[r].push(" ");
  }

  // BORDER COMPRESSION (matches qiskit's merge_lines + should_compress)
  // Qiskit merges adjacent border rows: └───┘ + ┌───┐ → ├───┤
  // But does NOT compress when ┴ is above ┬/╥ (measurement/control connectors)
  // Also does NOT compress quantum borders (─) with classical wire (═) or vice versa
  // And does NOT compress the last qubit bottom border with the classical wire row
  const lines = grid.map(row => row.join(""));
  const merged = [];
  for (let r = 0; r < lines.length; r++) {
    if (merged.length === 0) { merged.push(lines[r]); continue; }
    const top = merged[merged.length - 1];
    const bot = lines[r];
    // Check if we should compress (merge) these two lines
    let shouldCompress = true;
    for (let c = 0; c < Math.min(top.length, bot.length); c++) {
      const tc = top[c], bc = bot[c];
      // Don't compress ┴ above ┬/╥ (qiskit's should_compress rule)
      if (tc === "\u2534" && (bc === "\u252C" || bc === "\u2565")) { shouldCompress = false; break; }
      if (tc === "\u2568" && (bc === "\u252C" || bc === "\u2565")) { shouldCompress = false; break; }
      // Don't compress if alphanumeric chars overlap
      if ((tc >= "a" && tc <= "z" || tc >= "A" && tc <= "Z" || tc >= "0" && tc <= "9") && bc !== " ") { shouldCompress = false; break; }
      if ((bc >= "a" && bc <= "z" || bc >= "A" && bc <= "Z" || bc >= "0" && bc <= "9") && tc !== " ") { shouldCompress = false; break; }
      // Don't compress quantum wire char (─) with classical wire char (═)
      // These should stay on separate lines
      if (tc === "\u2500" && bc === "\u2550") { shouldCompress = false; break; }
      // Don't compress ┘/└/┤/├ (quantum box borders) with ═ (classical wire)
      if ((tc === "\u2518" || tc === "\u2514" || tc === "\u2524" || tc === "\u251C") && bc === "\u2550") { shouldCompress = false; break; }
    }
    if (shouldCompress) {
      merged[merged.length - 1] = _mergeLines(top, bot, "top");
    } else {
      merged.push(bot);
    }
  }

  return merged.join("\n");
}

// Merge two lines (qiskit's merge_lines algorithm)
function _mergeLines(top, bot, icod) {
  icod = icod || "top";
  let ret = "";
  const len = Math.max(top.length, bot.length);
  for (let c = 0; c < len; c++) {
    const tc = c < top.length ? top[c] : " ";
    const bc = c < bot.length ? bot[c] : " ";
    if (tc === bc) { ret += tc; }
    else if ((tc === "\u253C" || tc === "\u256A") && bc === " ") { ret += "\u2502"; } // ┼/╫ + space → │
    else if (tc === " ") { ret += bc; }
    else if ((tc === "\u252C" || tc === "\u2565") && (bc === " " || bc === "\u2551" || bc === "\u2502") && icod === "top") { ret += tc; } // ┬/╥ + space/║/│ → ┬/╥
    else if (tc === "\u252C" && bc === " " && icod === "bot") { ret += "\u2502"; } // ┬ + space → │
    else if (tc === "\u2565" && bc === " " && icod === "bot") { ret += "\u2551"; } // ╥ + space → ║
    else if ((tc === "\u252C" || tc === "\u2502") && bc === "\u2550") { ret += "\u256A"; } // ┬/│ + ═ → ╫
    else if ((tc === "\u252C" || tc === "\u2502") && bc === "\u2500") { ret += "\u253C"; } // ┬/│ + ─ → ┼
    else if ((tc === "\u2514" || tc === "\u2518" || tc === "\u2551" || tc === "\u2502" || tc === "\u2591") && bc === " " && icod === "top") { ret += tc; } // └┘║│░
    else if ((tc === "\u2500" || tc === "\u2550") && bc === " " && icod === "top") { ret += tc; } // ─═
    else if ((tc === "\u2500" || tc === "\u2550") && bc === " " && icod === "bot") { ret += bc; }
    else if ((tc === "\u2551" || tc === "\u2565") && bc === "\u2550") { ret += "\u256C"; } // ║/╥ + ═ → ╬
    else if ((tc === "\u2551" || tc === "\u2565") && bc === "\u2500") { ret += "\u256B"; } // ║/╥ + ─ → ╫
    else if ((tc === "\u2551" || tc === "\u256B" || tc === "\u256C") && bc === " ") { ret += "\u2551"; } // ║/╫/╬ + space → ║
    else if ((tc === "\u2502" || tc === "\u253C" || tc === "\u256A") && bc === " ") { ret += "\u2502"; } // │/┼/╫ + space → │
    else if (tc === "\u2514" && bc === "\u250C" && icod === "top") { ret += "\u251C"; } // └ + ┌ → ├
    else if (tc === "\u2518" && bc === "\u2510" && icod === "top") { ret += "\u2524"; } // ┘ + ┐ → ┤
    else if ((bc === "\u250C" || bc === "\u2510") && icod === "top") { ret += "\u252C"; } // ┌/┐ → ┬
    else if ((tc === "\u2518" || tc === "\u2514") && bc === "\u2500" && icod === "top") { ret += "\u2534"; } // ┘/└ + ─ → ┴
    else if (bc === " " && icod === "top") { ret += tc; }
    else { ret += bc; }
  }
  return ret;
}

// Qiskit uses mixed-case for parameterized gates: Rx, Ry, Rz, P, etc.
// Standard gates are uppercase: H, X, Y, Z, S, T, SX, etc.
function _gateLabel(op) {
  const name = op.name;
  if (name === "measure") return "M";
  if (name === "reset") return "|0>";
  if (name === "barrier") return "";
  if (name === "delay") return `D${op.params[0] || 0}`;

  // Mixed-case parameterized gates (match qiskit's naming)
  const mixedCase = {
    "rx": "Rx", "ry": "Ry", "rz": "Rz",
    "rxx": "Rxx", "ryy": "Ryy", "rzz": "Rzz", "rzx": "Rzx",
    "crx": "CRx", "cry": "CRy", "crz": "CRz",
    "p": "P", "cp": "P", "u1": "U1", "u2": "U2", "u3": "U3", "u": "U",
    "cu1": "CU1", "cu3": "CU3", "cu": "CU",
    "iswap": "Iswap", "dcx": "Dcx",
    "sdg": "Sdg", "tdg": "Tdg", "sxdg": "Sxdg",
  };

  const params = op.params.length
    ? `(${op.params.map(_fmtParam).join(",")})`
    : "";

  const baseName = mixedCase[name] || name.toUpperCase();
  return `${baseName}${params}`;
}

// For controlled gates, return the label to show on the target qubit.
// CX -> "X", CY -> "Y", CZ -> (special, no box), CH -> "H", CRX -> "Rx(0.5)"
function _controlledTargetLabel(name, op) {
  const paramMap = {
    "cx": "X", "cy": "Y", "cz": null, "ch": "H", "csx": "SX",
    "crx": "Rx", "cry": "Ry", "crz": "Rz",
    "cp": "P", "cu1": "U1", "cu3": "U3", "cu": "U",
  };
  const base = paramMap[name];
  if (base === null) return null; // CZ: no box
  if (!base) return name.toUpperCase();
  // For parameterized controlled gates, include params
  if (op.params.length > 0 && ["crx", "cry", "crz", "cp", "cu1", "cu3", "cu"].includes(name)) {
    return `${base}(${op.params.map(_fmtParam).join(",")})`;
  }
  return base;
}

function _fmtParam(p) {
  if (typeof p === "number") {
    // Qiskit uses ~4 decimal places for floats
    if (Number.isInteger(p)) return p.toString();
    // Round to 4 decimal places and strip trailing zeros
    const rounded = Math.round(p * 10000) / 10000;
    let str = rounded.toFixed(4);
    // Remove trailing zeros after decimal point
    str = str.replace(/\.?0+$/, "");
    return str;
  }
  if (p && typeof p.toString === "function") return p.toString();
  return String(p);
}

function _maxLabelWidth(nq, nc) {
  let maxW = 0;
  for (let i = 0; i < nq; i++) {
    maxW = Math.max(maxW, `q_${i}: `.length);
  }
  if (nc > 0) {
    maxW = Math.max(maxW, `c: ${nc}/`.length);
  }
  return maxW;
}

function _addLabels(grid, nq, nc, labelWidth) {
  for (let i = 0; i < nq; i++) {
    const wireRow = 2 * i + 1;
    const label = `q_${i}: `;
    _setRow(grid, wireRow, label, labelWidth);
    if (wireRow > 0) _setRow(grid, wireRow - 1, "", labelWidth);
    if (wireRow + 1 < grid.length) _setRow(grid, wireRow + 1, "", labelWidth);
  }
  if (nc > 0) {
    const cWireRow = 2 * nq + 1;
    const label = `c: ${nc}/`;
    _setRow(grid, cWireRow, label, labelWidth);
    _setRow(grid, cWireRow - 1, "", labelWidth);
    if (cWireRow + 1 < grid.length) _setRow(grid, cWireRow + 1, "", labelWidth);
  }
}

function _setRow(grid, row, text, width) {
  while (grid[row].length < width) {
    grid[row].push(text[grid[row].length] || " ");
  }
}

// Draw a single sub-column of instructions (all at the same offset)
function _drawColumn(grid, col, width, nq, nc, colOffset) {
  const measurements = col.filter(i => i.op.name === "measure");
  const otherInstrs = col.filter(i => i.op.name !== "measure" && i.op.name !== "barrier");

  // Collect single-qubit box info
  const boxQubits = new Set();
  const boxLabels = new Map();
  for (const instr of otherInstrs) {
    if (instr.op.numQubits === 1 && instr.op.name !== "reset") {
      boxQubits.add(instr.qubits[0]);
      boxLabels.set(instr.qubits[0], _gateLabel(instr.op));
    } else if (instr.op.name === "reset") {
      boxQubits.add(instr.qubits[0]);
      boxLabels.set(instr.qubits[0], "|0>");
    }
  }

  // Collect controlled gate info
  const controlLines = [];
  const targetBoxes = new Set();
  for (const instr of otherInstrs) {
    if (CONTROLLED_2Q.includes(instr.op.name) || instr.op.name === "ccx") {
      const centerCol = colOffset + Math.floor(width / 2);
      if (instr.op.name === "ccx") {
        controlLines.push({ from: Math.min(...instr.qubits), to: Math.max(...instr.qubits), centerCol });
        targetBoxes.add(instr.qubits[2]);
      } else {
        controlLines.push({ from: Math.min(instr.qubits[0], instr.qubits[1]), to: Math.max(instr.qubits[0], instr.qubits[1]), centerCol });
        if (instr.op.name !== "cz") targetBoxes.add(instr.qubits[1]);
      }
    }
  }

  // Draw single-qubit boxes first
  _drawBoxesWithSharedBorders(grid, boxQubits, boxLabels, width, colOffset, nq);

  // Draw multi-qubit gates
  for (const instr of otherInstrs) {
    if (instr.op.numQubits >= 2) {
      _drawInstruction(grid, instr, width, nq, nc, colOffset);
    }
  }

  // Draw measurements LAST (so they can merge with existing box borders)
  if (measurements.length > 0) {
    _drawMeasurementLayer(grid, measurements, width, nq, nc, colOffset);
  }

  // Draw control line crossings through single-qubit boxes
  for (const line of controlLines) {
    for (let q = line.from + 1; q < line.to; q++) {
      const wireRow = 2 * q + 1;
      if (boxQubits.has(q) && !targetBoxes.has(q)) {
        _setCell(grid, wireRow, line.centerCol, "\u253C");
      }
    }
  }
}

function _drawLayer(grid, layer, width, nq, nc, colOffset) {
  const isBarrierLayer = layer.length > 0 && layer.every(i => i.op.name === "barrier");
  if (isBarrierLayer) {
    _drawBarrierLayer(grid, layer, width, nq, nc, colOffset);
    return;
  }

  const measurements = layer.filter(i => i.op.name === "measure");
  const otherInstrs = layer.filter(i => i.op.name !== "measure" && i.op.name !== "barrier");

  // Collect single-qubit box info
  const boxQubits = new Set();
  const boxLabels = new Map();
  for (const instr of otherInstrs) {
    if (instr.op.numQubits === 1 && instr.op.name !== "reset") {
      boxQubits.add(instr.qubits[0]);
      boxLabels.set(instr.qubits[0], _gateLabel(instr.op));
    } else if (instr.op.name === "reset") {
      boxQubits.add(instr.qubits[0]);
      boxLabels.set(instr.qubits[0], "|0>");
    }
  }

  // Collect controlled gate info
  const controlLines = [];
  const targetBoxes = new Set();
  for (const instr of otherInstrs) {
    if (CONTROLLED_2Q.includes(instr.op.name) || instr.op.name === "ccx") {
      const centerCol = colOffset + Math.floor(width / 2);
      if (instr.op.name === "ccx") {
        controlLines.push({ from: Math.min(instr.qubits[0], instr.qubits[1], instr.qubits[2]), to: Math.max(instr.qubits[0], instr.qubits[1], instr.qubits[2]), centerCol });
        targetBoxes.add(instr.qubits[2]);
      } else {
        controlLines.push({ from: Math.min(instr.qubits[0], instr.qubits[1]), to: Math.max(instr.qubits[0], instr.qubits[1]), centerCol });
        if (instr.op.name !== "cz") targetBoxes.add(instr.qubits[1]);
      }
    }
  }

  // Draw single-qubit boxes with shared borders
  _drawBoxesWithSharedBorders(grid, boxQubits, boxLabels, width, colOffset, nq);

  // Draw multi-qubit gates
  for (const instr of otherInstrs) {
    if (instr.op.numQubits >= 2) {
      _drawInstruction(grid, instr, width, nq, nc, colOffset);
    }
  }

  // Draw control line crossings through single-qubit boxes
  for (const line of controlLines) {
    for (let q = line.from + 1; q < line.to; q++) {
      const wireRow = 2 * q + 1;
      if (boxQubits.has(q) && !targetBoxes.has(q)) {
        _setCell(grid, wireRow, line.centerCol, "\u253C"); // ┼
      }
    }
  }

  // Draw measurements
  if (measurements.length > 0) {
    _drawMeasurementLayer(grid, measurements, width, nq, nc, colOffset);
  }

  // Fill empty wire segments
  const touchedQubits = new Set();
  for (const instr of layer) {
    for (const q of instr.qubits) touchedQubits.add(q);
  }
  for (let i = 0; i < nq; i++) {
    if (!touchedQubits.has(i)) {
      const wireRow = 2 * i + 1;
      for (let c = 0; c < width; c++) {
        if (!_touchedCells.has(`${wireRow},${colOffset + c}`)) {
          _setCell(grid, wireRow, colOffset + c, "\u2500");
        }
      }
    }
  }

  if (nc > 0 && measurements.length === 0) {
    const cWireRow = 2 * nq + 1;
    for (let c = 0; c < width; c++) {
      if (!_touchedCells.has(`${cWireRow},${colOffset + c}`)) {
        _setCell(grid, cWireRow, colOffset + c, "\u2550");
      }
    }
  }
}

// Draw single-qubit boxes with shared borders between vertically adjacent boxes.
// Qiskit uses ├───┤ for shared borders, ┌───┐ for top of a group, └───┘ for bottom.
function _drawBoxesWithSharedBorders(grid, boxQubits, boxLabels, width, colOffset, nq) {
  if (boxQubits.size === 0) return;

  const innerW = Math.max(3, width - 2);
  const sortedQubits = Array.from(boxQubits).sort((a, b) => a - b);

  for (const q of sortedQubits) {
    const wireRow = 2 * q + 1;
    const topRow = wireRow - 1;
    const botRow = wireRow + 1;
    const label = boxLabels.get(q);
    const padL = Math.floor((innerW - label.length) / 2);
    const padR = innerW - label.length - padL;

    const hasBoxAbove = boxQubits.has(q - 1);
    const hasBoxBelow = boxQubits.has(q + 1);

    // Top border
    if (hasBoxAbove) {
      _setCell(grid, topRow, colOffset, "\u251C"); // ├
    } else {
      _setCell(grid, topRow, colOffset, "\u250C"); // ┌
    }
    for (let c = 0; c < innerW; c++) _setCell(grid, topRow, colOffset + 1 + c, "\u2500");
    if (hasBoxAbove) {
      _setCell(grid, topRow, colOffset + 1 + innerW, "\u2524"); // ┤
    } else {
      _setCell(grid, topRow, colOffset + 1 + innerW, "\u2510"); // ┐
    }

    // Wire row: ┤ label ├
    _setCell(grid, wireRow, colOffset, "\u2524"); // ┤
    for (let c = 0; c < padL; c++) _setCell(grid, wireRow, colOffset + 1 + c, " ");
    for (let c = 0; c < label.length; c++) _setCell(grid, wireRow, colOffset + 1 + padL + c, label[c]);
    for (let c = 0; c < padR; c++) _setCell(grid, wireRow, colOffset + 1 + padL + label.length + c, " ");
    _setCell(grid, wireRow, colOffset + 1 + innerW, "\u251C"); // ├

    // Bottom border
    if (hasBoxBelow) {
      _setCell(grid, botRow, colOffset, "\u251C"); // ├
    } else {
      _setCell(grid, botRow, colOffset, "\u2514"); // └
    }
    for (let c = 0; c < innerW; c++) _setCell(grid, botRow, colOffset + 1 + c, "\u2500");
    if (hasBoxBelow) {
      _setCell(grid, botRow, colOffset + 1 + innerW, "\u2524"); // ┤
    } else {
      _setCell(grid, botRow, colOffset + 1 + innerW, "\u2518"); // ┘
    }
  }
}

// Draw a layer of measurements.
// Each measurement gets its own column within the layer. Measurements on
// adjacent qubits can have their vertical lines pass through each other's
// boxes, using ┬ on top borders and ╫ on wire rows.
function _drawMeasurementLayer(grid, measurements, width, nq, nc, colOffset) {
  const cWireRow = 2 * nq + 1;
  const cLabelRow = 2 * nq + 2;

  const mWidth = 3; // ┤M├
  const totalMeasureWidth = measurements.length * mWidth;
  let mStartOffset = colOffset + Math.floor((width - totalMeasureWidth) / 2);
  if (mStartOffset < colOffset) mStartOffset = colOffset;

  // Assign each measurement a position
  const measurePositions = [];
  for (let mi = 0; mi < measurements.length; mi++) {
    const m = measurements[mi];
    const mOffset = mStartOffset + mi * mWidth;
    const centerCol = mOffset + 1;
    measurePositions.push({
      qubit: m.qubits[0],
      clbit: m.clbits[0],
      mOffset,
      centerCol,
      qWireRow: 2 * m.qubits[0] + 1,
    });
  }

  // Sort by qubit (top to bottom) for drawing
  const sorted = measurePositions.slice().sort((a, b) => a.qubit - b.qubit);

  // First: draw all vertical pass-through lines (║ and ╫)
  for (const mp of sorted) {
    for (let r = mp.qWireRow + 2; r < cWireRow; r++) {
      const isQuantumWire = (r % 2 === 1) && (r < 2 * nq + 1);
      if (isQuantumWire) {
        _setCell(grid, r, mp.centerCol, "\u256B"); // ╫
      } else {
        _setCell(grid, r, mp.centerCol, "\u2551"); // ║
      }
    }
  }

  // Second: draw measurement boxes, with pass-through markers on borders
  for (const mp of sorted) {
    const { qubit, clbit, mOffset, centerCol, qWireRow } = mp;

    // Check if any measurement above this one has a line at this column
    let lineFromAbove = false;
    for (const other of sorted) {
      if (other.qubit < qubit && other.centerCol === centerCol) {
        lineFromAbove = true;
        break;
      }
    }

    // Top border: ┌─┐ or ┌─┬─┐ (if line from above)
    // When the M box shares a row with a gate box bottom border (same layer),
    // merge them: ─+┌ → ┬, ─+┐ → ┬, └+┌ → ├, ┘+┐ → ┤
    const topRow = qWireRow - 1;
    // Left corner
    const existingLeft = grid[topRow] && grid[topRow][mOffset] || " ";
    if (existingLeft === "\u2500") {
      _setCell(grid, topRow, mOffset, "\u252C"); // ─ + ┌ → ┬
    } else if (existingLeft === "\u2514") {
      _setCell(grid, topRow, mOffset, "\u251C"); // └ + ┌ → ├
    } else if (existingLeft === "\u2518") {
      _setCell(grid, topRow, mOffset, "\u2524"); // ┘ + ┌ → ┤ (shouldn't happen)
    } else {
      _setCell(grid, topRow, mOffset, "\u250C"); // ┌
    }
    // Middle
    if (lineFromAbove) {
      _setCell(grid, topRow, centerCol, "\u252C"); // ┬
    } else {
      _setCell(grid, topRow, centerCol, "\u2500"); // ─
    }
    // Right corner
    const existingRight = grid[topRow] && grid[topRow][mOffset + 2] || " ";
    if (existingRight === "\u2500") {
      _setCell(grid, topRow, mOffset + 2, "\u252C"); // ─ + ┐ → ┬
    } else if (existingRight === "\u2518") {
      _setCell(grid, topRow, mOffset + 2, "\u2524"); // ┘ + ┐ → ┤
    } else if (existingRight === "\u2514") {
      _setCell(grid, topRow, mOffset + 2, "\u251C"); // └ + ┐ → ├ (shouldn't happen)
    } else {
      _setCell(grid, topRow, mOffset + 2, "\u2510"); // ┐
    }

    // Wire row: ┤M├
    _setCell(grid, qWireRow, mOffset, "\u2524");
    _setCell(grid, qWireRow, centerCol, "M");
    _setCell(grid, qWireRow, mOffset + 2, "\u251C");

    // Bottom border: └╥┘
    _setCell(grid, qWireRow + 1, mOffset, "\u2514");
    _setCell(grid, qWireRow + 1, centerCol, "\u2565"); // ╥
    _setCell(grid, qWireRow + 1, mOffset + 2, "\u2518");

    // ╩ on classical wire at target
    _setCell(grid, cWireRow, centerCol, "\u2569");

    // Bit number
    _setCell(grid, cLabelRow, centerCol, clbit.toString());
  }

  // Fill classical wire with ═
  for (let c = 0; c < width; c++) {
    if (!_touchedCells.has(`${cWireRow},${colOffset + c}`)) {
      _setCell(grid, cWireRow, colOffset + c, "\u2550");
    }
  }
}

function _drawBarrierLayer(grid, layer, width, nq, nc, colOffset) {
  const barrierQubits = new Set();
  for (const instr of layer) {
    for (const q of instr.qubits) barrierQubits.add(q);
  }
  for (let i = 0; i < nq; i++) {
    const wireRow = 2 * i + 1;
    if (barrierQubits.has(i)) {
      for (let c = 0; c < width; c++) {
        _setCell(grid, wireRow, colOffset + c, "\u2591"); // ░
      }
    } else {
      for (let c = 0; c < width; c++) {
        _setCell(grid, wireRow, colOffset + c, "\u2500"); // ─
      }
    }
  }
  if (nc > 0) {
    const cWireRow = 2 * nq + 1;
    for (let c = 0; c < width; c++) {
      _setCell(grid, cWireRow, colOffset + c, "\u2550"); // ═
    }
  }
}

function _drawInstruction(grid, instr, width, nq, nc, colOffset) {
  const op = instr.op;
  const name = op.name;

  if (name === "barrier") return;
  if (name === "measure") return; // handled by _drawMeasurementLayer

  if (name === "reset") {
    _drawBox(grid, instr.qubits[0], "|0>", width, colOffset);
    return;
  }

  if (op.numQubits === 1) {
    _drawBox(grid, instr.qubits[0], _gateLabel(op), width, colOffset);
    return;
  }

  if (op.numQubits === 2) {
    if (CONTROLLED_2Q.includes(name)) {
      const targetLabel = _controlledTargetLabel(name, op);
      _drawControlled(grid, instr.qubits[0], instr.qubits[1], targetLabel, width, nq, colOffset, name);
    } else {
      _drawSymmetric(grid, instr.qubits[0], instr.qubits[1], name, width, nq, colOffset);
    }
    return;
  }

  if (op.numQubits === 3) {
    if (name === "ccx") {
      _drawControlled(grid, instr.qubits[0], instr.qubits[2], "X", width, nq, colOffset, "x", instr.qubits[1]);
    } else if (name === "cswap") {
      _drawCSwap(grid, instr.qubits[0], instr.qubits[1], instr.qubits[2], width, nq, colOffset);
    } else {
      _drawControlled(grid, instr.qubits[0], instr.qubits[2], _gateLabel(op), width, nq, colOffset, name, instr.qubits[1]);
    }
    return;
  }
}

// Draw a single-qubit gate box: ┌───┐ / ┤ H ├ / └───┘
// Minimum interior width is 3.
function _drawBox(grid, qubit, label, width, colOffset) {
  const wireRow = 2 * qubit + 1;
  const topRow = wireRow - 1;
  const botRow = wireRow + 1;
  const innerW = Math.max(3, width - 2);
  const padL = Math.floor((innerW - label.length) / 2);
  const padR = innerW - label.length - padL;

  // Top: ┌───┐
  _setCell(grid, topRow, colOffset, "\u250C");
  for (let c = 0; c < innerW; c++) _setCell(grid, topRow, colOffset + 1 + c, "\u2500");
  _setCell(grid, topRow, colOffset + 1 + innerW, "\u2510");

  // Wire: ┤ label ├
  _setCell(grid, wireRow, colOffset, "\u2524");
  for (let c = 0; c < padL; c++) _setCell(grid, wireRow, colOffset + 1 + c, " ");
  for (let c = 0; c < label.length; c++) _setCell(grid, wireRow, colOffset + 1 + padL + c, label[c]);
  for (let c = 0; c < padR; c++) _setCell(grid, wireRow, colOffset + 1 + padL + label.length + c, " ");
  _setCell(grid, wireRow, colOffset + 1 + innerW, "\u251C");

  // Bottom: └───┘
  _setCell(grid, botRow, colOffset, "\u2514");
  for (let c = 0; c < innerW; c++) _setCell(grid, botRow, colOffset + 1 + c, "\u2500");
  _setCell(grid, botRow, colOffset + 1 + innerW, "\u2518");
}

// Draw a target box for a controlled gate.
// If lineFromAbove is true, the top border has ┴ at centerCol (┌─┴─┐).
// If lineFromAbove is false (line from below), the bottom border has ┬ (└─┬─┘).
function _drawTargetBox(grid, qubit, label, width, colOffset, lineFromAbove, centerCol) {
  const wireRow = 2 * qubit + 1;
  const topRow = wireRow - 1;
  const botRow = wireRow + 1;
  const innerW = Math.max(3, width - 2);
  const padL = Math.floor((innerW - label.length) / 2);
  const padR = innerW - label.length - padL;

  // Top: ┌─┴─┐ (if line from above) or ┌───┐
  _setCell(grid, topRow, colOffset, "\u250C");
  for (let c = 0; c < innerW; c++) {
    if (lineFromAbove && colOffset + 1 + c === centerCol) {
      _setCell(grid, topRow, colOffset + 1 + c, "\u2534"); // ┴
    } else {
      _setCell(grid, topRow, colOffset + 1 + c, "\u2500");
    }
  }
  _setCell(grid, topRow, colOffset + 1 + innerW, "\u2510");

  // Wire: ┤ label ├
  _setCell(grid, wireRow, colOffset, "\u2524");
  for (let c = 0; c < padL; c++) _setCell(grid, wireRow, colOffset + 1 + c, " ");
  for (let c = 0; c < label.length; c++) _setCell(grid, wireRow, colOffset + 1 + padL + c, label[c]);
  for (let c = 0; c < padR; c++) _setCell(grid, wireRow, colOffset + 1 + padL + label.length + c, " ");
  _setCell(grid, wireRow, colOffset + 1 + innerW, "\u251C");

  // Bottom: └───┘ (or └─┬─┘ if line from below)
  _setCell(grid, botRow, colOffset, "\u2514");
  for (let c = 0; c < innerW; c++) {
    if (!lineFromAbove && colOffset + 1 + c === centerCol) {
      _setCell(grid, botRow, colOffset + 1 + c, "\u252C"); // ┬
    } else {
      _setCell(grid, botRow, colOffset + 1 + c, "\u2500");
    }
  }
  _setCell(grid, botRow, colOffset + 1 + innerW, "\u2518");
}

// Draw a controlled gate.
// CX: ■ on control, ┤ X ├ on target, ┴ connecting
// CY: ■ on control, ┤ Y ├ on target
// CZ: ■ on both qubits (no box)
// CH: ■ on control, ┤ H ├ on target
// CRX/CRY/CRZ/CP: ■ on control, ┤ Rx(0.5) ├ on target
// CCX: ■ ■ on controls, ┤ X ├ on target
function _drawControlled(grid, control, target, label, width, nq, colOffset, gateName, extraControl = null) {
  const ctrlRow = 2 * control + 1;
  const centerCol = colOffset + Math.floor(width / 2);

  // Control dot(s): ■
  _setCell(grid, ctrlRow, centerCol, "\u25A0");
  if (extraControl !== null) {
    _setCell(grid, 2 * extraControl + 1, centerCol, "\u25A0");
  }

  // Vertical line connecting all controls to target
  // On border rows (even): use │
  // On wire rows (odd, not control/target): use ┼ (crossing marker)
  const allNodes = extraControl !== null ? [control, extraControl, target] : [control, target];
  allNodes.sort((a, b) => a - b);
  for (let i = 0; i < allNodes.length - 1; i++) {
    for (let r = 2 * allNodes[i] + 2; r <= 2 * allNodes[i + 1]; r++) {
      if (r % 2 === 1) {
        // Wire row — use ┼ (crossing) unless it's the target's wire
        const qIdx = (r - 1) / 2;
        if (qIdx !== target) {
          _setCell(grid, r, centerCol, "\u253C"); // ┼
        }
      } else {
        _setCell(grid, r, centerCol, "\u2502"); // │
      }
    }
  }

  if (gateName === "cz") {
    // CZ: ■ on target too (no box)
    _setCell(grid, 2 * target + 1, centerCol, "\u25A0");
  } else {
    // Boxed target — draw with connection to control line
    _drawTargetBox(grid, target, label, width, colOffset, control < target, centerCol);
  }
}

// Draw symmetric 2-qubit gate.
// SWAP: X on both qubits, │ connecting (qiskit uses capital X, not ×)
// iSWAP/DCX: boxed on both
// RXX/RYY/RZZ/RZX: boxed on both with label
function _drawSymmetric(grid, q0, q1, gateName, width, nq, colOffset) {
  const centerCol = colOffset + Math.floor(width / 2);

  if (gateName === "swap") {
    _setCell(grid, 2 * q0 + 1, centerCol, "X");
    _setCell(grid, 2 * q1 + 1, centerCol, "X");
    const lo = Math.min(q0, q1);
    const hi = Math.max(q0, q1);
    for (let r = 2 * lo + 2; r <= 2 * hi; r++) {
      if (r % 2 === 1) {
        _setCell(grid, r, centerCol, "\u253C"); // ┼
      } else {
        _setCell(grid, r, centerCol, "\u2502"); // │
      }
    }
  } else {
    // Other symmetric gates: draw box on both qubits.
    // iswap/dcx use a per-qubit box (matching qiskit's text-drawer output
    // for two-qubit symmetric gates without a dedicated glyph).
    const label = _symmetricLabel(gateName);
    _drawBox(grid, q0, label, width, colOffset);
    _drawBox(grid, q1, label, width, colOffset);
  }
}

function _symmetricLabel(name) {
  const map = {
    "iswap": "Iswap",
    "dcx": "Dcx",
    "rxx": "Rxx",
    "ryy": "Ryy",
    "rzz": "ZZ",  // qiskit shows RZZ as ZZ
    "rzx": "Rzx",
  };
  return map[name] || name.toUpperCase();
}

// Draw CSWAP: ■ on control, X on both targets, │ connecting
function _drawCSwap(grid, control, t1, t2, width, nq, colOffset) {
  const centerCol = colOffset + Math.floor(width / 2);

  // Control dot
  _setCell(grid, 2 * control + 1, centerCol, "\u25A0");

  // X on both targets (qiskit uses X for SWAP markers)
  _setCell(grid, 2 * t1 + 1, centerCol, "X");
  _setCell(grid, 2 * t2 + 1, centerCol, "X");

  // Vertical lines connecting all
  const allQs = [control, t1, t2].sort((a, b) => a - b);
  for (let i = 0; i < allQs.length - 1; i++) {
    for (let r = 2 * allQs[i] + 2; r <= 2 * allQs[i + 1]; r++) {
      if (r % 2 === 1) {
        _setCell(grid, r, centerCol, "\u253C"); // ┼
      } else {
        _setCell(grid, r, centerCol, "\u2502"); // │
      }
    }
  }
}

function _fillUntouchedWires(grid, nq, nc, startCol, endCol) {
  for (let i = 0; i < nq; i++) {
    const wireRow = 2 * i + 1;
    for (let c = startCol; c < endCol; c++) {
      if (!_touchedCells.has(`${wireRow},${c}`)) {
        grid[wireRow][c] = "\u2500"; // ─
      }
    }
  }
  if (nc > 0) {
    const cWireRow = 2 * nq + 1;
    for (let c = startCol; c < endCol; c++) {
      if (!_touchedCells.has(`${cWireRow},${c}`)) {
        grid[cWireRow][c] = "\u2550"; // ═
      }
    }
  }
}

// LaTeX source exporter

function _drawLatexSource(circuit) {
  const nq = circuit.numQubits;
  const nc = circuit.numClbits;
  const lines = [];
  lines.push("\\documentclass[border=2pt]{standalone}");
  lines.push("\\usepackage[braket, qm]{qcircuit}");
  lines.push("\\begin{document}");
  lines.push("\\Qcircuit @C=1.0em @R=0.7em {");
  const header = [];
  for (let i = 0; i < nq; i++) header.push(`\\lstick{\\ket{q_${i}}}`);
  for (let i = 0; i < nc; i++) header.push(`\\lstick{c_{i}}`);
  header.push("\\qw");
  lines.push("  " + header.join(" & ") + " \\\\");
  for (const ci of circuit.data) {
    const row = new Array(nq + nc).fill("\\qw");
    const name = ci.operation.name;
    if (name === "barrier") {
      for (const q of ci.qubits) row[circuit._qubit_index.get(q)] = "\\barrier{0}";
    } else if (name === "measure") {
      const qi = circuit._qubit_index.get(ci.qubits[0]);
      const ci2 = circuit._clbit_index.get(ci.clbits[0]);
      row[qi] = "\\meter";
      row[nq + ci2] = "\\cw";
    } else if (ci.operation.numQubits === 1) {
      const qi = circuit._qubit_index.get(ci.qubits[0]);
      const p = ci.operation.params.length ? `(${ci.operation.params.map(_fmtParam).join(",")})` : "";
      row[qi] = `\\gate{${name.toUpperCase()}${p}}`;
    } else if (ci.operation.numQubits === 2) {
      const q0 = circuit._qubit_index.get(ci.qubits[0]);
      const q1 = circuit._qubit_index.get(ci.qubits[1]);
      if (CONTROLLED_2Q.includes(name)) {
        row[q0] = "\\ctrl";
        row[q1] = `\\gate{${name.toUpperCase()}}`;
      } else {
        row[q0] = `\\gate{${name.toUpperCase()}}`;
        row[q1] = `\\gate{${name.toUpperCase()}}`;
      }
    } else if (ci.operation.numQubits === 3) {
      const q0 = circuit._qubit_index.get(ci.qubits[0]);
      const q1 = circuit._qubit_index.get(ci.qubits[1]);
      const q2 = circuit._qubit_index.get(ci.qubits[2]);
      row[q0] = "\\ctrl";
      row[q1] = "\\ctrl";
      row[q2] = `\\gate{${name.toUpperCase()}}`;
    }
    lines.push("  " + row.join(" & ") + " \\\\");
  }
  lines.push("}");
  lines.push("\\end{document}");
  return lines.join("\n");
}

// HTML table exporter

function _drawHtml(circuit) {
  const rows = [];
  const nq = circuit.numQubits;
  const nc = circuit.numClbits;
  const cols = [];
  for (const ci of circuit.data) {
    const col = new Array(nq + nc).fill("");
    const name = ci.operation.name;
    if (name === "barrier") {
      for (const q of ci.qubits) col[circuit._qubit_index.get(q)] = "\u2502";
    } else if (name === "measure") {
      col[circuit._qubit_index.get(ci.qubits[0])] = "M";
      col[nq + circuit._clbit_index.get(ci.clbits[0])] = "\u2190";
    } else if (ci.operation.numQubits === 1) {
      col[circuit._qubit_index.get(ci.qubits[0])] = name.toUpperCase();
    } else if (ci.operation.numQubits === 2) {
      if (CONTROLLED_2Q.includes(name)) {
        col[circuit._qubit_index.get(ci.qubits[0])] = "\u25A0";
        col[circuit._qubit_index.get(ci.qubits[1])] = name.toUpperCase();
      } else {
        col[circuit._qubit_index.get(ci.qubits[0])] = name.toUpperCase();
        col[circuit._qubit_index.get(ci.qubits[1])] = name.toUpperCase();
      }
    } else if (ci.operation.numQubits === 3) {
      col[circuit._qubit_index.get(ci.qubits[0])] = "\u25A0";
      col[circuit._qubit_index.get(ci.qubits[1])] = "\u25A0";
      col[circuit._qubit_index.get(ci.qubits[2])] = name.toUpperCase();
    }
    cols.push(col);
  }
  rows.push("<table class=\"qiskit-circuit\" style=\"border-collapse:collapse;font-family:monospace\">");
  for (let r = 0; r < nq + nc; r++) {
    rows.push("<tr>");
    const label = r < nq ? `q_${r}` : `c_${r - nq}`;
    rows.push(`<td style="padding:4px 8px;border:1px solid #ccc;font-weight:bold">${label}</td>`);
    for (const col of cols) {
      const cell = col[r] || "";
      rows.push(`<td style="padding:4px 8px;border:1px solid #eee;text-align:center">${cell}</td>`);
    }
    rows.push("</tr>");
  }
  rows.push("</table>");
  return rows.join("");
}

// Plot functions

export function plot_histogram(data, kwargs = {}) {
  const counts = data && data.toDict ? data.toDict() : data;
  const total = Object.values(counts).reduce((a, b) => a + b, 0);
  const normalized = {};
  for (const k in counts) normalized[k] = counts[k] / total;
  return _makePlaceholder("histogram", { counts, normalized, kwargs });
}

export function plot_state_city(statevector, kwargs = {}) {
  const dim = statevector.numQubits ? (1 << statevector.numQubits) : statevector.size;
  const data = statevector.data || statevector;
  const realPart = new Array(dim);
  const imagPart = new Array(dim);
  for (let i = 0; i < dim; i++) {
    realPart[i] = new Array(dim);
    imagPart[i] = new Array(dim);
    for (let j = 0; j < dim; j++) {
      const a = data.get ? data.get(i) : data[i];
      const b = data.get ? data.get(j) : data[j];
      const conjB = b.conjugate ? b.conjugate() : { re: b.re, im: -b.im };
      const prod = (a.re !== undefined)
        ? { re: a.re * conjB.re - a.im * conjB.im, im: a.re * conjB.im + a.im * conjB.re }
        : a * b;
      realPart[i][j] = prod.re;
      imagPart[i][j] = prod.im;
    }
  }
  return _makePlaceholder("state_city", { realPart, imagPart, numQubits: statevector.numQubits, kwargs });
}

export function plot_bloch_vector(blochVector, kwargs = {}) {
  return _makePlaceholder("bloch", { vector: blochVector, kwargs });
}

export function plot_state_hinton(statevector, kwargs = {}) {
  return plot_state_city(statevector, kwargs);
}

export function plot_state_qsphere(statevector, kwargs = {}) {
  return _makePlaceholder("qsphere", { statevector: statevector, kwargs });
}

function _makePlaceholder(kind, data) {
  if (typeof document !== "undefined") {
    const div = document.createElement("div");
    div.textContent = "Placeholder";
    return div;
  }
  return { kind, data };
}

_setDrawHook(drawCircuit);
