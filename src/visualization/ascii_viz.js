// ASCII Histogram
export function render_histogram_ascii(counts, options = {}) {
  const maxBarWidth = options.maxBarWidth || 40;
  const data = counts && counts.toDict ? counts.toDict() : counts;
  const entries = Object.entries(data).sort((a, b) => b[1] - a[1]);
  const maxCount = Math.max(...entries.map(e => e[1]));
  const total = entries.reduce((s, e) => s + e[1], 0);
  const barChar = options.barChar || "█";
  const emptyChar = options.emptyChar || " ";

  const labelWidth = Math.max(...entries.map(e => e[0].length));
  const countWidth = String(maxCount).length;
  const lines = [];

  for (const [key, count] of entries) {
    const barLen = Math.round((count / maxCount) * maxBarWidth);
    const bar = barChar.repeat(barLen);
    const pct = ((count / total) * 100).toFixed(1);
    const paddedKey = key.padEnd(labelWidth);
    const paddedCount = String(count).padStart(countWidth);
    lines.push(`${paddedKey} │${bar} ${paddedCount} (${pct}%)`);
  }

  // Add axis
  const axisLine = " ".repeat(labelWidth) + "└" + "─".repeat(maxBarWidth + countWidth + 8);
  lines.push(axisLine);

  return lines.join("\n");
}

// ASCII Bloch Sphere
export function render_bloch_ascii(vector, options = {}) {
  // vector: [x, y, z] Bloch vector
  const [x, y, z] = vector;
  const width = options.width || 24;
  const height = options.height || 12;
  const radius = Math.min(width, height * 2) / 2 - 2;

  // Project Bloch sphere to 2D (orthographic projection, looking from +y axis)
  // Screen x = Bloch x, Screen y = -Bloch z (up is +z)
  const px = Math.round(radius * x + width / 2);
  const py = Math.round(-radius * z + height / 2);

  const grid = [];
  for (let r = 0; r < height; r++) {
    grid.push(new Array(width).fill(" "));
  }

  // Draw circle outline
  for (let theta = 0; theta < 2 * Math.PI; theta += 0.05) {
    const cx = Math.round(radius * Math.cos(theta) + width / 2);
    const cy = Math.round(radius * Math.sin(theta) * 0.5 + height / 2);
    if (cy >= 0 && cy < height && cx >= 0 && cx < width) {
      grid[cy][cx] = "·";
    }
  }

  // Draw axes
  // X axis (horizontal)
  for (let i = -radius; i <= radius; i++) {
    const cx = Math.round(i + width / 2);
    const cy = Math.round(height / 2);
    if (cy >= 0 && cy < height && cx >= 0 && cx < width && grid[cy][cx] === " ") {
      grid[cy][cx] = "─";
    }
  }
  // Z axis (vertical)
  for (let i = -radius / 2; i <= radius / 2; i++) {
    const cx = Math.round(width / 2);
    const cy = Math.round(i + height / 2);
    if (cy >= 0 && cy < height && cx >= 0 && cx < width && grid[cy][cx] === " ") {
      grid[cy][cx] = "│";
    }
  }

  // Draw the vector
  // Draw line from center to (px, py)
  const cx0 = Math.round(width / 2);
  const cy0 = Math.round(height / 2);
  const steps = Math.max(Math.abs(px - cx0), Math.abs(py - cy0), 1);
  for (let s = 0; s <= steps; s++) {
    const ix = Math.round(cx0 + (px - cx0) * s / steps);
    const iy = Math.round(cy0 + (py - cy0) * s / steps);
    if (iy >= 0 && iy < height && ix >= 0 && ix < width) {
      grid[iy][ix] = "•";
    }
  }

  // Mark the endpoint
  if (py >= 0 && py < height && px >= 0 && px < width) {
    grid[py][px] = "◉";
  }

  // Label axes
  _setGrid(grid, 0, width - 1, "x");
  _setGrid(grid, 0, 0, "-x");
  _setGrid(grid, 0, Math.floor(width / 2), "z");

  return grid.map(row => row.join("")).join("\n");
}

// ASCII State City (3D-ish bar chart for density matrix)
export function render_state_city_ascii(statevector, options = {}) {
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

  const lines = [];
  lines.push("=== State City (Real part of ρ) ===");
  lines.push(_renderMatrix3D(realPart, options));
  lines.push("");
  lines.push("=== State City (Imaginary part of ρ) ===");
  lines.push(_renderMatrix3D(imagPart, options));
  return lines.join("\n");
}

function _renderMatrix3D(matrix, options = {}) {
  const maxBarHeight = options.maxBarHeight || 8;
  const dim = matrix.length;
  const maxVal = Math.max(...matrix.flat().map(Math.abs));
  if (maxVal < 1e-15) return "(all zeros)";

  const lines = [];
  const barChar = "█";

  // Render as a 2D grid of bars (height proportional to value)
  // Each cell is a column of characters
  for (let row = 0; row < dim; row++) {
    let line = "";
    for (let col = 0; col < dim; col++) {
      const val = matrix[row][col];
      const height = Math.round(Math.abs(val) / maxVal * maxBarHeight);
      const sign = val >= 0 ? "+" : "-";
      // Show as a single char representing the bar
      if (height === 0) {
        line += "  ·  ";
      } else {
        line += sign + barChar.repeat(Math.min(height, 4)) + " ";
      }
    }
    lines.push(line);
  }

  // Also show numeric values
  lines.push("");
  lines.push("Values:");
  for (let row = 0; row < dim; row++) {
    let line = "";
    for (let col = 0; col < dim; col++) {
      line += matrix[row][col].toFixed(3).padStart(8) + " ";
    }
    lines.push(line);
  }

  return lines.join("\n");
}

export function render_qsphere_ascii(statevector, options = {}) {
  const n = statevector.numQubits;
  const dim = 1 << n;
  const probs = statevector.probabilities();
  const lines = [];
  lines.push("=== QSphere ===");

  // Group by Hamming weight
  const groups = new Array(n + 1);
  for (let i = 0; i <= n; i++) groups[i] = [];
  for (let i = 0; i < dim; i++) {
    const weight = _hammingWeight(i);
    groups[weight].push({ idx: i, prob: probs[i] });
  }

  // Render as concentric rings
  for (let w = 0; w <= n; w++) {
    const states = groups[w];
    if (states.length === 0) continue;
    const totalProb = states.reduce((s, st) => s + st.prob, 0);
    const ringRadius = w;
    lines.push(`  Ring ${w} (|${"1".repeat(w)}${"0".repeat(n - w)}⟩): total prob = ${totalProb.toFixed(4)}`);
    for (const st of states) {
      const bits = st.idx.toString(2).padStart(n, "0");
      const barLen = Math.round(st.prob * 30);
      lines.push(`    |${bits}⟩ ${"█".repeat(barLen)} ${st.prob.toFixed(4)}`);
    }
  }

  return lines.join("\n");
}

function _hammingWeight(n) {
  let count = 0;
  while (n) { count += n & 1; n >>= 1; }
  return count;
}

function _setGrid(grid, row, col, char) {
  if (row >= 0 && row < grid.length && col >= 0 && col < grid[0].length) {
    if (char.length === 1) {
      grid[row][col] = char;
    } else {
      // Multi-char: place each char
      for (let i = 0; i < char.length && col + i < grid[0].length; i++) {
        grid[row][col + i] = char[i];
      }
    }
  }
}
