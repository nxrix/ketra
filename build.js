import { build } from "esbuild";
import fs from "fs";

const index = fs.readFileSync("src/index.js","utf8");
const banner = "/* Copyright nxrix, 2026 */";

const common = {
  entryPoints: ["src/index.js"],
  bundle: true,
  format: "esm",
  banner: { js: banner }
};

await build({
  ...common,
  minify: false,
  outfile: "dist/ketra.js",
});

await build({
  ...common,
  minify: true,
  outfile: "dist/ketra.min.js",
});
