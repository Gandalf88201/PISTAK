#!/usr/bin/env node
/**
 * Esegue l'algoritmo ORIGINALE del sito (src/lib/crystal/*.ts) sui CIF di `examples/` e salva i risultati
 * in tests/fixtures/site_expected.json, che `test_parity.py` confronta con il porting Python.
 *
 *   SITE_DIR=/percorso/SITOMIO node tests/parity/dump_site.js
 *
 * Serve la cartella del sito con `node_modules/typescript` (basta quella di Next.js): i file .ts vengono
 * compilati in una cartella temporanea, senza toccare il sito.
 */
const fs = require("fs");
const os = require("os");
const path = require("path");

const SITE = process.env.SITE_DIR;
if (!SITE) throw new Error("imposta SITE_DIR con il percorso della cartella del sito");
const ts = require(path.join(SITE, "node_modules", "typescript"));
const tmp = fs.mkdtempSync(path.join(os.tmpdir(), "pistak-site-"));
const srcDir = path.join(SITE, "src", "lib", "crystal");
for (const f of fs.readdirSync(srcDir).filter((f) => f.endsWith(".ts"))) {
  const out = ts.transpileModule(fs.readFileSync(path.join(srcDir, f), "utf8"), {
    compilerOptions: { module: ts.ModuleKind.CommonJS, target: ts.ScriptTarget.ES2020 },
  });
  fs.writeFileSync(path.join(tmp, f.replace(/\.ts$/, ".js")), out.outputText);
}
const { parseCif } = require(path.join(tmp, "cif.js"));
const { knownElement } = require(path.join(tmp, "radii.js"));
const { buildStructure } = require(path.join(tmp, "structure.js"));
const { analyse } = require(path.join(tmp, "pistack.js"));

const root = path.resolve(__dirname, "..", "..");
const result = {};
for (const f of fs.readdirSync(path.join(root, "examples")).filter((f) => f.endsWith(".cif") && !f.startsWith(".")).sort()) {
  // (i file "._*" sono metadati di macOS su volumi exFAT, non CIF)
  const text = fs.readFileSync(path.join(root, "examples", f), "utf8");
  try {
    const s = buildStructure(parseCif(text, knownElement));
    const a = analyse(s);
    result[f] = {
      cellAtoms: s.cellAtoms.length,
      molecules: s.molecules.map((m, i) => ({
        atoms: m.atoms.length, bonds: m.bonds.length, extended: m.extended,
        core: a.cores[i] ? { kind: a.cores[i].kind, atoms: a.cores[i].atoms.length, rms: a.cores[i].rms } : null,
      })),
      pairs: a.pairs.map((p) => ({
        kind: p.kind, interplanar: p.interplanar, longSlip: p.longSlip, latSlip: p.latSlip, slip: p.slip,
        dihedral: p.dihedral, centroidDistance: p.centroidDistance, contact: p.contact, overlap: p.overlap, count: p.count,
      })),
    };
  } catch (e) {
    result[f] = { error: e.code || String(e) };
  }
}
fs.mkdirSync(path.join(root, "tests", "fixtures"), { recursive: true });
fs.writeFileSync(path.join(root, "tests", "fixtures", "site_expected.json"), JSON.stringify(result, null, 1));
console.log("scritto tests/fixtures/site_expected.json:", Object.keys(result).join(", "));
fs.rmSync(tmp, { recursive: true, force: true });
