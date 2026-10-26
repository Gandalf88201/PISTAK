"""Istantanee dei risultati: ogni modifica al calcolo che cambia un numero deve essere voluta.

Per ogni CIF di `examples/` si salva in `tests/snapshots/examples/<nome>.json` tutto ciò che l'analisi produce con
le soglie del sito (`Limits()`): molecole, nuclei π, coppie con tutti i parametri predefiniti (builtin ed extra,
non i plugin), coppie di anelli e colonne di stack. Il test confronta i risultati attuali con quelli salvati.

    pytest -m snapshot                       confronto
    pytest -m snapshot --update-snapshots    riscrive le istantanee (dopo un cambiamento voluto: guarda il diff!)

Strutture reali che non si possono distribuire (per esempio dal CSD): metti i CIF in una cartella e indicala con
PISTAK_SNAPSHOT_DIR; le loro istantanee vanno in `tests/snapshots/local/`, che non entra nel repository.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest

from pistak.core.analysis import Analysis, analyse_file
from pistak.core.config import Limits
from pistak.core.ringpairs import ring_pairs
from pistak.params import compute_all, specs

ROOT = Path(__file__).resolve().parents[1]
SNAP_DIR = ROOT / "tests" / "snapshots"
REL_TOL, ABS_TOL = 1e-7, 1e-6


def _cases() -> list[tuple[str, Path, Path]]:
    cases = [(f"examples/{p.name}", p, SNAP_DIR / "examples" / f"{p.stem}.json")
             for p in sorted((ROOT / "examples").glob("*.cif")) if not p.name.startswith(".")]
    extra = os.environ.get("PISTAK_SNAPSHOT_DIR")
    if extra:
        for p in sorted(Path(extra).expanduser().glob("*.cif")):
            if not p.name.startswith("."):
                cases.append((f"local/{p.name}", p, SNAP_DIR / "local" / f"{p.stem}.json"))
    return cases


def _num(x) -> float | None:
    x = float(x)
    if not math.isfinite(x):
        return None  # NaN e ±inf diventano null (il confronto li tratta come uguali tra loro)
    return float(f"{x:.10g}")


def _vec(v) -> list[float | None]:
    return [_num(x) for x in v]


def _pair_sort_key(p: dict) -> tuple:
    # a→b e b→a con lo stesso interplanare possono scambiarsi di posto per rumore numerico (~1e-15):
    # l'ordine si fissa con valori arrotondati e indici
    v = p["params"]
    return (p["kind"], round(v["interplanar"] or 0, 4), round(v["contact"] or 0, 4), p["a"], p["b"], p["shift"])


def snapshot(a: Analysis) -> dict:
    """Tutto il risultato di un'analisi in forma JSON (solo tipi semplici, numeri a 10 cifre significative)."""
    keys = sorted(s.key for s in specs() if s.source in ("builtin", "extra"))
    molecules = []
    for m, core in zip(a.structure.molecules, a.cores):
        molecules.append({
            "formula": "".join(f"{e}{n}" for e, n in m.formula()), "atoms": len(m), "bonds": len(m.bonds),
            "extended": m.extended, "centroid": _vec(m.centroid),
            "core": None if core is None else {
                "kind": core.kind, "atoms": [int(i) for i in core.atoms], "rms": _num(core.rms),
                "centroid": _vec(core.centroid), "normal": _vec(core.normal), "long": _vec(core.long),
            },
        })
    pairs = []
    for p in a.pairs:
        ctx = a.pair_context(p)
        pairs.append({
            "a": p.a, "b": p.b, "shift": list(p.shift), "kind": p.kind, "count": p.count,
            "members": sorted([m[0], m[1], list(m[2])] for m in p.members),
            "params": {k: _num(p.get(k)) for k in keys},
            "ring_pairs": [{"rings": [r.ring_a, r.ring_b], "sizes": [r.size_a, r.size_b], "cg": _num(r.cg_distance),
                            "alpha": _num(r.alpha), "beta": _num(r.beta), "gamma": _num(r.gamma),
                            "slippage": _num(r.slippage), "closest": _num(r.closest)} for r in ring_pairs(ctx)],
        })
    pairs.sort(key=_pair_sort_key)
    stacks = [{"rank": c.rank, "molecules": sorted(c.molecules), "translations": [list(t) for t in c.translations],
               "repeat": None if c.repeat is None else _num(c.repeat), "edges": len(c.edges),
               "mean_interplanar": _num(c.mean("interplanar"))} for c in a.stacks]
    return {
        "space_group": a.cif.space_group, "cell_atoms": a.structure.n_cell_atoms,
        "skipped_disorder": a.cif.skipped_disorder, "molecules": molecules, "pairs": pairs, "stacks": stacks,
    }


def _diff(exp, got, path: str = "$") -> list[str]:
    """Differenze tra due istantanee, con il percorso del valore; i numeri si confrontano con tolleranza."""
    if isinstance(exp, dict) and isinstance(got, dict):
        out = [f"{path}.{k}: manca nel risultato attuale" for k in exp.keys() - got.keys()]
        out += [f"{path}.{k}: nuovo nel risultato attuale" for k in got.keys() - exp.keys()]
        for k in sorted(exp.keys() & got.keys()):
            out += _diff(exp[k], got[k], f"{path}.{k}")
        return out
    if isinstance(exp, list) and isinstance(got, list):
        if len(exp) != len(got):
            return [f"{path}: {len(exp)} elementi salvati, {len(got)} attuali"]
        return [d for i, (x, y) in enumerate(zip(exp, got)) for d in _diff(x, y, f"{path}[{i}]")]
    if isinstance(exp, float) or isinstance(got, float):
        if isinstance(exp, (int, float)) and isinstance(got, (int, float)) and not isinstance(exp, bool) \
                and math.isclose(exp, got, rel_tol=REL_TOL, abs_tol=ABS_TOL):
            return []
        return [f"{path}: salvato {exp!r}, attuale {got!r}"]
    return [] if exp == got else [f"{path}: salvato {exp!r}, attuale {got!r}"]


@pytest.mark.parametrize("label,cif,snap", _cases(), ids=[c[0] for c in _cases()])
def test_risultati_invariati(label, cif, snap, request):
    got = snapshot(compute_all(analyse_file(cif, Limits())))
    if request.config.getoption("--update-snapshots"):
        snap.parent.mkdir(parents=True, exist_ok=True)
        snap.write_text(json.dumps(got, indent=1, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
        return
    if not snap.exists():
        pytest.fail(f"manca l'istantanea {snap.relative_to(ROOT)}: creala con  pytest -m snapshot --update-snapshots")
    diffs = _diff(json.loads(snap.read_text(encoding="utf-8")), got)
    assert not diffs, f"{label}: {len(diffs)} differenze rispetto all'istantanea\n" + "\n".join(diffs[:40])


def test_il_confronto_vede_le_differenze():
    base = {"x": 3.4, "v": [1.0, None], "k": "stack", "n": 2}
    assert _diff(base, dict(base)) == []
    assert _diff(base, {**base, "x": 3.4 + 1e-9}) == []
    assert _diff(base, {**base, "x": 3.41})
    assert _diff(base, {**base, "v": [1.0, 2.0]})
    assert _diff(base, {**base, "k": "parallel"})
    assert _diff(base, {**base, "n": 3})
    assert _diff(base, {k: v for k, v in base.items() if k != "n"})
