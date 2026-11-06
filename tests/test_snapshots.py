"""Result snapshots: every change to the computation that changes a number must be intended.

For every CIF in `examples/` and `tests/data/cod/` everything the analysis produces is saved (molecules, π cores,
pairs with all the default parameters — builtin and extra, not plugins —, ring pairs, stack columns, structure
quantities, quick-check alerts and the molecular geometry: how many measurements per kind and the list of van der
Waals contacts) in two modes:

    tests/snapshots/site/<name>.json        Limits.site(): the website conventions (must stay as they always were)
    tests/snapshots/canonical/<name>.json   Limits(): the defaults (unique pairs, symmetry equivalence)

The test compares the current results with the saved ones.

    pytest -m snapshot                       comparison
    pytest -m snapshot --update-snapshots    rewrite the snapshots (after an intended change: look at the diff!)

Real structures that cannot be redistributed (for example from the CSD): put the CIFs in a folder and point
PISTAK_SNAPSHOT_DIR to it; their snapshots go to `tests/snapshots/local/<mode>/`, which stays out of the repository.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import pytest

from pistak.core.analysis import Analysis, analyse_file
from pistak.core.config import Limits
from pistak.core.intramolecular import molecule_geometry
from pistak.core.ringpairs import ring_pairs
from pistak.params import compute_all, compute_structure, specs, structure_specs

ROOT = Path(__file__).resolve().parents[1]
SNAP_DIR = ROOT / "tests" / "snapshots"
REL_TOL, ABS_TOL = 1e-7, 1e-6


MODES = {"site": Limits.site, "canonical": Limits}


def _cases() -> list[tuple[str, Path, Path, str]]:
    cifs = [("", p) for d in (ROOT / "examples", ROOT / "tests" / "data" / "cod")
            for p in sorted(d.glob("*.cif")) if not p.name.startswith(".")]
    extra = os.environ.get("PISTAK_SNAPSHOT_DIR")
    if extra:
        cifs += [("local/", p) for p in sorted(Path(extra).expanduser().glob("*.cif")) if not p.name.startswith(".")]
    return [(f"{where}{mode}/{p.name}", p, SNAP_DIR / f"{where}{mode}" / f"{p.stem}.json", mode)
            for mode in MODES for where, p in cifs]


def _num(x) -> float | None:
    x = float(x)
    if not math.isfinite(x):
        return None  # NaN and ±inf become null (the comparison treats them as equal to each other)
    return float(f"{x:.10g}")


def _vec(v) -> list[float | None]:
    return [_num(x) for x in v]


def _pair_sort_key(p: dict) -> tuple:
    # a→b and b→a with the same interplanar distance can swap places because of numerical noise (~1e-15):
    # the order is fixed with rounded values and indices
    v = p["params"]
    return (p["kind"], round(v["interplanar"] or 0, 4), round(v["contact"] or 0, 4), p["a"], p["b"], p["shift"])


def snapshot(a: Analysis) -> dict:
    """The whole result of an analysis as JSON (simple types only, numbers to 10 significant digits)."""
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
    structure = {s.key: _num(a.structure_values.get(s.key, math.nan)) for s in structure_specs() if s.source == "builtin"}
    alerts = sorted(([al.level, al.code, {k: _num(v) if isinstance(v, float) else v for k, v in al.values.items()}]
                     for al in a.alerts or []), key=lambda x: (x[0], x[1], json.dumps(x[2], sort_keys=True)))
    geometry = molecule_geometry(a.structure, a.cif)
    counts: dict[str, int] = {}
    for g in geometry:
        counts[g.kind] = counts.get(g.kind, 0) + 1
    contacts = [{"atoms": list(g.atoms), "codes": list(g.codes), "d": _num(g.value), "vdw": _num(g.vdw)}
                for g in geometry if g.kind == "contact"]
    return {
        "space_group": a.cif.space_group, "cell_atoms": a.structure.n_cell_atoms,
        "skipped_disorder": a.cif.skipped_disorder, "molecules": molecules, "pairs": pairs, "stacks": stacks,
        "structure": structure, "alerts": alerts, "geometry": {"counts": counts, "contacts": contacts},
    }


def _diff(exp, got, path: str = "$") -> list[str]:
    """Differences between two snapshots, with the path of each value; numbers are compared with a tolerance."""
    if isinstance(exp, dict) and isinstance(got, dict):
        out = [f"{path}.{k}: missing from the current result" for k in exp.keys() - got.keys()]
        out += [f"{path}.{k}: new in the current result" for k in got.keys() - exp.keys()]
        for k in sorted(exp.keys() & got.keys()):
            out += _diff(exp[k], got[k], f"{path}.{k}")
        return out
    if isinstance(exp, list) and isinstance(got, list):
        if len(exp) != len(got):
            return [f"{path}: {len(exp)} items saved, {len(got)} now"]
        return [d for i, (x, y) in enumerate(zip(exp, got)) for d in _diff(x, y, f"{path}[{i}]")]
    if isinstance(exp, float) or isinstance(got, float):
        if isinstance(exp, (int, float)) and isinstance(got, (int, float)) and not isinstance(exp, bool) \
                and math.isclose(exp, got, rel_tol=REL_TOL, abs_tol=ABS_TOL):
            return []
        return [f"{path}: saved {exp!r}, now {got!r}"]
    return [] if exp == got else [f"{path}: saved {exp!r}, now {got!r}"]


@pytest.mark.parametrize("label,cif,snap,mode", _cases(), ids=[c[0] for c in _cases()])
def test_results_unchanged(label, cif, snap, mode, request):
    got = snapshot(compute_structure(compute_all(analyse_file(cif, MODES[mode]()))))
    if request.config.getoption("--update-snapshots"):
        snap.parent.mkdir(parents=True, exist_ok=True)
        snap.write_text(json.dumps(got, indent=1, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
        return
    if not snap.exists():
        pytest.fail(f"snapshot {snap.relative_to(ROOT)} is missing: create it with  pytest -m snapshot --update-snapshots")
    diffs = _diff(json.loads(snap.read_text(encoding="utf-8")), got)
    assert not diffs, f"{label}: {len(diffs)} differences from the snapshot\n" + "\n".join(diffs[:40])


def test_comparison_sees_differences():
    base = {"x": 3.4, "v": [1.0, None], "k": "stack", "n": 2}
    assert _diff(base, dict(base)) == []
    assert _diff(base, {**base, "x": 3.4 + 1e-9}) == []
    assert _diff(base, {**base, "x": 3.41})
    assert _diff(base, {**base, "v": [1.0, 2.0]})
    assert _diff(base, {**base, "k": "parallel"})
    assert _diff(base, {**base, "n": 3})
    assert _diff(base, {k: v for k, v in base.items() if k != "n"})
