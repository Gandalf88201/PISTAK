"""Esportazione dei risultati: tabella CSV, JSON completo, coppia in formato XYZ."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from pistak.core.analysis import Analysis
from pistak.core.pairs import Pair
from pistak.i18n import tf, tr
from pistak.params import get_spec


def table(analysis: Analysis, keys: list[str] | None = None) -> tuple[list[str], list[list]]:
    """Intestazioni e righe della tabella delle coppie: indice, tipo, relazione e i parametri scelti
    (per impostazione tutti quelli del registro)."""
    keys = keys if keys is not None else list(analysis.param_keys)
    headers = ["n", tr("tipo"), tr("relazione")] + [k for k in keys]
    rows = []
    for p in analysis.pairs:
        rows.append([p.index + 1, p.kind, p.relation()] + [p.get(k) for k in keys])
    return headers, rows


def _clean(v):
    return None if isinstance(v, float) and not math.isfinite(v) else v


def write_csv(analysis: Analysis, path: str | Path, keys: list[str] | None = None) -> Path:
    path = Path(path)
    headers, rows = table(analysis, keys)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in rows:
            w.writerow(["" if _clean(v) is None else (f"{v:.5g}" if isinstance(v, float) else v) for v in r])
    return path


def summary_dict(analysis: Analysis) -> dict:
    c = analysis.cif.cell
    return {
        "file": analysis.name,
        "block": analysis.cif.name,
        "space_group": analysis.cif.space_group,
        "cell": {"a": c.a, "b": c.b, "c": c.c, "alpha": c.alpha, "beta": c.beta, "gamma": c.gamma},
        "cell_atoms": analysis.structure.n_cell_atoms,
        "molecules": [{"formula": g.text, "count": g.count, "with_pi_core": g.with_core} for g in analysis.molecule_groups()],
        "has_extended_network": analysis.has_extended,
        "skipped_disorder_atoms": analysis.cif.skipped_disorder,
        "limits": analysis.limits.to_dict(),
    }


def pair_dict(p: Pair, keys: list[str]) -> dict:
    return {
        "n": p.index + 1, "kind": p.kind, "molecule_a": p.a + 1, "molecule_b": p.b + 1, "shift": list(p.shift),
        "equivalent_pairs": p.count, "parameters": {k: _clean(p.get(k)) for k in keys},
    }


def write_json(analysis: Analysis, path: str | Path, keys: list[str] | None = None) -> Path:
    path = Path(path)
    keys = keys if keys is not None else list(analysis.param_keys)
    units = {k: (get_spec(k).unit if get_spec(k) else "") for k in keys}
    data = {"summary": summary_dict(analysis), "units": units, "pairs": [pair_dict(p, keys) for p in analysis.pairs]}
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def pair_xyz(analysis: Analysis, pair: Pair) -> str:
    """Le due molecole della coppia, in coordinate cartesiane (Å), nel formato XYZ."""
    ctx = analysis.pair_context(pair)
    atoms = [(e, p) for e, p in zip(ctx.A.el, ctx.pos_a)] + [(e, p) for e, p in zip(ctx.B.el, ctx.pos_b)]
    head = tf("{name} coppia {n} {kind} interplanare={ip:.3f} slittamento_L={ls:.3f} slittamento_W={ws:.3f}", name=analysis.name,
              n=pair.index + 1, kind=pair.kind, ip=pair.interplanar, ls=pair.long_slip, ws=pair.lat_slip)
    body = "\n".join(f"{e:<2} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}" for e, p in atoms)
    return f"{len(atoms)}\n{head}\n{body}\n"


def write_pair_xyz(analysis: Analysis, pair: Pair, path: str | Path) -> Path:
    path = Path(path)
    path.write_text(pair_xyz(analysis, pair), encoding="utf-8")
    return path
