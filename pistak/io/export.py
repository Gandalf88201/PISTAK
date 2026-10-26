"""Esportazione dei risultati: tabella CSV, JSON completo, coppia in formato XYZ.

I file scritti qui hanno intestazioni e chiavi fisse (inglesi, mai tradotte): uno script che li legge funziona
qualunque sia la lingua dell'interfaccia. I nomi leggibili dei parametri stanno nel registro (`pistak params`).
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from pistak.core.analysis import Analysis
from pistak.core.pairs import Pair
from pistak.params import get_spec

#: Formato del JSON scritto da `write_json`: cambia `JSON_VERSION` se cambia la struttura.
JSON_FORMAT, JSON_VERSION = "pistak-analysis", 2


def table(analysis: Analysis, keys: list[str] | None = None) -> tuple[list[str], list[list]]:
    """Intestazioni e righe della tabella delle coppie: indice, tipo, le due molecole (da 1), la traslazione in celle
    della seconda e i parametri scelti (per impostazione tutti quelli del registro, per chiave)."""
    keys = keys if keys is not None else list(analysis.param_keys)
    headers = ["n", "kind", "molecule_a", "molecule_b", "shift", *keys]
    rows = []
    for p in analysis.pairs:
        rows.append([p.index + 1, p.kind, p.a + 1, p.b + 1, " ".join(map(str, p.shift))] + [p.get(k) for k in keys])
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
        "symops_source": analysis.cif.symops_source,
        "cell": {"a": c.a, "b": c.b, "c": c.c, "alpha": c.alpha, "beta": c.beta, "gamma": c.gamma},
        "cell_atoms": analysis.structure.n_cell_atoms,
        "molecules": [{"formula": g.text, "count": g.count, "with_pi_core": g.with_core, "pi_fragments": g.fragments}
                      for g in analysis.molecule_groups()],
        "has_extended_network": analysis.has_extended,
        "skipped_disorder_atoms": analysis.cif.skipped_disorder,
        "disorder_policy": analysis.cif.disorder,
        "limits": analysis.limits.to_dict(),
    }


def pair_dict(p: Pair, keys: list[str]) -> dict:
    return {
        "n": p.index + 1, "kind": p.kind, "molecule_a": p.a + 1, "molecule_b": p.b + 1, "shift": list(p.shift),
        "fragment_a": p.fa + 1, "fragment_b": p.fb + 1,
        "equivalent_pairs": p.count, "parameters": {k: _clean(p.get(k)) for k in keys},
    }


def fragments_list(analysis: Analysis) -> list[dict]:
    """I frammenti π di ogni molecola (da 1), con atomi, planarità, isotropia e origine dell'asse lungo."""
    out = []
    for m, frags in enumerate(analysis.fragments or [[c] if c else [] for c in analysis.cores]):
        for k, f in enumerate(frags):
            out.append({"molecule": m + 1, "fragment": k + 1, "kind": f.kind, "atoms": len(f.atoms),
                        "rms": _clean(float(f.rms)), "aspect": _clean(float(f.plane.aspect)), "axis": f.axis_source,
                        "labels": [analysis.structure.molecules[m].labels[i] for i in f.atoms]})
    return out


def contacts_list(analysis: Analysis) -> list[dict]:
    """Il riepilogo per molecola: un elemento per contatto tra molecole, con le sue coppie di frammenti (numeri n)."""
    return [{"molecule_a": c.a + 1, "molecule_b": c.b + 1, "shift": list(c.shift), "intramolecular": c.intramolecular,
             "count": c.count, "contact": _clean(c.contact), "pairs": [p.index + 1 for p in c.pairs],
             "n_stack": len(c.stacks), "best": c.best.index + 1} for c in analysis.molecule_contacts()]


def write_json(analysis: Analysis, path: str | Path, keys: list[str] | None = None) -> Path:
    path = Path(path)
    keys = keys if keys is not None else list(analysis.param_keys)
    units = {k: (get_spec(k).unit if get_spec(k) else "") for k in keys}
    from pistak import __version__

    data = {"format": JSON_FORMAT, "version": JSON_VERSION, "pistak": __version__, "summary": summary_dict(analysis),
            "units": units, "fragments": fragments_list(analysis), "pairs": [pair_dict(p, keys) for p in analysis.pairs],
            "molecule_contacts": contacts_list(analysis)}
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return path


def pair_xyz(analysis: Analysis, pair: Pair) -> str:
    """Le due molecole della coppia, in coordinate cartesiane (Å), nel formato XYZ."""
    ctx = analysis.pair_context(pair)
    atoms = [(e, p) for e, p in zip(ctx.A.el, ctx.pos_a)] + [(e, p) for e, p in zip(ctx.B.el, ctx.pos_b)]
    head = (f"{analysis.name} pair={pair.index + 1} kind={pair.kind} interplanar={pair.interplanar:.3f} "
            f"long_slip={pair.long_slip:.3f} lat_slip={pair.lat_slip:.3f}")
    body = "\n".join(f"{e:<2} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}" for e, p in atoms)
    return f"{len(atoms)}\n{head}\n{body}\n"


def write_pair_xyz(analysis: Analysis, pair: Pair, path: str | Path) -> Path:
    path = Path(path)
    path.write_text(pair_xyz(analysis, pair), encoding="utf-8")
    return path
