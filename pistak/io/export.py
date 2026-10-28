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
from pistak.params import get_spec, get_structure_spec

#: Formato del JSON scritto da `write_json`: cambia `JSON_VERSION` se cambia la struttura.
JSON_FORMAT, JSON_VERSION = "pistak-analysis", 4


def table(analysis: Analysis, keys: list[str] | None = None) -> tuple[list[str], list[list]]:
    """Intestazioni e righe della tabella delle coppie: indice, tipo, le due molecole (da 1), la traslazione in celle
    della seconda e i parametri scelti (per impostazione tutti quelli del registro, per chiave)."""
    keys = keys if keys is not None else list(analysis.param_keys)
    headers = ["n", "kind", "molecule_a", "molecule_b", "shift", *keys]
    u = analysis.uncertainty
    su_keys = [k for k in keys if u is not None and k in u.keys]
    if u is not None:  # incertezze dal Monte Carlo: una colonna <chiave>_su per parametro, più la stabilità del tipo
        headers += [f"{k}_su" for k in su_keys] + ["kind_agreement"]
    rows = []
    for p in analysis.pairs:
        row = [p.index + 1, p.kind, p.a + 1, p.b + 1, " ".join(map(str, p.shift))] + [p.get(k) for k in keys]
        if u is not None:
            row += [u.get(p, k) for k in su_keys] + [u.kind_agreement.get(p.index, math.nan)]
        rows.append(row)
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
        "uncertainty": None if analysis.uncertainty is None else {
            "method": "monte_carlo", "samples": analysis.uncertainty.samples, "seed": analysis.uncertainty.seed,
            "include_cell": analysis.uncertainty.include_cell, "has_su": analysis.uncertainty.has_su},
        "structure": {k: _clean(v) for k, v in analysis.structure_values.items()},
        "alerts": alerts_list(analysis),
    }


def alerts_list(analysis: Analysis) -> list[dict]:
    """Gli avvisi dei controlli rapidi: codice, livello e valori (chiavi fisse; il testo lo fa l'interfaccia)."""
    return [{"code": a.code, "level": a.level, "values": {k: _clean(v) for k, v in a.values.items()}}
            for a in analysis.alerts or []]


def write_alerts_csv(analysis: Analysis, path: str | Path) -> Path:
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["level", "code", "values"])
        for a in analysis.alerts or []:
            w.writerow([a.level, a.code, "; ".join(f"{k}={v}" for k, v in a.values.items())])
    return path


def write_geometry_csv(analysis: Analysis, path: str | Path, include_h: bool = True) -> Path:
    """Legami, angoli, torsioni e contatti di van der Waals delle molecole indipendenti (`pistak.core.intramolecular`),
    con le s.u. Le colonne symop_i danno l'operazione che genera l'atomo i dall'atomo del CIF (vuota: quello del CIF);
    vdw_sum e delta (distanza − somma dei raggi) solo per i contatti."""
    from pistak.core.intramolecular import molecule_geometry

    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["kind", "molecule", "atom_1", "atom_2", "atom_3", "atom_4", "symop_1", "symop_2", "symop_3", "symop_4",
                    "value", "su", "vdw_sum", "delta"])
        for g in molecule_geometry(analysis.structure, analysis.cif, include_h=include_h):
            atoms = list(g.atoms) + [""] * (4 - len(g.atoms))
            ops = list(g.codes or [""] * len(g.atoms)) + [""] * (4 - len(g.atoms))
            contact = g.kind == "contact"
            w.writerow([g.kind, g.molecule + 1, *atoms, *ops, f"{g.value:.5g}", "" if not math.isfinite(g.su) else f"{g.su:.2g}",
                        f"{g.vdw:.2f}" if contact else "", f"{g.delta:.3f}" if contact else ""])
    return path


def pair_dict(p: Pair, keys: list[str], uncertainty=None) -> dict:
    d = {
        "n": p.index + 1, "kind": p.kind, "molecule_a": p.a + 1, "molecule_b": p.b + 1, "shift": list(p.shift),
        "fragment_a": p.fa + 1, "fragment_b": p.fb + 1,
        "equivalent_pairs": p.count, "parameters": {k: _clean(p.get(k)) for k in keys},
    }
    if uncertainty is not None:
        d["su"] = {k: _clean(uncertainty.get(p, k)) for k in keys if k in uncertainty.keys}
        d["kind_agreement"] = _clean(uncertainty.kind_agreement.get(p.index, math.nan))
    return d


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


def molecules_list(analysis: Analysis) -> list[dict]:
    """Le molecole della cella (da 1, come `molecule_a`/`molecule_b` delle coppie): formula, atomi, se sono una rete
    estesa e se sono quelle indipendenti per simmetria su cui si misura la geometria."""
    from pistak.core.intramolecular import unique_molecules

    unique = set(unique_molecules(analysis.structure))
    return [{"n": k + 1, "formula": " ".join(f"{e}{c}" for e, c in m.formula()), "atoms": len(m), "extended": m.extended,
             "independent": k in unique, "centroid": [_clean(float(x)) for x in m.centroid], "labels": list(m.labels)}
            for k, m in enumerate(analysis.structure.molecules)]


def stacks_list(analysis: Analysis) -> list[dict]:
    """Le colonne di stack (aggregati, colonne, strati, reti 3D) con molecole, traslazioni e passo."""
    return [{"kind": c.kind_key, "rank": c.rank, "molecules": [m + 1 for m in c.molecules],
             "translations": [list(t) for t in c.translations], "repeat": None if c.repeat is None else _clean(c.repeat),
             "pairs": sorted({e[3].index + 1 for e in c.edges}), "mean_interplanar": _clean(c.mean("interplanar"))}
            for c in analysis.stacks]


def ring_pairs_list(analysis: Analysis) -> list[dict]:
    """Le coppie di anelli (stile Janiak) di ogni coppia di molecole, con il numero `pair` della coppia."""
    from dataclasses import asdict

    from pistak.core.ringpairs import ring_pairs

    out = []
    for p in analysis.pairs:
        for r in ring_pairs(analysis.pair_context(p)):
            d = {k: (list(v) if isinstance(v, tuple) else _clean(v)) for k, v in asdict(r).items()}
            out.append({"pair": p.index + 1, **d})
    return out


def geometry_dict(analysis: Analysis, include_h: bool = True) -> dict:
    """Legami, angoli, torsioni e contatti di van der Waals delle molecole indipendenti, come la scheda Geometria:
    per ogni misura gli atomi, l'operazione che genera ciascuno (vuota per l'atomo del CIF), valore e s.u.; per i
    contatti anche somma dei raggi e differenza."""
    from pistak.core.intramolecular import molecule_geometry

    out: dict[str, list] = {"bonds": [], "angles": [], "torsions": [], "contacts": []}
    plural = {"bond": "bonds", "angle": "angles", "torsion": "torsions", "contact": "contacts"}
    for g in molecule_geometry(analysis.structure, analysis.cif, include_h=include_h):
        d = {"molecule": g.molecule + 1, "atoms": list(g.atoms), "symops": list(g.codes or [""] * len(g.atoms)),
             "value": _clean(g.value), "su": _clean(g.su)}
        if g.kind == "contact":
            d |= {"molecule_b": g.places[1][0] + 1 if g.places else None, "vdw_sum": _clean(g.vdw), "delta": _clean(g.delta)}
        out[plural[g.kind]].append(d)
    return out


def write_json(analysis: Analysis, path: str | Path, keys: list[str] | None = None) -> Path:
    """Tutto il risultato: riepilogo (cella, soglie, grandezze della struttura, avvisi), molecole, operazioni di
    simmetria, frammenti π, coppie con i parametri (e le σ se calcolate), contatti tra molecole, colonne di stack,
    coppie di anelli e geometria delle molecole indipendenti con i contatti di van der Waals."""
    path = Path(path)
    keys = keys if keys is not None else list(analysis.param_keys)
    units = {k: (get_spec(k).unit if get_spec(k) else "") for k in keys}
    structure_units = {k: (s.unit if (s := get_structure_spec(k)) else "") for k in analysis.structure_values}
    from pistak import __version__

    data = {"format": JSON_FORMAT, "version": JSON_VERSION, "pistak": __version__, "summary": summary_dict(analysis),
            "units": units, "structure_units": structure_units,
            "symops": [op.text for op in analysis.cif.symops], "molecule_list": molecules_list(analysis),
            "fragments": fragments_list(analysis), "pairs": [pair_dict(p, keys, analysis.uncertainty) for p in analysis.pairs],
            "molecule_contacts": contacts_list(analysis), "stacks": stacks_list(analysis), "ring_pairs": ring_pairs_list(analysis),
            "geometry": geometry_dict(analysis)}
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
