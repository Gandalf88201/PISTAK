"""Esportazione dei risultati: tabella CSV, JSON completo, coppia in formato XYZ.

I file scritti qui hanno intestazioni e chiavi fisse (inglesi, mai tradotte): uno script che li legge funziona
qualunque sia la lingua dell'interfaccia. I nomi leggibili dei parametri stanno nel registro (`pistak params`).
"""
from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

from pistak.core.analysis import Analysis
from pistak.core.pairs import Pair
from pistak.params import get_spec, get_structure_spec

#: Formato del JSON scritto da `write_json`: cambia `JSON_VERSION` se cambia la struttura.
JSON_FORMAT, JSON_VERSION = "pistak-analysis", 5


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
            "include_cell": analysis.uncertainty.include_cell, "has_su": analysis.uncertainty.has_su,
            "covariance": analysis.uncertainty.covariance, "reported_value": "nominal"},
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
        "equivalent_pairs": p.count, "polymer": p.polymer, "bonded": p.bonded, "bridged": p.bridged, "clamped": p.clamped,
        "long_axis_a": p.axis_a, "long_axis_b": p.axis_b, "parameters": {k: _clean(p.get(k)) for k in keys},
    }
    if uncertainty is not None:
        d["su"] = {k: _clean(uncertainty.get(p, k)) for k in keys if k in uncertainty.keys}
        d["kind_agreement"] = _clean(uncertainty.kind_agreement.get(p.index, math.nan))
        if p.index in uncertainty.incomplete:  # σ da meno campioni: quanti hanno dato un valore
            d["su_samples"] = dict(uncertainty.incomplete[p.index])
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


def _xyz_text(atoms: list, head: str) -> str:
    body = "\n".join(f"{e:<2} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}" for e, p in atoms)
    return f"{len(atoms)}\n{head}\n{body}\n"


def same_molecule(pair: Pair) -> bool:
    """Vero se la coppia è tra due frammenti π della stessa molecola, nella stessa cella: le «due molecole» sono la
    stessa, una sopra l'altra (con una traslazione, come nelle catene di un polimero, sono due copie distinte)."""
    return pair.a == pair.b and not any(pair.shift)


def pair_xyz(analysis: Analysis, pair: Pair) -> str:
    """Le due molecole della coppia, in coordinate cartesiane (Å), nel formato XYZ. Se i due frammenti sono della
    stessa molecola (`same_molecule`), la molecola compare una volta sola, non due volte sovrapposta."""
    ctx = analysis.pair_context(pair)
    atoms = [(e, p) for e, p in zip(ctx.A.el, ctx.pos_a)]
    if not same_molecule(pair):
        atoms += [(e, p) for e, p in zip(ctx.B.el, ctx.pos_b)]
    head = (f"{analysis.name} pair={pair.index + 1} kind={pair.kind} interplanar={pair.interplanar:.3f} "
            f"long_slip={pair.long_slip:.3f} lat_slip={pair.lat_slip:.3f}")
    if same_molecule(pair):
        head += " same_molecule"
    return _xyz_text(atoms, head)


def write_pair_xyz(analysis: Analysis, pair: Pair, path: str | Path) -> Path:
    path = Path(path)
    path.write_text(pair_xyz(analysis, pair), encoding="utf-8")
    return path


#: Il tipo di coppia nei nomi dei file XYZ: fisso e in inglese, come le chiavi degli altri file esportati.
XYZ_KIND_NAMES = {"stack": "pistaked", "parallel": "parallel", "tilted": "tilted"}


def xyz_stem(analysis: Analysis) -> str:
    """Il nome della struttura nei nomi dei file: quello del CIF senza estensione, con soli caratteri sicuri."""
    return re.sub(r"[^\w.-]+", "_", Path(analysis.name).stem) or "pistak"


@dataclass
class Dimer:
    """Due molecole vicine e tutte le coppie di frammenti π che le legano, nell'ordine della tabella. Una molecola con
    più frammenti dà una coppia per ogni coppia di frammenti vicini, di tipo anche diverso: tutte descrivono lo stesso
    dimero. La prima (il tipo più interessante, poi la distanza più corta) dà nome e coordinate."""

    pairs: list[Pair]
    #: «<tipo>-<n>», come nel nome del file (pistaked-1, parallel-3, …): la colonna «Dimero» della tabella.
    label: str = ""
    name: str = ""

    @property
    def kind(self) -> str:
        return self.pairs[0].kind


def _dimer_image(a: int, b: int, shift) -> tuple:
    """Lo stesso dimero visto da A o da B: (a, b, t) e (b, a, −t) sono le stesse due molecole."""
    t = tuple(int(v) for v in shift)
    return min((a, b, t), (b, a, tuple(-v for v in t)))


def dimers(analysis: Analysis) -> list[Dimer]:
    """Le coppie raccolte per dimero di molecole. Due coppie sono dello stesso dimero se una delle loro copie
    equivalenti per simmetria (`Pair.members`) unisce le stesse due molecole con la stessa traslazione. Le coppie tra
    frammenti della stessa molecola (`same_molecule`) non sono dimeri e restano fuori. Ordine: quello della prima
    coppia di ogni dimero nella tabella; il nome è «<struttura>-<tipo>-<n>.xyz» con n che conta i dimeri di quel tipo."""
    pairs = [p for p in analysis.pairs if not same_molecule(p)]
    parent = list(range(len(pairs)))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owner: dict[tuple, int] = {}
    for i, p in enumerate(pairs):
        for m in p.members or [p.key]:
            img = _dimer_image(m[0], m[1], m[2])
            if img in owner:
                ra, rb = root(i), root(owner[img])
                parent[max(ra, rb)] = min(ra, rb)  # la radice resta la coppia che viene prima nella tabella
            else:
                owner[img] = i
    groups: dict[int, list[Pair]] = {}
    for i, p in enumerate(pairs):
        groups.setdefault(root(i), []).append(p)
    out = [Dimer(g) for _r, g in sorted(groups.items(), key=lambda kv: kv[1][0].index)]
    stem, seen = xyz_stem(analysis), dict.fromkeys(XYZ_KIND_NAMES, 0)
    for d in out:
        seen[d.kind] = seen.get(d.kind, 0) + 1
        d.label = f"{XYZ_KIND_NAMES.get(d.kind, d.kind)}-{seen[d.kind]}"
        d.name = f"{stem}-{d.label}.xyz"
    return out


def dimer_of_pairs(analysis: Analysis) -> dict[int, tuple[int, Dimer]]:
    """Per ogni coppia (indice nella tabella) il suo dimero e la posizione del dimero in `dimers`; le coppie tra
    frammenti della stessa molecola non ci sono. Serve alla tabella e allo scatter per legare righe e punti ai file."""
    return {p.index: (k, d) for k, d in enumerate(dimers(analysis)) for p in d.pairs}


def dimer_xyz(analysis: Analysis, dimer: Dimer) -> str:
    """Il dimero in formato XYZ, con le coordinate della sua prima coppia. La riga di commento ha i valori della prima
    coppia (`pair=`, come l'XYZ di una coppia) e l'elenco di tutte le sue coppie: numeri nella tabella (`pairs=`), tipi
    (`kinds=`) e frammenti, numerati da 1 come nel JSON (`fragments=`)."""
    text = pair_xyz(analysis, dimer.pairs[0])
    n, head, body = text.split("\n", 2)
    head += (f" pairs={','.join(str(p.index + 1) for p in dimer.pairs)}"
             f" kinds={','.join(p.kind for p in dimer.pairs)}"
             f" fragments={','.join(f'{p.fa + 1}-{p.fb + 1}' for p in dimer.pairs)}")
    return f"{n}\n{head}\n{body}"


def write_dimers_xyz(analysis: Analysis, folder: str | Path) -> list[Path]:
    """Un file XYZ per ogni dimero distinto (`dimers`), nella cartella `folder` (creata se manca). Restituisce i file
    scritti, nell'ordine dei dimeri."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    out = []
    for d in dimers(analysis):
        path = folder / d.name
        path.write_text(dimer_xyz(analysis, d), encoding="utf-8")
        out.append(path)
    return out


def pair_files(analysis: Analysis) -> list[tuple[Pair, str]]:
    """Ogni coppia della tabella con il nome del suo file XYZ, «<struttura>-<tipo>-<n>.xyz», n che conta le coppie di
    quel tipo nell'ordine della tabella. Gli stessi nomi dei dimeri ma un'altra numerazione: vanno in un'altra
    cartella."""
    stem, seen = xyz_stem(analysis), dict.fromkeys(XYZ_KIND_NAMES, 0)
    out = []
    for p in analysis.pairs:
        seen[p.kind] = seen.get(p.kind, 0) + 1
        out.append((p, f"{stem}-{XYZ_KIND_NAMES.get(p.kind, p.kind)}-{seen[p.kind]}.xyz"))
    return out


def write_pairs_xyz(analysis: Analysis, folder: str | Path) -> list[Path]:
    """Un file XYZ per ogni coppia della tabella (`pair_files`), nella cartella `folder` (creata se manca). Le coppie
    dello stesso dimero danno file con le stesse due molecole: la riga di commento dice il dimero (`dimer=`); due
    frammenti della stessa molecola danno la molecola da sola (`same_molecule`)."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    dimer_of = dimer_of_pairs(analysis)
    out = []
    for p, name in pair_files(analysis):
        n, head, body = pair_xyz(analysis, p).split("\n", 2)
        if p.index in dimer_of:
            head += f" dimer={dimer_of[p.index][1].label}"
        path = folder / name
        path.write_text(f"{n}\n{head}\n{body}", encoding="utf-8")
        out.append(path)
    return out


@dataclass
class XyzCounts:
    """Quanti file dà ciascuna esportazione XYZ e perché differiscono: il testo dell'avviso lo fa l'interfaccia."""

    pairs: int  # coppie della tabella: un file ciascuna
    dimers: int  # dimeri distinti: un file ciascuno
    repeated: int  # file per coppia che ripetono un dimero già scritto da un'altra coppia
    same_molecule: int  # coppie tra due frammenti della stessa molecola
    example: Dimer | None  # il dimero con più coppie, per l'esempio dell'avviso (None se nessun dimero ne ha due)


def xyz_counts(analysis: Analysis) -> XyzCounts:
    ds = dimers(analysis)
    same = sum(1 for p in analysis.pairs if same_molecule(p))
    big = max(ds, key=lambda d: len(d.pairs), default=None)
    return XyzCounts(len(analysis.pairs), len(ds), sum(len(d.pairs) - 1 for d in ds), same,
                     big if big is not None and len(big.pairs) > 1 else None)
