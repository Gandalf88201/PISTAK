"""Checklist, sezione 1: lettura dei CIF di PISTAK contro gemmi (oracolo esterno) sul campione COD di `scarica_cod.py`,
più i casi limite di `cif_su` ed `element_of` controllati a mano.

Per ogni CIF: cella, numero di siti, coordinate, occupanze ed elemento di ogni sito (per etichetta), operazioni di
simmetria (come insiemi di R, t modulo il reticolo). PISTAK legge con `disorder="all"` per avere tutti i siti.

    .venv/bin/python validation/check_cif.py
"""
from __future__ import annotations

import math
import sys
from collections import Counter
from pathlib import Path

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from pistak.core.cif import cif_su, element_of, parse_cif  # noqa: E402
from pistak.core.errors import CifError  # noqa: E402


def op_key(r, t) -> tuple:
    return tuple(int(round(x)) for x in np.ravel(r)) + tuple(int(round(24 * (x % 1))) % 24 for x in t)


def gemmi_symops(block: gemmi.cif.Block) -> set | None:
    for tag in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"):
        col = block.find_values(tag)
        if len(col):
            out = set()
            for s in col:
                op = gemmi.Op(gemmi.cif.as_string(s))
                out.add(op_key(np.array(op.rot) / 24, np.array(op.tran) / 24))
            return out
    return None


def compare(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    problems: list[str] = []
    try:
        mine = parse_cif(text, disorder="all")
    except CifError as e:
        return [f"PISTAK non legge: {e}"]
    doc = gemmi.cif.read_string(text)
    block = next(b for b in doc if b.find_value("_cell_length_a") and len(b.find_values("_atom_site_fract_x")))
    if block.name != mine.name:
        problems.append(f"blocco diverso: {mine.name} / {block.name}")
    ss = gemmi.make_small_structure_from_block(block)
    ref_cell = ss.cell.parameters
    my_cell = (mine.cell.a, mine.cell.b, mine.cell.c, mine.cell.alpha, mine.cell.beta, mine.cell.gamma)
    if max(abs(x - y) for x, y in zip(my_cell, ref_cell)) > 1e-6:
        problems.append(f"cella {my_cell} / {ref_cell}")
    ref_sites = {s.label: s for s in ss.sites}
    my_sites = {a.label: a for a in mine.atoms}
    if len(ss.sites) != len(mine.atoms):
        only_ref = sorted(set(ref_sites) - set(my_sites))
        only_mine = sorted(set(my_sites) - set(ref_sites))
        problems.append(f"siti {len(mine.atoms)} / {len(ss.sites)}; solo gemmi {only_ref[:6]}, solo PISTAK {only_mine[:6]}")
    for lab, a in my_sites.items():
        s = ref_sites.get(lab)
        if s is None:
            continue
        if np.abs(a.frac - np.array(s.fract.tolist())).max() > 1e-9:
            problems.append(f"coordinate di {lab}")
        if abs(a.occ - s.occ) > 1e-9:
            problems.append(f"occupanza di {lab}: {a.occ} / {s.occ}")
        ref_el = s.element.name
        if ref_el == "D":
            ref_el = "H"
        if a.el != ref_el:
            problems.append(f"elemento di {lab}: {a.el} / {ref_el} (type_symbol «{s.type_symbol}»)")
    ref_ops = gemmi_symops(block)
    if mine.symops_source == "cif":
        my_ops = {op_key(o.r, o.t) for o in mine.symops}
        if ref_ops is not None and my_ops != ref_ops:
            problems.append(f"operazioni {len(my_ops)} / {len(ref_ops)}")
    return problems


def check_sample() -> None:
    files = sorted((HERE / "cod_cache").glob("[0-9]*.cif"))
    print(f"== parse_cif contro gemmi su {len(files)} CIF della COD ==")
    kinds: Counter = Counter()
    bad = 0
    for f in files:
        problems = compare(f)
        if problems:
            bad += 1
            print(f"  {f.stem}: " + "; ".join(problems[:4]) + (f" (+{len(problems) - 4})" if len(problems) > 4 else ""))
            kinds.update(p.split(" ")[0] for p in problems)
    print(f"  {len(files) - bad} uguali, {bad} con differenze; per tipo: {dict(kinds)}")


def check_su() -> None:
    print("\n== cif_su: incertezze tra parentesi, calcolate a mano ==")
    cases = {"1.234(5)": 0.005, "1.2345(12)": 0.0012, "123(4)": 4.0, "0.5": 0.0, "-0.0123(4)": 0.0004, ".5(1)": 0.1,
             "1.5(12)": 1.2, "10.3528(16)": 0.0016, "0.12345(100)": 0.001, "1.2e-3(2)": 0.2e-3, "?": 0.0, ".": 0.0}
    bad = 0
    for v, expected in cases.items():
        got = cif_su(v)
        ok = math.isclose(got, expected, rel_tol=1e-12, abs_tol=1e-15)
        bad += not ok
        print(f"  {v:14s} → {got:<10g} atteso {expected:<10g} {'' if ok else 'DIVERSO'}")
    print(f"  {len(cases)} casi, {bad} diversi")


def check_elements() -> None:
    print("\n== element_of: etichette ambigue ==")
    print("  attesi: type_symbol se c'è; senza, etichette tutte maiuscole a una lettera tranne CL e BR (CA1 è un C-alfa,")
    print("  CD1 un C-delta); vuoto (sito scartato) per centroidi Cg e picchi Q")
    cases = [  # (type_symbol, etichetta, atteso)
        (None, "C1", "C"), (None, "Cl1", "Cl"), (None, "CL1", "Cl"), (None, "C1A", "C"), (None, "N1'", "N"),
        (None, "Co1", "Co"), (None, "CO1", "C"), (None, "CA1", "C"), (None, "CD1", "C"), (None, "HG1", "H"), (None, "D1", "H"), (None, "H1A", "H"), (None, "Br2", "Br"),
        (None, "BR2", "Br"), (None, "Cg1", ""), (None, "Q1", ""), (None, "O1W", "O"), (None, "Fe1", "Fe"),
        ("Fe2+", "Fe1", "Fe"), ("O-2", "O1", "O"), ("D", "D1", "H"), ("Cl", "CL1", "Cl"), ("C", "Co1", "C"),
        ("Cl1-", "Cl1", "Cl"), ("Q", "Q1", ""),
    ]
    for sym, lab, expected in cases:
        got = element_of(sym, lab)
        flag = "" if got == expected else "DIVERSO"
        print(f"  type_symbol {sym!s:6s} etichetta {lab:5s} → {got or '(vuoto)':7s} atteso {expected or '(vuoto)':7s} {flag}")


def expected_major(block: gemmi.cif.Block) -> tuple[set, list[str]]:
    """Regola `major` ricalcolata dalle colonne lette da gemmi: per ogni assemblaggio si tiene il gruppo (preso in
    valore assoluto: −1 e 1 sono lo stesso gruppo su un elemento di simmetria) con l'occupanza media più alta, a
    parità il numero più basso; gli atomi senza gruppo si tengono sempre."""
    tab = block.find("_atom_site_", ["label", "?disorder_group", "?disorder_assembly", "?occupancy", "?calc_flag"])
    rows = []
    for r in tab:
        lab = gemmi.cif.as_string(r[0])
        grp = gemmi.cif.as_string(r[1]) if r.has(1) else "."
        asm = gemmi.cif.as_string(r[2]) if r.has(2) else "."
        occ = gemmi.cif.as_number(r[3]) if r.has(3) and r[3] not in ("?", ".") else 1.0
        flag = gemmi.cif.as_string(r[4]) if r.has(4) else ""
        if flag.lower() == "dum":
            continue
        rows.append((lab, grp if grp not in ("", "?") else ".", asm if asm not in ("", "?") else ".", occ))
    by_group: dict = {}
    for _lab, g, a, o in rows:
        if g != ".":
            by_group.setdefault(a, {}).setdefault(g.lstrip("-"), []).append(o)
    best = {a: min(gs, key=lambda g: (-sum(gs[g]) / len(gs[g]), int(g) if g.isdigit() else 10**6, g)) for a, gs in by_group.items()}
    keep = {lab for lab, g, a, _ in rows if g == "." or best[a] == g.lstrip("-")}
    notes = [f"assemblaggio {a}: " + ", ".join(f"gruppo {g} occ. media {sum(v) / len(v):.3f} ({len(v)} atomi)" for g, v in sorted(gs.items()))
             + f" → {best[a]}" for a, gs in by_group.items()]
    return keep, notes


def check_disorder(show: int = 3) -> None:
    files = sorted((HERE / "cod_cache").glob("[0-9]*.cif"))
    print("\n== disordine, regola major: PISTAK contro la regola ricalcolata ==")
    n = bad = shown = 0
    for f in files:
        text = f.read_text(encoding="utf-8", errors="replace")
        block = next(b for b in gemmi.cif.read_string(text) if len(b.find_values("_atom_site_fract_x")))
        groups = [gemmi.cif.as_string(v) for v in block.find_values("_atom_site_disorder_group")]
        if not any(g not in (".", "?", "") for g in groups):
            continue
        n += 1
        keep, notes = expected_major(block)
        mine = parse_cif(text, disorder="major")
        got = {a.label for a in mine.atoms}
        neg = any(g.startswith("-") for g in groups)
        if got != keep & {a.label for a in parse_cif(text, disorder="all").atoms}:
            bad += 1
            print(f"  {f.stem}: DIVERSO; solo PISTAK {sorted(got - keep)[:5]}, solo regola {sorted(keep - got)[:5]}")
        if shown < show or (neg and shown < show + 2):
            shown += 1
            print(f"  {f.stem}{' (gruppi negativi)' if neg else ''}: tenuti {len(got)} atomi, scartati {mine.skipped_disorder}")
            for line in notes:
                print(f"      {line}")
    print(f"  {n} CIF con disordine, {bad} diversi")


if __name__ == "__main__":
    which = sys.argv[1:] or ["sample", "su", "elements", "disorder"]
    if "su" in which:
        check_su()
    if "elements" in which:
        check_elements()
    if "sample" in which:
        check_sample()
    if "disorder" in which:
        check_disorder()
