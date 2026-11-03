"""Checklist, sezioni 1 (crystal.py), 2 (struttura e molecole) e 3 (anelli) sul campione COD di `scarica_cod.py`.

    cell      densità e contenuto della cella contro i valori del CIF (_exptl_crystal_density_diffrn, Z × formula)
    moiety    formule delle molecole ricostruite contro _chemical_formula_moiety
    oblique   la stessa struttura in una cella obliqua equivalente (det = 1): stesse molecole, stesse coppie
    rings     anelli di `small_rings` contro una base minima dei cicli (networkx)
    radii     raggi covalenti e di van der Waals contro le tabelle di gemmi (confronto, non oracolo: fonti non dichiarate)

    .venv/bin/python validation/check_struttura.py [cell moiety oblique rings radii]
"""
from __future__ import annotations

import math
import re
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from pistak.core.analysis import analyse_text  # noqa: E402
from pistak.core.cif import parse_cif  # noqa: E402
from pistak.core.config import Limits  # noqa: E402
from pistak.core.crystal import cell_contents, cell_volume, cif_value, density, parse_formula, reported_formula  # noqa: E402
from pistak.core.structure import build_structure  # noqa: E402

FILES = sorted((HERE / "cod_cache").glob("[0-9]*.cif"))


def text_of(f: Path) -> str:
    return f.read_text(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------------------------------------
def check_cell() -> None:
    print("== densità e contenuto della cella contro il CIF ==")
    rel: list[float] = []
    bad_d, bad_z, n_z = [], [], 0
    for f in FILES:
        cif = parse_cif(text_of(f), disorder="major")
        contents = cell_contents(cif)
        d_cif = cif_value(cif, "density_cif")
        if math.isfinite(d_cif):
            r = density(contents, cell_volume(cif)) / d_cif - 1
            rel.append(r)
            if abs(r) > 0.005:
                bad_d.append((f.stem, round(100 * r, 2)))
        z, formula = cif_value(cif, "z"), reported_formula(cif)
        if math.isfinite(z) and formula:
            n_z += 1
            ratios = {e: contents.get(e, 0.0) / (z * n) for e, n in formula.items() if e != "H"}
            if any(abs(x - 1) > 0.01 for x in ratios.values()):
                bad_z.append((f.stem, {e: round(x, 3) for e, x in ratios.items() if abs(x - 1) > 0.01}))
    a = np.abs(rel)
    print(f"  densità: {len(rel)} CIF, scarto mediano {100 * np.median(a):.3f} %, 95° percentile {100 * np.percentile(a, 95):.3f} %, "
          f"oltre 0,5 %: {len(bad_d)}")
    for x in bad_d:
        print(f"    {x[0]}: {x[1]} %")
    print(f"  contenuto della cella / (Z × formula), atomi pesanti: {n_z} CIF, diversi oltre 1 %: {len(bad_z)}")
    for x in bad_z[:15]:
        print(f"    {x[0]}: {x[1]}")


# ---------------------------------------------------------------------------------------------------------
_CHARGE = re.compile(r"^\d*[+-]\d*$")


def moiety_formulas(text: str) -> list[tuple[int, Counter]]:
    """'C20 H24 N4 2+, 2(Cl -), 0.5(C H2 Cl2)' → [(1, C20H24N4), (2, Cl), (0.5, CH2Cl2)] (senza cariche)."""
    out = []
    for part in re.split(r",(?![^(]*\))", text):
        part = part.strip()
        m = re.match(r"^(\d*\.?\d*)\s*\((.*)\)\s*([+-]?\d*[+-]?)$", part)
        mult = 1.0
        if m:
            mult = float(m.group(1)) if m.group(1) else 1.0
            part = m.group(2)
        tokens = [t for t in part.split() if not _CHARGE.match(t)]
        counts = Counter({e: n for e, n in parse_formula(" ".join(tokens)).items()})
        if counts:
            out.append((mult, counts))
    return out


def key_of(c: Counter) -> str:
    return " ".join(f"{e}{c[e]:g}" for e in sorted(c, key=lambda e: (e != "C", e != "H", e)))


def check_moiety() -> None:
    print("\n== formule delle molecole contro _chemical_formula_moiety ==")
    n = ok = only_h = 0
    rows = []
    for f in FILES:
        text = text_of(f)
        cif = parse_cif(text, disorder="major")
        moiety = cif.items.get("_chemical_formula_moiety")
        if not moiety or moiety in ("?", "."):
            continue
        n += 1
        s = build_structure(cif, Limits())
        mine = {key_of(Counter(dict(m.formula()))) for m in s.molecules if not m.extended}
        ext = sum(m.extended for m in s.molecules)
        ref = {key_of(c) for _, c in moiety_formulas(moiety)}
        if mine == ref:
            ok += 1
            continue
        strip_h = lambda keys: {re.sub(r"\bH[\d.]+\s?", "", k).strip() for k in keys}  # noqa: E731
        if strip_h(mine) == strip_h(ref):
            only_h += 1
            kind = "solo H"
        else:
            kind = "diverse"
        rows.append((f.stem, kind, sorted(mine), sorted(ref), ext))
    print(f"  {n} CIF con la formula delle unità: {ok} uguali, {only_h} diverse solo negli H, {len(rows) - only_h} diverse")
    for stem, kind, mine, ref, ext in rows:
        if kind == "diverse":
            print(f"    {stem}: PISTAK {mine[:4]}{f' (+{len(mine) - 4})' if len(mine) > 4 else ''}"
                  f"{f', {ext} reti estese' if ext else ''}; CIF {ref[:4]}")


# ---------------------------------------------------------------------------------------------------------
def _frac_str(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def transformed_cif(text: str, T: np.ndarray) -> str:
    """Lo stesso cristallo nella cella a' = T·a (righe di T: nuovi vettori in funzione dei vecchi, det = 1). Coordinate
    f' = f·T⁻¹, operazioni R' = T⁻ᵀ R Tᵀ, t' = T⁻ᵀ t. Si scrive un CIF minimo: cella, operazioni, atomi con s.u. tolte."""
    cif = parse_cif(text, disorder="major")
    from pistak.core.structure import cell_vectors

    M = cell_vectors(cif.cell)
    M2 = T @ M
    a, b, c = (np.linalg.norm(v) for v in M2)
    ang = lambda u, v: math.degrees(math.acos(u @ v / np.linalg.norm(u) / np.linalg.norm(v)))  # noqa: E731
    Ti = np.linalg.inv(T)
    TiT = Ti.T
    ops = []
    for op in cif.symops:
        R = np.rint(TiT @ op.r @ T.T).astype(int)
        t = TiT @ op.t
        comps = []
        for row, ti in zip(R, t):
            s = "".join(f"{'+' if v > 0 else '-'}{abs(v) if abs(v) != 1 else ''}{x}" for v, x in zip(row, "xyz") if v)
            f = Fraction(float(ti)).limit_denominator(24)
            comps.append((s.lstrip("+") or "0") + (f"+{_frac_str(f)}" if f > 0 else f"-{_frac_str(-f)}" if f < 0 else ""))
        ops.append("'" + ",".join(comps) + "'")
    atoms = [f"{at.label} {at.el} {x:.8f} {y:.8f} {z:.8f} {at.occ:g}" for at in cif.atoms for x, y, z in [at.frac @ Ti]]
    return "\n".join([
        "data_obliqua", f"_cell_length_a {a:.8f}", f"_cell_length_b {b:.8f}", f"_cell_length_c {c:.8f}",
        f"_cell_angle_alpha {ang(M2[1], M2[2]):.8f}", f"_cell_angle_beta {ang(M2[0], M2[2]):.8f}",
        f"_cell_angle_gamma {ang(M2[0], M2[1]):.8f}", "loop_", "_symmetry_equiv_pos_as_xyz", *ops,
        "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
        "_atom_site_fract_z", "_atom_site_occupancy", *atoms, ""])


def summary(an) -> tuple:
    mols = Counter(key_of(Counter(dict(m.formula()))) for m in an.structure.molecules if not m.extended)
    pairs = Counter((p.kind, round(p.interplanar, 3), round(p.slip, 3), round(p.contact, 3), p.count) for p in an.pairs)
    return mols, pairs


def check_oblique(n: int = 40) -> None:
    print("\n== stessa struttura in una cella obliqua equivalente (a' = a, b' = b + 2a, c' = c + a + b) ==")
    T = np.array([[1, 0, 0], [2, 1, 0], [1, 1, 1]])
    same = diff = 0
    for f in FILES[:n]:
        text = text_of(f)
        lim = Limits()
        try:
            base = analyse_text(text, f.stem, lim)
        except Exception as e:  # noqa: BLE001
            print(f"  {f.stem}: non analizzabile ({type(e).__name__})")
            continue
        obl = analyse_text(transformed_cif(text, T), f.stem, lim)
        (m0, p0), (m1, p1) = summary(base), summary(obl)
        if m0 == m1 and p0 == p1:
            same += 1
        else:
            diff += 1
            print(f"  {f.stem}: molecole {sum(m0.values())}→{sum(m1.values())}, coppie {len(base.pairs)}→{len(obl.pairs)}"
                  f"{'' if m0 == m1 else ', formule diverse'}")
    print(f"  {same + diff} strutture: {same} identiche, {diff} diverse")


# ---------------------------------------------------------------------------------------------------------
def check_rings(n: int = 200) -> None:
    import networkx as nx

    from pistak.core.cores import small_rings

    print("\n== small_rings contro la base minima dei cicli di networkx (atomi pesanti, anelli fino a 8) ==")
    mols = missing = extra = 0
    examples = []
    for f in FILES[:n]:
        s = build_structure(parse_cif(text_of(f), disorder="major"), Limits())
        seen = set()
        for m in s.molecules:
            key = key_of(Counter(dict(m.formula())))
            if m.extended or key in seen:
                continue
            seen.add(key)
            heavy = m.heavy
            g = nx.Graph()
            g.add_edges_from((int(p), int(q)) for p, q in m.bonds if heavy[p] and heavy[q])
            ref = {frozenset(c) for c in nx.minimum_cycle_basis(g) if len(c) <= 8}
            mine = {frozenset(r) for r in small_rings(len(m), m.bonds, lambda i, h=heavy: h[i], 8)}
            mols += 1
            lost = ref - mine
            missing += len(lost)
            extra += len(mine - ref)
            if lost:
                examples.append((f.stem, key, [len(r) for r in lost]))
    print(f"  {mols} molecole distinte: anelli della base minima non trovati {missing}; anelli in più (cicli più corti di "
          f"qualche legame fuori dalla base, non errori) {extra}")
    for x in examples[:10]:
        print(f"    {x[0]} {x[1]}: mancano anelli di {x[2]} atomi")


def check_radii() -> None:
    import gemmi

    from pistak.core.elements import COVALENT_RADII, VDW_RADII

    print("\n== raggi: PISTAK contro le tabelle di gemmi (scarti oltre 0,02 Å) ==")
    for name, table, attr in (("covalenti", COVALENT_RADII, "covalent_r"), ("van der Waals", VDW_RADII, "vdw_r")):
        d = [(e, r, getattr(gemmi.Element(e), attr)) for e, r in table.items()]
        off = [(e, r, round(g, 3)) for e, r, g in d if abs(r - g) > 0.02]
        print(f"  {name}: {len(d)} elementi, {len(off)} diversi: {off}")


CHECKS = {"cell": check_cell, "moiety": check_moiety, "oblique": check_oblique, "rings": check_rings, "radii": check_radii}

if __name__ == "__main__":
    for name in sys.argv[1:] or CHECKS:
        CHECKS[name]()
