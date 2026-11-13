"""Checklist, sezioni 1 (crystal.py), 2 (struttura e molecole) e 3 (anelli) sul campione COD di `scarica_cod.py`.

    cell      densità e contenuto della cella contro i valori del CIF (_exptl_crystal_density_diffrn, Z × formula)
    moiety    formule delle molecole ricostruite contro _chemical_formula_moiety
    oblique   la stessa struttura in una cella obliqua equivalente (det = 1): stesse molecole, stesse coppie
    rings     anelli di `small_rings` contro una base minima dei cicli (networkx)
    rdkit     legami e anelli contro RDKit: legami rifatti dalle sole coordinate (DetermineConnectivity), anelli dalla
              SSSR simmetrizzata, sia sui legami di PISTAK sia su quelli di RDKit
    radii     raggi covalenti e di van der Waals contro le tabelle di gemmi (confronto, non oracolo: fonti non dichiarate)
    hthresh   voce 23: le soglie del controllo sugli idrogeni (`TRIPLE_CX`, `SP_ANGLE`) contro le distanze C=O dei
              carbonili e gli angoli dei C a due vicini misurati con numpy sulle molecole del campione

    .venv/bin/python validation/check_struttura.py [cell moiety oblique rings rdkit radii hthresh]
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
        if math.isfinite(z) and z > 0 and formula and all(n > 0 for n in formula.values()):  # Z = 0 nel CIF: non confrontabile
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


def _with_groups(text: str) -> Counter:
    """Formula con gruppi tra parentesi e moltiplicatore: '[Sn (C6 H5) (C H3)2]' → Sn C8 H11 (i due metili contano due
    volte; `parse_formula` di PISTAK legge solo _chemical_formula_sum, che le parentesi non le ha)."""
    total: Counter = Counter()
    for inner, mult in re.findall(r"\(([^()]*)\)\s*(\d*)", text):
        for e, n in parse_formula(inner).items():
            total[e] += n * (int(mult) if mult else 1)
    rest = re.sub(r"\(([^()]*)\)\s*\d*", " ", text).replace("[", " ").replace("]", " ")
    for e, n in parse_formula(rest).items():
        total[e] += n
    return total


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
        counts = _with_groups(" ".join(tokens))
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


def transformed_cif(text: str, T: np.ndarray, shift: np.ndarray | None = None) -> str:
    """Lo stesso cristallo nella cella a' = T·a (righe di T: nuovi vettori in funzione dei vecchi, det = 1) e con
    l'origine spostata di −`shift` (coordinate nuove). Coordinate f' = f·T⁻¹ + s, operazioni R' = T⁻ᵀ R Tᵀ,
    t' = T⁻ᵀ t + (I − R') s. Si scrive un CIF minimo: cella, operazioni, atomi con s.u. tolte."""
    origin = np.zeros(3) if shift is None else np.asarray(shift, float)
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
        t = TiT @ op.t + (np.eye(3) - R) @ origin
        comps = []
        for row, ti in zip(R, t):
            s = "".join(f"{'+' if v > 0 else '-'}{abs(v) if abs(v) != 1 else ''}{x}" for v, x in zip(row, "xyz") if v)
            f = Fraction(float(ti) % 1.0 if abs(ti) > 1e-9 else 0.0).limit_denominator(1000000)
            comps.append((s.lstrip("+") or "0") + (f"+{_frac_str(f)}" if f > 0 else f"-{_frac_str(-f)}" if f < 0 else ""))
        ops.append("'" + ",".join(comps) + "'")
    atoms = [f"{at.label} {at.el} {x:.8f} {y:.8f} {z:.8f} {at.occ:g}" for at in cif.atoms for x, y, z in [at.frac @ Ti + origin]]
    return "\n".join([
        "data_obliqua", f"_cell_length_a {a:.8f}", f"_cell_length_b {b:.8f}", f"_cell_length_c {c:.8f}",
        f"_cell_angle_alpha {ang(M2[1], M2[2]):.8f}", f"_cell_angle_beta {ang(M2[0], M2[2]):.8f}",
        f"_cell_angle_gamma {ang(M2[0], M2[1]):.8f}", "loop_", "_symmetry_equiv_pos_as_xyz", *ops,
        "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
        "_atom_site_fract_z", "_atom_site_occupancy", *atoms, ""])


def summary(an) -> tuple:
    mols = Counter(key_of(Counter(dict(m.formula()))) for m in an.structure.molecules if not m.extended)
    pairs = sorted((p.kind, p.count, p.interplanar, p.slip, p.contact) for p in an.pairs)
    return mols, pairs


def same_pairs(p0: list, p1: list, tol: float = 1e-6) -> bool:
    """Stesse coppie: tipo e molteplicità uguali, distanze entro `tol` (l'arrotondamento non conta)."""
    if len(p0) != len(p1):
        return False
    free = list(p1)
    for x in p0:  # abbinamento uno a uno: coppie quasi uguali possono scambiarsi di posto nell'ordinamento
        hit = next((k for k, y in enumerate(free)
                    if x[:2] == y[:2] and max(abs(u - v) for u, v in zip(x[2:], y[2:])) < tol), None)
        if hit is None:
            return False
        free.pop(hit)
    return True


def check_oblique(n: int = 40, shift: np.ndarray | None = None) -> None:
    if shift is None:
        print("\n== stessa struttura in una cella obliqua equivalente (a' = a, b' = b + 2a, c' = c + a + b) ==")
        T = np.array([[1, 0, 0], [2, 1, 0], [1, 1, 1]])
    else:
        print(f"\n== stessa struttura con l'origine spostata (coordinate + {shift.tolist()}) ==")
        T = np.eye(3, dtype=int)
    same = diff = 0
    for f in FILES[:n]:
        text = text_of(f)
        lim = Limits()
        try:
            base = analyse_text(text, f.stem, lim)
        except Exception as e:  # noqa: BLE001
            print(f"  {f.stem}: non analizzabile ({type(e).__name__})")
            continue
        obl = analyse_text(transformed_cif(text, T, shift), f.stem, lim)
        (m0, p0), (m1, p1) = summary(base), summary(obl)
        if m0 == m1 and same_pairs(p0, p1):
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


def check_rdkit(n: int = 200) -> None:
    """Oracolo indipendente per legami e anelli: RDKit non usa né i raggi né i legami di PISTAK."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdDetermineBonds
    from rdkit.Geometry import Point3D

    from pistak.core.cores import small_rings

    RDLogger.DisableLog("rdApp.*")
    print("\n== legami e anelli contro RDKit (molecole distinte; anelli fino a 8 atomi, solo atomi pesanti) ==")
    mols = organic = failed = 0
    bonds = Counter()  # (gruppo, uguali/diversi)
    rings = Counter()  # (gruppo, confronto, esito): confronto «algoritmo» = SSSR di RDKit sui legami di PISTAK,
    #                   «catena» = SSSR di RDKit sui legami rifatti da RDKit; tutti contro small_rings di PISTAK
    examples: list[str] = []

    def rdkit_mol(m, bonds=None):
        rw = Chem.RWMol()
        for e in m.el:
            rw.AddAtom(Chem.Atom(e))
        if bonds is not None:
            for p, q in bonds:
                rw.AddBond(int(p), int(q), Chem.BondType.SINGLE)
        conf = Chem.Conformer(len(m.el))
        for i, x in enumerate(m.pos):
            conf.SetAtomPosition(i, Point3D(*map(float, x)))
        rw.AddConformer(conf, assignId=True)
        return rw

    def rings_of(mol, heavy):
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
        return {frozenset(r) for r in Chem.GetSymmSSSR(mol) if len(r) <= 8 and all(heavy[i] for i in r)}

    for f in FILES[:n]:
        s = build_structure(parse_cif(text_of(f), disorder="major"), Limits())
        seen = set()
        for m in s.molecules:
            key = key_of(Counter(dict(m.formula())))
            if m.extended or key in seen:
                continue
            seen.add(key)
            mols += 1
            heavy = m.heavy
            is_organic = set(m.el) <= {"C", "H", "N", "O", "F", "S", "Cl", "Br", "I", "P", "B", "Si", "Se"}
            organic += is_organic
            mine_bonds = {frozenset((int(p), int(q))) for p, q in m.bonds}
            try:
                ref = rdkit_mol(m)
                rdDetermineBonds.DetermineConnectivity(ref)
            except Exception:  # noqa: BLE001  (RDKit rifiuta la molecola: contata a parte)
                failed += 1
                continue
            ref_bonds = {frozenset((b.GetBeginAtomIdx(), b.GetEndAtomIdx())) for b in ref.GetBonds()}
            same = mine_bonds == ref_bonds
            group = "organiche" if is_organic else "con metalli"
            bonds[group, same] += 1
            mine = {frozenset(r) for r in small_rings(len(m), m.bonds, lambda i, h=heavy: h[i], 8)}
            algo = rings_of(rdkit_mol(m, m.bonds), heavy)
            full = rings_of(ref, heavy)
            for name, ref_rings in (("algoritmo", algo), ("catena", full)):
                rings[group, name, "uguali"] += mine == ref_rings
                rings[group, name, "solo RDKit"] += len(ref_rings - mine)
                rings[group, name, "solo PISTAK"] += len(mine - ref_rings)
                if mine != ref_rings and is_organic:
                    extra = sorted(len(r) for r in mine - ref_rings)
                    lost = sorted(len(r) for r in ref_rings - mine)
                    bond_note = "" if same else (f"; legami solo PISTAK {len(mine_bonds - ref_bonds)}, "
                                                 f"solo RDKit {len(ref_bonds - mine_bonds)}")
                    examples.append(f"{f.stem} {key} [{name}]: anelli solo PISTAK {extra}, solo RDKit {lost}{bond_note}")
    print(f"  {mols} molecole distinte ({organic} organiche, solo non metalli), {failed} rifiutate da RDKit")
    for group, n_group in (("organiche", organic), ("con metalli", mols - failed - organic)):
        print(f"  {group} ({n_group}): legami uguali in {bonds[group, True]}, diversi in {bonds[group, False]}")
        for name, what in (("algoritmo", "sui legami di PISTAK"), ("catena", "sui legami di RDKit")):
            print(f"    anelli {what}: molecole uguali {rings[group, name, 'uguali']}; anelli solo RDKit "
                  f"{rings[group, name, 'solo RDKit']}, solo PISTAK {rings[group, name, 'solo PISTAK']}")
    print("  differenze nelle organiche:")
    for x in examples:
        print("    " + x)


def check_radii() -> None:
    import gemmi

    from pistak.core.elements import COVALENT_RADII, VDW_RADII

    print("\n== raggi: PISTAK contro le tabelle di gemmi (scarti oltre 0,02 Å) ==")
    for name, table, attr in (("covalenti", COVALENT_RADII, "covalent_r"), ("van der Waals", VDW_RADII, "vdw_r")):
        d = [(e, r, getattr(gemmi.Element(e), attr)) for e, r in table.items()]
        off = [(e, r, round(g, 3)) for e, r, g in d if abs(r - g) > 0.02]
        print(f"  {name}: {len(d)} elementi, {len(off)} diversi: {off}")


def check_hthresh() -> None:
    from pistak.core.validation import SP_ANGLE, TRIPLE_CX

    print(f"\n== soglie degli idrogeni: TRIPLE_CX = {TRIPLE_CX} Å, SP_ANGLE = {SP_ANGLE}° (numpy sulle molecole) ==")
    co, sp2, bare = [], [], []
    for f in FILES:
        try:
            a = analyse_text(f.read_text(encoding="utf-8", errors="replace"), f.stem, Limits())
        except Exception:  # noqa: BLE001  (CIF che PISTAK rifiuta: fuori dal confronto)
            continue
        for m in a.structure.molecules:
            nb = [[] for _ in m.el]
            for i, j in m.bonds:
                nb[i].append(j)
                nb[j].append(i)
            for k, e in enumerate(m.el):
                if e != "C":
                    continue
                heavy = [x for x in nb[k] if m.el[x] != "H"]
                n_h = len(nb[k]) - len(heavy)
                if len(heavy) == 3 and n_h == 0:  # carbonile: C a tre vicini con un O terminale
                    co += [float(np.linalg.norm(m.pos[x] - m.pos[k])) for x in heavy if m.el[x] == "O" and len(nb[x]) == 1]
                if len(heavy) == 2 and n_h <= 1:
                    u, v = m.pos[heavy[0]] - m.pos[k], m.pos[heavy[1]] - m.pos[k]
                    ang = math.degrees(math.acos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
                    (sp2 if n_h else bare).append(ang)
    co, sp2, bare = np.array(co), np.array(sp2), np.array(bare)
    print(f"  C=O: {len(co)} carbonili, da {co.min():.3f} a {co.max():.3f} Å (mediana {np.median(co):.3f}); "
          f"sotto TRIPLE_CX: {int((co <= TRIPLE_CX).sum())}")
    print(f"  C sp² con un H: {len(sp2)}, angolo massimo {sp2.max():.1f}°; sopra SP_ANGLE: {int((sp2 >= SP_ANGLE).sum())}")
    print(f"  C con due vicini e nessun H: {len(bare)}, mediana {np.median(bare):.1f}° (nitrili, alchini); "
          f"sotto SP_ANGLE: {int((bare < SP_ANGLE).sum())} (H mancanti o C legati a metalli)")


CHECKS = {"cell": check_cell, "moiety": check_moiety, "oblique": check_oblique,
          "origin": lambda: check_oblique(shift=np.array([0.137, 0.291, 0.443])), "rings": check_rings, "rdkit": check_rdkit, "radii": check_radii,
          "hthresh": check_hthresh}

if __name__ == "__main__":
    # "rdkit=500", "rings=500": the check on the first N structures of the sample (default 200, as in the first rounds)
    for arg in sys.argv[1:] or CHECKS:
        name, _, n = arg.partition("=")
        CHECKS[name](int(n)) if n else CHECKS[name]()
