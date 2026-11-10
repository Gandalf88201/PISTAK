"""Checklist, sezione 6: incertezze (Monte Carlo) e geometria intramolecolare.

    mc         σ Monte Carlo contro la propagazione lineare (oracolo in numpy) su tre cristalli sintetici
    converge   σ in funzione del numero di campioni N su tre coppie reali, e scarto media − valore nominale
    axis       σ degli slittamenti L/W per un nucleo isotropo (benzene) e uno allungato (naftalene)
    threshold  kind_agreement di una coppia costruita sulla soglia dell'angolo tra i piani (atteso ≈ 0,5)
    geom       legami, angoli, torsioni e codici di simmetria contro le tabelle _geom_ dei CIF COD (letti con gemmi)

L'oracolo non usa `pistak.core`: i cristalli sintetici sono scritti qui come testo CIF, la geometria della coppia
(piani ai minimi quadrati, normale media, slittamento) e la propagazione lineare sono ricalcolate con numpy, le tabelle
_geom_ si leggono con gemmi.

    .venv/bin/python validation/check_incertezze.py [mc converge axis threshold geom]
"""
from __future__ import annotations

import math
import re
import sys
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from pistak.core.analysis import analyse_text  # noqa: E402
from pistak.core.config import Limits  # noqa: E402
from pistak.core.uncertainty import monte_carlo  # noqa: E402

COD = sorted((HERE / "cod_cache").glob("[0-9]*.cif"))
KEYS = ("interplanar", "centroid_distance", "dihedral", "slip")


# ---- cristalli sintetici scritti a mano ---------------------------------------------------------------------------
def ring(n_atoms: int = 6, r: float = 1.39) -> np.ndarray:
    t = np.arange(n_atoms) * 2 * np.pi / n_atoms
    return np.c_[r * np.cos(t), r * np.sin(t), np.zeros(n_atoms)]


NAPHTHALENE_HALF = np.array([[0.0, 0.71, 0], [1.23, 1.40, 0], [2.44, 0.70, 0], [2.44, -0.70, 0], [1.23, -1.40, 0]])
NAPHTHALENE = np.r_[NAPHTHALENE_HALF, -NAPHTHALENE_HALF]


def rot(axis: str, deg: float) -> np.ndarray:
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    i, j = {"x": (1, 2), "y": (2, 0), "z": (0, 1)}[axis]
    m = np.eye(3)
    m[i, i] = m[j, j] = c
    m[i, j], m[j, i] = -s, s
    return m


def cif_text(cell: tuple[float, float, float], ops: list[str], frac: np.ndarray, su: np.ndarray) -> str:
    """CIF minimo, cella ortorombica; coordinate con s.u. nella notazione x(u) sull'ultima cifra (6 decimali)."""
    lines = ["data_sintetico", f"_cell_length_a {cell[0]}", f"_cell_length_b {cell[1]}", f"_cell_length_c {cell[2]}",
             "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90",
             "loop_", "_symmetry_equiv_pos_as_xyz", *(f"'{o}'" for o in ops),
             "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
             "_atom_site_fract_z"]
    for k, (f, u) in enumerate(zip(frac, su)):
        coords = " ".join(f"{x:.6f}({round(e * 1e6)})" for x, e in zip(f, u))
        lines.append(f"C{k + 1} C {coords}")
    return "\n".join(lines) + "\n"


# ---- oracolo: geometria della coppia e propagazione lineare -------------------------------------------------------
def plane(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    c = p.mean(axis=0)
    return c, np.linalg.svd(p - c)[2][-1]


def pair_geometry(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    ca, na = plane(a)
    cb, nb = plane(b)
    if na @ nb < 0:
        nb = -nb
    n = (na + nb) / np.linalg.norm(na + nb)
    d = cb - ca
    z = d @ n
    return {"interplanar": abs(z), "centroid_distance": float(np.linalg.norm(d)),
            "dihedral": math.degrees(math.acos(min(1.0, abs(na @ nb)))), "slip": float(np.linalg.norm(d - z * n))}


def linear_sigma(f, x0: np.ndarray, su: np.ndarray, step: float = 1e-7) -> dict[str, float]:
    """σ² = Σ (∂f/∂x_k)² σ_k², coordinate indipendenti (la stessa ipotesi del Monte Carlo di PISTAK)."""
    var = dict.fromkeys(f(x0), 0.0)
    for idx in np.ndindex(x0.shape):
        if su[idx] == 0:
            continue
        xp, xm = x0.copy(), x0.copy()
        xp[idx] += step
        xm[idx] -= step
        fp, fm = f(xp), f(xm)
        for k in var:
            var[k] += ((fp[k] - fm[k]) / (2 * step) * su[idx]) ** 2
    return {k: math.sqrt(v) for k, v in var.items()}


def case_general_p1(rng):
    """P1: due benzeni indipendenti, B inclinato di 6° e spostato; cella larga, una coppia sola."""
    cell = np.array([25.0, 25.0, 25.0])
    a = ring() + [7.5, 7.5, 7.5]
    b = ring() @ rot("x", 6).T + [7.5 + 1.2, 7.5 + 0.5, 7.5 + 3.5]
    frac = np.r_[a, b] / cell
    su = rng.uniform(1e-4, 3e-4, frac.shape)
    return ("P1, due benzeni indipendenti", cif_text(tuple(cell), ["x,y,z"], frac, su), frac, su,
            lambda x: pair_geometry(x[:6] * cell, x[6:] * cell))


def case_inversion(rng):
    """P-1: un benzene in posizione generale e la sua immagine per inversione (stesse s.u., spostamenti opposti)."""
    cell = np.array([12.0, 14.0, 16.0])
    a = ring() @ rot("y", 25).T @ rot("x", 10).T + [6.0 + 1.0, 7.0 + 0.6, 8.0 + 1.75]
    frac = a / cell
    su = rng.uniform(1e-4, 3e-4, frac.shape)

    def f(x):
        p = x * cell
        return pair_geometry(p, (1.0 - x) * cell)  # -x+1,-y+1,-z+1: inversione nel centro (½, ½, ½)

    return ("P-1, benzene e immagine per inversione", cif_text(tuple(cell), ["x,y,z", "-x,-y,-z"], frac, su), frac, su, f)


def case_half_naphthalene(rng):
    """P-1, Z' = ½: naftalene sul centro di inversione nell'origine, coppia con la copia traslata di a."""
    cell = np.array([4.0, 12.0, 14.0])
    orient = rot("z", 30) @ rot("y", 90)  # normale del naftalene (asse z) portata a 30° da a nel piano ab
    half = NAPHTHALENE_HALF @ orient.T
    frac = half / cell
    su = rng.uniform(1e-4, 3e-4, frac.shape)

    def f(x):
        mol = np.r_[x, -x] * cell
        return pair_geometry(mol, mol + [cell[0], 0, 0])

    return ("P-1, naftalene sul centro di inversione (Z' = ½)", cif_text(tuple(cell), ["x,y,z", "-x,-y,-z"], frac, su),
            frac, su, f)


def find_pair(an, nominal: dict[str, float]):
    best = min(an.pairs, key=lambda p: abs(p.centroid_distance - nominal["centroid_distance"]))
    assert abs(best.centroid_distance - nominal["centroid_distance"]) < 1e-4, "coppia non trovata"
    return best


def check_mc(samples: int = 4000) -> None:
    print(f"== σ Monte Carlo (N = {samples}) contro propagazione lineare indipendente ==")
    rng = np.random.default_rng(7)
    tol = 3 / math.sqrt(2 * (samples - 1))  # 3 volte l'errore statistico relativo di una deviazione standard
    worst = 0.0
    for case in (case_general_p1, case_inversion, case_half_naphthalene):
        name, text, frac, su, f = case(rng)
        an = analyse_text(text, "sintetico", Limits())
        nominal = f(frac)
        p = find_pair(an, nominal)
        lin = linear_sigma(f, frac, su)
        u = monte_carlo(an, samples=samples, seed=1, include_cell=False)
        print(f"  {name}: coppia {p.kind}, {len(an.pairs)} coppie nella struttura")
        for k in KEYS:
            mc, nom = u.get(p, k), getattr(p, k)
            if k == "dihedral" and nominal[k] < 1e-4:
                # piani paralleli per simmetria (inversione, traslazione) in ogni campione: |θ| ha una cuspide in 0 e
                # la derivata numerica non ha senso; il Monte Carlo deve dare σ ≈ 0
                ok = mc < 1e-4
                worst = max(worst, 0.0 if ok else 1.0)
                print(f"    {k:18s} nominale {nom:8.4f}: piani paralleli per simmetria, σ MC {mc:.1e} "
                      f"({'≈ 0, OK' if ok else 'NON NULLA'})")
                continue
            if lin[k] < 1e-9:
                print(f"    {k:18s} nominale {nom:8.4f} (oracolo {nominal[k]:8.4f})  σ lineare 0, σ MC {mc:.2e}")
                continue
            rel = mc / lin[k] - 1
            worst = max(worst, abs(rel))
            print(f"    {k:18s} nominale {nom:8.4f} (oracolo {nominal[k]:8.4f})  σ lineare {lin[k]:.5f}  "
                  f"σ MC {mc:.5f}  scarto {100 * rel:+.1f} %")
    print(f"  scarto relativo massimo {100 * worst:.1f} % (errore statistico atteso ±{100 * tol / 3:.1f} %, soglia "
          f"±{100 * tol:.1f} %): {'OK' if worst < tol else 'FUORI'}")


# ---- convergenza e bias su coppie reali ---------------------------------------------------------------------------
#: Coppie della convergenza (revisione 3, N8): (CIF, etichette del frammento A, del frammento B, Cg···Cg nominale in Å),
#: le stesse di `paper/scripts/mc_data.py` e `theta_kurtosis.py`. 2014747: 4,5-diazafluoren-9-one coordinato al Cu
#: sopra la stessa molecola di solvente (θ = 12,2°); 2210998: due nuclei C26H16FN3O, θ = 0 per simmetria; 4032966: due
#: benzeni di due pillar[5]areni, θ = 0 per simmetria. Sostituisce la coppia acido glutarico/piridina di 7212873, che
#: dopo la voce 70 (ripiego solo sp²/sp) non è più una coppia π.
_FLUORENONE_A = ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "N1", "N2")
_FLUORENONE_B = ("C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "N3", "N4")
_PYRIDOQ = ("C10A", "C3", "C3A", "C4", "C4A", "C5", "C5A", "C6", "C7", "C8", "C9", "C9A", "C9B", "N1", "N10", "N2")
_BENZENE = ("C29", "C30", "C31", "C32", "C33", "C34")
CONVERGE = [("2014747", _FLUORENONE_A, _FLUORENONE_B, 3.742),
            ("2210998", _PYRIDOQ, _PYRIDOQ, 5.455),
            ("4032966", _BENZENE, _BENZENE, 3.852)]


def pick_pair(an, frag_a, frag_b, cg: float, tol: float = 0.005):
    """La coppia π-stack intermolecolare i cui frammenti hanno queste etichette del CIF (in un ordine o nell'altro) e
    il cui Cg···Cg dista meno di `tol` da `cg`. Deve essercene esattamente una: altrimenti è un errore, non un ripiego
    silenzioso (la stessa regola di `paper/scripts/mc_data.pick_pair`)."""
    want = sorted([tuple(sorted(frag_a)), tuple(sorted(frag_b))])
    hits = []
    s = an.structure
    for p in an.pairs:
        if p.intramolecular:
            continue
        la = tuple(sorted(s.molecules[p.a].labels[i] for i in an.fragments[p.a][p.fa].atoms))
        lb = tuple(sorted(s.molecules[p.b].labels[i] for i in an.fragments[p.b][p.fb].atoms))
        if sorted([la, lb]) == want and abs(p.centroid_distance - cg) < tol:
            hits.append(p)
    if len(hits) != 1:
        raise RuntimeError(f"{len(hits)} coppie corrispondono a {frag_a} / {frag_b} con Cg···Cg {cg} Å")
    if hits[0].kind != "stack":
        raise RuntimeError(f"la coppia {hits[0].index} è {hits[0].kind}, non un π-stack")
    return hits[0]


def check_converge(ns=(50, 100, 200, 500, 1000), seeds: int = 5) -> None:
    print("\n== σ in funzione di N (media e dispersione relativa su", seeds, "semi) ==")
    for stem, frag_a, frag_b, cg in CONVERGE:
        text = (HERE / "cod_cache" / f"{stem}.cif").read_text(encoding="utf-8", errors="replace")
        an = analyse_text(text, stem, Limits())
        p = pick_pair(an, frag_a, frag_b, cg)
        t0 = time.time()
        print(f"  {stem}, coppia {p.index} ({p.kind}, interplanare {p.interplanar:.3f} Å, angolo {p.dihedral:.2f}°), "
              f"{len(an.pairs)} coppie nella struttura")
        for n in ns:
            sig = {k: [] for k in KEYS}
            bias = {k: [] for k in KEYS}
            for seed in range(seeds):
                u = monte_carlo(an, samples=n, seed=seed)
                for k in KEYS:
                    sig[k].append(u.get(p, k))
                    bias[k].append(u.mean[p.index][k] - getattr(p, k))
            row = "  ".join(f"{k} {np.mean(v):.4f} ±{100 * np.std(v, ddof=1) / np.mean(v):4.1f}%" for k, v in sig.items())
            print(f"    N = {n:4d}: {row}")
        b = "  ".join(f"{k} {np.mean(v) / np.mean(sig[k]):+.2f}σ" for k, v in bias.items())
        print(f"    media − nominale a N = {ns[-1]}, in unità di σ: {b}   ({time.time() - t0:.0f} s)")
        print("    dispersione relativa attesa di una σ stimata con N campioni: 1/√(2(N−1)) = "
              + ", ".join(f"{100 / math.sqrt(2 * (n - 1)):.1f}% (N={n})" for n in ns))


# ---- asse lungo di un nucleo isotropo ------------------------------------------------------------------------------
def check_axis(samples: int = 1000) -> None:
    print("\n== σ degli slittamenti lungo L e W: benzene (asse lungo non definito) contro naftalene ==")
    cell = np.array([25.0, 25.0, 25.0])
    for name, mol in (("benzene", ring()), ("naftalene", NAPHTHALENE)):
        a = mol + [7.5, 7.5, 7.5]
        b = mol + [7.5 + 1.5, 7.5 + 0.6, 7.5 + 3.5]
        frac = np.r_[a, b] / cell
        su = np.full(frac.shape, 2e-4)
        an = analyse_text(cif_text(tuple(cell), ["x,y,z"], frac, su), name, Limits())
        p = an.pairs[0]
        u = monte_carlo(an, samples=samples, seed=2, include_cell=False)
        print(f"  {name:9s}: slittamento {p.slip:.3f} ± {u.get(p, 'slip'):.3f} Å; L {p.long_slip:.3f} ± "
              f"{u.get(p, 'long_slip'):.3f} Å; W {p.lat_slip:.3f} ± {u.get(p, 'lat_slip'):.3f} Å")


# ---- coppia sulla soglia -------------------------------------------------------------------------------------------
def check_threshold(samples: int = 2000) -> None:
    lim = Limits()
    print(f"\n== kind_agreement di una coppia con angolo tra i piani = max_dihedral ({lim.max_dihedral}°) ==")
    rng = np.random.default_rng(5)
    cell = np.array([25.0, 25.0, 25.0])
    a = ring() + [7.5, 7.5, 7.5]
    b = ring() @ rot("x", lim.max_dihedral).T + [7.5 + 0.4, 7.5, 7.5 + 3.5]
    frac = np.round(np.r_[a, b] / cell, 6)
    su = rng.uniform(1e-4, 3e-4, frac.shape)
    an = analyse_text(cif_text(tuple(cell), ["x,y,z"], frac, su), "soglia", lim)
    p = an.pairs[0]
    u = monte_carlo(an, samples=samples, seed=4, include_cell=False)
    k = u.kind_agreement[p.index]
    err = math.sqrt(0.25 / samples)
    print(f"  angolo nominale {p.dihedral:.4f}° ({p.kind}), σ {u.get(p, 'dihedral'):.3f}°: kind_agreement {k:.3f} "
          f"(atteso 0,5 ± {err:.3f}, se la distribuzione dell'angolo è simmetrica attorno al nominale)")


# ---- tabelle _geom_ dei CIF ----------------------------------------------------------------------------------------
_TERM = re.compile(r"([+-]?)(\d*/?\d*)\*?([xyz]?)")


def parse_op(text: str) -> tuple[np.ndarray, np.ndarray]:
    """'-x+1/2,y,z' → (R, t), parser indipendente da quello di PISTAK."""
    r, t = np.zeros((3, 3)), np.zeros(3)
    for i, comp in enumerate(text.replace(" ", "").lower().split(",")):
        for sign, num, var in _TERM.findall(comp):
            if not num and not var:
                continue
            v = float(Fraction(num)) if num else 1.0
            v = -v if sign == "-" else v
            if var:
                r[i, "xyz".index(var)] += v
            else:
                t[i] += v
    return r, t


def compose(a, b):
    return a[0] @ b[0], a[0] @ b[1] + a[1]


def inverse(a):
    ri = np.linalg.inv(a[0])
    return ri, -ri @ a[1]


IDENT = (np.eye(3), np.zeros(3))


def cif_ops(block) -> list[tuple[np.ndarray, np.ndarray]]:
    for tag in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"):
        col = block.find_values(tag)
        if col:
            return [parse_op(gemmi_str(x)) for x in col]
    return []


def gemmi_str(x: str) -> str:
    import gemmi

    return gemmi.cif.as_string(x)


def site_op(code: str, ops) -> tuple[np.ndarray, np.ndarray] | None:
    code = gemmi_str(code)  # «.» (nessuna operazione) diventa ""
    if code in ("", ".", "?", "1_555", "1"):
        return IDENT if code != "?" else None
    m = re.fullmatch(r"(\d+)(?:_(\d)(\d)(\d))?", code)
    if not m or int(m.group(1)) > len(ops):
        return None
    r, t = ops[int(m.group(1)) - 1]
    shift = np.array([int(m.group(k)) - 5 for k in (2, 3, 4)]) if m.group(2) else np.zeros(3)
    return r, t + shift


def same_op(a, b) -> bool:
    return bool(np.allclose(a[0], b[0]) and np.allclose(a[1], b[1], atol=1e-6))


def value_su(text: str) -> tuple[float, float, int]:
    text = gemmi_str(text)
    m = re.fullmatch(r"(-?\d*\.?(\d*))(?:\((\d+)\))?", text)
    if not m or not re.search(r"\d", m.group(1)):
        return math.nan, math.nan, 0
    dec = len(m.group(2))
    su = int(m.group(3)) * 10.0 ** -dec if m.group(3) else math.nan
    return float(m.group(1)), su, dec


def measure(kind: str, p: np.ndarray) -> float:
    """Distanza, angolo, torsione (IUPAC: positiva se, guardando lungo 2→3, il legame 1–2 ruota in senso orario
    per eclissare 3–4)."""
    if kind == "bond":
        return float(np.linalg.norm(p[1] - p[0]))
    if kind == "angle":
        u, v = p[0] - p[1], p[2] - p[1]
        return math.degrees(math.acos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
    b1, b2, b3 = p[1] - p[0], p[2] - p[1], p[3] - p[2]
    return math.degrees(math.atan2(np.linalg.norm(b2) * b1 @ np.cross(b2, b3), np.cross(b1, b2) @ np.cross(b2, b3)))


def vals_text(x: str) -> str:
    return gemmi_str(x)


def check_geom() -> None:
    import gemmi

    from pistak.core.intramolecular import molecule_geometry

    print("\n== legami, angoli, torsioni contro le tabelle _geom_ dei CIF COD ==")
    stats = {k: {"d": [], "half": 0, "n": 0, "ratio": [], "worst": []} for k in ("bond", "angle", "torsion")}
    recomputed, recomputed_bad = 0, []
    sign_bad, code_ok, code_bad, n_cif, examples = 0, 0, 0, 0, []
    n_lab = {"bond": 2, "angle": 3, "torsion": 4}
    value_tag = {"bond": "_geom_bond_distance", "angle": "_geom_angle", "torsion": "_geom_torsion"}
    for path in COD:
        block = gemmi.cif.read_file(str(path)).sole_block()
        if not block.find_values("_geom_bond_distance"):
            continue
        try:
            an = analyse_text(path.read_text(encoding="utf-8", errors="replace"), path.stem, Limits())
        except Exception:  # noqa: BLE001  (CIF che PISTAK non analizza: contati nel registro, sezione 1)
            continue
        n_cif += 1
        ops = cif_ops(block)
        mine: dict[tuple, list] = {}
        geom = molecule_geometry(an.structure, an.cif, kinds=("bond", "angle", "torsion"))
        # codici di simmetria: la misura rifatta dalle coordinate del CIF con le operazioni scritte da PISTAK
        cell = gemmi.UnitCell(*(block.find_value(f"_cell_{t}") and gemmi.cif.as_number(block.find_value(f"_cell_{t}"))
                                for t in ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")))
        frac = {gemmi.cif.as_string(r[0]): np.array([gemmi.cif.as_number(r[i]) for i in (1, 2, 3)])
                for r in block.find("_atom_site_", ["label", "fract_x", "fract_y", "fract_z"])}
        for g in geom:
            if not any(g.codes) or not all(lab in frac for lab in g.atoms):
                continue
            pts = []
            for lab, code in zip(g.atoms, g.codes):
                r, t = parse_op(code) if code else IDENT
                f = r @ frac[lab] + t
                pts.append(np.array(cell.orthogonalize(gemmi.Fractional(*f)).tolist()))
            v = measure(g.kind, np.array(pts))
            recomputed += 1
            dv = (v - g.value + 180) % 360 - 180 if g.kind == "torsion" else v - g.value
            if abs(dv) > 1e-5:  # 1e-5: acos vicino a 180° perde ~1e-6°
                recomputed_bad.append((path.stem, g.kind, g.atoms, g.codes, g.value, v))
        for g in geom:
            mine.setdefault((g.kind, g.atoms), []).append(g)
            if g.kind != "torsion":
                mine.setdefault((g.kind, g.atoms[::-1]), []).append(g)
        for kind, n in n_lab.items():
            pre = value_tag[kind] if kind != "bond" else "_geom_bond"
            cols = [f"{pre}_atom_site_label_{i + 1}" for i in range(n)] + [value_tag[kind]]
            cols += [f"?{pre}_site_symmetry_{i + 1}" for i in range(n)]
            table = block.find("", cols)
            for row in table:
                key = tuple(row.str(i) for i in range(n))
                value, su, dec = value_su(row[n])
                if not math.isfinite(value):
                    continue
                syms = [row[n + 1 + i] if row.has(n + 1 + i) else "." for i in range(n)]
                codes = [site_op(x, ops) for x in syms]
                if any(c is None for c in codes):
                    continue
                cands = mine.get((kind, key), [])
                if not cands:
                    continue
                # tra le misure di PISTAK con queste etichette, quella con la stessa relazione di simmetria tra gli atomi
                want = [compose(inverse(codes[0]), c) for c in codes]
                pick, improper = None, False
                for g in cands:
                    gc = [parse_op(c) if c else IDENT for c in (g.codes or [""] * len(g.atoms))]
                    if kind != "torsion" and g.atoms != key:
                        gc = gc[::-1]
                    rel = [compose(inverse(gc[0]), c) for c in gc]
                    if all(same_op(a, b) for a, b in zip(rel, want)):
                        pick = g
                        # le due copie differiscono per G = O_PISTAK ∘ O_CIF⁻¹: se G è impropria la torsione cambia segno
                        improper = bool(np.linalg.det(compose(gc[0], inverse(codes[0]))[0]) < 0)
                        break
                has_sym = any(not same_op(c, IDENT) for c in want)
                if pick is None:
                    if has_sym:
                        tol = 3 * su if math.isfinite(su) and su > 0 else 0.5 * 10.0 ** -dec + 0.01
                        # torsioni: la copia equivalente per un'operazione impropria ha il segno opposto
                        if any(abs(abs(g.value) - abs(value)) <= tol if kind == "torsion" else abs(g.value - value) <= tol
                               for g in cands):
                            code_bad += 1  # stessa misura su un'altra copia equivalente
                        else:
                            examples.append((path.stem, kind, key, [gemmi_str(x) for x in syms]))
                    continue
                if has_sym:
                    code_ok += 1
                st = stats[kind]
                d = pick.value - value
                if kind == "torsion":
                    if improper:
                        value = -value
                        d = pick.value - value
                    d = (d + 180) % 360 - 180
                    if abs(value) > 5 and abs(abs(value) - 180) > 5 and np.sign(pick.value) != np.sign(value):
                        sign_bad += 1
                if not math.isfinite(su):
                    continue  # misure con idrogeni in posizione calcolata: senza s.u., coordinate scritte con meno cifre
                st["n"] += 1
                st["d"].append(abs(d))
                st["half"] += abs(d) <= 0.5 * 10.0 ** -dec + 1e-9
                st["worst"].append((abs(d) / su if su > 0 else 0.0, path.stem, key, vals_text(row[n]), pick.value))
                if math.isfinite(pick.su) and su > 0:
                    st["ratio"].append(pick.su / su)
    print(f"  {n_cif} CIF con le tabelle _geom_ (solo misure con s.u. nel CIF):")
    lim = {"bond": 0.001, "angle": 0.1, "torsion": 0.1}
    for kind, st in stats.items():
        if not st["n"]:
            continue
        d = np.array(st["d"])
        r = np.array(st["ratio"])
        print(f"    {kind:8s} {st['n']:6d} misure: |Δ| mediano {np.median(d):.4f}, 99° perc. {np.percentile(d, 99):.4f}, "
              f"max {d.max():.4f}; entro mezza unità dell'ultima cifra {100 * st['half'] / st['n']:.1f} %, entro "
              f"{lim[kind]} {100 * np.mean(d <= lim[kind] + 1e-9):.1f} %")
        z = np.array([w[0] for w in st["worst"]])
        print(f"             |Δ|/s.u.(CIF): oltre 1 {100 * np.mean(z > 1):.1f} %, oltre 3 {100 * np.mean(z > 3):.2f} % "
              f"({int(np.sum(z > 3))} misure)")
        for w in sorted(st["worst"], reverse=True)[:3]:
            print(f"             {w[0]:5.1f} s.u.: {w[1]} {'–'.join(w[2])} CIF {w[3]} PISTAK {w[4]:.4f}")
        print(f"             s.u. PISTAK / s.u. CIF: mediana {np.median(r):.2f}, 5–95 % {np.percentile(r, 5):.2f}–"
              f"{np.percentile(r, 95):.2f}")
    print(f"  codici di simmetria: {recomputed} misure con atomi trasformati rifatte dalle coordinate del CIF con le "
          f"operazioni scritte da PISTAK, {len(recomputed_bad)} diverse")
    for x in recomputed_bad[:5]:
        print(f"    diversa: {x}")
    print(f"  misure del CIF con atomi trasformati: {code_ok} abbinate con la stessa relazione di simmetria, {code_bad} "
          f"con lo stesso valore su un'altra copia equivalente (PISTAK sceglie quella con meno atomi trasformati), "
          f"{len(examples)} senza corrispondente")
    print("    (esaminate a mano il 5/11/2026: tutte contengono un legame a un atomo trasformato che il programma di "
          "affinamento considera tale e PISTAK no, come Cu–Cl a 3,0 Å, Na–O, Ag–O, o un atomo di disordine)")
    by = Counter((x[0], x[1]) for x in examples)
    print("    senza corrispondente, per CIF e tipo: " + ", ".join(f"{k[0]} {k[1]} {v}" for k, v in by.most_common(10)))
    print(f"  torsioni con segno opposto (|τ| fra 5° e 175°): {sign_bad}")


CHECKS = {"mc": check_mc, "converge": check_converge, "axis": check_axis, "threshold": check_threshold, "geom": check_geom}

if __name__ == "__main__":
    for name in sys.argv[1:] or CHECKS:
        CHECKS[name]()
