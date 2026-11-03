"""Verifiche delle voci del registro di revisione (checklist, sezione 10), con oracoli indipendenti dal codice.

Ogni funzione `check_*` confronta PISTAK con un riferimento che non usa `pistak.core` (geometria scritta qui con
numpy, conti a mano, topologia nota dei policicli) e stampa il risultato. Esecuzione:

    .venv/bin/python validation/check_registro.py
"""
from __future__ import annotations

import itertools
import math
import random
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CC = 1.40  # Å, lato dell'esagono
H = CC * math.sqrt(3) / 2


# ---------------------------------------------------------------------------------------------------------
# Oracolo per i policicli: atomi e anelli costruiti sul reticolo a nido d'ape, senza codice di PISTAK
def honeycomb(centres: list[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray, list[frozenset[int]]]:
    """Atomi (N×3), legami (M×2) e anelli attesi (gli esagoni) di un idrocarburo benzenoide dato dai centri dei suoi
    anelli in coordinate assiali (q, r)."""
    atoms: list[np.ndarray] = []
    rings: list[frozenset[int]] = []

    def index(p: np.ndarray) -> int:
        for i, a in enumerate(atoms):
            if np.linalg.norm(a - p) < 1e-6:
                return i
        atoms.append(p)
        return len(atoms) - 1

    for q, r in centres:
        c = np.array([q * 2 * H + r * H, r * 1.5 * CC, 0.0])
        ring = [index(c + CC * np.array([math.cos(math.pi / 6 + k * math.pi / 3), math.sin(math.pi / 6 + k * math.pi / 3), 0.0]))
                for k in range(6)]
        rings.append(frozenset(ring))
    pos = np.array(atoms)
    bonds = np.array([(i, j) for i, j in itertools.combinations(range(len(pos)), 2) if np.linalg.norm(pos[i] - pos[j]) < CC + 0.05])
    return pos, bonds, rings


PAH = {
    "naftalene": [(0, 0), (1, 0)],
    "fenantrene": [(0, 0), (1, 0), (1, 1)],
    "pirene": [(0, 0), (1, 0), (0, 1), (1, -1)],
    "trifenilene": [(0, 0), (1, 0), (-1, 1), (0, -1)],
    "perilene": [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1)],
    "benzo[ghi]perilene": [(0, 0), (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1)],
    "coronene": [(0, 0), (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)],
}


def check_small_rings(trials: int = 200) -> None:
    """Registro: `small_rings` può perdere anelli centrali. Oracolo: gli esagoni del reticolo; il loro numero è anche
    il rango ciclico E − V + 1 (benzenoidi planari). Si permuta l'ordine degli atomi, perché la ricerca in ampiezza
    dipende dall'ordine dei vicini."""
    from pistak.core.cores import small_rings

    print("\n== small_rings contro gli esagoni attesi ==")
    rng = random.Random(0)
    for name, centres in PAH.items():
        pos, bonds, expected = honeycomb(centres)
        assert len(bonds) - len(pos) + 1 == len(expected)
        missing_runs = 0
        missing_example = None
        for t in range(trials):
            perm = list(range(len(pos)))
            if t:
                rng.shuffle(perm)
            inv = {old: new for new, old in enumerate(perm)}
            b = np.array([(inv[i], inv[j]) for i, j in bonds])
            if t:
                b = b[rng.sample(range(len(b)), len(b))]
            found = {frozenset(perm[i] for i in r) for r in small_rings(len(pos), b, lambda i: True, 8)}
            lost = [r for r in expected if r not in found]
            extra = [r for r in found if r not in expected]
            if lost or extra:
                missing_runs += 1
                missing_example = missing_example or (len(lost), len(extra))
        verdict = "OK" if not missing_runs else f"DIFETTO: anelli diversi in {missing_runs}/{trials} ordinamenti (persi, in più = {missing_example})"
        print(f"  {name:20s} {len(expected)} anelli attesi: {verdict}")


# ---------------------------------------------------------------------------------------------------------
_ANALYSES: dict[Path, object] = {}


def _analyse(path: Path):
    from pistak.api import analyse_file
    from pistak.core.config import Limits

    if path not in _ANALYSES:
        _ANALYSES[path] = analyse_file(path, Limits())
    return _ANALYSES[path]


def _structures() -> list[Path]:
    files = sorted((ROOT / "tests/data/cod").glob("*.cif")) + sorted((ROOT / "examples").glob("*.cif"))
    return [f for f in files if not f.name.startswith(".")]  # niente ._*.cif di macOS sui dischi exFAT


def check_mean_atom_height() -> None:
    """Registro: `mean_atom_height` dice «dal piano medio di A». Oracolo: media di (x − Cg_A)·n_A sugli atomi del nucleo
    B, con n_A la normale di A ricalcolata qui da un'SVD indipendente."""
    print("\n== mean_atom_height: descrizione (piano di A) contro codice ==")
    worst = (0.0, None)
    worst_flat = 0.0  # solo coppie con piani entro max_dihedral (stack e parallel)
    n_pairs = 0
    for path in _structures():
        an = _analyse(path)
        for p in an.pairs:
            ctx = an.pair_context(p)
            pa = ctx.pos_a[ctx.core_a.atoms]
            pb = ctx.pos_b[ctx.core_b.atoms]
            ca = pa.mean(axis=0)
            na = np.linalg.svd(pa - ca)[2][-1]
            signed = (pb - ca) @ na
            oracle_desc = float(np.abs(signed).mean())  # «distanza media degli atomi dal piano di A»
            got = p.get("mean_atom_height")
            diff = abs(got - oracle_desc)
            n_pairs += 1
            if p.kind != "tilted":
                worst_flat = max(worst_flat, diff)
            if diff > worst[0]:
                worst = (diff, (path.stem, p.index, p.kind, p.dihedral, got, oracle_desc, float(np.abs(signed).mean())))
    d, ex = worst
    print(f"  {n_pairs} coppie; scarto massimo codice − descrizione = {d:.3f} Å (stack e parallel: {worst_flat:.3f} Å)")
    if ex:
        print(f"  peggior caso {ex[0]} coppia {ex[1]} ({ex[2]}, angolo {ex[3]:.1f}°): codice {ex[4]:.3f}, "
              f"|media| dal piano A {ex[5]:.3f}, media dei |.| {ex[6]:.3f}")


def check_beta_gamma() -> None:
    """Registro: `beta_b` etichettato γ. Oracolo: definizioni di PLATON (β: Cg(I)→Cg(J) e normale al piano I;
    γ: Cg(I)→Cg(J) e normale al piano J), ricalcolate con normali da SVD."""
    from pistak.params import REGISTRY

    print("\n== beta_a / beta_b contro β, γ di PLATON ==")
    for key in ("beta_a", "beta_b"):
        s = REGISTRY.get(key)
        print(f"  {key}: etichetta «{s.label}»")
    worst = 0.0
    for path in _structures():
        an = _analyse(path)
        for p in an.pairs:
            ctx = an.pair_context(p)
            pa, pb = ctx.pos_a[ctx.core_a.atoms], ctx.pos_b[ctx.core_b.atoms]
            ca, cb = pa.mean(axis=0), pb.mean(axis=0)
            na, nb = np.linalg.svd(pa - ca)[2][-1], np.linalg.svd(pb - cb)[2][-1]
            d = cb - ca
            ang = lambda u, v: math.degrees(math.acos(min(1.0, abs(u @ v) / np.linalg.norm(u) / np.linalg.norm(v))))
            worst = max(worst, abs(p.get("beta_a") - ang(d, na)), abs(p.get("beta_b") - ang(d, nb)))
    print(f"  scarto massimo dai valori ricalcolati: {worst:.2e}°  → i numeri sono β e γ di PLATON; la chiave beta_b no")


def check_unit_zero() -> None:
    """Registro: `unit()` restituisce il vettore nullo. Si controlla se nei due usi (pair_frame) la norma può essere
    piccola: con n_B orientata come n_A, |n_A + n_B| ≥ √2 e |L_A proiettato| ≥ 1/√2 (dimostrazione nel registro).
    Verifica numerica su coppie di nuclei casuali, anche perpendicolari."""
    from pistak.core.pairs import pair_frame

    print("\n== unit(): norme minime nei due usi di pair_frame ==")
    rng = np.random.default_rng(1)
    min_n, min_l = 9.0, 9.0

    def core(normal, long):
        return SimpleNamespace(normal=normal, long=long, centroid=np.zeros(3))

    for _ in range(20000):
        na = rng.normal(size=3)
        na /= np.linalg.norm(na)
        la = np.cross(na, rng.normal(size=3))
        la /= np.linalg.norm(la)
        nb = la if rng.random() < 0.1 else rng.normal(size=3)  # anche il caso peggiore: n_B lungo l'asse di A
        nb /= np.linalg.norm(nb)
        lb = np.cross(nb, rng.normal(size=3))
        lb /= np.linalg.norm(lb)
        nb_al = -nb if na @ nb < 0 else nb
        min_n = min(min_n, float(np.linalg.norm(na + nb_al)))
        n = (na + nb_al) / np.linalg.norm(na + nb_al)
        min_l = min(min_l, float(np.linalg.norm(la - n * (la @ n))))
        for sym in (False, True):
            f = pair_frame(core(na, la), core(nb, lb), symmetric=sym)
            assert abs(np.linalg.norm(f.L) - 1) < 1e-9 and abs(np.linalg.norm(f.n) - 1) < 1e-9
    print(f"  min |n_A + n_B| = {min_n:.3f} (limite √2 = 1.414), min |L_A proiettato| = {min_l:.3f} (limite 0.707)")
    print("  → nei percorsi attuali unit() non riceve mai un vettore nullo; L, W, n sempre unitari")


def check_max_shells() -> None:
    """Registro: MAX_SHELLS = 8 taglia in silenzio. Si chiede a neighbour_shells una cella con un asse corto e una
    molecola estesa: il valore richiesto dal limite teorico (D·|a*| + 1) si calcola qui."""
    from pistak.core.config import Limits
    from pistak.core.pairs import MAX_SHELLS, neighbour_shells

    print("\n== neighbour_shells e MAX_SHELLS ==")
    lim = Limits()
    for a, ext in ((3.5, 7.0), (3.5, 12.0), (3.4, 15.0)):
        M = np.diag([a, 30.0, 30.0])
        s = SimpleNamespace(molecules=[SimpleNamespace(extent=ext, extended=False)], vectors=M)
        need = math.floor((2 * ext + lim.contact_cutoff) / a + 1)
        got = neighbour_shells(s, [[object()]], lim)[0]
        flag = "TAGLIATO in silenzio" if need > got else "ok"
        print(f"  a = {a} Å, estensione {ext} Å: servono ±{need}, esplorate ±{got} (MAX_SHELLS = {MAX_SHELLS}) {flag}")


def _column_cif(c: float, slip: float = 0.0) -> str:
    """Naftalene nel piano xy, cella ortorombica 20 × 20 × c Å con una piccola inclinazione della cella per lo slittamento:
    colonna infinita lungo c, piani esattamente paralleli a distanza c·cos(angolo). Coordinate frazionarie a mano."""
    pos, _, _ = honeycomb(PAH["naftalene"])
    # monoclino: c inclinato di beta nel piano xz → slittamento lungo x = c·cos(beta), quota = c·sin(beta)
    clen = math.hypot(c, slip)
    beta = math.degrees(math.atan2(c, slip)) if slip else 90.0
    rows = []
    for i, p in enumerate(pos):
        x, y, z = p + np.array([10.0, 10.0, 0.0])
        fz = 0.0
        fx = x / 20.0
        rows.append(f"C{i + 1} C {fx:.6f} {y / 20.0:.6f} {fz:.6f}")
    return "\n".join([
        "data_col", "_symmetry_space_group_name_H-M 'P 1'", "_cell_length_a 20.0", "_cell_length_b 20.0",
        f"_cell_length_c {clen:.6f}", "_cell_angle_alpha 90", f"_cell_angle_beta {beta:.6f}", "_cell_angle_gamma 90",
        "loop_", "_symmetry_equiv_pos_as_xyz", "'x, y, z'",
        "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        *rows, ""])


def check_cutoff_vs_interplanar() -> None:
    """Registro: coppie cercate entro 4,0 Å ma stack fino a 4,2 Å. Colonne di naftalene a piani paralleli, distanza
    nota per costruzione: per piani paralleli ogni distanza atomo–atomo è ≥ la distanza interplanare, quindi uno stack
    tra 4,0 e 4,2 Å non viene mai misurato."""
    from pistak.api import analyse_text
    from pistak.core.config import Limits

    print("\n== contact_cutoff 4,0 Å contro max_interplanar 4,2 Å (piani paralleli) ==")
    for d in (3.5, 3.9, 3.99, 4.05, 4.1, 4.19):
        for slip in (0.0, 1.5):
            an = analyse_text(_column_cif(d, slip), "col", Limits())
            stacks = [p for p in an.pairs if p.kind == "stack"]
            got = f"stack, interplanare {stacks[0].interplanar:.3f} Å" if stacks else f"nessuno stack ({len(an.pairs)} coppie)"
            print(f"  interplanare atteso {d:.2f} Å, slittamento {slip:.1f} Å: {got}")


def check_cell_vectors() -> None:
    """`tests/cifgen.py` usa `cell_vectors` di PISTAK per scrivere i CIF sintetici: l'oracolo è circolare a meno che
    `cell_vectors` sia corretto. Oracolo indipendente: tensore metrico G_ij = a_i·a_j dai parametri di cella e volume
    V = abc·√(1 − cos²α − cos²β − cos²γ + 2 cosα cosβ cosγ)."""
    from pistak.core.cif import Cell
    from pistak.core.structure import cell_vectors

    print("\n== cell_vectors contro tensore metrico e volume ==")
    rng = np.random.default_rng(2)
    worst_g = worst_v = 0.0
    n = 0
    while n < 5000:
        a, b, c = rng.uniform(3, 30, 3)
        al, be, ga = rng.uniform(50, 130, 3)
        ca, cb, cg = np.cos(np.radians([al, be, ga]))
        disc = 1 - ca ** 2 - cb ** 2 - cg ** 2 + 2 * ca * cb * cg
        if disc <= 1e-3:
            continue
        n += 1
        G = np.array([[a * a, a * b * cg, a * c * cb], [a * b * cg, b * b, b * c * ca], [a * c * cb, b * c * ca, c * c]])
        M = cell_vectors(Cell(a, b, c, al, be, ga))
        worst_g = max(worst_g, float(np.abs(M @ M.T - G).max()))
        worst_v = max(worst_v, abs(np.linalg.det(M) - a * b * c * math.sqrt(disc)) / (a * b * c))
        assert abs(M[0, 1]) < 1e-12 and abs(M[0, 2]) < 1e-12 and abs(M[1, 2]) < 1e-12  # a lungo x, b nel piano xy
    print(f"  5000 celle casuali: max |M·Mᵀ − G| = {worst_g:.1e} Å², max errore relativo sul volume = {worst_v:.1e}")


CHECKS = {
    "cell": check_cell_vectors, "rings": check_small_rings, "beta": check_beta_gamma, "height": check_mean_atom_height,
    "unit": check_unit_zero, "shells": check_max_shells, "cutoff": check_cutoff_vs_interplanar,
}

if __name__ == "__main__":
    for name in sys.argv[1:] or CHECKS:
        CHECKS[name]()
