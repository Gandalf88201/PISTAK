"""Checklist, sezione 7: i file per i calcoli quantistici e periodici rimessi in coordinate e confrontati con il CIF.

    periodic  POSCAR, pw.x, CP2K, CRYSTAL: letti qui (parser propri), confrontati con la cella espansa da gemmi
              (stessi atomi e posizioni, densità); con gli H normalizzati, le distanze X–H da neutroni
    dimer     il dimero dell'input ORCA: ogni atomo è un'immagine per simmetria di un sito del CIF (gemmi), con una
              traslazione di reticolo intera

I programmi non sono installati qui: questo controlla che i file descrivano la struttura giusta, non che i programmi
li accettino (le parole chiave sono state confrontate a mano con i manuali, vedi REGISTRO.md).

    .venv/bin/python validation/check_export.py [periodic dimer]
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

import gemmi
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from pistak.core.analysis import analyse_text  # noqa: E402
from pistak.core.config import Limits  # noqa: E402
from pistak.io import periodic, qc  # noqa: E402

STEMS = (HERE / "platon" / "elenco.txt").read_text().split()[:8]
SYMBOL = {gemmi.Element(z).name: z for z in range(1, 119)}


def oracle(stem: str):
    st = gemmi.read_small_structure(str(HERE / "cod_cache" / f"{stem}.cif"))
    st.setup_cell_images()
    sites = st.get_all_unit_cell_sites()
    el = [s.element.name for s in sites]
    frac = np.array([[s.fract.x, s.fract.y, s.fract.z] for s in sites]) % 1.0
    return st, el, frac


# ---- parser indipendenti dei formati ------------------------------------------------------------------------------
def read_poscar(text: str):
    lines = text.splitlines()
    scale = float(lines[1])
    vec = np.array([[float(x) for x in lines[i].split()] for i in (2, 3, 4)]) * scale
    species, counts = lines[5].split(), [int(x) for x in lines[6].split()]
    el = [e for e, n in zip(species, counts) for _ in range(n)]
    frac = np.array([[float(x) for x in lines[8 + k].split()[:3]] for k in range(len(el))])
    return vec, el, frac


def read_qe(text: str):
    lines = text.splitlines()
    i = next(k for k, line in enumerate(lines) if line.startswith("CELL_PARAMETERS"))
    vec = np.array([[float(x) for x in lines[i + j].split()] for j in (1, 2, 3)])
    i = next(k for k, line in enumerate(lines) if line.startswith("ATOMIC_POSITIONS"))
    el, frac = [], []
    for line in lines[i + 1:]:
        if line.startswith("K_POINTS"):
            break
        p = line.split()
        el.append(p[0])
        frac.append([float(x) for x in p[1:4]])
    return vec, el, np.array(frac)


def read_cp2k(text: str):
    block = text.split("&CELL", 1)[1].split("&END CELL", 1)[0]  # solo il blocco della cella: "C x y z" è anche un atomo
    cell = dict(re.findall(r"^\s+([ABC])\s+(\S+\s+\S+\s+\S+)\s*$", block, re.M))
    vec = np.array([[float(x) for x in cell[k].split()] for k in "ABC"])
    block = text.split("SCALED", 1)[1].split("&END COORD", 1)[0]
    el, frac = [], []
    for line in block.strip().splitlines():
        p = line.split()
        el.append(p[0])
        frac.append([float(x) for x in p[1:4]])
    return vec, el, np.array(frac)


def read_crystal(text: str):
    lines = text.splitlines()
    a, b, c, al, be, ga = (float(x) for x in lines[4].split())
    uc = gemmi.UnitCell(a, b, c, al, be, ga)
    vec = np.array(uc.orth.mat.tolist()).T  # righe = vettori di cella (gemmi: a lungo x, come CRYSTAL)
    n = int(lines[5])
    names = {z: e for e, z in SYMBOL.items()}
    el, frac = [], []
    for line in lines[6:6 + n]:
        p = line.split()
        el.append(names[int(p[0])])
        frac.append([float(x) for x in p[1:4]])
    return vec, el, np.array(frac)


READERS = {"vasp": ("POSCAR", read_poscar), "qe": (".pw.in", read_qe), "cp2k": (".inp", read_cp2k),
           "crystal": (".d12", read_crystal)}


def match(vec, el, frac, ref_el, ref_frac) -> float:
    """Massima distanza (Å) tra ogni sito del riferimento e l'atomo esportato dello stesso elemento più vicino
    (immagine minima)."""
    worst = 0.0
    for e in set(ref_el):
        a = frac[[k for k, x in enumerate(el) if x == e]]
        b = ref_frac[[k for k, x in enumerate(ref_el) if x == e]]
        d = (b[:, None, :] - a[None, :, :] + 0.5) % 1.0 - 0.5
        dist = np.linalg.norm(d @ vec, axis=2).min(axis=1)
        worst = max(worst, float(dist.max()))
    return worst


def check_periodic() -> None:
    print("== file periodici contro la cella espansa da gemmi (H come nel CIF) ==")
    for stem in STEMS:
        st, ref_el, ref_frac = oracle(stem)
        if any(s.occ < 0.99 for s in st.sites):
            continue
        text = (HERE / "cod_cache" / f"{stem}.cif").read_text(encoding="utf-8", errors="replace")
        a = analyse_text(text, stem, Limits())
        row = []
        for prog, (suffix, reader) in READERS.items():
            files = periodic.periodic_files(a, periodic.PeriodicOptions(prog, normalize_h=False))
            name = next(k for k in files if k.endswith(suffix))
            vec, el, frac = reader(files[name])
            same = Counter(el) == Counter(ref_el)
            vol = abs(np.linalg.det(vec))
            dvol = vol / st.cell.volume - 1
            row.append(f"{prog} {'=' if same else '≠'} {match(vec, el, frac, ref_el, ref_frac):.1e} Å {100 * dvol:+.0e}%")
        print(f"  {stem}: {len(ref_el)} atomi; " + "; ".join(row))
    print("  (= stessa composizione; distanza massima sito del CIF → atomo esportato; volume rispetto al CIF)")
    print("\n== con gli H normalizzati: distanze X–H nel POSCAR ==")
    for stem in STEMS[:4]:
        text = (HERE / "cod_cache" / f"{stem}.cif").read_text(encoding="utf-8", errors="replace")
        a = analyse_text(text, stem, Limits())
        vec, el, frac = read_poscar(periodic.periodic_files(a, periodic.PeriodicOptions("vasp"))["POSCAR"])
        st, ref_el, ref_frac = oracle(stem)
        heavy = [k for k, e in enumerate(el) if e != "H"]
        hv = [k for k, e in enumerate(ref_el) if e != "H"]
        moved = match(vec, [el[k] for k in heavy], frac[heavy], [ref_el[k] for k in hv], ref_frac[hv])
        xh = Counter()
        for h in (k for k, e in enumerate(el) if e == "H"):
            d = (frac[heavy] - frac[h] + 0.5) % 1.0 - 0.5
            r = np.linalg.norm(d @ vec, axis=1)
            x = heavy[int(r.argmin())]
            xh[(el[x], round(float(r.min()), 3))] += 1
        print(f"  {stem}: atomi pesanti spostati al massimo di {moved:.1e} Å; X–H: "
              + ", ".join(f"{e}–H {d} ×{n}" for (e, d), n in sorted(xh.items())))


def check_dimer() -> None:
    print("\n== dimero ORCA: ogni atomo è un sito del CIF portato da un'operazione di simmetria + traslazione ==")
    for stem in STEMS:
        text = (HERE / "cod_cache" / f"{stem}.cif").read_text(encoding="utf-8", errors="replace")
        a = analyse_text(text, stem, Limits())
        pairs = qc.exportable(a.pairs)
        if not pairs:
            continue
        inp = qc.dimer_input(a, pairs[0], qc.QCOptions("orca", normalize_h=False, counterpoise=False))
        block = inp.split("* xyz", 1)[1].split("\n*", 1)[0].splitlines()[1:]
        el = [line.split()[0] for line in block]
        cart = np.array([[float(x) for x in line.split()[1:4]] for line in block])
        st, ref_el, ref_frac = oracle(stem)
        frac = np.array([st.cell.fractionalize(gemmi.Position(*p)).tolist() for p in cart])
        worst = 0.0
        for e, f in zip(el, frac):
            cand = ref_frac[[k for k, x in enumerate(ref_el) if x == e]]
            d = (cand - f + 0.5) % 1.0 - 0.5
            dist = np.linalg.norm(np.array([st.cell.orthogonalize(gemmi.Fractional(*v)).tolist() for v in d]), axis=1)
            worst = max(worst, float(dist.min()))
        print(f"  {stem}: coppia {pairs[0].index + 1}, {len(el)} atomi, distanza massima da un sito del cristallo "
              f"{worst:.1e} Å")


CHECKS = {"periodic": check_periodic, "dimer": check_dimer}

if __name__ == "__main__":
    for name in sys.argv[1:] or CHECKS:
        CHECKS[name]()
