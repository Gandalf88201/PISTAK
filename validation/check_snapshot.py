"""Checklist, sezione 8: le istantanee COD dei test (`tests/snapshots/canonical/*.json`) rifatte a parte.

Le istantanee sono state scritte dal programma stesso: prima di considerarle «vere» si ricalcola ogni coppia con un
calcolo indipendente in numpy. Per ogni coppia salvata si prendono le posizioni degli atomi dei due nuclei π (quali
atomi formano un nucleo è una scelta di PISTAK, verificata nella sezione 3) e si rifanno piani ai minimi quadrati,
normale media, distanza interplanare, slittamento, angolo tra i piani, distanza tra i baricentri e contatto più corto
tra atomi pesanti delle due molecole.

    .venv/bin/python validation/check_snapshot.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from pistak.core.analysis import analyse_file  # noqa: E402
from pistak.core.config import Limits  # noqa: E402

SNAP = ROOT / "tests" / "snapshots" / "canonical"
COD = ROOT / "tests" / "data" / "cod"


def plane(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    c = p.mean(axis=0)
    return c, np.linalg.svd(p - c)[2][-1]


def geometry(a: np.ndarray, b: np.ndarray, mol_a: np.ndarray, mol_b: np.ndarray) -> dict[str, float]:
    ca, na = plane(a)
    cb, nb = plane(b)
    if na @ nb < 0:
        nb = -nb
    n = (na + nb) / np.linalg.norm(na + nb)
    d = cb - ca
    z = d @ n
    contact = float(np.min(np.linalg.norm(mol_a[:, None, :] - mol_b[None, :, :], axis=2)))
    return {"interplanar": abs(z), "slip": float(np.linalg.norm(d - z * n)),
            "dihedral": math.degrees(math.acos(min(1.0, abs(na @ nb)))), "centroid_distance": float(np.linalg.norm(d)),
            "contact": contact}


def main() -> None:
    print("== istantanee COD contro un calcolo indipendente (numpy) ==")
    worst_all: dict[str, float] = {}
    for snap in sorted(SNAP.glob("[0-9]*.json")):
        data = json.loads(snap.read_text())
        an = analyse_file(COD / f"{snap.stem}.cif", Limits())  # le soglie predefinite, come le istantanee
        s = an.structure
        frags = an.fragments or [[c] if c else [] for c in an.cores]
        worst = {}
        n = 0
        for sp in data["pairs"]:
            if sp["a"] == sp["b"] and not any(sp["shift"]):
                continue  # coppie intramolecolari: il contatto è tra frammenti, non tra molecole
            par = sp["params"]
            cands = [q for q in an.pairs if q.a == sp["a"] and q.b == sp["b"] and list(q.shift) == list(sp["shift"])]
            # l'istantanea non salva gli indici dei frammenti: tra le coppie con le stesse molecole e traslazione (più
            # frammenti per molecola, es. i quattro isoindoli di H₂Pc) si sceglie quella più vicina in Cg···Cg e d⊥
            p = min(cands, key=lambda q: abs(q.centroid_distance - par["centroid_distance"]) + abs(q.interplanar - par["interplanar"]))
            A, B = s.molecules[p.a], s.molecules[p.b]
            v = s.shift_vector(p.shift)
            fa, fb = frags[p.a][p.fa], frags[p.b][p.fb]
            # per le coppie di polimero PISTAK misura il contatto tra i due frammenti, non tra le unità ripetute (voce 47)
            ca, cb = (A.pos[fa.atoms], B.pos[fb.atoms] + v) if p.polymer else (A.pos[A.heavy], B.pos[B.heavy] + v)
            g = geometry(A.pos[fa.atoms], B.pos[fb.atoms] + v, ca, cb)
            n += 1
            for k, val in g.items():
                worst[k] = max(worst.get(k, 0.0), abs(val - par[k]))
        for k, x in worst.items():
            worst_all[k] = max(worst_all.get(k, 0.0), x)
        print(f"  {snap.stem}: {n} coppie; scarto massimo " + ", ".join(f"{k} {x:.1e}" for k, x in worst.items()))
    print("  in tutto: " + ", ".join(f"{k} {x:.1e}" for k, x in worst_all.items()))


if __name__ == "__main__":
    main()
