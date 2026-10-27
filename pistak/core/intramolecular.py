"""Geometria delle molecole: lunghezze di legame, angoli di legame e angoli di torsione, con le loro incertezze.

Si misurano solo le molecole indipendenti per simmetria: per ogni insieme di molecole della cella legate da
un'operazione del gruppo si prende quella che contiene più atomi dell'unità asimmetrica così come sono nel CIF.
Gli atomi generati da un'altra operazione (molecola su un elemento di simmetria) hanno l'etichetta con l'apice
(C1′); le misure uguali per simmetria interna della molecola (C1–C2 e C1′–C2′) compaiono una volta sola.

Incertezze: si propagano le s.u. delle coordinate frazionarie del CIF attraverso il gradiente della misura, tenendo
conto che due atomi copie dello stesso atomo dell'unità asimmetrica si muovono insieme. Come in `uncertainty.py`
le coordinate sono considerate indipendenti tra loro e la cella esatta; se un atomo coinvolto non ha s.u. (idrogeni
in posizione calcolata) l'incertezza è NaN.
"""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

from pistak.core.cif import Cif
from pistak.core.structure import Structure

KINDS = ("bond", "angle", "torsion")
PRIME = "′"
_STEP = 1e-4  # Å, per le differenze finite


@dataclass(frozen=True)
class GeometryItem:
    kind: str  # bond | angle | torsion
    #: Indice della molecola in `Structure.molecules`.
    molecule: int
    #: Etichette degli atomi, nell'ordine della misura (per gli angoli l'atomo centrale è il secondo).
    atoms: tuple[str, ...]
    value: float  # Å o gradi
    su: float  # stessa unità; NaN se non si può stimare


def distance(p: np.ndarray) -> float:
    return float(np.linalg.norm(p[1] - p[0]))


def angle(p: np.ndarray) -> float:
    u, v = p[0] - p[1], p[2] - p[1]
    c = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def torsion(p: np.ndarray) -> float:
    """Angolo di torsione 1-2-3-4 con segno (convenzione IUPAC/Klyne–Prelog, −180…180°)."""
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1n, v), w), np.dot(v, w))))


MEASURES: dict[str, Callable[[np.ndarray], float]] = {"bond": distance, "angle": angle, "torsion": torsion}


def _identity_index(s: Structure) -> int:
    for k, op in enumerate(s.symops):
        if np.allclose(op.r, np.eye(3)) and np.allclose(op.t - np.round(op.t), 0):
            return k
    return -1


def unique_molecules(s: Structure) -> list[int]:
    """Una molecola per ogni insieme di molecole equivalenti per simmetria (quella con più atomi dell'unità
    asimmetrica nella loro posizione del CIF), in ordine di indice."""
    n = len(s.molecules)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for row in s.molecule_images():
        for m, img in enumerate(row):
            if img is not None:
                parent[find(m)] = find(img[0])
    ident = _identity_index(s)
    orbits: dict[int, list[int]] = {}
    for m in range(n):
        orbits.setdefault(find(m), []).append(m)

    def own(m: int) -> int:
        ci = s.molecules[m].cell_index
        return int(np.count_nonzero(s.cell_origin[ci, 1] == ident)) if len(s.cell_origin) else 0

    return sorted(max(ms, key=lambda m: (own(m), -m)) for ms in orbits.values())


class _Propagator:
    """s.u. di una misura dalle s.u. delle coordinate degli atomi dell'unità asimmetrica."""

    def __init__(self, s: Structure, cif: Cif):
        self.s = s
        self.su = [a.su for a in cif.atoms]

    def __call__(self, fn: Callable[[np.ndarray], float], pos: np.ndarray, cell_atoms: np.ndarray) -> float:
        if not len(self.s.cell_origin):
            return math.nan
        jac: dict[int, np.ndarray] = {}
        for i, ci in enumerate(cell_atoms):
            ia, k = (int(x) for x in self.s.cell_origin[ci])
            if not np.any(self.su[ia]):
                return math.nan
            g = np.zeros(3)
            for d in range(3):
                plus, minus = pos.copy(), pos.copy()
                plus[i, d] += _STEP
                minus[i, d] -= _STEP
                g[d] = _wrapped(fn, plus, minus) / (2 * _STEP)
            # posizione cartesiana = Vᵀ·(R·f + t): derivata rispetto alle coordinate frazionarie dell'atomo del CIF
            jac[ia] = jac.get(ia, np.zeros(3)) + g @ self.s.vectors.T @ self.s.symops[k].r
        return float(math.sqrt(sum(float(np.sum((j * self.su[ia]) ** 2)) for ia, j in jac.items())))


def _wrapped(fn, plus: np.ndarray, minus: np.ndarray) -> float:
    d = fn(plus) - fn(minus)
    return (d + 180.0) % 360.0 - 180.0 if fn is torsion else d


def molecule_geometry(s: Structure, cif: Cif, molecules: list[int] | None = None, include_h: bool = True,
                      kinds: tuple[str, ...] = KINDS) -> list[GeometryItem]:
    """Legami, angoli e torsioni delle molecole indicate (per impostazione quelle indipendenti per simmetria)."""
    ident = _identity_index(s)
    prop = _Propagator(s, cif)
    out: list[GeometryItem] = []
    for mi in unique_molecules(s) if molecules is None else molecules:
        m = s.molecules[mi]
        ops = s.cell_origin[m.cell_index, 1] if len(s.cell_origin) else np.full(len(m), ident)
        names = [lab + ("" if op == ident else PRIME) for lab, op in zip(m.labels, ops)]
        nb: list[list[int]] = [[] for _ in m.el]
        for i, j in m.bonds:
            nb[i].append(int(j))
            nb[j].append(int(i))
        for row in nb:
            row.sort()
        groups: dict[str, list[tuple[int, ...]]] = {"bond": [], "angle": [], "torsion": []}
        groups["bond"] = [(int(i), int(j)) for i, j in m.bonds]
        groups["angle"] = [(i, j, k) for j in range(len(m.el)) for x, i in enumerate(nb[j]) for k in nb[j][x + 1:]]
        groups["torsion"] = [(i, j, k, t) for j, k in groups["bond"] for i in nb[j] if i != k
                             for t in nb[k] if t not in (j, i)]
        seen: set[tuple[str, ...]] = set()
        for kind in kinds:
            fn = MEASURES[kind]
            for idx in groups[kind]:
                if not include_h and any(m.el[i] == "H" for i in idx):
                    continue
                base = tuple(m.labels[i] for i in idx)
                key = (kind, *min(base, base[::-1]))
                if key in seen:
                    continue
                seen.add(key)
                if base[::-1] < base and kind != "torsion":  # ordine di lettura stabile (la torsione no: ha un verso)
                    idx = idx[::-1]
                p = m.pos[list(idx)]
                out.append(GeometryItem(kind, mi, tuple(names[i] for i in idx), fn(p), prop(fn, p, m.cell_index[list(idx)])))
    return out
