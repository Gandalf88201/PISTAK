"""Coppie di molecole vicine e geometria del loro π-stacking.

Per ogni coppia, con d il vettore tra i baricentri dei nuclei e n la normale media:
- distanza interplanare = |d·n|;
- slittamento = componente di d nel piano, s = d − (d·n)n;
- slittamento longitudinale = |s·L| lungo l'asse lungo L della prima molecola,
  slittamento laterale = |s·W| lungo l'asse corto W = n × L.
Una coppia è un π-stack se i piani sono quasi paralleli, abbastanza vicini e i nuclei,
proiettati sul piano, si sovrappongono almeno in parte.

Porting di `measurePair`/`neighbourPairs` di `src/lib/crystal/pistack.ts`.
"""
from __future__ import annotations

import itertools
from collections.abc import Iterator
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial.distance import cdist

from pistak.core.config import Limits
from pistak.core.cores import Core
from pistak.core.geometry import convex_hull, points_inside, unit
from pistak.core.structure import Structure

KINDS = ("stack", "parallel", "tilted")  # in ordine di interesse
Shift = tuple[int, int, int]


@dataclass
class Frame:
    """Sistema di riferimento della coppia: origine nel baricentro del primo nucleo,
    assi L, W nel piano medio, n normale."""

    origin: np.ndarray
    L: np.ndarray
    W: np.ndarray
    n: np.ndarray

    def project(self, p: np.ndarray) -> np.ndarray:
        """Coordinate (L, W, n) di uno o più punti 3D."""
        return (np.asarray(p) - self.origin) @ np.column_stack([self.L, self.W, self.n])


def pair_frame(ca: Core, cb: Core) -> Frame:
    nb = -cb.normal if np.dot(ca.normal, cb.normal) < 0 else cb.normal
    n = unit(ca.normal + nb)
    L = unit(ca.long - n * np.dot(ca.long, n))
    return Frame(ca.centroid, L, np.cross(n, L), n)


@dataclass
class Pair:
    """Una coppia di molecole vicine. Gli indici si riferiscono a `Structure.molecules`."""

    a: int
    b: int
    shift: Shift  # traslazione (in celle) della seconda molecola
    kind: str  # stack | parallel | tilted
    interplanar: float
    long_slip: float
    lat_slip: float
    slip: float
    dihedral: float
    centroid_distance: float
    contact: float
    #: Frazione degli atomi del nucleo di una molecola che cade, in proiezione, dentro il contorno dell'altra.
    overlap: float
    #: Quante coppie equivalenti (stessi valori) sono state raccolte in questa.
    count: int = 1
    #: Tutte le coppie raccolte (a, b, shift), la prima è questa.
    members: list[tuple[int, int, Shift]] = field(default_factory=list)
    #: Valori di tutti i parametri del registro (`pistak.params`), per chiave; riempiti da `compute_all`.
    values: dict[str, float] = field(default_factory=dict)
    index: int = -1

    def get(self, key: str, default: float = float("nan")) -> float:
        """Valore di un parametro del registro per questa coppia."""
        return self.values.get(key, default)


def measure_pair(s: Structure, cores: list[Core | None], a: int, b: int, shift: Shift, limits: Limits) -> Pair | None:
    ca, cb0 = cores[a], cores[b]
    if ca is None or cb0 is None:
        return None
    v = s.shift_vector(shift)
    cb = cb0.moved(v)
    A, B = s.molecules[a], s.molecules[b]
    pa, pb = A.pos, B.pos + v
    ha, hb = A.heavy, B.heavy
    contact = float(cdist(pa[ha], pb[hb]).min()) if ha.any() and hb.any() else float("inf")
    f = pair_frame(ca, cb)
    d = cb.centroid - ca.centroid
    z = float(np.dot(d, f.n))
    sv = d - f.n * z
    dihedral = float(np.degrees(np.arccos(min(1.0, abs(np.dot(ca.normal, cb.normal))))))
    qa = f.project(pa[ca.atoms])[:, :2]
    qb = f.project(pb[cb.atoms])[:, :2]
    hull_a, hull_b = convex_hull(qa), convex_hull(qb)
    overlap = max(points_inside(hull_a, qb).mean(), points_inside(hull_b, qa).mean())
    interplanar = abs(z)
    if dihedral > limits.max_dihedral:
        kind = "tilted"
    elif interplanar <= limits.max_interplanar and overlap > 0:
        kind = "stack"
    else:
        kind = "parallel"
    return Pair(
        a=a, b=b, shift=shift, kind=kind, interplanar=interplanar,
        long_slip=abs(float(np.dot(sv, f.L))), lat_slip=abs(float(np.dot(sv, f.W))), slip=float(np.linalg.norm(sv)),
        dihedral=dihedral, centroid_distance=float(np.linalg.norm(d)), contact=contact, overlap=float(overlap),
    )


def iter_neighbour_pairs(s: Structure, cores: list[Core | None], limits: Limits) -> Iterator[Pair]:
    """Tutte le coppie di molecole vicine (contatto tra atomi pesanti entro `contact_cutoff`), una per traslazione,
    senza raccogliere quelle equivalenti."""
    R = limits.neighbour_shells
    shifts = np.array(list(itertools.product(range(-R, R + 1), repeat=3)), dtype=int)
    shift_vec = shifts @ s.vectors
    extent = [m.extent for m in s.molecules]
    usable = [i for i, m in enumerate(s.molecules) if cores[i] is not None and not m.extended]
    for a in usable:
        for b in usable:
            # scarto veloce: baricentri troppo lontani perché le molecole si tocchino
            gap = np.linalg.norm(s.molecules[b].centroid + shift_vec - s.molecules[a].centroid, axis=1)
            near = gap <= extent[a] + extent[b] + limits.contact_cutoff
            for k in np.flatnonzero(near):
                sh = tuple(int(x) for x in shifts[k])
                if a == b and sh == (0, 0, 0):
                    continue
                p = measure_pair(s, cores, a, b, sh, limits)
                if p is not None and p.contact <= limits.contact_cutoff:
                    yield p


StackEdge = tuple[int, int, Shift, Pair]


def collect_pairs(s: Structure, cores: list[Core | None], limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
    """Coppie di molecole vicine, con le coppie dagli stessi valori (legate dalla simmetria) raccolte in una sola,
    e, a parte, tutti gli archi `stack` non raccolti (servono a ricostruire le colonne)."""
    found: dict[str, Pair] = {}
    stack_edges: list[StackEdge] = []
    for p in iter_neighbour_pairs(s, cores, limits):
        if p.kind == "stack":
            stack_edges.append((p.a, p.b, p.shift, p))
        key = "|".join([p.kind, f"{p.interplanar:.2f}", f"{p.long_slip:.2f}", f"{p.lat_slip:.2f}",
                        f"{p.dihedral:.1f}", f"{p.contact:.2f}"])
        prev = found.get(key)
        if prev is not None:
            prev.count += 1
            prev.members.append((p.a, p.b, p.shift))
        else:
            p.members = [(p.a, p.b, p.shift)]
            found[key] = p
    order = {k: i for i, k in enumerate(KINDS)}
    pairs = sorted(found.values(), key=lambda p: (order[p.kind], p.interplanar, p.contact))
    for i, p in enumerate(pairs):
        p.index = i
    return pairs, stack_edges


def neighbour_pairs(s: Structure, cores: list[Core | None], limits: Limits) -> list[Pair]:
    return collect_pairs(s, cores, limits)[0]
