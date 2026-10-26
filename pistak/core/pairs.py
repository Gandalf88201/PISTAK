"""Coppie di molecole vicine e geometria del loro π-stacking.

Per ogni coppia, con d il vettore tra i baricentri dei nuclei e n la normale media:
- distanza interplanare = |d·n|;
- slittamento = componente di d nel piano, s = d − (d·n)n;
- slittamento longitudinale = |s·L| lungo l'asse lungo L, slittamento laterale = |s·W| lungo W = n × L.
Una coppia è un π-stack se i piani sono quasi paralleli, abbastanza vicini e i nuclei, proiettati sul piano,
si sovrappongono quanto chiede il criterio di stack (`Limits.overlap_measure`, `min_overlap`, `max_slip`).

Due modi di raccogliere le coppie (`Limits.pair_mode`):

* ``site`` (porting di `measurePair`/`neighbourPairs` di `src/lib/crystal/pistack.ts`): ogni contatto compare due
  volte, a→b e b→a; l'asse L è quello della prima molecola, quindi i due versi possono dare slittamenti diversi;
  le coppie «equivalenti» si riconoscono dai valori arrotondati; si esplorano sempre ±`neighbour_shells` celle.
* ``canonical``: ogni contatto compare una volta (a < b, oppure a = b con la traslazione di segno negativo);
  L è la bisettrice degli assi lunghi delle due molecole, quindi i valori non dipendono da quale si chiama A;
  due contatti sono equivalenti se un'operazione del gruppo spaziale porta l'uno nell'altro; le celle esplorate
  bastano sempre a trovare tutti i contatti (almeno ±`neighbour_shells`). `count` = contatti per cella.
"""
from __future__ import annotations

import itertools
import math
from collections.abc import Iterator
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial.distance import cdist

from pistak.core.config import Limits
from pistak.core.cores import Core
from pistak.core.geometry import clip_convex, convex_hull, points_inside, polygon_area, unit
from pistak.core.structure import Structure

KINDS = ("stack", "parallel", "tilted")  # in ordine di interesse
Shift = tuple[int, int, int]
#: Limite alle celle esplorate per asse in modalità canonical (celle molto piccole o molto oblique).
MAX_SHELLS = 8


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


def pair_frame(ca: Core, cb: Core, symmetric: bool = False) -> Frame:
    """n = normale media; L = asse lungo di A proiettato sul piano medio, oppure (``symmetric``) la bisettrice degli
    assi lunghi di A e B, che non cambia scambiando le due molecole."""
    nb = -cb.normal if np.dot(ca.normal, cb.normal) < 0 else cb.normal
    n = unit(ca.normal + nb)
    la = ca.long - n * np.dot(ca.long, n)
    if symmetric:
        lb = cb.long - n * np.dot(cb.long, n)
        if np.dot(la, lb) < 0:
            lb = -lb
        if np.linalg.norm(la + lb) > 1e-6:
            la = la + lb
    L = unit(la)
    return Frame(ca.centroid, L, np.cross(n, L), n)


def frame_for(ca: Core, cb: Core, limits: Limits) -> Frame:
    return pair_frame(ca, cb, symmetric=limits.pair_mode != "site")


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
    #: Quante coppie equivalenti sono state raccolte in questa (canonical: contatti per cella).
    count: int = 1
    #: Tutte le coppie raccolte (a, b, shift), la prima è questa.
    members: list[tuple[int, int, Shift]] = field(default_factory=list)
    #: Valori di tutti i parametri del registro (`pistak.params`), per chiave; riempiti da `compute_all`.
    values: dict[str, float] = field(default_factory=dict)
    index: int = -1

    def get(self, key: str, default: float = float("nan")) -> float:
        """Valore di un parametro del registro per questa coppia."""
        return self.values.get(key, default)


def area_overlap(hull_a: np.ndarray, hull_b: np.ndarray) -> float:
    """Area comune dei due contorni divisa per l'area del più piccolo (0–1)."""
    if len(hull_a) < 3 or len(hull_b) < 3:
        return 0.0
    smaller = min(polygon_area(hull_a), polygon_area(hull_b))
    return polygon_area(clip_convex(hull_a, hull_b)) / smaller if smaller > 0 else 0.0


def classify(dihedral: float, interplanar: float, overlap: float, slip: float, limits: Limits) -> str:
    """stack | parallel | tilted secondo le soglie. `overlap` è già misurata come chiede `limits.overlap_measure`."""
    if dihedral > limits.max_dihedral:
        return "tilted"
    if interplanar <= limits.max_interplanar and overlap > limits.min_overlap and (limits.max_slip <= 0 or slip <= limits.max_slip):
        return "stack"
    return "parallel"


def measure_pair(s: Structure, cores: list[Core | None], a: int, b: int, shift: Shift, limits: Limits,
                 max_contact: float | None = None) -> Pair | None:
    """Geometria della coppia (a nella cella 0, b traslata di `shift`). None se manca un nucleo o se il contatto
    più corto tra atomi pesanti supera `max_contact`."""
    ca, cb0 = cores[a], cores[b]
    if ca is None or cb0 is None:
        return None
    v = s.shift_vector(shift)
    A, B = s.molecules[a], s.molecules[b]
    pa, pb = A.pos, B.pos + v
    ha, hb = A.heavy, B.heavy
    contact = float(cdist(pa[ha], pb[hb]).min()) if ha.any() and hb.any() else float("inf")
    if max_contact is not None and contact > max_contact:
        return None
    cb = cb0.moved(v)
    f = frame_for(ca, cb, limits)
    d = cb.centroid - ca.centroid
    z = float(np.dot(d, f.n))
    sv = d - f.n * z
    dihedral = float(np.degrees(np.arccos(min(1.0, abs(np.dot(ca.normal, cb.normal))))))
    qa = f.project(pa[ca.atoms])[:, :2]
    qb = f.project(pb[cb.atoms])[:, :2]
    hull_a, hull_b = convex_hull(qa), convex_hull(qb)
    overlap = max(points_inside(hull_a, qb).mean(), points_inside(hull_b, qa).mean())
    interplanar = abs(z)
    slip = float(np.linalg.norm(sv))
    measured = overlap if limits.overlap_measure != "area" else area_overlap(hull_a, hull_b)
    kind = classify(dihedral, interplanar, float(measured), slip, limits)
    return Pair(
        a=a, b=b, shift=shift, kind=kind, interplanar=interplanar,
        long_slip=abs(float(np.dot(sv, f.L))), lat_slip=abs(float(np.dot(sv, f.W))), slip=slip,
        dihedral=dihedral, centroid_distance=float(np.linalg.norm(d)), contact=contact, overlap=float(overlap),
    )


def neighbour_shells(s: Structure, cores: list[Core | None], limits: Limits) -> tuple[int, int, int]:
    """Celle da esplorare (±) lungo a, b, c. Nel modo del sito sempre `neighbour_shells`; altrimenti quante ne
    servono perché due molecole qualunque a contatto siano trovate, ma almeno `neighbour_shells`.

    Due molecole si toccano solo se |Δbaricentri| ≤ D = E_a + E_b + contact_cutoff (E = estensione); la
    componente frazionaria i di un vettore g è al più |g|·|a*_i| e i baricentri stanno nella cella [0, 1), quindi
    basta esplorare |k_i| ≤ D·|a*_i| + 1."""
    R = limits.neighbour_shells
    if limits.pair_mode == "site":
        return (R, R, R)
    ext = [m.extent for m, c in zip(s.molecules, cores) if c is not None and not m.extended]
    if not ext:
        return (R, R, R)
    D = 2 * max(ext) + limits.contact_cutoff
    recip = np.linalg.norm(np.linalg.inv(s.vectors), axis=0)  # |a*_i|: colonne di M⁻¹
    need = np.floor(D * recip + 1).astype(int)
    return tuple(int(min(MAX_SHELLS, max(R, n))) for n in need)  # type: ignore[return-value]


def _neg(sh: Shift) -> Shift:
    return (-sh[0], -sh[1], -sh[2])


def _canon(a: int, b: int, sh: Shift) -> tuple[int, int, Shift]:
    """Lo stesso contatto scritto in un solo modo: (a, b, s) e (b, a, −s) sono la stessa coppia."""
    return min((a, b, sh), (b, a, _neg(sh)))


def iter_neighbour_pairs(s: Structure, cores: list[Core | None], limits: Limits) -> Iterator[Pair]:
    """Tutte le coppie di molecole vicine (contatto tra atomi pesanti entro `contact_cutoff`), una per traslazione,
    senza raccogliere quelle equivalenti. In modalità site ogni contatto compare nei due versi, altrimenti in uno."""
    na, nb, nc = neighbour_shells(s, cores, limits)
    shifts = np.array(list(itertools.product(range(-na, na + 1), range(-nb, nb + 1), range(-nc, nc + 1))), dtype=int)
    shift_vec = shifts @ s.vectors
    extent = [m.extent for m in s.molecules]
    usable = [i for i, m in enumerate(s.molecules) if cores[i] is not None and not m.extended]
    one_way = limits.pair_mode != "site"
    for a in usable:
        for b in usable:
            if one_way and b < a:
                continue
            # scarto veloce: baricentri troppo lontani perché le molecole si tocchino
            gap = np.linalg.norm(s.molecules[b].centroid + shift_vec - s.molecules[a].centroid, axis=1)
            near = gap <= extent[a] + extent[b] + limits.contact_cutoff
            for k in np.flatnonzero(near):
                sh = tuple(int(x) for x in shifts[k])
                if a == b and sh == (0, 0, 0):
                    continue
                if one_way and _canon(a, b, sh) != (a, b, sh):
                    continue
                p = measure_pair(s, cores, a, b, sh, limits, max_contact=limits.contact_cutoff)
                if p is not None:
                    yield p


StackEdge = tuple[int, int, Shift, Pair]


def collect_pairs(s: Structure, cores: list[Core | None], limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
    """Coppie di molecole vicine, con quelle equivalenti raccolte in una sola, e, a parte, tutti gli archi `stack`
    (nei due versi: servono a ricostruire le colonne)."""
    if limits.pair_mode == "site":
        pairs, stack_edges = _collect_site(s, cores, limits)
    else:
        pairs, stack_edges = _collect_canonical(s, cores, limits)
    order = {k: i for i, k in enumerate(KINDS)}
    pairs.sort(key=lambda p: (order[p.kind], p.interplanar, p.contact))
    for i, p in enumerate(pairs):
        p.index = i
    return pairs, stack_edges


def _collect_site(s: Structure, cores: list[Core | None], limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
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
    return list(found.values()), stack_edges


def symmetry_key(s: Structure, a: int, b: int, sh: Shift) -> tuple[int, int, Shift]:
    """Rappresentante dell'orbita del contatto sotto il gruppo spaziale: la più piccola delle sue immagini (scritte
    con `_canon`). Due contatti sono equivalenti per simmetria se e solo se hanno la stessa chiave; le operazioni che
    non mappano le molecole (`Structure.molecule_images` = None) non si usano, così non si uniscono mai contatti
    diversi."""
    best = _canon(a, b, sh)
    for k, row in enumerate(s.molecule_images()):
        ia, ib = row[a], row[b]
        if ia is None or ib is None:
            continue
        (a2, ta), (b2, tb) = ia, ib
        rs = s.rotation(k) @ np.array(sh)
        new = (tb[0] + int(rs[0]) - ta[0], tb[1] + int(rs[1]) - ta[1], tb[2] + int(rs[2]) - ta[2])
        if a2 == b2 and new == (0, 0, 0):
            continue
        best = min(best, _canon(a2, b2, new))
    return best


def _collect_canonical(s: Structure, cores: list[Core | None], limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
    classes: dict[tuple, list[Pair]] = {}
    stack_edges: list[StackEdge] = []
    for p in iter_neighbour_pairs(s, cores, limits):
        if p.kind == "stack":
            stack_edges += [(p.a, p.b, p.shift, p), (p.b, p.a, _neg(p.shift), p)]
        classes.setdefault(symmetry_key(s, p.a, p.b, p.shift), []).append(p)
    pairs = []
    for key, members in classes.items():
        members.sort(key=lambda p: (p.a, p.b, p.shift))
        rep = next((p for p in members if (p.a, p.b, p.shift) == key), members[0])
        rep.count = len(members)
        rep.members = [(rep.a, rep.b, rep.shift)] + [(p.a, p.b, p.shift) for p in members if p is not rep]
        pairs.append(rep)
    return pairs, stack_edges


def class_spread(s: Structure, cores: list[Core | None], limits: Limits, pair: Pair) -> float:
    """Differenza massima di distanza interplanare e slittamenti tra i membri della classe di `pair` (Å): ~0 se
    l'equivalenza per simmetria è esatta. Serve ai controlli."""
    vals = []
    for a, b, sh in pair.members:
        q = measure_pair(s, cores, a, b, sh, limits)
        if q is not None:
            vals.append((q.interplanar, q.long_slip, q.lat_slip, q.contact))
    if len(vals) < 2:
        return 0.0
    arr = np.array(vals)
    return float(np.max(arr.max(axis=0) - arr.min(axis=0))) if math.isfinite(arr.sum()) else float("nan")


def neighbour_pairs(s: Structure, cores: list[Core | None], limits: Limits) -> list[Pair]:
    return collect_pairs(s, cores, limits)[0]
