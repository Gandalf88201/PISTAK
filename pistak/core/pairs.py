"""Coppie di frammenti π vicini e geometria del loro π-stacking.

Per ogni coppia, con d il vettore tra i baricentri dei due nuclei (frammenti π) e n la normale media:
- distanza interplanare = |d·n|;
- slittamento = componente di d nel piano, s = d − (d·n)n;
- slittamento longitudinale = |s·L| lungo l'asse lungo L, slittamento laterale = |s·W| lungo W = n × L.
Una coppia è un π-stack se i piani sono quasi paralleli, abbastanza vicini e i nuclei, proiettati sul piano,
si sovrappongono quanto chiede il criterio di stack (`Limits.overlap_measure`, `min_overlap`, `max_slip`).

Ogni molecola ha uno o più frammenti π (`cores.pi_fragments`, `Limits.fragments`): una coppia è tra il frammento
`fa` della molecola `a` e il frammento `fb` della molecola `b` (traslata di `shift`). Tra due molecole a contatto si
misurano il frammento principale di ciascuna (0–0) e le altre combinazioni di frammenti abbastanza vicini; con più
frammenti nella stessa molecola anche le coppie intramolecolari (a = b, shift = 0).

Due modi di raccogliere le coppie (`Limits.pair_mode`):

* ``site`` (porting di `measurePair`/`neighbourPairs` di `src/lib/crystal/pistack.ts`): ogni contatto compare due
  volte, a→b e b→a; l'asse L è quello della prima molecola, quindi i due versi possono dare slittamenti diversi;
  le coppie «equivalenti» si riconoscono dai valori arrotondati; si esplorano sempre ±`neighbour_shells` celle.
* ``canonical``: ogni contatto compare una volta (a < b, oppure a = b con la traslazione di segno negativo);
  L è la bisettrice degli assi lunghi dei due frammenti, quindi i valori non dipendono da quale si chiama A;
  due contatti sono equivalenti se un'operazione del gruppo spaziale porta l'uno nell'altro (molecole e frammenti);
  le celle esplorate bastano sempre a trovare tutti i contatti (almeno ±`neighbour_shells`). `count` = contatti
  per cella.
"""
from __future__ import annotations

import itertools
import math
from collections.abc import Iterator, Sequence
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
#: Due frammenti secondari formano una coppia se i loro atomi sono entro `contact_cutoff` + questo margine (Å).
FRAGMENT_MARGIN = 1.0
#: (a, b, shift, fa, fb): un contatto tra il frammento fa della molecola a e il frammento fb di b traslata di shift.
Member = tuple[int, int, Shift, int, int]
Fragments = list[list[Core]]


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
    """Una coppia di frammenti π vicini. `a`, `b` sono indici in `Structure.molecules`, `fa`, `fb` i frammenti."""

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
    #: Contatto più corto tra atomi pesanti delle due molecole (intramolecolare: tra i due frammenti).
    contact: float
    #: Frazione degli atomi del nucleo di una molecola che cade, in proiezione, dentro il contorno dell'altra.
    overlap: float
    #: Quante coppie equivalenti sono state raccolte in questa (canonical: contatti per cella).
    count: int = 1
    #: Tutte le coppie raccolte (a, b, shift, fa, fb), la prima è questa.
    members: list[Member] = field(default_factory=list)
    #: Valori di tutti i parametri del registro (`pistak.params`), per chiave; riempiti da `compute_all`.
    values: dict[str, float] = field(default_factory=dict)
    index: int = -1
    fa: int = 0
    fb: int = 0
    #: Vero se almeno una delle due molecole ha più di un frammento π (allora conta quali sono fa e fb).
    fragmented: bool = False

    def get(self, key: str, default: float = float("nan")) -> float:
        """Valore di un parametro del registro per questa coppia."""
        return self.values.get(key, default)

    @property
    def intramolecular(self) -> bool:
        return self.a == self.b and self.shift == (0, 0, 0)

    @property
    def key(self) -> Member:
        return (self.a, self.b, self.shift, self.fa, self.fb)


def as_fragments(cores: Sequence) -> Fragments:
    """Accetta sia l'elenco dei nuclei principali (`Core | None` per molecola) sia quello dei frammenti."""
    return [list(c) if isinstance(c, (list, tuple)) else ([c] if c is not None else []) for c in cores]


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


def _heavy_contact(s: Structure, a: int, b: int, v: np.ndarray) -> float:
    A, B = s.molecules[a], s.molecules[b]
    ha, hb = A.heavy, B.heavy
    return float(cdist(A.pos[ha], B.pos[hb] + v).min()) if ha.any() and hb.any() else float("inf")


def _fragment_distance(s: Structure, ca: Core, cb: Core, a: int, b: int, v: np.ndarray) -> float:
    return float(cdist(s.molecules[a].pos[ca.atoms], s.molecules[b].pos[cb.atoms] + v).min())


def _measure(s: Structure, ca: Core, cb0: Core, a: int, b: int, shift: Shift, limits: Limits, contact: float,
             fa: int, fb: int, fragmented: bool) -> Pair:
    v = s.shift_vector(shift)
    pa, pb = s.molecules[a].pos, s.molecules[b].pos + v
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
        fa=fa, fb=fb, fragmented=fragmented,
    )


def measure_pair(s: Structure, cores: Sequence, a: int, b: int, shift: Shift, limits: Limits,
                 max_contact: float | None = None, fa: int = 0, fb: int = 0) -> Pair | None:
    """Geometria della coppia (frammento fa di a nella cella 0, frammento fb di b traslata di `shift`). None se manca
    un frammento o se il contatto più corto tra atomi pesanti supera `max_contact`. `cores` può essere l'elenco dei
    nuclei principali o quello dei frammenti di ogni molecola."""
    frags = as_fragments(cores)
    if fa >= len(frags[a]) or fb >= len(frags[b]):
        return None
    ca, cb0 = frags[a][fa], frags[b][fb]
    v = s.shift_vector(shift)
    if a == b and shift == (0, 0, 0):
        contact = _fragment_distance(s, ca, cb0, a, b, v)
    else:
        contact = _heavy_contact(s, a, b, v)
    if max_contact is not None and contact > max_contact:
        return None
    return _measure(s, ca, cb0, a, b, shift, limits, contact, fa, fb, len(frags[a]) > 1 or len(frags[b]) > 1)


def remeasure(s: Structure, frags: Fragments, pair: Pair, limits: Limits) -> Pair:
    """La stessa coppia (stesse molecole, frammenti e traslazione) misurata su un'altra geometria: serve al Monte Carlo
    delle incertezze, dove gli atomi si spostano di poco ma la topologia resta quella dell'analisi."""
    ca, cb0 = frags[pair.a][pair.fa], frags[pair.b][pair.fb]
    v = s.shift_vector(pair.shift)
    contact = _fragment_distance(s, ca, cb0, pair.a, pair.b, v) if pair.intramolecular else _heavy_contact(s, pair.a, pair.b, v)
    q = _measure(s, ca, cb0, pair.a, pair.b, pair.shift, limits, contact, pair.fa, pair.fb, pair.fragmented)
    q.index, q.count, q.members = pair.index, pair.count, pair.members
    return q


def neighbour_shells(s: Structure, cores: Sequence, limits: Limits) -> tuple[int, int, int]:
    """Celle da esplorare (±) lungo a, b, c. Nel modo del sito sempre `neighbour_shells`; altrimenti quante ne
    servono perché due molecole qualunque a contatto siano trovate, ma almeno `neighbour_shells`.

    Due molecole si toccano solo se |Δbaricentri| ≤ D = E_a + E_b + contact_cutoff (E = estensione); la
    componente frazionaria i di un vettore g è al più |g|·|a*_i| e i baricentri stanno nella cella [0, 1), quindi
    basta esplorare |k_i| ≤ D·|a*_i| + 1."""
    R = limits.neighbour_shells
    if limits.pair_mode == "site":
        return (R, R, R)
    frags = as_fragments(cores)
    ext = [m.extent for m, f in zip(s.molecules, frags) if f and not m.extended]
    if not ext:
        return (R, R, R)
    D = 2 * max(ext) + limits.contact_cutoff
    recip = np.linalg.norm(np.linalg.inv(s.vectors), axis=0)  # |a*_i|: colonne di M⁻¹
    need = np.floor(D * recip + 1).astype(int)
    return tuple(int(min(MAX_SHELLS, max(R, n))) for n in need)  # type: ignore[return-value]


def _neg(sh: Shift) -> Shift:
    return (-sh[0], -sh[1], -sh[2])


def _canon(a: int, b: int, sh: Shift) -> tuple[int, int, Shift]:
    """Lo stesso contatto tra molecole scritto in un solo modo: (a, b, s) e (b, a, −s) sono la stessa coppia."""
    return min((a, b, sh), (b, a, _neg(sh)))


def _canon5(a: int, b: int, sh: Shift, fa: int, fb: int) -> Member:
    """Come `_canon`, con i frammenti: (a, b, s, fa, fb) e (b, a, −s, fb, fa) sono lo stesso contatto."""
    return min((a, b, sh, fa, fb), (b, a, _neg(sh), fb, fa))


def iter_neighbour_pairs(s: Structure, cores: Sequence, limits: Limits) -> Iterator[Pair]:
    """Tutte le coppie di frammenti vicini, una per traslazione, senza raccogliere quelle equivalenti. Le molecole
    devono toccarsi (atomi pesanti entro `contact_cutoff`); oltre alla coppia dei frammenti principali (0–0) si
    misurano quelle di frammenti entro `contact_cutoff + FRAGMENT_MARGIN`. In modalità site ogni contatto compare nei
    due versi, altrimenti in uno."""
    frags = as_fragments(cores)
    na, nb, nc = neighbour_shells(s, frags, limits)
    shifts = np.array(list(itertools.product(range(-na, na + 1), range(-nb, nb + 1), range(-nc, nc + 1))), dtype=int)
    shift_vec = shifts @ s.vectors
    extent = [m.extent for m in s.molecules]
    usable = [i for i, m in enumerate(s.molecules) if frags[i] and not m.extended]
    one_way = limits.pair_mode != "site"
    reach = limits.contact_cutoff + FRAGMENT_MARGIN
    for a in usable:
        for b in usable:
            if one_way and b < a:
                continue
            # scarto veloce: baricentri troppo lontani perché le molecole si tocchino
            gap = np.linalg.norm(s.molecules[b].centroid + shift_vec - s.molecules[a].centroid, axis=1)
            near = gap <= extent[a] + extent[b] + limits.contact_cutoff
            fragmented = len(frags[a]) > 1 or len(frags[b]) > 1
            for k in np.flatnonzero(near):
                sh = tuple(int(x) for x in shifts[k])
                if a == b and sh == (0, 0, 0):
                    continue
                if one_way and _canon(a, b, sh) != (a, b, sh):
                    continue
                v = s.shift_vector(sh)
                contact = _heavy_contact(s, a, b, v)
                if contact > limits.contact_cutoff:
                    continue
                for fa, fb in itertools.product(range(len(frags[a])), range(len(frags[b]))):
                    if (fa or fb) and _fragment_distance(s, frags[a][fa], frags[b][fb], a, b, v) > reach:
                        continue
                    yield _measure(s, frags[a][fa], frags[b][fb], a, b, sh, limits, contact, fa, fb, fragmented)
    # coppie di frammenti della stessa molecola (molecole con più sistemi π, p. es. ripiegate su se stesse)
    zero = np.zeros(3)
    for a in usable:
        for fa, fb in itertools.permutations(range(len(frags[a])), 2):
            if one_way and fb < fa:
                continue
            dist = _fragment_distance(s, frags[a][fa], frags[a][fb], a, a, zero)
            if dist <= reach:
                yield _measure(s, frags[a][fa], frags[a][fb], a, a, (0, 0, 0), limits, dist, fa, fb, True)


StackEdge = tuple[int, int, Shift, Pair]


def collect_pairs(s: Structure, cores: Sequence, limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
    """Coppie di frammenti vicini, con quelle equivalenti raccolte in una sola, e, a parte, tutti gli archi `stack`
    tra molecole diverse (nei due versi: servono a ricostruire le colonne)."""
    frags = as_fragments(cores)
    if limits.pair_mode == "site":
        pairs, stack_edges = _collect_site(s, frags, limits)
    else:
        pairs, stack_edges = _collect_canonical(s, frags, limits)
    order = {k: i for i, k in enumerate(KINDS)}
    pairs.sort(key=lambda p: (order[p.kind], p.interplanar, p.contact))
    for i, p in enumerate(pairs):
        p.index = i
    return pairs, stack_edges


def _collect_site(s: Structure, frags: Fragments, limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
    found: dict[str, Pair] = {}
    stack_edges: list[StackEdge] = []
    for p in iter_neighbour_pairs(s, frags, limits):
        if p.kind == "stack" and not p.intramolecular:
            stack_edges.append((p.a, p.b, p.shift, p))
        key = "|".join([p.kind, f"{p.interplanar:.2f}", f"{p.long_slip:.2f}", f"{p.lat_slip:.2f}",
                        f"{p.dihedral:.1f}", f"{p.contact:.2f}"])
        prev = found.get(key)
        if prev is not None:
            prev.count += 1
            prev.members.append(p.key)
        else:
            p.members = [p.key]
            found[key] = p
    return list(found.values()), stack_edges


def fragment_maps(s: Structure, frags: Fragments) -> list[list[list[int | None] | None]]:
    """Per ogni operazione e molecola: in quale frammento della molecola immagine va ciascun frammento (None se non
    corrisponde a un frammento)."""
    sets = [[frozenset(int(i) for i in f.atoms) for f in fs] for fs in frags]
    out = []
    for row_img, row_map in zip(s.molecule_images(), s.atom_maps()):
        row = []
        for m, (img, amap) in enumerate(zip(row_img, row_map)):
            if img is None or amap is None:
                row.append(None)
                continue
            target = {fs: j for j, fs in enumerate(sets[img[0]])}
            row.append([target.get(frozenset(int(amap[i]) for i in fs)) for fs in sets[m]])
        out.append(row)
    return out


def symmetry_key(s: Structure, a: int, b: int, sh: Shift, fa: int = 0, fb: int = 0,
                 fmaps: list[list[list[int | None] | None]] | None = None) -> Member:
    """Rappresentante dell'orbita del contatto sotto il gruppo spaziale: la più piccola delle sue immagini (scritte
    con `_canon5`). Due contatti sono equivalenti per simmetria se e solo se hanno la stessa chiave; le operazioni che
    non mappano molecole e frammenti non si usano, così non si uniscono mai contatti diversi."""
    best = _canon5(a, b, sh, fa, fb)
    for k, row in enumerate(s.molecule_images()):
        ia, ib = row[a], row[b]
        if ia is None or ib is None:
            continue
        fa2, fb2 = fa, fb
        if fmaps is not None:
            ma, mb = fmaps[k][a], fmaps[k][b]
            if ma is None or mb is None or ma[fa] is None or mb[fb] is None:
                continue
            fa2, fb2 = ma[fa], mb[fb]
        (a2, ta), (b2, tb) = ia, ib
        rs = s.rotation(k) @ np.array(sh)
        new = (tb[0] + int(rs[0]) - ta[0], tb[1] + int(rs[1]) - ta[1], tb[2] + int(rs[2]) - ta[2])
        if a2 == b2 and new == (0, 0, 0) and fa2 == fb2:
            continue
        best = min(best, _canon5(a2, b2, new, fa2, fb2))
    return best


def _collect_canonical(s: Structure, frags: Fragments, limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
    classes: dict[tuple, list[Pair]] = {}
    stack_edges: list[StackEdge] = []
    fmaps = fragment_maps(s, frags) if any(len(f) > 1 for f in frags) else None
    for p in iter_neighbour_pairs(s, frags, limits):
        if p.kind == "stack" and not p.intramolecular:
            stack_edges += [(p.a, p.b, p.shift, p), (p.b, p.a, _neg(p.shift), p)]
        classes.setdefault(symmetry_key(s, p.a, p.b, p.shift, p.fa, p.fb, fmaps), []).append(p)
    pairs = []
    for key, members in classes.items():
        members.sort(key=lambda p: p.key)
        rep = next((p for p in members if p.key == key), members[0])
        rep.count = len(members)
        rep.members = [rep.key] + [p.key for p in members if p is not rep]
        pairs.append(rep)
    return pairs, stack_edges


def molecule_key(s: Structure, pair: Pair, limits: Limits) -> tuple:
    """Chiave del contatto tra le due molecole della coppia (senza i frammenti): le coppie di frammenti con la stessa
    chiave appartengono allo stesso contatto tra molecole, a meno della simmetria in modalità canonical."""
    if pair.intramolecular:
        if limits.pair_mode == "site":
            return ("intra", pair.a)
        images = [row[pair.a][0] for row in s.molecule_images() if row[pair.a] is not None]
        return ("intra", min([pair.a, *images]))
    if limits.pair_mode == "site":
        return ("inter", *_canon(pair.a, pair.b, pair.shift))
    return ("inter", *symmetry_key(s, pair.a, pair.b, pair.shift)[:3])


def class_spread(s: Structure, cores: Sequence, limits: Limits, pair: Pair) -> float:
    """Differenza massima di distanza interplanare, slittamenti e contatto tra i membri della classe di `pair` (Å):
    ~0 se l'equivalenza per simmetria è esatta. Serve ai controlli."""
    vals = []
    for a, b, sh, *f in pair.members:
        fa, fb = f if f else (0, 0)
        q = measure_pair(s, cores, a, b, sh, limits, fa=fa, fb=fb)
        if q is not None:
            vals.append((q.interplanar, q.long_slip, q.lat_slip, q.contact))
    if len(vals) < 2:
        return 0.0
    arr = np.array(vals)
    return float(np.max(arr.max(axis=0) - arr.min(axis=0))) if math.isfinite(arr.sum()) else float("nan")


def neighbour_pairs(s: Structure, cores: Sequence, limits: Limits) -> list[Pair]:
    return collect_pairs(s, cores, limits)[0]
