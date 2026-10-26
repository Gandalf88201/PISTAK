"""Dal contenuto di un CIF alle molecole del cristallo: si applicano le operazioni di simmetria,
si eliminano i doppioni sulla cella, si riconoscono i legami (anche attraverso i bordi della cella)
e si ricompongono le molecole intere in coordinate cartesiane (Å).

Porting di `src/lib/crystal/structure.ts`; i cicli annidati del sito sono qui vettorizzati con numpy/scipy.
Convenzione: una matrice di cella ha i vettori a, b, c come RIGHE, quindi cart = frac @ M.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from pistak.core.cif import Cell, Cif
from pistak.core.config import Limits
from pistak.core.elements import covalent_radius


class StructureError(Exception):
    def __init__(self, code: str, detail: str = ""):
        super().__init__(f"{code}: {detail}" if detail else code)
        self.code = code
        self.detail = detail


def cell_vectors(c: Cell) -> np.ndarray:
    """Vettori di cella in coordinate cartesiane, con a lungo x e b nel piano xy (convenzione abituale)."""
    ca, cb, cg = np.cos(np.radians([c.alpha, c.beta, c.gamma]))
    sg = np.sin(np.radians(c.gamma))
    cx = c.c * cb
    cy = c.c * (ca - cb * cg) / sg
    cz = np.sqrt(max(0.0, c.c ** 2 - cx ** 2 - cy ** 2))
    return np.array([[c.a, 0.0, 0.0], [c.b * cg, c.b * sg, 0.0], [cx, cy, cz]])


def to_cart(m: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Coordinate frazionarie → cartesiane (accetta un punto o un array N×3)."""
    return np.asarray(f) @ m


def to_frac(m: np.ndarray, p: np.ndarray) -> np.ndarray:
    return np.asarray(p) @ np.linalg.inv(m)


@dataclass
class Molecule:
    el: list[str]
    labels: list[str]
    pos: np.ndarray  # N×3, Å
    bonds: np.ndarray  # M×2, indici in `el`
    centroid: np.ndarray
    #: Vero se la rete di legami non si chiude: polimero o solido esteso, non una molecola.
    extended: bool
    #: Per ogni atomo: indice dell'atomo nella cella (`Structure.cell_el`) e cella di appartenenza.
    cell_index: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    cell_shift: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), dtype=int))

    def __len__(self) -> int:
        return len(self.el)

    @property
    def heavy(self) -> np.ndarray:
        """Maschera booleana degli atomi che non sono idrogeno."""
        return np.array([e != "H" for e in self.el])

    @property
    def extent(self) -> float:
        """Distanza massima di un atomo dal baricentro: serve a scartare in fretta le molecole lontane."""
        return float(np.linalg.norm(self.pos - self.centroid, axis=1).max(initial=0.0))

    def formula(self) -> list[tuple[str, int]]:
        """Formula in ordine di Hill: C, H, poi gli altri in ordine alfabetico."""
        n: dict[str, int] = {}
        for e in self.el:
            n[e] = n.get(e, 0) + 1
        order = sorted(n, key=lambda e: (e != "C", e != "H", e))
        return [(e, n[e]) for e in order]


@dataclass
class Structure:
    vectors: np.ndarray
    cell_el: list[str]
    cell_labels: list[str]
    cell_frac: np.ndarray
    molecules: list[Molecule]

    @property
    def n_cell_atoms(self) -> int:
        return len(self.cell_el)

    def shift_vector(self, shift) -> np.ndarray:
        return to_cart(self.vectors, np.asarray(shift, dtype=float))


def _wrap(x: np.ndarray) -> np.ndarray:
    return x - np.floor(x + 1e-9)


def fill_cell(cif: Cif, vectors: np.ndarray, limits: Limits) -> tuple[list[str], list[str], np.ndarray]:
    """Tutti gli atomi della cella: ogni operazione di simmetria applicata a ogni atomo, senza doppioni."""
    els: list[str] = []
    labels: list[str] = []
    buf = np.zeros((16, 3))  # posizioni frazionarie accettate finora (cresce al bisogno)
    by_el: dict[str, list[int]] = {}
    for at in cif.atoms:
        for op in cif.symops:
            f = _wrap(op.apply(at.frac))
            idx = by_el.get(at.el)
            if idx:
                d = f - buf[idx]
                d -= np.round(d)
                if (np.linalg.norm(d @ vectors, axis=1) < limits.same_site).any():
                    continue
            if len(els) == len(buf):
                buf = np.vstack([buf, np.zeros_like(buf)])
            buf[len(els)] = f
            by_el.setdefault(at.el, []).append(len(els))
            els.append(at.el)
            labels.append(at.label)
            if len(els) > limits.max_cell_atoms:
                raise StructureError("tooLarge", str(limits.max_cell_atoms))
    return els, labels, buf[:len(els)].copy()


_SHIFTS = np.array([(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)], dtype=int)


def cell_bonds(els: list[str], pos: np.ndarray, vectors: np.ndarray, limits: Limits) -> list[list[tuple[int, tuple[int, int, int]]]]:
    """Legami tra atomi della cella, anche verso le celle vicine: per ogni atomo, lista di (vicino, traslazione in celle)."""
    n = len(els)
    r = np.array([covalent_radius(e) for e in els])
    is_h = np.array([e == "H" for e in els])
    max_cut = 2 * r.max() + limits.bond_tolerance
    tree = cKDTree(pos)
    rows: list[np.ndarray] = []
    for s in _SHIFTS:
        v = s @ vectors
        other = cKDTree(pos + v)
        # coppie (i, j) con |pos[j] + v − pos[i]| < max_cut
        sm = tree.sparse_distance_matrix(other, max_cut, output_type="ndarray")
        if not len(sm):
            continue
        i, j, d = sm["i"], sm["j"], sm["v"]
        keep = (d < r[i] + r[j] + limits.bond_tolerance) & (d * d >= 0.16)
        keep &= ~(is_h[i] & is_h[j])  # due idrogeni non sono mai legati tra loro in un cristallo molecolare
        if not s.any():
            keep &= i != j
        i, j = i[keep], j[keep]
        rows.append(np.column_stack([i, j, np.broadcast_to(s, (len(i), 3))]))
    edges: list[list[tuple[int, tuple[int, int, int]]]] = [[] for _ in range(n)]
    if rows:
        allrows = np.concatenate(rows)
        order = np.lexsort((allrows[:, 4], allrows[:, 3], allrows[:, 2], allrows[:, 1], allrows[:, 0]))
        for i, j, sx, sy, sz in allrows[order]:
            edges[i].append((int(j), (int(sx), int(sy), int(sz))))
    return edges


def build_structure(cif: Cif, limits: Limits | None = None) -> Structure:
    """Ricompone le molecole intere partendo dagli atomi della cella e dai loro legami."""
    limits = limits or Limits()
    vectors = cell_vectors(cif.cell)
    els, labels, frac = fill_cell(cif, vectors, limits)
    cart = to_cart(vectors, frac)
    edges = cell_bonds(els, cart, vectors, limits)
    seen: dict[int, tuple[int, int, int]] = {}  # atomo della cella → cella in cui è stato messo
    molecules: list[Molecule] = []
    inv = np.linalg.inv(vectors)
    for start in range(len(els)):
        if start in seen:
            continue
        members: list[tuple[int, tuple[int, int, int]]] = [(start, (0, 0, 0))]
        index = {start: 0}
        bonds: list[tuple[int, int]] = []
        extended = False
        seen[start] = (0, 0, 0)
        k = 0
        while k < len(members):
            i, shift = members[k]
            for j, es in edges[i]:
                s = (shift[0] + es[0], shift[1] + es[1], shift[2] + es[2])
                prev = seen.get(j)
                if prev is None:
                    seen[j] = s
                    index[j] = len(members)
                    members.append((j, s))
                    bonds.append((k, len(members) - 1))
                elif prev != s:
                    extended = True  # lo stesso atomo raggiunto in due celle diverse: la rete continua all'infinito
                elif index[j] > k:
                    bonds.append((k, index[j]))
            k += 1
        ci = np.array([m[0] for m in members])
        cs = np.array([m[1] for m in members], dtype=int)
        pos = to_cart(vectors, frac[ci] + cs)
        # la molecola si sposta di celle intere finché il suo baricentro non cade nella cella di partenza
        c = pos.mean(axis=0)
        fc = np.floor(c @ inv + 1e-6)
        back = -(fc @ vectors)
        pos = pos + back
        molecules.append(Molecule(
            el=[els[i] for i in ci], labels=[labels[i] for i in ci], pos=pos,
            bonds=_dedupe(bonds), centroid=c + back, extended=extended,
            cell_index=ci, cell_shift=cs - fc.astype(int),
        ))
    return Structure(vectors, els, labels, frac, molecules)


def _dedupe(b: list[tuple[int, int]]) -> np.ndarray:
    if not b:
        return np.zeros((0, 2), dtype=int)
    arr = np.sort(np.array(b, dtype=int), axis=1)
    return np.unique(arr, axis=0)
