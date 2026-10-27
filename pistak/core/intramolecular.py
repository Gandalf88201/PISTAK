"""Geometria delle molecole: lunghezze di legame, angoli di legame, angoli di torsione e contatti intermolecolari più
corti della somma dei raggi di van der Waals, con le loro incertezze.

Si misurano solo le molecole indipendenti per simmetria: per ogni insieme di molecole della cella legate da
un'operazione del gruppo si prende quella che contiene più atomi dell'unità asimmetrica così come sono nel CIF.
Ogni atomo porta l'operazione di simmetria che lo genera dall'atomo del CIF (`GeometryItem.codes`, "" se è
l'atomo così com'è); `symmetry_codes` le numera come nei CIF di SHELXL (C2#1, con #1 = -x+1,-y,-z). Le misure
uguali per simmetria interna della molecola (C1–C2 e C1#1–C2#1) compaiono una volta sola.

Contatti: coppie di atomi di due molecole diverse (o di due copie della stessa) più vicini della somma dei raggi
di van der Waals (`elements.VDW_RADII`); il primo atomo è sempre di una molecola indipendente.

Incertezze: si propagano le s.u. delle coordinate frazionarie del CIF attraverso il gradiente della misura, tenendo
conto che due atomi copie dello stesso atomo dell'unità asimmetrica si muovono insieme. Come in `uncertainty.py`
le coordinate sono considerate indipendenti tra loro e la cella esatta; se un atomo coinvolto non ha s.u. (idrogeni
in posizione calcolata) l'incertezza è NaN.
"""
from __future__ import annotations

import math
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
from scipy.spatial import cKDTree

from pistak.core.cif import Cif
from pistak.core.elements import vdw_radius
from pistak.core.structure import Structure

KINDS = ("bond", "angle", "torsion", "contact")
_STEP = 1e-4  # Å, per le differenze finite


@dataclass(frozen=True)
class GeometryItem:
    kind: str  # bond | angle | torsion | contact
    #: Indice della molecola in `Structure.molecules` (per i contatti, quella del primo atomo).
    molecule: int
    #: Etichette degli atomi del CIF, nell'ordine della misura (per gli angoli l'atomo centrale è il secondo).
    atoms: tuple[str, ...]
    value: float  # Å o gradi
    su: float  # stessa unità; NaN se non si può stimare
    #: Per ogni atomo l'operazione che lo genera dall'atomo del CIF («-x+1,-y,-z»), "" se è quello del CIF.
    codes: tuple[str, ...] = ()
    #: Solo per i contatti: somma dei raggi di van der Waals (Å).
    vdw: float = math.nan

    @property
    def delta(self) -> float:
        """Per i contatti: distanza meno somma dei raggi di van der Waals (negativa = più corto)."""
        return self.value - self.vdw


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


def symop_text(r: np.ndarray, t: np.ndarray) -> str:
    """Un'operazione di simmetria come nei CIF: «-x+1,-y,z+1/2»."""
    rows = []
    for i in range(3):
        s = ""
        for c, v in zip(np.rint(r[i]).astype(int), "xyz"):
            if c:
                s += ("+" if c > 0 else "-") + ("" if abs(c) == 1 else str(abs(c))) + v
        f = Fraction(float(t[i])).limit_denominator(12)
        if f:
            s += f"{'+' if f > 0 else '-'}{abs(f)}"
        rows.append(s.lstrip("+") or "0")
    return ",".join(rows)


def symmetry_codes(items: list[GeometryItem]) -> tuple[list[tuple[str, ...]], list[tuple[str, str]]]:
    """Le etichette con i codici di simmetria numerati in ordine di comparsa (C2#1) e la legenda [("#1", op), …]."""
    num: dict[str, str] = {}
    labels = []
    for g in items:
        row = []
        for lab, op in zip(g.atoms, g.codes or [""] * len(g.atoms)):
            if op:
                num.setdefault(op, f"#{len(num) + 1}")
                lab += num[op]
            row.append(lab)
        labels.append(tuple(row))
    return labels, [(n, op) for op, n in num.items()]


class _Codes:
    """L'operazione che porta ogni atomo del CIF nella posizione in cui sta in una molecola (o accanto a essa)."""

    def __init__(self, s: Structure, cif: Cif):
        self.s, self.cif = s, cif
        self.ident = _identity_index(s)
        self.inv = np.linalg.inv(s.vectors)

    def ops(self, cell_atoms: np.ndarray, pos: np.ndarray) -> list[tuple[int, np.ndarray]]:
        """(operazione, traslazione di reticolo) per atomi della cella portati nelle posizioni cartesiane `pos`."""
        out = []
        for ci, p in zip(cell_atoms, pos):
            ia, k = (int(x) for x in self.s.cell_origin[ci])
            base = self.s.symops[k].apply(self.cif.atoms[ia].frac)
            out.append((k, np.rint(p @ self.inv - base).astype(int)))
        return out

    def reference(self, ops: list[tuple[int, np.ndarray]]) -> np.ndarray:
        """La traslazione più comune tra gli atomi non trasformati: la si toglie a tutti, così la molecola indipendente
        ha gli atomi del CIF senza codice anche se la cella la mette in un'altra cella."""
        own = Counter(tuple(n) for k, n in ops if k == self.ident)
        return np.array(own.most_common(1)[0][0]) if own else np.zeros(3, dtype=int)

    def text(self, op: tuple[int, np.ndarray], ref: np.ndarray) -> str:
        k, n = op
        n = n - ref
        if k == self.ident and not n.any():
            return ""
        so = self.s.symops[k]
        return symop_text(so.r, so.t + n)


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
    """Legami, angoli, torsioni e contatti di van der Waals delle molecole indicate (per impostazione quelle
    indipendenti per simmetria)."""
    prop = _Propagator(s, cif)
    sym = _Codes(s, cif) if len(s.cell_origin) else None
    mols = unique_molecules(s) if molecules is None else molecules
    out: list[GeometryItem] = []
    for mi in mols:
        m = s.molecules[mi]
        if sym is not None:
            ops = sym.ops(m.cell_index, m.pos)
            ref = sym.reference(ops)
            codes = [sym.text(op, ref) for op in ops]
        else:
            codes = [""] * len(m)
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
        # misure uguali per simmetria interna (stesse etichette): una sola, quella con meno atomi trasformati
        best: dict[tuple[str, ...], tuple[int, ...]] = {}
        for kind in (k for k in kinds if k in MEASURES):
            for idx in groups[kind]:
                if not include_h and any(m.el[i] == "H" for i in idx):
                    continue
                base = tuple(m.labels[i] for i in idx)
                key = (kind, *min(base, base[::-1]))
                n_coded = sum(bool(codes[i]) for i in idx)
                if key not in best or n_coded < sum(bool(codes[i]) for i in best[key]):
                    best[key] = idx
        for (kind, *_), idx in best.items():
            if [m.labels[i] for i in idx[::-1]] < [m.labels[i] for i in idx] and kind != "torsion":
                idx = idx[::-1]  # ordine di lettura stabile (la torsione no: ha un verso)
            fn = MEASURES[kind]
            p = m.pos[list(idx)]
            out.append(GeometryItem(kind, mi, tuple(m.labels[i] for i in idx), fn(p),
                                    prop(fn, p, m.cell_index[list(idx)]), tuple(codes[i] for i in idx)))
    if "contact" in kinds:
        out += vdw_contacts(s, cif, mols, include_h, prop, sym)
    return out


def vdw_contacts(s: Structure, cif: Cif, molecules: list[int] | None = None, include_h: bool = True,
                 prop: _Propagator | None = None, sym: _Codes | None = None) -> list[GeometryItem]:
    """Contatti tra atomi di molecole diverse più corti della somma dei raggi di van der Waals, dal più corto in
    proporzione; ogni contatto equivalente per simmetria compare una volta sola."""
    prop = prop or _Propagator(s, cif)
    if sym is None and len(s.cell_origin):
        sym = _Codes(s, cif)
    cart = s.cell_frac @ s.vectors
    radii = np.array([vdw_radius(e) for e in s.cell_el])
    cutoff = 2 * float(radii.max()) if len(radii) else 0.0
    inv = np.linalg.inv(s.vectors)
    lengths = np.linalg.norm(s.vectors, axis=1)
    out: list[GeometryItem] = []
    seen: set[tuple] = set()
    for mi in unique_molecules(s) if molecules is None else molecules:
        m = s.molecules[mi]
        if m.extended or not len(m):
            continue
        own = {(int(ci), tuple(int(x) for x in sh)) for ci, sh in zip(m.cell_index, m.cell_shift)}
        frac = m.pos @ inv
        reach = np.ceil(cutoff / lengths).astype(int)
        lo, hi = np.floor(frac.min(axis=0)).astype(int) - reach, np.floor(frac.max(axis=0)).astype(int) + reach
        shifts = np.array([(x, y, z) for x in range(lo[0], hi[0] + 1) for y in range(lo[1], hi[1] + 1)
                           for z in range(lo[2], hi[2] + 1)], dtype=int)
        pts = (cart[None, :, :] + (shifts @ s.vectors)[:, None, :]).reshape(-1, 3)
        tree = cKDTree(pts)
        n = len(cart)
        ref = None
        if sym is not None:
            ref = sym.reference(sym.ops(m.cell_index, m.pos))
        for a, near in enumerate(tree.query_ball_point(m.pos, cutoff)):
            for q in near:
                j, sh = q % n, shifts[q // n]
                if (j, tuple(int(x) for x in sh)) in own:
                    continue
                if not include_h and "H" in (m.el[a], s.cell_el[j]):
                    continue
                p = np.array([m.pos[a], pts[q]])
                d = distance(p)
                lim = float(vdw_radius(m.el[a]) + radii[j])
                if d >= lim or d < 1e-3:
                    continue
                key = (*sorted((m.labels[a], s.cell_labels[j])), round(d, 3))
                if key in seen:
                    continue
                seen.add(key)
                cells = np.array([m.cell_index[a], j])
                codes = ("", "")
                if sym is not None:
                    ops = sym.ops(cells, p)
                    codes = (sym.text(ops[0], ref), sym.text(ops[1], ref))
                out.append(GeometryItem("contact", mi, (m.labels[a], s.cell_labels[j]), d, prop(distance, p, cells),
                                        codes, lim))
    out.sort(key=lambda g: (g.delta, g.atoms))
    return out
