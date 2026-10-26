"""Anelli e nuclei π di una molecola.

Per ogni molecola si cerca il nucleo π: il sistema di anelli condensati planari più grande
(in mancanza, tutta la molecola senza idrogeni, se è planare). Del nucleo si calcolano il piano medio
(minimi quadrati), il baricentro e l'asse lungo (direzione di massima estensione nel piano).

Porting di `piCore`/`smallRings` di `src/lib/crystal/pistack.ts`.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, replace

import numpy as np

from pistak.core.config import Limits
from pistak.core.geometry import Plane, fit_plane
from pistak.core.structure import Molecule


@dataclass(frozen=True)
class Core:
    """Nucleo π: piano medio e atomi (indici nella molecola) che lo compongono."""

    plane: Plane
    atoms: np.ndarray
    kind: str  # "rings" | "molecule"
    #: Da dove viene l'asse lungo: "pca" (direzione di massima estensione) | "atoms" (definito dall'utente).
    axis_source: str = "pca"

    @property
    def centroid(self) -> np.ndarray:
        return self.plane.centroid

    @property
    def normal(self) -> np.ndarray:
        return self.plane.normal

    @property
    def long(self) -> np.ndarray:
        return self.plane.long

    @property
    def rms(self) -> float:
        return self.plane.rms

    def moved(self, v: np.ndarray) -> Core:
        """Lo stesso nucleo traslato di `v` (la molecola vicina in un'altra cella)."""
        return replace(self, plane=replace(self.plane, centroid=self.plane.centroid + v))

    def with_long(self, direction: np.ndarray) -> Core | None:
        """Lo stesso nucleo con l'asse lungo lungo `direction` (proiettata sul piano); None se è quasi normale al piano."""
        n = self.plane.normal
        v = np.asarray(direction, dtype=float) - n * np.dot(direction, n)
        norm = np.linalg.norm(v)
        if norm < 1e-3 * max(1.0, float(np.linalg.norm(direction))):
            return None
        return replace(self, plane=replace(self.plane, long=v / norm), axis_source="atoms")


def small_rings(n: int, bonds: np.ndarray, allowed, max_size: int = 8) -> list[list[int]]:
    """Anelli fino a `max_size` atomi: per ogni legame, il ciclo più corto che lo contiene."""
    adj: list[list[int]] = [[] for _ in range(n)]
    ok = [(int(p), int(q)) for p, q in bonds if allowed(p) and allowed(q)]
    for p, q in ok:
        adj[p].append(q)
        adj[q].append(p)
    rings: dict[tuple[int, ...], list[int]] = {}
    for p, q in ok:
        # cammino più corto da p a q senza usare il legame p–q
        prev = {p: -1}
        queue = deque([p])
        while queue and q not in prev:
            u = queue.popleft()
            for w in adj[u]:
                if w in prev or (u == p and w == q):
                    continue
                prev[w] = u
                queue.append(w)
        if q not in prev:
            continue
        ring = []
        u = q
        while u != -1:
            ring.append(u)
            u = prev[u]
        if len(ring) < 3 or len(ring) > max_size:
            continue
        rings[tuple(sorted(ring))] = ring
    return list(rings.values())


def planar_rings(m: Molecule, limits: Limits) -> list[tuple[list[int], Plane]]:
    """Gli anelli planari della molecola (senza idrogeni), ciascuno con il suo piano medio."""
    heavy = m.heavy
    out = []
    for r in small_rings(len(m), m.bonds, lambda i: heavy[i], limits.max_ring_size):
        pl = fit_plane(m.pos[r])
        if pl.rms < limits.ring_rms:
            out.append((r, pl))
    return out


def pi_core(m: Molecule, limits: Limits) -> Core | None:
    """Nucleo π principale di una molecola, o None se non ha parti planari adatte."""
    frags = pi_fragments(m, limits.with_values(fragments="largest"))
    return frags[0] if frags else None


def pi_fragments(m: Molecule, limits: Limits) -> list[Core]:
    """I frammenti π della molecola: i sistemi di anelli planari condensati (che condividono almeno un legame), dal
    più grande. Con `limits.fragments == "largest"` solo il primo, il nucleo π del sito; con "all" anche gli altri
    con almeno `min_fragment_atoms` atomi. Senza anelli adatti: tutta la molecola senza idrogeni, se è planare."""
    # anelli condensati: condividono almeno un legame (due atomi)
    groups: list[set[int]] = []
    for r, _ in planar_rings(m, limits):
        touching = [g for g in groups if sum(1 for i in r if i in g) >= 2]
        merged = set(r)
        for g in touching:
            merged |= g
            groups.remove(g)
        groups.append(merged)
    candidates = [sorted(g) for g in groups]
    candidates = [g for g in candidates if fit_plane(m.pos[g]).rms < limits.molecule_rms]
    candidates.sort(key=len, reverse=True)  # ordinamento stabile come nel sito
    if candidates:
        keep = candidates[:1] if limits.fragments == "largest" else \
            candidates[:1] + [g for g in candidates[1:] if len(g) >= limits.min_fragment_atoms]
        return [Core(fit_plane(m.pos[g]), np.array(g), "rings") for g in keep]
    allh = np.flatnonzero(m.heavy)
    if len(allh) < 4:
        return []
    pl = fit_plane(m.pos[allh])
    return [Core(pl, allh, "molecule")] if pl.rms < limits.molecule_rms else []


def axis_atoms(m: Molecule, spec: str) -> tuple[int, int] | None:
    """Gli atomi che definiscono l'asse lungo secondo `spec` ("C1,C9": etichette del CIF), in questa molecola: tra
    tutti gli atomi con quelle etichette (una molecola su un elemento di simmetria ne ha più copie, "C1,C1" vale) la
    coppia più lontana. None se la molecola non li contiene."""
    parts = [p.strip() for p in spec.replace(";", ",").split(",") if p.strip()]
    if len(parts) != 2:
        return None
    ia = [i for i, lab in enumerate(m.labels) if lab == parts[0]]
    ib = [i for i, lab in enumerate(m.labels) if lab == parts[1]]
    pairs = [(i, j) for i in ia for j in ib if i != j]
    if not pairs:
        return None
    return max(pairs, key=lambda ij: (float(np.linalg.norm(m.pos[ij[1]] - m.pos[ij[0]])), -ij[0], -ij[1]))


def apply_axis(m: Molecule, frags: list[Core], spec: str) -> list[Core]:
    """Asse lungo definito dall'utente: si applica ai frammenti che contengono entrambi gli atomi, o, se nessuno li
    contiene, al frammento principale. Il verso va dal primo atomo al secondo."""
    ij = axis_atoms(m, spec) if spec else None
    if ij is None or not frags:
        return frags
    i, j = ij
    direction = m.pos[j] - m.pos[i]
    targets = [k for k, f in enumerate(frags) if i in f.atoms and j in f.atoms] or [0]
    out = list(frags)
    for k in targets:
        out[k] = frags[k].with_long(direction) or frags[k]
    return out
