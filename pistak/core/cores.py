"""Anelli e nuclei π di una molecola.

Per ogni molecola si cerca il nucleo π: il sistema di anelli condensati planari più grande
(in mancanza, tutta la molecola senza idrogeni, se è planare). Del nucleo si calcolano il piano medio
(minimi quadrati), il baricentro e l'asse lungo (direzione di massima estensione nel piano).

Porting di `piCore`/`smallRings` di `src/lib/crystal/pistack.ts`.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

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
        p = self.plane
        return Core(Plane(p.centroid + v, p.normal, p.long, p.rms), self.atoms, self.kind)


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
    """Nucleo π di una molecola, o None se non ha parti planari adatte."""
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
        best = candidates[0]
        return Core(fit_plane(m.pos[best]), np.array(best), "rings")
    allh = np.flatnonzero(m.heavy)
    if len(allh) < 4:
        return None
    pl = fit_plane(m.pos[allh])
    return Core(pl, allh, "molecule") if pl.rms < limits.molecule_rms else None
