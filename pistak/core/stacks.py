"""Colonne di π-stack: molecole collegate da coppie `stack` nel cristallo periodico.

Le coppie `stack` sono gli archi di un grafo sulle molecole della cella; percorrendolo si assegna a ogni molecola
la cella in cui si trova. Se si torna a una molecola già vista in una cella diversa, la differenza è una
traslazione del reticolo lungo la quale la colonna continua. Il numero di traslazioni indipendenti (rango) dice
che cosa si è trovato:

    rango 0   aggregato finito (dimero, trimero…)
    rango 1   colonna infinita lungo un asse, il caso classico dell'impilamento
    rango 2   strato di molecole impilate in due direzioni
    rango 3   rete tridimensionale di contatti di stack
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import reduce
from math import gcd

import numpy as np

from pistak.core.geometry import angle_deg
from pistak.core.pairs import Pair, Shift
from pistak.core.structure import Structure

RANK_NAMES = {0: "aggregato finito", 1: "colonna", 2: "strato", 3: "rete 3D"}  # testi italiani: pistak.labels.rank_name
RANK_KEYS = {0: "finite", 1: "column", 2: "layer", 3: "network"}  # identificatori stabili (non tradotti)


@dataclass
class StackColumn:
    molecules: list[int]  # indici in Structure.molecules
    rank: int
    #: Cella in cui sta ciascuna molecola rispetto alla prima (percorso ad albero del grafo).
    positions: dict[int, Shift]
    #: Traslazioni del reticolo (in celle) lungo cui la colonna continua, indipendenti e ridotte.
    translations: list[Shift] = field(default_factory=list)
    axis: np.ndarray | None = None  # rango 1: vettore di traslazione in Å (passo della colonna)
    edges: list[tuple[int, int, Shift, Pair]] = field(default_factory=list)

    @property
    def kind_key(self) -> str:
        """Identificatore stabile del tipo: finite | column | layer | network."""
        return RANK_KEYS[self.rank]

    @property
    def repeat(self) -> float | None:
        return float(np.linalg.norm(self.axis)) if self.axis is not None else None

    def mean(self, attr: str) -> float:
        return float(np.mean([getattr(p, attr) for *_, p in self.edges])) if self.edges else float("nan")

    def axis_normal_angle(self, normal: np.ndarray) -> float:
        """Angolo (°) tra l'asse della colonna e la normale ai piani: 0° = molecole sovrapposte, sopra 0° = slittate."""
        if self.axis is None:
            return float("nan")
        a = angle_deg(self.axis, normal)
        return a if a <= 90 else 180 - a


def _primitive(t: tuple[int, int, int]) -> tuple[int, int, int]:
    g = reduce(gcd, (abs(x) for x in t)) or 1
    return (t[0] // g, t[1] // g, t[2] // g)


def _reduce_translations(ts: list[tuple[int, int, int]], rank: int) -> list[Shift]:
    """Per il rango 1 una sola traslazione primitiva (con il passo giusto); altrimenti le più corte indipendenti."""
    if not ts:
        return []
    if rank == 1:
        d = _primitive(ts[0])
        k = reduce(gcd, (abs(int(np.dot(t, d) / np.dot(d, d))) for t in ts))
        return [(d[0] * k, d[1] * k, d[2] * k)]
    chosen: list[tuple[int, int, int]] = []
    for t in sorted(set(ts), key=lambda t: sum(x * x for x in t)):
        if np.linalg.matrix_rank(np.array(chosen + [t])) > len(chosen):
            chosen.append(t)
    return chosen


def find_stacks(s: Structure, edges: list[tuple[int, int, Shift, Pair]]) -> list[StackColumn]:
    """Raggruppa le molecole in colonne a partire dagli archi `stack` (a, b, traslazione di b, coppia)."""
    adj: dict[int, list[tuple[int, Shift, Pair]]] = {}
    for a, b, sh, p in edges:
        adj.setdefault(a, []).append((b, sh, p))
    columns: list[StackColumn] = []
    seen: set[int] = set()
    for start in sorted(adj):
        if start in seen:
            continue
        pos: dict[int, Shift] = {start: (0, 0, 0)}
        order = [start]
        trans: list[tuple[int, int, int]] = []
        comp_edges: list[tuple[int, int, Shift, Pair]] = []
        k = 0
        while k < len(order):
            a = order[k]
            k += 1
            for b, sh, p in adj.get(a, []):
                comp_edges.append((a, b, sh, p))
                want = (pos[a][0] + sh[0], pos[a][1] + sh[1], pos[a][2] + sh[2])
                if b not in pos:
                    pos[b] = want
                    order.append(b)
                elif pos[b] != want:
                    t = (want[0] - pos[b][0], want[1] - pos[b][1], want[2] - pos[b][2])
                    trans.append(t)
        seen.update(order)
        rank = int(np.linalg.matrix_rank(np.array(trans))) if trans else 0
        reduced = _reduce_translations(trans, rank)
        col = StackColumn(order, rank, pos, reduced, edges=comp_edges)
        if rank == 1:
            col.axis = s.shift_vector(reduced[0])
        columns.append(col)
    columns.sort(key=lambda c: (-c.rank if c.rank else 9, -len(c.molecules)))
    return columns
