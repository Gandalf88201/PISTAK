"""Colonne di π-stack: molecole collegate da coppie `stack` nel cristallo periodico.

Le coppie `stack` sono gli archi di un grafo sulle molecole della cella; percorrendolo si assegna a ogni molecola
la cella in cui si trova. Se si torna a una molecola già vista in una cella diversa, la differenza è una
traslazione del reticolo lungo la quale la colonna continua. Il numero di traslazioni indipendenti (rango) dice
che cosa si è trovato:

    rango 0   aggregato finito (dimero, trimero…)
    rango 1   colonna infinita lungo un asse, il caso classico dell'impilamento
    rango 2   strato di molecole impilate in due direzioni
    rango 3   rete tridimensionale di contatti di stack

I nodi del grafo sono le molecole discrete e, per le reti estese (polimeri di coordinazione, catene covalenti), i
singoli frammenti π (revisione 5, N1). L'unità ripetuta di una rete è un taglio arbitrario: se un anello passa
all'unità vicina, gli archi che lo toccano si spostano di una traslazione e gli altri archi dell'unità no, quindi le
traslazioni dei cicli, e il rango, dipenderebbero dal taglio. Un frammento è invece un oggetto finito: spostarlo di
una traslazione sposta tutti i suoi archi insieme e non cambia le traslazioni dei cicli. Il rango si riferisce agli
archi di stack soltanto, indipendentemente dai legami covalenti o di coordinazione tra i nodi.
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


#: Nodo del grafo di stack: (molecola, frammento) per un frammento di una rete estesa, (molecola, -1) per una molecola.
Node = tuple[int, int]


@dataclass
class StackColumn:
    molecules: list[int]  # indici in Structure.molecules, nell'ordine in cui la visita li incontra, senza ripetizioni
    rank: int
    #: Cella in cui sta ciascuna molecola rispetto alla prima (percorso ad albero del grafo); per una rete estesa la
    #: cella del primo dei suoi frammenti incontrato.
    positions: dict[int, Shift]
    #: Traslazioni del reticolo (in celle) lungo cui la colonna continua, indipendenti e ridotte.
    translations: list[Shift] = field(default_factory=list)
    axis: np.ndarray | None = None  # rango 1: vettore di traslazione in Å (passo della colonna)
    edges: list[tuple] = field(default_factory=list)  # archi `pairs.StackEdge`: (a, b, shift, coppia, fa, fb)
    #: I nodi del grafo (`Node`) e la cella di ciascuno.
    nodes: list[Node] = field(default_factory=list)
    node_positions: dict[Node, Shift] = field(default_factory=dict)

    @property
    def kind_key(self) -> str:
        """Identificatore stabile del tipo: finite | column | layer | network."""
        return RANK_KEYS[self.rank]

    @property
    def repeat(self) -> float | None:
        return float(np.linalg.norm(self.axis)) if self.axis is not None else None

    def mean(self, attr: str) -> float:
        return float(np.mean([getattr(e[3], attr) for e in self.edges])) if self.edges else float("nan")

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


def stack_node(s: Structure, m: int, f: int) -> Node:
    """Il nodo del grafo di stack per il frammento f della molecola m: la molecola, o il frammento se è una rete."""
    return (m, f) if s.molecules[m].extended else (m, -1)


def find_stacks(s: Structure, edges: list[tuple]) -> list[StackColumn]:
    """Raggruppa le molecole (e i frammenti delle reti estese) in colonne a partire dagli archi `stack`
    (a, b, traslazione di b, coppia, fa, fb)."""
    adj: dict[Node, list[tuple[Node, Shift, Pair, tuple]]] = {}
    for e in edges:
        a, b, sh, p = e[:4]
        fa, fb = (e[4], e[5]) if len(e) >= 6 else (p.fa, p.fb)
        u, w = stack_node(s, a, fa), stack_node(s, b, fb)
        if u == w and not any(sh):
            continue  # due frammenti della stessa molecola discreta: nessun arco
        adj.setdefault(u, []).append((w, sh, p, e))
    columns: list[StackColumn] = []
    seen: set[Node] = set()
    for start in sorted(adj):
        if start in seen:
            continue
        pos: dict[Node, Shift] = {start: (0, 0, 0)}
        order = [start]
        trans: list[tuple[int, int, int]] = []
        comp_edges: list[tuple] = []
        k = 0
        while k < len(order):
            u = order[k]
            k += 1
            for w, sh, p, e in adj.get(u, []):
                comp_edges.append(e)
                want = (pos[u][0] + sh[0], pos[u][1] + sh[1], pos[u][2] + sh[2])
                if w not in pos:
                    pos[w] = want
                    order.append(w)
                elif pos[w] != want:
                    t = (want[0] - pos[w][0], want[1] - pos[w][1], want[2] - pos[w][2])
                    trans.append(t)
        seen.update(order)
        rank = int(np.linalg.matrix_rank(np.array(trans))) if trans else 0
        reduced = _reduce_translations(trans, rank)
        mols: list[int] = []
        mpos: dict[int, Shift] = {}
        for node in order:
            if node[0] not in mpos:
                mols.append(node[0])
                mpos[node[0]] = pos[node]
        col = StackColumn(mols, rank, mpos, reduced, edges=comp_edges, nodes=order, node_positions=pos)
        if rank == 1:
            col.axis = s.shift_vector(reduced[0])
        columns.append(col)
    columns.sort(key=lambda c: (-c.rank if c.rank else 9, -len(c.nodes)))
    return columns
