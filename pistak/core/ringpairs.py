"""Parametri anello–anello (stile Janiak) per una coppia di molecole.

Il nucleo π di una coppia è un solo piano medio; qui si guardano invece tutti gli anelli planari delle due
molecole, uno contro l'altro: distanza tra i centroidi Cg···Cg, angolo α tra i piani degli anelli,
angoli β e γ tra il vettore Cg···Cg e le normali, distanze perpendicolari e slittamento.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.distance import cdist

from pistak.core.analysis import PairContext
from pistak.core.cores import planar_rings
from pistak.core.geometry import acute_angle_deg


@dataclass(frozen=True)
class RingPair:
    ring_a: int
    ring_b: int
    size_a: int
    size_b: int
    atoms_a: tuple[int, ...]  # indici nella molecola A
    atoms_b: tuple[int, ...]  # indici nella molecola B
    cg_distance: float  # Å
    alpha: float  # ° tra i piani
    beta: float  # ° tra Cg···Cg e la normale dell'anello A
    gamma: float  # ° tra Cg···Cg e la normale dell'anello B
    perp_a: float  # Å: distanza del centroide B dal piano A
    perp_b: float  # Å: distanza del centroide A dal piano B
    slippage: float  # Å: spostamento laterale dei centroidi
    closest: float  # Å: contatto più corto tra atomi dei due anelli


HEADERS = {
    "ring_a": "Anello A", "ring_b": "Anello B", "size_a": "Atomi A", "size_b": "Atomi B", "cg_distance": "Cg···Cg (Å)",
    "alpha": "α (°)", "beta": "β (°)", "gamma": "γ (°)", "perp_a": "⟂ A (Å)", "perp_b": "⟂ B (Å)",
    "slippage": "Slittamento (Å)", "closest": "Contatto min (Å)",
}


def ring_pairs(ctx: PairContext, max_cg: float | None = None) -> list[RingPair]:
    """Tutte le coppie di anelli planari (uno per molecola) con Cg···Cg entro `max_cg`, dalla più vicina."""
    max_cg = max_cg if max_cg is not None else ctx.limits.ring_pair_max_distance
    ra = planar_rings(ctx.A, ctx.limits)
    rb = planar_rings(ctx.B, ctx.limits)
    out: list[RingPair] = []
    for ia, (atoms_a, pa) in enumerate(ra):
        for ib, (atoms_b, pb) in enumerate(rb):
            d = (pb.centroid + ctx.shift_vector) - pa.centroid
            cg = float(np.linalg.norm(d))
            if cg > max_cg or cg == 0:
                continue
            perp_a = abs(float(np.dot(d, pa.normal)))
            perp_b = abs(float(np.dot(d, pb.normal)))
            closest = float(cdist(ctx.pos_a[atoms_a], ctx.pos_b[atoms_b]).min())
            out.append(RingPair(
                ia, ib, len(atoms_a), len(atoms_b), tuple(atoms_a), tuple(atoms_b), cg,
                acute_angle_deg(pa.normal, pb.normal), acute_angle_deg(d, pa.normal), acute_angle_deg(d, pb.normal),
                perp_a, perp_b, float(np.sqrt(max(0.0, cg ** 2 - perp_a ** 2))), closest,
            ))
    return sorted(out, key=lambda r: r.cg_distance)
