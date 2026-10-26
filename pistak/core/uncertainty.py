"""Incertezze dei parametri, propagate con un Monte Carlo dalle incertezze standard (s.u.) del CIF.

Ogni campione sposta le coordinate frazionarie degli atomi dell'unità asimmetrica di quantità casuali gaussiane con
deviazione standard pari alla loro s.u. (e, se richiesto, i parametri di cella); lo spostamento passa a tutti gli
atomi della cella attraverso la parte rotazionale dell'operazione che li genera, quindi le copie per simmetria si
muovono in modo coerente. Si ricalcolano piani e assi dei frammenti π (con lo stesso verso di quelli nominali) e si
rimisura ogni coppia, con la stessa topologia dell'analisi (molecole, frammenti, traslazioni). La σ di un parametro è
la deviazione standard dei suoi valori sui campioni.

Approssimazioni, le stesse dei programmi di cristallografia che non hanno la matrice di covarianza completa: le
coordinate sono indipendenti tra loro (nessuna correlazione), le s.u. degli atomi su posizioni speciali non sono
vincolate, gli atomi senza s.u. (spesso gli idrogeni) restano fermi.
"""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

import numpy as np

from pistak.core.cif import Cell
from pistak.core.cores import Core, axis_atoms
from pistak.core.geometry import Plane, fit_plane
from pistak.core.pairs import Pair, remeasure
from pistak.core.structure import cell_vectors

if TYPE_CHECKING:
    from pistak.core.analysis import Analysis

#: Parametri della geometria della coppia, già nei campi di `Pair`: quelli del Monte Carlo senza registro.
GEOMETRY_KEYS = ("interplanar", "long_slip", "lat_slip", "slip", "dihedral", "centroid_distance", "contact", "overlap")
#: Parametri angolari periodici (−180…180°): la dispersione si misura attorno al valore nominale.
PERIODIC = {"slip_direction"}

Evaluator = Callable[["Analysis", Pair], dict[str, float]]


def geometry_values(_analysis: Analysis, pair: Pair) -> dict[str, float]:
    return {k: float(getattr(pair, k)) for k in GEOMETRY_KEYS}


@dataclass
class Uncertainty:
    """Risultato del Monte Carlo: per ogni coppia (per `Pair.index`) e parametro, σ e media dei campioni."""

    samples: int
    seed: int
    include_cell: bool
    #: Vero se il CIF ha almeno una s.u. diversa da zero (altrimenti tutte le σ sono 0).
    has_su: bool
    keys: list[str] = field(default_factory=list)
    sigma: dict[int, dict[str, float]] = field(default_factory=dict)
    mean: dict[int, dict[str, float]] = field(default_factory=dict)
    #: Frazione dei campioni in cui la coppia ha lo stesso tipo (stack/parallel/tilted) di quella nominale.
    kind_agreement: dict[int, float] = field(default_factory=dict)

    def get(self, pair: Pair | int, key: str) -> float:
        i = pair if isinstance(pair, int) else pair.index
        return self.sigma.get(i, {}).get(key, float("nan"))


def has_su(analysis: Analysis, include_cell: bool = True) -> bool:
    atoms = any(np.any(a.su > 0) for a in analysis.cif.atoms)
    return atoms or (include_cell and any(x > 0 for x in analysis.cif.cell_su))


def perturbed(analysis: Analysis, rng: np.random.Generator, include_cell: bool = True, scale: float = 1.0) -> Analysis:
    """Una copia dell'analisi con gli atomi spostati secondo le s.u. (moltiplicate per `scale`); stesse coppie."""
    s, cif = analysis.structure, analysis.cif
    su = np.array([a.su for a in cif.atoms]).reshape(-1, 3)
    delta = rng.standard_normal(su.shape) * su * scale
    vectors = s.vectors
    if include_cell and any(x > 0 for x in cif.cell_su):
        c = cif.cell
        vals = np.array([c.a, c.b, c.c, c.alpha, c.beta, c.gamma]) + rng.standard_normal(6) * np.array(cif.cell_su) * scale
        vectors = cell_vectors(Cell(*map(float, vals)))
    rot = np.array([op.r for op in s.symops])
    d_cell = np.einsum("jab,jb->ja", rot[s.cell_origin[:, 1]], delta[s.cell_origin[:, 0]]) if len(s.cell_origin) else 0
    mols = []
    for m in s.molecules:
        frac = s.cell_frac[m.cell_index] + m.cell_shift + (d_cell[m.cell_index] if len(s.cell_origin) else 0)
        pos = frac @ vectors
        mols.append(replace(m, pos=pos, centroid=pos.mean(axis=0)))
    s2 = replace(s, vectors=vectors, molecules=mols, _images={}, _atom_maps={})
    spec = analysis.limits.long_axis_atoms
    frags = []
    for m_old, m_new, fs in zip(s.molecules, mols, analysis.fragments or [[c] if c else [] for c in analysis.cores]):
        ij = axis_atoms(m_old, spec) if spec else None
        frags.append([_refit(f, m_new.pos, ij) for f in fs])
    return replace(analysis, structure=s2, fragments=frags, cores=[f[0] if f else None for f in frags], uncertainty=None)


def _refit(f: Core, pos: np.ndarray, ij: tuple[int, int] | None) -> Core:
    pl = fit_plane(pos[f.atoms])
    normal = pl.normal if np.dot(pl.normal, f.normal) >= 0 else -pl.normal
    long = pl.long if np.dot(pl.long, f.long) >= 0 else -pl.long
    core = Core(Plane(pl.centroid, normal, long, pl.rms, pl.aspect), f.atoms, f.kind, f.axis_source)
    if f.axis_source == "atoms" and ij is not None:
        core = core.with_long(pos[ij[1]] - pos[ij[0]]) or core
    return core


def _spread(values: np.ndarray, nominal: float, periodic: bool) -> tuple[float, float]:
    ok = values[np.isfinite(values)]
    if len(ok) < 2:
        return float("nan"), float("nan")
    if periodic:
        ok = nominal + (ok - nominal + 180.0) % 360.0 - 180.0
    return float(np.std(ok, ddof=1)), float(np.mean(ok))


def monte_carlo(analysis: Analysis, samples: int = 100, seed: int = 0, evaluate: Evaluator = geometry_values,
                include_cell: bool = True, progress: Callable[[int, int], None] | None = None,
                periodic: set[str] = PERIODIC) -> Uncertainty:
    """σ di tutti i valori di `evaluate` per ogni coppia dell'analisi, su `samples` geometrie perturbate."""
    rng = np.random.default_rng(seed)
    nominal = {p.index: evaluate(analysis, p) for p in analysis.pairs}
    keys = list(next(iter(nominal.values()), {}))
    acc = {p.index: {k: np.full(samples, np.nan) for k in keys} for p in analysis.pairs}
    same = {p.index: 0 for p in analysis.pairs}
    for n in range(samples):
        a2 = perturbed(analysis, rng, include_cell)
        for p in analysis.pairs:
            q = remeasure(a2.structure, a2.fragments, p, a2.limits)
            same[p.index] += q.kind == p.kind
            for k, v in evaluate(a2, q).items():
                if k in acc[p.index]:
                    acc[p.index][k][n] = v
        if progress:
            progress(n + 1, samples)
    out = Uncertainty(samples, seed, include_cell, has_su(analysis, include_cell), keys)
    for i, by_key in acc.items():
        out.sigma[i], out.mean[i] = {}, {}
        for k, vals in by_key.items():
            nom = nominal[i].get(k, math.nan)
            out.sigma[i][k], out.mean[i][k] = _spread(vals, nom, k in periodic)
        out.kind_agreement[i] = same[i] / samples if samples else math.nan
    return out
