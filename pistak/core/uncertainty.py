"""Incertezze dei parametri, propagate con un Monte Carlo dalle incertezze standard (s.u.) del CIF.

Ogni campione sposta le coordinate frazionarie degli atomi dell'unità asimmetrica di quantità casuali gaussiane con
deviazione standard pari alla loro s.u. (e, se richiesto, i parametri di cella); lo spostamento passa a tutti gli
atomi della cella attraverso la parte rotazionale dell'operazione che li genera, quindi le copie per simmetria si
muovono in modo coerente. Si ricalcolano piani e assi dei frammenti π (con lo stesso verso di quelli nominali) e si
rimisura ogni coppia, con la stessa topologia dell'analisi (molecole, frammenti, traslazioni). La σ di un parametro è
la deviazione standard dei suoi valori sui campioni.

Approssimazioni, le stesse dei programmi di cristallografia che non hanno la matrice di covarianza completa: le
coordinate sono indipendenti tra loro (nessuna correlazione), gli atomi senza s.u. (spesso gli idrogeni) restano
fermi. Per un atomo su una posizione speciale si estrae una deviata per ogni parametro libero e i vincoli del sito danno
le altre coordinate (`site_constraints`): ogni parametro libero ha la sua s.u. del CIF.
"""
from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from typing import TYPE_CHECKING

import numpy as np

from pistak.core.cif import Cell
from pistak.core.cores import FALLBACK_SOURCES, Core, axis_atoms, fallback_axis
from pistak.core.geometry import Plane, fit_plane
from pistak.core.pairs import Pair, remeasure
from pistak.core.structure import cell_vectors

if TYPE_CHECKING:
    from pistak.core.analysis import Analysis

#: Parametri della geometria della coppia, già nei campi di `Pair`: quelli del Monte Carlo senza registro.
GEOMETRY_KEYS = ("interplanar", "long_slip", "lat_slip", "slip", "dihedral", "centroid_distance", "contact", "overlap",
                 "long_slip_a", "lat_slip_a", "long_slip_b", "lat_slip_b")
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
    #: Campioni con un valore finito, per coppia e parametro: la σ è calcolata solo su questi.
    valid: dict[int, dict[str, int]] = field(default_factory=dict)
    #: Parametri con un valore nominale ma mancanti (NaN) in parte dei campioni, per coppia: {chiave: campioni validi}.
    #: Non si scartano in silenzio: la σ viene da meno campioni, e chi fallisce può essere proprio la geometria estrema.
    incomplete: dict[int, dict[str, int]] = field(default_factory=dict)
    #: Intervallo dei percentili 2,5 e 97,5 dei campioni (95 % centrale), per coppia e parametro: per le grandezze
    #: limitate (angolo tra i piani vicino a 0°, sovrapposizione vicina a 0 o 1) la media ± σ non lo descrive bene.
    interval: dict[int, dict[str, tuple[float, float]]] = field(default_factory=dict)
    #: Vero se gli spostamenti vengono da una matrice di covarianza data dall'utente e non dalle s.u. indipendenti.
    covariance: bool = False

    def get(self, pair: Pair | int, key: str) -> float:
        i = pair if isinstance(pair, int) else pair.index
        return self.sigma.get(i, {}).get(key, float("nan"))


def has_su(analysis: Analysis, include_cell: bool = True) -> bool:
    atoms = any(np.any(a.su > 0) for a in analysis.cif.atoms)
    return atoms or (include_cell and any(x > 0 for x in analysis.cif.cell_su))


_PROJECTORS: dict[int, tuple[object, dict[int, np.ndarray]]] = {}


def site_projectors(cif) -> dict[int, np.ndarray]:
    """Per ogni atomo del CIF su una posizione speciale, la matrice (frazionaria) che proietta uno spostamento sul
    sottospazio permesso dalla simmetria del sito: la media delle parti rotazionali delle operazioni che mandano il
    sito in sé (a meno di traslazioni). Gli atomi in posizione generale non compaiono."""
    hit = _PROJECTORS.get(id(cif))
    if hit is not None and hit[0] is cif:
        return hit[1]
    out: dict[int, np.ndarray] = {}
    for i, at in enumerate(cif.atoms):
        rots = []
        for op in cif.symops:
            d = op.apply(at.frac) - at.frac
            if np.all(np.abs(d - np.round(d)) < 1e-4):
                rots.append(op.r)
        if len(rots) > 1:
            out[i] = np.mean(rots, axis=0)
    _PROJECTORS[id(cif)] = (cif, out)
    return out


def constraint_matrix(proj: np.ndarray, su: np.ndarray) -> np.ndarray:
    """La matrice C (3×3, frazionaria) che dà lo spostamento di un atomo su una posizione speciale dalle deviate delle
    sue coordinate: δ = C·δ₀, con δ₀ le deviate indipendenti (s.u. del CIF) di x, y, z. Le colonne delle coordinate
    libere sono una base del sottospazio V permesso dalla simmetria del sito (l'immagine di `proj`) con 1 sulla
    coordinata stessa e 0 sulle altre libere; le colonne delle coordinate dipendenti o fissate sono nulle.

    Così ogni parametro libero riceve una sola deviata con la sua s.u. e i vincoli fanno il resto: per (x, x, z)
    δx = δy = δ₀x e δz = δ₀z, per (x, 2x, z) δy = 2δ₀x, per (x, −x, z) δy = −δ₀x, per uno specchio z = 0 δz = 0.
    La proiezione `proj` da sola (media di R_g sulle operazioni del sito) mediava due deviate indipendenti: per (x, x, z)
    lo spostamento lungo x = y aveva σ/√2 invece della s.u. del parametro raffinato (revisione 3, N2). Le coordinate
    libere si scelgono in ordine x, y, z preferendo quelle con una s.u. nel CIF (un programma può scrivere la s.u. solo
    sulla coordinata raffinata). Anche C è un proiettore su V (C·v = v per ogni v in V), ma lungo le coordinate
    dipendenti invece che ortogonale."""
    n = proj.shape[0]  # 3 per le coordinate, 6 per i parametri di cella (`cell_constraints`)
    u, sv, _ = np.linalg.svd(proj)
    k = int(np.sum(sv > 1e-6))
    if k == 0:
        return np.zeros((n, n))
    if k == n:
        return np.eye(n)
    basis = u[:, :k]
    free: list[int] = []
    for j in sorted(range(n), key=lambda j: (not su[j] > 0, j)):
        if np.linalg.matrix_rank(basis[free + [j], :], tol=1e-6) == len(free) + 1:
            free.append(j)
        if len(free) == k:
            break
    out = np.zeros((n, n))
    out[:, free] = basis @ np.linalg.inv(basis[free, :])
    out[np.abs(out) < 1e-12] = 0.0
    return out


_CONSTRAINTS: dict[int, tuple[object, dict[int, np.ndarray]]] = {}


def site_constraints(cif) -> dict[int, np.ndarray]:
    """Per ogni atomo del CIF su una posizione speciale, la matrice `constraint_matrix` (dalle s.u. del CIF)."""
    hit = _CONSTRAINTS.get(id(cif))
    if hit is not None and hit[0] is cif:
        return hit[1]
    out = {i: constraint_matrix(p, np.asarray(cif.atoms[i].su, dtype=float).reshape(3))
           for i, p in site_projectors(cif).items()}
    _CONSTRAINTS[id(cif)] = (cif, out)
    return out


def _metric(p: np.ndarray) -> np.ndarray:
    a, b, c = p[:3]
    ca, cb, cg = np.cos(np.radians(p[3:]))
    return np.array([[a * a, a * b * cg, a * c * cb], [a * b * cg, b * b, b * c * ca], [a * c * cb, b * c * ca, c * c]])


def _cell_of(g: np.ndarray) -> np.ndarray:
    a, b, c = np.sqrt(np.diag(g))
    return np.array([a, b, c, *np.degrees(np.arccos(np.clip([g[1, 2] / (b * c), g[0, 2] / (a * c), g[0, 1] / (a * b)],
                                                           -1.0, 1.0)))])


_CELL: dict[int, tuple[object, np.ndarray | None]] = {}
#: Tolleranza sulla metrica nominale (frazione del più grande a², b², c²): oltre, la cella del CIF non ha la simmetria
#: delle operazioni (cella o operazioni sbagliate) e i vincoli non si possono ricavare.
CELL_METRIC_TOLERANCE = 5e-4


def cell_constraints(cif) -> np.ndarray | None:
    """La matrice dei vincoli (6×6, `constraint_matrix`) dei parametri di cella a, b, c, α, β, γ, ricavata dalle
    operazioni di simmetria e non dal sistema cristallino: la metrica G deve restare invariata, G = Rᵀ G R per la parte
    rotazionale R di ogni operazione. La media di Rᵀ G R sulle operazioni proietta una metrica qualunque su quelle
    permesse; il suo jacobiano nei parametri di cella è il proiettore sui parametri liberi, da cui `constraint_matrix`
    sceglie le colonne libere (prima quelle con una s.u. nel CIF). Tetragonale: δb = δa, α, β, γ fermi; esagonale e
    trigonale: δb = δa, γ fermo a 120°; cubico: δb = δc = δa; monoclinico: α e γ fermi; triclino: tutto libero; vale
    in qualunque impostazione (assi unici, celle romboedriche). None se la metrica nominale non è invariante entro
    `CELL_METRIC_TOLERANCE` (allora il Monte Carlo perturba i sei parametri indipendenti, come prima, e i controlli
    rapidi danno l'avviso `cellConstraints`) (revisione 3)."""
    hit = _CELL.get(id(cif))
    if hit is not None and hit[0] is cif:
        return hit[1]
    c = cif.cell
    p0 = np.array([c.a, c.b, c.c, c.alpha, c.beta, c.gamma], dtype=float)
    rots = {tuple(np.asarray(op.r, dtype=float).round(6).ravel()) for op in cif.symops}
    rots = [np.array(r).reshape(3, 3) for r in rots]

    def project(g: np.ndarray) -> np.ndarray:
        return np.mean([r.T @ g @ r for r in rots], axis=0)

    g0 = _metric(p0)
    out: np.ndarray | None
    if not rots or np.max(np.abs(project(g0) - g0)) > CELL_METRIC_TOLERANCE * np.max(np.diag(g0)):
        out = None
    else:
        jac = np.zeros((6, 6))
        for j in range(6):
            h = 1e-6 * max(1.0, abs(p0[j]))
            dp = np.zeros(6)
            dp[j] = h
            jac[:, j] = (_cell_of(project(_metric(p0 + dp))) - _cell_of(project(_metric(p0 - dp)))) / (2 * h)
        jac[np.abs(jac) < 1e-6] = 0.0
        out = constraint_matrix(jac, np.asarray(cif.cell_su, dtype=float))
        out[np.abs(out - np.round(out)) < 1e-6] = np.round(out[np.abs(out - np.round(out)) < 1e-6])
    _CELL[id(cif)] = (cif, out)
    return out


def covariance_factor(cov: np.ndarray) -> np.ndarray:
    """Il fattore F (3N×3N) con F·Fᵀ = `cov`: F·z, con z deviate normali standard indipendenti, ha covarianza `cov`.
    Dalla decomposizione agli autovalori, non da Cholesky: la matrice di un raffinamento con vincoli (posizioni
    speciali, atomi fissati) è singolare; gli autovalori negativi per arrotondamento si azzerano."""
    cov = np.asarray(cov, dtype=float)
    values, vectors = np.linalg.eigh((cov + cov.T) / 2)
    return vectors * np.sqrt(np.clip(values, 0.0, None))


def read_covariance(path, cif) -> np.ndarray:
    """Legge la matrice di covarianza delle coordinate frazionarie degli atomi del CIF (revisione 4, M4.5): 3N×3N,
    righe e colonne nell'ordine x₁, y₁, z₁, x₂, … degli atomi di `cif.atoms` (quelli tenuti dalla regola di
    disordine, nell'ordine del CIF), in unità frazionarie al quadrato. File `.npy` di NumPy o testo con 3N righe di 3N
    numeri. Il CIF non la contiene: va esportata dal programma di raffinamento. ValueError se la forma non è quella
    attesa, se non è simmetrica o se la diagonale è negativa."""
    from pathlib import Path

    path = Path(path)
    cov = np.load(path) if path.suffix == ".npy" else np.loadtxt(path)
    n = 3 * len(cif.atoms)
    if cov.shape != (n, n):
        raise ValueError(f"covariance matrix {cov.shape}: expected {n}×{n} (3 × {len(cif.atoms)} atoms of the CIF block)")
    if not np.allclose(cov, cov.T, rtol=1e-6, atol=1e-14) or (np.diag(cov) < 0).any():
        raise ValueError("covariance matrix: not symmetric or negative variances")
    return cov


def perturbed(analysis: Analysis, rng: np.random.Generator, include_cell: bool = True, scale: float = 1.0,
              factor: np.ndarray | None = None) -> Analysis:
    """Una copia dell'analisi con gli atomi spostati secondo le s.u. (moltiplicate per `scale`); stesse coppie. Con
    `factor` (`covariance_factor` della matrice di covarianza) gli spostamenti sono correlati come nel raffinamento; i
    vincoli delle posizioni speciali sono allora già nella matrice e non si applicano di nuovo."""
    s, cif = analysis.structure, analysis.cif
    su = np.array([a.su for a in cif.atoms]).reshape(-1, 3)
    if factor is not None:
        delta = (factor @ rng.standard_normal(factor.shape[1])).reshape(-1, 3) * scale
    else:
        delta = rng.standard_normal(su.shape) * su * scale
        # atomi su posizioni speciali: una deviata per parametro libero, le altre coordinate dai vincoli del sito, così
        # x = y o z = ¼ restano soddisfatti (revisione 2, voce 61) e ogni parametro libero ha la sua s.u. (revisione 3, N2)
        for i, c in site_constraints(cif).items():
            delta[i] = c @ delta[i]
    vectors = s.vectors
    if include_cell and any(x > 0 for x in cif.cell_su):
        c = cif.cell
        d = rng.standard_normal(6) * np.array(cif.cell_su) * scale
        # una deviata per parametro di cella libero, gli altri dai vincoli della simmetria (a = b, γ = 120°…): le
        # coppie equivalenti per simmetria restano equivalenti in ogni campione (revisione 3)
        cons = cell_constraints(cif)
        if cons is not None:
            d = cons @ d
        vals = np.array([c.a, c.b, c.c, c.alpha, c.beta, c.gamma]) + d
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
        frags.append([_refit(f, m_new, ij, analysis.limits) for f in fs])
    return replace(analysis, structure=s2, fragments=frags, cores=[f[0] if f else None for f in frags], uncertainty=None,
                   cache=None)


def _refit(f: Core, m, ij: tuple[int, int] | None, limits=None) -> Core:
    """Lo stesso frammento sulla geometria perturbata `m`: piano e asse ricalcolati con il verso di quelli nominali.
    L'asse definito da due atomi e quello di ripiego (molecola o sostituenti) si ricalcolano dalle nuove posizioni; se in
    un campione il ripiego non è più definito (isotropia oltre la soglia) resta la direzione nominale."""
    pos = m.pos
    pl = fit_plane(pos[f.atoms])
    normal = pl.normal if np.dot(pl.normal, f.normal) >= 0 else -pl.normal
    long = pl.long if np.dot(pl.long, f.long) >= 0 else -pl.long
    core = Core(Plane(pl.centroid, normal, long, pl.rms, pl.aspect), f.atoms, f.kind, f.axis_source, f.rings)
    if f.axis_source == "atoms" and ij is not None:
        core = core.with_long(pos[ij[1]] - pos[ij[0]]) or core
    elif f.axis_source in FALLBACK_SOURCES:
        hit = fallback_axis(m, core, limits.with_values(axis_fallback=f.axis_source, max_axis_isotropy=1.0)) if limits else None
        v = hit[0] if hit is not None else f.long
        core = core.with_long(v if np.dot(v, f.long) >= 0 else -v, source=f.axis_source) or core
    elif f.axis_source == "atoms":
        core = core.with_long(f.long) or core
    return core


def _spread(values: np.ndarray, nominal: float, periodic: bool) -> tuple[float, float, tuple[float, float]]:
    ok = values[np.isfinite(values)]
    if len(ok) < 2:
        return float("nan"), float("nan"), (float("nan"), float("nan"))
    if periodic:
        ok = nominal + (ok - nominal + 180.0) % 360.0 - 180.0
    lo, hi = np.percentile(ok, [2.5, 97.5])
    return float(np.std(ok, ddof=1)), float(np.mean(ok)), (float(lo), float(hi))


def monte_carlo(analysis: Analysis, samples: int = 100, seed: int = 0, evaluate: Evaluator = geometry_values,
                include_cell: bool = True, progress: Callable[[int, int], None] | None = None,
                periodic: set[str] = PERIODIC, covariance: np.ndarray | None = None) -> Uncertainty:
    """σ di tutti i valori di `evaluate` per ogni coppia dell'analisi, su `samples` geometrie perturbate.

    Il valore riportato di ogni descrittore resta quello nominale (dalle coordinate del CIF); `Uncertainty.mean` è la
    media dei campioni, che per le grandezze non negative vicine a zero (θ, |s|) è spostata verso l'alto e non va
    usata come stima. Il generatore è `numpy.random.default_rng(seed)`: stesso seme, stessi risultati. `covariance`:
    matrice di covarianza delle coordinate frazionarie (`read_covariance`) al posto delle s.u. indipendenti."""
    rng = np.random.default_rng(seed)
    factor = covariance_factor(covariance) if covariance is not None else None
    if factor is not None and factor.shape[0] != 3 * len(analysis.cif.atoms):
        raise ValueError(f"covariance matrix {covariance.shape}: expected 3 × {len(analysis.cif.atoms)} atoms")
    nominal = {p.index: evaluate(analysis, p) for p in analysis.pairs}
    keys = list(next(iter(nominal.values()), {}))
    acc = {p.index: {k: np.full(samples, np.nan) for k in keys} for p in analysis.pairs}
    same = {p.index: 0 for p in analysis.pairs}
    for n in range(samples):
        a2 = perturbed(analysis, rng, include_cell, factor=factor)
        for p in analysis.pairs:
            q = remeasure(a2.structure, a2.fragments, p, a2.limits)
            same[p.index] += q.kind == p.kind
            for k, v in evaluate(a2, q).items():
                if k in acc[p.index]:
                    acc[p.index][k][n] = v
        if progress:
            progress(n + 1, samples)
    out = Uncertainty(samples, seed, include_cell, has_su(analysis, include_cell) or factor is not None, keys)
    out.covariance = factor is not None
    for i, by_key in acc.items():
        out.sigma[i], out.mean[i], out.valid[i], out.interval[i] = {}, {}, {}, {}
        for k, vals in by_key.items():
            nom = nominal[i].get(k, math.nan)
            out.sigma[i][k], out.mean[i][k], out.interval[i][k] = _spread(vals, nom, k in periodic)
            n_ok = out.valid[i][k] = int(np.count_nonzero(np.isfinite(vals)))
            if math.isfinite(nom) and n_ok < samples:
                out.incomplete.setdefault(i, {})[k] = n_ok
        out.kind_agreement[i] = same[i] / samples if samples else math.nan
    return out
