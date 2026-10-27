"""Perché un'analisi non ha coppie: il motivo, con i numeri che lo spiegano, per dirlo all'utente invece di mostrare
viste vuote. Come per gli avvisi dei controlli, qui ci sono solo un codice stabile e i valori; il testo lo costruisce
`pistak.labels.no_pairs_text`.

Codici, nell'ordine in cui si controllano:

    extendedOnly       solo reti di legami estese (polimeri, solidi): nessuna molecola separata
    noContacts         ci sono nuclei π, ma nessuna coppia entro il contatto massimo
    noRings            nessun anello e molecola non planare nel complesso
    nonPlanarRings     anelli trovati, nessuno planare, e sono saturi (carboni con quattro legami)
    nearlyPlanarRings  anelli trovati, nessuno planare, ma almeno uno è insaturo: forse un aromatico distorto
    coreNotPlanar      anelli planari, ma i sistemi che formano non sono planari nel complesso

Un anello è saturo se uno dei suoi carboni ha quattro vicini (idrogeni compresi). Senza idrogeni nel modello non si
può dire: allora conta lo scarto dal piano, e un anello sotto il doppio della soglia vale come «quasi planare».
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from pistak.core.cores import small_rings
from pistak.core.geometry import fit_plane
from pistak.core.intramolecular import unique_molecules

if TYPE_CHECKING:
    from pistak.core.analysis import Analysis


@dataclass(frozen=True)
class NoPairsReason:
    code: str
    values: dict = field(default_factory=dict)


@dataclass(frozen=True)
class RingInfo:
    molecule: int
    labels: tuple[str, ...]
    rms: float  # Å, scarto dal piano medio
    #: True se ha un carbonio con quattro vicini; None se la molecola non ha idrogeni (non si può dire).
    saturated: bool | None


def ring_planarity(analysis: Analysis) -> list[RingInfo]:
    """Gli anelli (senza idrogeni) delle molecole indipendenti per simmetria, con planarità e saturazione."""
    s, lim = analysis.structure, analysis.limits
    out = []
    for mi in unique_molecules(s):
        m = s.molecules[mi]
        if m.extended:
            continue
        degree = np.bincount(m.bonds.ravel(), minlength=len(m)) if len(m.bonds) else np.zeros(len(m), dtype=int)
        has_h = "H" in m.el
        heavy = m.heavy
        for r in small_rings(len(m), m.bonds, lambda i, h=heavy: h[i], lim.max_ring_size):
            sat = any(m.el[i] == "C" and degree[i] >= 4 for i in r) if has_h else None
            out.append(RingInfo(mi, tuple(m.labels[i] for i in r), float(fit_plane(m.pos[r]).rms), sat))
    return out


def why_no_pairs(analysis: Analysis) -> NoPairsReason | None:
    """Il motivo per cui l'analisi non ha coppie; None se ne ha."""
    if analysis.pairs:
        return None
    s, lim = analysis.structure, analysis.limits
    molecules = [m for m in s.molecules if not m.extended]
    if not molecules:
        return NoPairsReason("extendedOnly")
    with_core = sum(c is not None for c in analysis.cores)
    if with_core:
        return NoPairsReason("noContacts", {"n": with_core, "limit": lim.contact_cutoff})
    rings = ring_planarity(analysis)
    if not rings:
        rms = min((float(fit_plane(m.pos[np.flatnonzero(m.heavy)]).rms) for m in molecules if m.heavy.sum() >= 4),
                  default=float("nan"))
        return NoPairsReason("noRings", {"size": lim.max_ring_size, "rms": f"{rms:.2f}", "limit": f"{lim.molecule_rms:.2f}"})
    values = [r.rms for r in rings]
    if min(values) >= lim.ring_rms:
        # forse un aromatico distorto: un anello insaturo, o (senza idrogeni) uno sotto il doppio della soglia
        nearly = [r.rms for r in rings if r.saturated is False or (r.saturated is None and r.rms < 2 * lim.ring_rms)]
        shown = nearly or values
        return NoPairsReason("nearlyPlanarRings" if nearly else "nonPlanarRings",
                             {"n": len(shown), "min": f"{min(shown):.3f}", "max": f"{max(shown):.3f}",
                              "limit": f"{lim.ring_rms:.2f}"})
    return NoPairsReason("coreNotPlanar", {"limit": f"{lim.molecule_rms:.2f}", "planar": sum(r < lim.ring_rms for r in values)})
