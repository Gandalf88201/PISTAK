"""Punto d'ingresso per chi usa PISTAK come libreria: analisi completa, con tutti i parametri del registro."""
from __future__ import annotations

from pathlib import Path

from pistak.core import analysis as _core
from pistak.core.analysis import Analysis
from pistak.core.config import Limits, load_limits
from pistak.params import REGISTRY, compute_all, load_plugins


def ensure_plugins(force: bool = False) -> list[str]:
    """Carica i plugin la prima volta (o sempre, con force=True). Ritorna i file caricati."""
    if REGISTRY.plugins_loaded and not force:
        return []
    return load_plugins()


def finish(analysis: Analysis) -> Analysis:
    """Calcola i parametri del registro per tutte le coppie di un'analisi già fatta."""
    ensure_plugins()
    return compute_all(analysis)


def analyse_text(text: str, name: str = "", limits: Limits | None = None, block: str | int | None = None) -> Analysis:
    return finish(_core.analyse_text(text, name, limits or load_limits(), block))


def analyse_file(path: str | Path, limits: Limits | None = None, block: str | int | None = None) -> Analysis:
    return finish(_core.analyse_file(path, limits or load_limits(), block))


def reanalyse(analysis: Analysis, limits: Limits) -> Analysis:
    """L'analisi con altre soglie, rifacendo solo gli stadi che le soglie cambiate toccano (`pistak.core.stages`) e
    riusando i parametri geometrici già calcolati. Se le soglie sono le stesse ritorna `analysis`."""
    new = _core.reanalyse(analysis, limits)
    return new if new is analysis else finish(new)


def uncertainty(analysis: Analysis, samples: int = 100, seed: int = 0, keys=None, include_cell: bool = True, progress=None):
    """σ dei parametri dalle s.u. del CIF con un Monte Carlo (vedi `pistak.core.uncertainty`); il risultato va anche
    in `analysis.uncertainty`."""
    from pistak.params.uncertainty import uncertainty as _mc

    ensure_plugins()
    return _mc(analysis, samples, seed, keys, include_cell, progress)
