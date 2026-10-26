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


def analyse_text(text: str, name: str = "", limits: Limits | None = None) -> Analysis:
    return finish(_core.analyse_text(text, name, limits or load_limits()))


def analyse_file(path: str | Path, limits: Limits | None = None) -> Analysis:
    return finish(_core.analyse_file(path, limits or load_limits()))
