"""Il Monte Carlo delle incertezze (`pistak.core.uncertainty`) con i parametri del registro: σ per ogni parametro,
compresi quelli dei plugin."""
from __future__ import annotations

import math
from collections.abc import Callable, Iterable

from pistak.core.analysis import Analysis
from pistak.core.pairs import Pair
from pistak.core.uncertainty import Uncertainty, monte_carlo
from pistak.params.registry import REGISTRY, ParamSpec, Registry


def default_keys(registry: Registry = REGISTRY) -> list[str]:
    """I parametri per cui ha senso una σ: quelli con decimali (esclusi conteggi e numeri di frammento)."""
    return [s.key for s in registry.specs() if s.digits > 0]


def registry_evaluator(keys: Iterable[str] | None = None, registry: Registry = REGISTRY) -> Callable[[Analysis, Pair], dict[str, float]]:
    wanted = set(keys) if keys is not None else set(default_keys(registry))
    specs: list[ParamSpec] = [s for s in registry.specs() if s.key in wanted]

    def evaluate(analysis: Analysis, pair: Pair) -> dict[str, float]:
        ctx = analysis.pair_context(pair)
        out = {}
        for s in specs:
            try:
                out[s.key] = float(s.fn(ctx))
            except Exception:  # noqa: BLE001  (un plugin che fallisce su un campione: quel valore manca)
                out[s.key] = math.nan
        return out

    return evaluate


def uncertainty(analysis: Analysis, samples: int = 100, seed: int = 0, keys: Iterable[str] | None = None,
                include_cell: bool = True, progress: Callable[[int, int], None] | None = None,
                registry: Registry = REGISTRY) -> Uncertainty:
    """Calcola le σ di tutti i parametri (o di `keys`) e le salva in `analysis.uncertainty`."""
    result = monte_carlo(analysis, samples, seed, registry_evaluator(keys, registry), include_cell, progress)
    analysis.uncertainty = result
    return result
