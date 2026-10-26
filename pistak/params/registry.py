"""Registro dei parametri misurabili su una coppia di molecole.

Un parametro è una funzione che riceve un `PairContext` e restituisce un numero. Registrarla basta perché
compaia in tabella, negli assi dello scatter, negli istogrammi e nelle esportazioni::

    from pistak.params import parameter

    @parameter("mio_param", "Il mio parametro", unit="Å", group="plugin", digits=3)
    def mio_param(ctx):
        return float(abs(ctx.d @ ctx.frame.n))

`ctx` ha: .pair, .A, .B (molecole), .core_a, .core_b (nuclei), .frame (L, W, n), .d (vettore tra i
baricentri), .pos_a, .pos_b, .proj_a, .proj_b (atomi nel sistema della coppia), .limits, .analysis.
"""
from __future__ import annotations

import math
import traceback
from collections.abc import Callable
from dataclasses import dataclass

from pistak.core.analysis import Analysis, PairContext
from pistak.i18n import tr

ParamFn = Callable[[PairContext], float]


@dataclass(frozen=True)
class ParamSpec:
    key: str
    label: str
    fn: ParamFn
    unit: str = ""
    description: str = ""
    digits: int = 2
    group: str = "extra"
    source: str = "builtin"  # "builtin" | "extra" | nome del file plugin

    @property
    def name(self) -> str:
        """L'etichetta nella lingua corrente."""
        return tr(self.label)

    @property
    def info(self) -> str:
        """La descrizione nella lingua corrente."""
        return tr(self.description) if self.description else ""

    @property
    def group_name(self) -> str:
        return tr(self.group)

    @property
    def title(self) -> str:
        return f"{self.name} ({self.unit})" if self.unit else self.name


_REGISTRY: dict[str, ParamSpec] = {}
#: Errori incontrati calcolando i parametri: chiave → ultimo messaggio. Visibili nella GUI.
errors: dict[str, str] = {}
_current_source = "extra"


def set_source(name: str) -> str:
    """Imposta il nome della sorgente (file) a cui attribuire i parametri registrati da ora. Ritorna il precedente."""
    global _current_source
    prev, _current_source = _current_source, name
    return prev


def parameter(key: str, label: str, unit: str = "", description: str = "", digits: int = 2, group: str = "extra",
              source: str | None = None):
    """Decoratore: registra la funzione come parametro `key`. Se la chiave esiste già viene sostituita."""

    def deco(fn: ParamFn) -> ParamFn:
        _REGISTRY[key] = ParamSpec(key, label, fn, unit, description or (fn.__doc__ or "").strip(), digits, group,
                                   source or _current_source)
        return fn

    return deco


def specs() -> list[ParamSpec]:
    return list(_REGISTRY.values())


def get_spec(key: str) -> ParamSpec | None:
    return _REGISTRY.get(key)


def remove_source(source: str) -> None:
    for k in [k for k, s in _REGISTRY.items() if s.source == source]:
        del _REGISTRY[k]


def sources() -> list[str]:
    return sorted({s.source for s in _REGISTRY.values()})


def compute_all(analysis: Analysis) -> Analysis:
    """Riempie `Pair.values` di ogni coppia con tutti i parametri registrati. Un parametro che solleva
    un'eccezione vale NaN per le coppie in cui fallisce e l'errore resta in `errors`."""
    keys = [s.key for s in _REGISTRY.values()]
    for k in keys:
        errors.pop(k, None)
    for pair in analysis.pairs:
        ctx = analysis.pair_context(pair)
        pair.values = {}
        for s in _REGISTRY.values():
            try:
                v = float(s.fn(ctx))
            except Exception as e:  # un plugin difettoso non deve fermare l'analisi
                v = math.nan
                errors[s.key] = f"{type(e).__name__}: {e}\n{traceback.format_exc(limit=3)}"
            pair.values[s.key] = v
    analysis.param_keys = keys
    return analysis
