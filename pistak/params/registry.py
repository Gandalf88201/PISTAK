"""Registro dei parametri misurabili su una coppia di molecole.

Un parametro è una funzione che riceve un `PairContext` e restituisce un numero. Registrarla basta perché
compaia in tabella, negli assi dello scatter, negli istogrammi e nelle esportazioni::

    from pistak.params import parameter

    @parameter("mio_param", "Il mio parametro", unit="Å", group="plugin", digits=3)
    def mio_param(ctx):
        return float(abs(ctx.d @ ctx.frame.n))

`ctx` ha: .pair, .A, .B (molecole), .core_a, .core_b (nuclei), .frame (L, W, n), .d (vettore tra i
baricentri), .pos_a, .pos_b, .proj_a, .proj_b (atomi nel sistema della coppia), .limits, .analysis.

Con ``geometric=True`` si promette che il valore dipende solo dalla geometria della coppia (molecole, frammenti,
traslazione), non dal suo tipo, dalle altre coppie o dalle soglie del criterio di stack: allora, quando cambiano solo
quelle soglie (`pistak.core.stages`, stadio `pairs`), il valore si riusa invece di ricalcolarlo. Per impostazione è
False: i parametri si ricalcolano sempre.
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
    #: Il valore dipende solo dalla geometria della coppia: si riusa quando cambiano le soglie del criterio di stack.
    geometric: bool = False

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


class Registry:
    """Insieme di parametri, per chiave. Di solito si usa quello predefinito, `REGISTRY` (lo stesso che riempiono
    `@parameter`, i parametri builtin/extra e i plugin); un `Registry()` nuovo serve a calcolare un insieme
    diverso (per esempio solo alcuni parametri) senza toccare quello globale."""

    def __init__(self) -> None:
        self._specs: dict[str, ParamSpec] = {}
        self._source = "extra"
        #: Problemi incontrati caricando i plugin in questo registro: nome file → messaggio.
        self.load_errors: dict[str, str] = {}
        self.plugins_loaded = False

    # ---- registrazione -------------------------------------------------------------------
    def parameter(self, key: str, label: str, unit: str = "", description: str = "", digits: int = 2,
                  group: str = "extra", source: str | None = None, geometric: bool = False):
        """Decoratore: registra la funzione come parametro `key`. Se la chiave esiste già viene sostituita."""

        def deco(fn: ParamFn) -> ParamFn:
            self.add(ParamSpec(key, label, fn, unit, description or (fn.__doc__ or "").strip(), digits, group,
                               source or self._source, geometric))
            return fn

        return deco

    def add(self, spec: ParamSpec) -> None:
        self._specs[spec.key] = spec

    def set_source(self, name: str) -> str:
        """Imposta il nome della sorgente (file) a cui attribuire i parametri registrati da ora. Ritorna il precedente."""
        prev, self._source = self._source, name
        return prev

    def remove_source(self, source: str) -> None:
        for k in [k for k, s in self._specs.items() if s.source == source]:
            del self._specs[k]

    # ---- consultazione -------------------------------------------------------------------
    def specs(self) -> list[ParamSpec]:
        return list(self._specs.values())

    def get(self, key: str) -> ParamSpec | None:
        return self._specs.get(key)

    def keys(self) -> list[str]:
        return list(self._specs)

    def sources(self) -> list[str]:
        return sorted({s.source for s in self._specs.values()})

    def __contains__(self, key: str) -> bool:
        return key in self._specs

    def __len__(self) -> int:
        return len(self._specs)

    # ---- calcolo -------------------------------------------------------------------------
    def compute_all(self, analysis: Analysis) -> Analysis:
        """Riempie `Pair.values` di ogni coppia con tutti i parametri del registro. Un parametro che solleva
        un'eccezione vale NaN per le coppie in cui fallisce; l'errore resta in `analysis.param_errors`. I valori dei
        parametri `geometric` già calcolati per la stessa coppia misurata (`analysis.cache`) si riusano."""
        specs = self.specs()
        analysis.param_errors = {}
        cache = analysis.cache.values if analysis.cache is not None else None
        for pair in analysis.pairs:
            ctx = analysis.pair_context(pair)
            known = cache.setdefault(pair.key, {}) if cache is not None else {}
            pair.values = {}
            for s in specs:
                if s.geometric and s.key in known:
                    pair.values[s.key] = known[s.key]
                    continue
                try:
                    v = float(s.fn(ctx))
                except Exception as e:  # un plugin difettoso non deve fermare l'analisi
                    v = math.nan
                    analysis.param_errors[s.key] = f"{type(e).__name__}: {e}\n{traceback.format_exc(limit=3)}"
                else:
                    if s.geometric:  # gli errori no: si ripresentano (e si segnalano) a ogni calcolo
                        known[s.key] = v
                pair.values[s.key] = v
        analysis.param_keys = [s.key for s in specs]
        return analysis


#: Il registro predefinito: lo usano `@parameter`, la GUI, il terminale e gli export.
REGISTRY = Registry()

# Funzioni di modulo che lavorano sul registro predefinito (l'interfaccia usata da plugin e versioni precedenti).
parameter = REGISTRY.parameter
set_source = REGISTRY.set_source
remove_source = REGISTRY.remove_source
specs = REGISTRY.specs
get_spec = REGISTRY.get
sources = REGISTRY.sources


def compute_all(analysis: Analysis, registry: Registry | None = None) -> Analysis:
    """Calcola tutti i parametri di `registry` (per impostazione quello predefinito) sulle coppie di `analysis`."""
    return (registry or REGISTRY).compute_all(analysis)
