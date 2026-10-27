"""Registro delle grandezze dell'intera struttura (non di una coppia): densità, Z′, valori dell'affinamento, avvisi
dei controlli rapidi…

Funziona come quello dei parametri delle coppie: una funzione riceve l'`Analysis` e restituisce un numero; registrarla
basta perché compaia nel Dettaglio (scheda Struttura), nel riepilogo di `analyze`, nel JSON (`summary.structure`), nel
CSV di `batch` e nella scheda Strutture del confronto::

    from pistak.params import structure_parameter

    @structure_parameter("n_molecole", "Molecole nella cella", group="plugin", digits=0)
    def n_molecole(analysis):
        return len(analysis.structure.molecules)

I valori stanno in `analysis.structure_values` (chiave → numero, NaN se il calcolo fallisce: l'errore va in
`analysis.param_errors` con la chiave preceduta da «structure:»). Gli avvisi dei controlli rapidi
(`pistak.core.validation`) si calcolano una volta e stanno in `analysis.alerts`.
"""
from __future__ import annotations

import math
import traceback
from collections.abc import Callable
from dataclasses import dataclass

from pistak.core.analysis import Analysis
from pistak.core.crystal import cell_contents, cell_volume, cif_value, density, z_prime
from pistak.core.validation import counts, validate
from pistak.i18n import tr

StructureFn = Callable[[Analysis], float]


@dataclass(frozen=True)
class StructureSpec:
    key: str
    label: str
    fn: StructureFn
    unit: str = ""
    description: str = ""
    digits: int = 2
    group: str = "cristallo"
    source: str = "builtin"

    @property
    def name(self) -> str:
        return tr(self.label)

    @property
    def info(self) -> str:
        return tr(self.description) if self.description else ""

    @property
    def group_name(self) -> str:
        return tr(self.group)

    @property
    def title(self) -> str:
        return f"{self.name} ({self.unit})" if self.unit else self.name


class StructureRegistry:
    def __init__(self) -> None:
        self._specs: dict[str, StructureSpec] = {}
        self._source = "builtin"

    def parameter(self, key: str, label: str, unit: str = "", description: str = "", digits: int = 2,
                  group: str = "cristallo", source: str | None = None):
        def deco(fn: StructureFn) -> StructureFn:
            self._specs[key] = StructureSpec(key, label, fn, unit, description or (fn.__doc__ or "").strip(), digits, group,
                                             source or self._source)
            return fn

        return deco

    def set_source(self, name: str) -> str:
        prev, self._source = self._source, name
        return prev

    def remove_source(self, source: str) -> None:
        for k in [k for k, s in self._specs.items() if s.source == source]:
            del self._specs[k]

    def specs(self) -> list[StructureSpec]:
        return list(self._specs.values())

    def get(self, key: str) -> StructureSpec | None:
        return self._specs.get(key)

    def keys(self) -> list[str]:
        return list(self._specs)

    def compute_all(self, analysis: Analysis) -> Analysis:
        """Avvisi dei controlli (se non ci sono già) e tutte le grandezze del registro."""
        if analysis.alerts is None:
            analysis.alerts = validate(analysis)
        analysis.structure_values = {}
        for s in self.specs():
            try:
                v = float(s.fn(analysis))
            except Exception as e:  # noqa: BLE001  (un plugin difettoso non deve fermare l'analisi)
                v = math.nan
                analysis.param_errors[f"structure:{s.key}"] = f"{type(e).__name__}: {e}\n{traceback.format_exc(limit=3)}"
            analysis.structure_values[s.key] = v
        return analysis


#: Il registro predefinito delle grandezze di struttura.
STRUCTURE_REGISTRY = StructureRegistry()
structure_parameter = STRUCTURE_REGISTRY.parameter
structure_specs = STRUCTURE_REGISTRY.specs
get_structure_spec = STRUCTURE_REGISTRY.get


def compute_structure(analysis: Analysis, registry: StructureRegistry | None = None) -> Analysis:
    return (registry or STRUCTURE_REGISTRY).compute_all(analysis)


# ---- grandezze predefinite -------------------------------------------------------------------------------------
def _cif(key: str, label: str, unit: str, description: str, digits: int, group: str) -> None:
    structure_parameter(key, label, unit, description, digits, group)(lambda a, k=key: cif_value(a.cif, k))


structure_parameter("volume", "Volume della cella", "Å³", "Dai parametri di cella.", 1)(lambda a: cell_volume(a.cif))
_cif("z", "Z", "", "Unità formula nella cella (_cell_formula_units_Z).", 0, "cristallo")
structure_parameter("z_prime", "Z′", "", "Z diviso il numero di posizioni generali del gruppo spaziale: quante unità "
                    "formula indipendenti ci sono nell'unità asimmetrica.", 2)(lambda a: z_prime(a.cif))
structure_parameter("density_calc", "Densità calcolata", "g/cm³", "Dal contenuto della cella del modello (tutte le "
                    "componenti disordinate, pesate con l'occupanza) e dal volume.", 3)(
    lambda a: density(cell_contents(a.cif), cell_volume(a.cif)))
_cif("density_cif", "Densità dichiarata", "g/cm³", "_exptl_crystal_density_diffrn del CIF.", 3, "cristallo")
_cif("temperature", "Temperatura", "K", "Temperatura della misura (_diffrn_ambient_temperature).", 0, "cristallo")
_cif("r1", "R1", "", "_refine_ls_R_factor_gt: accordo tra modello e riflessi osservati.", 4, "affinamento")
_cif("wr2", "wR2", "", "_refine_ls_wR_factor_ref.", 4, "affinamento")
_cif("goof", "GoF", "", "_refine_ls_goodness_of_fit_ref: vicino a 1 per un modello e pesi corretti.", 3, "affinamento")
_cif("rho_max", "Densità residua massima", "e/Å³", "_refine_diff_density_max.", 2, "affinamento")
_cif("rho_min", "Densità residua minima", "e/Å³", "_refine_diff_density_min.", 2, "affinamento")
for _lev in "ABC":
    structure_parameter(f"alerts_{_lev.lower()}", f"Avvisi {_lev}", "", f"Avvisi di livello {_lev} dei controlli rapidi.", 0,
                        "controlli")(lambda a, lev=_lev: counts(a.alerts or [])[lev])
