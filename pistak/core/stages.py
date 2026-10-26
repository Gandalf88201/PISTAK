"""Stadi dell'analisi e soglie da cui dipende ciascuno: cambiando le soglie si rifà solo da dove serve.

La catena di `pistak.core.analysis` è divisa in stadi, ognuno dipendente solo dai precedenti e da alcune soglie:

    cif         lettura del blocco data_ (regola di disordine)
    structure   cella → atomi simmetrizzati → legami → molecole
    fragments   anelli e frammenti π, asse lungo
    contacts    coppie di frammenti vicini e la loro geometria (la parte costosa), equivalenza per simmetria
    pairs       classificazione stack/parallela/inclinata, raccolta delle coppie equivalenti, colonne, parametri

`first_stage(old, new)` dice il primo stadio che una modifica delle soglie rende da rifare; `analysis.reanalyse`
riparte da lì e riusa il resto (`StageCache`). Muovere le soglie del criterio di stack (angolo, distanza
interplanare, sovrapposizione, slittamento) tocca solo `pairs`: nessuna geometria si rimisura.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields

from pistak.core.config import Limits
from pistak.core.pairs import Member, Pair

STAGES = ("cif", "structure", "fragments", "contacts", "pairs")
#: Le soglie lette da ciascuno stadio. Ogni campo di `Limits` deve stare in uno stadio (lo controlla un test).
STAGE_FIELDS: dict[str, frozenset[str]] = {
    "cif": frozenset({"disorder"}),
    "structure": frozenset({"bond_tolerance", "same_site", "max_cell_atoms"}),
    "fragments": frozenset({"ring_rms", "molecule_rms", "max_ring_size", "fragments", "min_fragment_atoms", "long_axis_atoms"}),
    "contacts": frozenset({"contact_cutoff", "neighbour_shells", "pair_mode"}),
    "pairs": frozenset({"max_dihedral", "max_interplanar", "min_overlap", "max_slip", "overlap_measure",
                        "ring_pair_max_distance"}),
}


def stage_of(name: str) -> str:
    """Lo stadio che legge la soglia `name`."""
    for stage in STAGES:
        if name in STAGE_FIELDS[stage]:
            return stage
    raise KeyError(name)


def changed_fields(old: Limits, new: Limits) -> list[str]:
    return [f.name for f in fields(Limits) if getattr(old, f.name) != getattr(new, f.name)]


def first_stage(old: Limits, new: Limits) -> str | None:
    """Il primo stadio da rifare passando dalle soglie `old` a `new` (None se sono uguali)."""
    changed = {stage_of(name) for name in changed_fields(old, new)}
    return next((s for s in STAGES if s in changed), None)


def is_after(stage: str, other: str) -> bool:
    return STAGES.index(stage) > STAGES.index(other)


@dataclass
class StageCache:
    """Ciò che un'analisi conserva per essere rifatta in parte. Le misure non cambiano con le soglie dello stadio
    `pairs`, quindi si riusano; ogni rianalisi riparte da copie fresche delle coppie."""

    #: Il testo del CIF e il blocco analizzato (servono solo per rifare lo stadio `cif`).
    text: str | None = None
    block: str | int | None = None
    #: Tutte le coppie misurate, una per contatto, con il tipo secondo le soglie con cui sono state misurate.
    contacts: list[Pair] = field(default_factory=list)
    #: Chiave di simmetria di ciascuna (modalità canonical; None nel modo del sito, che raccoglie per valori).
    keys: list[Member] | None = None
    #: Valori dei parametri «geometrici» del registro per coppia misurata (`Pair.key` → chiave → valore).
    values: dict[Member, dict[str, float]] = field(default_factory=dict)

    def derived(self) -> StageCache:
        """La cache per un'analisi che riusa le stesse misure (copia dei valori: le due analisi restano indipendenti)."""
        return StageCache(self.text, self.block, self.contacts, self.keys, dict(self.values))
