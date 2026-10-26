"""Catena completa: testo CIF → cella → molecole → frammenti π → coppie vicine → colonne e contatti tra molecole.

La catena è divisa in stadi (`pistak.core.stages`): `reanalyse` rifà con nuove soglie solo gli stadi che quelle
soglie toccano e riusa le misure conservate nell'analisi (`Analysis.cache`)."""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from pistak.core.cif import Cif, parse_cif
from pistak.core.config import Limits, load_limits
from pistak.core.cores import Core, apply_axis, pi_fragments
from pistak.core.errors import InputError
from pistak.core.pairs import (
    Frame,
    Pair,
    Shift,
    _canon,
    as_fragments,
    contact_keys,
    frame_for,
    group_pairs,
    iter_neighbour_pairs,
    molecule_key,
)
from pistak.core.stacks import StackColumn, find_stacks
from pistak.core.stages import STAGES, StageCache, first_stage
from pistak.core.structure import Molecule, Structure, build_structure

if TYPE_CHECKING:
    from pistak.core.uncertainty import Uncertainty

MAX_BYTES = 50 * 1024 * 1024


@dataclass
class MoleculeGroup:
    formula: list[tuple[str, int]]
    count: int
    with_core: int
    #: Frammenti π per molecola (il massimo tra le molecole del gruppo).
    fragments: int = 0

    @property
    def text(self) -> str:
        return "".join(f"{e}{n if n > 1 else ''}" for e, n in self.formula)


@dataclass
class MoleculeContact:
    """Un contatto tra due molecole (o, intramolecolare, tra frammenti della stessa molecola) con tutte le sue coppie
    di frammenti. In modalità canonical i contatti equivalenti per simmetria sono uno solo (`count` per cella)."""

    a: int
    b: int
    shift: Shift
    pairs: list[Pair]
    count: int
    contact: float
    intramolecular: bool = False

    @property
    def stacks(self) -> list[Pair]:
        return [p for p in self.pairs if p.kind == "stack"]

    @property
    def best(self) -> Pair:
        """La coppia di frammenti più interessante: il primo stack, altrimenti la prima coppia (ordine della tabella)."""
        return (self.stacks or self.pairs)[0]


@dataclass
class Analysis:
    name: str
    cif: Cif
    structure: Structure
    cores: list[Core | None]
    pairs: list[Pair]
    limits: Limits
    #: Chiavi dei parametri presenti in `Pair.values`, nell'ordine di calcolo.
    param_keys: list[str] = field(default_factory=list)
    #: Errori dei parametri che hanno fallito su qualche coppia: chiave → ultimo messaggio (riempito da `compute_all`).
    param_errors: dict[str, str] = field(default_factory=dict)
    #: Colonne di π-stack (molecole impilate lungo una direzione del cristallo).
    stacks: list[StackColumn] = field(default_factory=list)
    #: Frammenti π di ogni molecola (il primo è `cores[i]`); vuoto se l'analisi è stata fatta senza frammenti.
    fragments: list[list[Core]] = field(default_factory=list)
    #: Incertezze dal Monte Carlo (`pistak.api.uncertainty`), se calcolate.
    uncertainty: Uncertainty | None = None
    #: Risultati intermedi per rifare l'analisi solo dagli stadi toccati da nuove soglie (`reanalyse`).
    cache: StageCache | None = field(default=None, repr=False, compare=False)

    def fragment(self, molecule: int, index: int = 0) -> Core | None:
        if self.fragments:
            fs = self.fragments[molecule]
            return fs[index] if index < len(fs) else None
        return self.cores[molecule] if index == 0 else None

    def molecule_contacts(self) -> list[MoleculeContact]:
        """Le coppie di frammenti raccolte per contatto tra molecole, nell'ordine della prima coppia di ciascuno."""
        groups: dict[tuple, list[Pair]] = {}
        for p in self.pairs:
            groups.setdefault(molecule_key(self.structure, p, self.limits), []).append(p)
        out = []
        for key, pairs in groups.items():
            first = pairs[0]
            intra = key[0] == "intra"
            if self.limits.pair_mode == "site":
                count = max(p.count for p in pairs)
            elif intra:  # quante molecole (equivalenti) della cella hanno questo contatto interno
                count = len({a for p in pairs for a, *_ in p.members})
            else:  # contatti tra molecole per cella: i membri delle coppie di frammenti, senza i frammenti
                count = len({_canon(a, b, sh) for p in pairs for a, b, sh, *_ in p.members})
            out.append(MoleculeContact(first.a, first.b, first.shift, pairs, count, min(p.contact for p in pairs), intra))
        return out

    def molecule_groups(self) -> list[MoleculeGroup]:
        """Molecole uguali raggruppate per formula (escluse quelle estese)."""
        groups: dict[str, MoleculeGroup] = {}
        for i, m in enumerate(self.structure.molecules):
            if m.extended:
                continue
            f = m.formula()
            key = "".join(f"{e}{n}" for e, n in f)
            g = groups.setdefault(key, MoleculeGroup(f, 0, 0))
            g.count += 1
            g.with_core += self.cores[i] is not None
            g.fragments = max(g.fragments, len(self.fragments[i]) if self.fragments else int(self.cores[i] is not None))
        return list(groups.values())

    @property
    def has_extended(self) -> bool:
        return any(m.extended for m in self.structure.molecules)

    def column_for(self, pair: Pair) -> StackColumn | None:
        """La colonna di π-stack a cui appartiene la coppia (solo per le coppie `stack`)."""
        if pair.kind != "stack":
            return None
        return next((c for c in self.stacks if any(pair is e[3] for e in c.edges)), None)

    def pair_context(self, pair: Pair) -> PairContext:
        return PairContext(self, pair)


class PairContext:
    """Tutto ciò che serve per misurare qualcosa su una coppia, calcolato una volta sola.

    I parametri del registro (`pistak.params`) ricevono un `PairContext`. Le proprietà sono "pigre":
    costano solo se qualcuno le usa.
    """

    def __init__(self, analysis: Analysis, pair: Pair):
        self.analysis = analysis
        self.pair = pair
        self.limits = analysis.limits
        self.structure = analysis.structure
        self.A: Molecule = analysis.structure.molecules[pair.a]
        self.B: Molecule = analysis.structure.molecules[pair.b]
        self.shift_vector = analysis.structure.shift_vector(pair.shift)
        #: I frammenti π della coppia (per le molecole con un solo frammento, il loro nucleo π).
        self.core_a: Core = analysis.fragment(pair.a, pair.fa)  # type: ignore[assignment]
        self.core_b: Core = analysis.fragment(pair.b, pair.fb).moved(self.shift_vector)  # type: ignore[union-attr]

    @cached_property
    def pos_a(self) -> np.ndarray:
        return self.A.pos

    @cached_property
    def pos_b(self) -> np.ndarray:
        """Posizioni degli atomi della seconda molecola, già traslata nella cella giusta."""
        return self.B.pos + self.shift_vector

    @cached_property
    def frame(self) -> Frame:
        return frame_for(self.core_a, self.core_b, self.limits)

    @cached_property
    def d(self) -> np.ndarray:
        """Vettore tra i baricentri dei due nuclei."""
        return self.core_b.centroid - self.core_a.centroid

    @cached_property
    def slip_vector(self) -> np.ndarray:
        return self.d - self.frame.n * np.dot(self.d, self.frame.n)

    @cached_property
    def proj_a(self) -> np.ndarray:
        """Atomi della prima molecola nel sistema (L, W, n) della coppia."""
        return self.frame.project(self.pos_a)

    @cached_property
    def proj_b(self) -> np.ndarray:
        return self.frame.project(self.pos_b)


def _fragments(structure: Structure, limits: Limits) -> list[list[Core]]:
    return [[] if m.extended else apply_axis(m, pi_fragments(m, limits), limits.long_axis_atoms) for m in structure.molecules]


def _contacts(structure: Structure, fragments: list[list[Core]], limits: Limits, text: str | None,
              block: str | int | None) -> StageCache:
    """Stadio `contacts`: misura tutte le coppie di frammenti vicini (la parte costosa) e le loro chiavi di simmetria."""
    contacts = list(iter_neighbour_pairs(structure, fragments, limits))
    return StageCache(text, block, contacts, contact_keys(structure, fragments, contacts, limits))


def _assemble(name: str, cif: Cif, structure: Structure, fragments: list[list[Core]], limits: Limits,
              cache: StageCache) -> Analysis:
    """Stadio `pairs`: classificazione, raccolta delle coppie equivalenti e colonne, dalle misure della cache."""
    pairs, stack_edges = group_pairs(cache.contacts, cache.keys, limits)
    return Analysis(name, cif, structure, [f[0] if f else None for f in fragments], pairs, limits,
                    stacks=find_stacks(structure, stack_edges), fragments=fragments, cache=cache)


def analyse_cif(cif: Cif, name: str = "", limits: Limits | None = None, text: str | None = None,
                block: str | int | None = None) -> Analysis:
    """Analisi di un CIF già letto. `text` e `block` (da dove viene) servono solo a `reanalyse`, per cambiare poi la
    regola di disordine."""
    limits = limits or load_limits()
    structure = build_structure(cif, limits)
    fragments = _fragments(structure, limits)
    return _assemble(name or cif.name, cif, structure, fragments, limits, _contacts(structure, fragments, limits, text, block))


def analyse_text(text: str, name: str = "", limits: Limits | None = None, block: str | int | None = None) -> Analysis:
    """Analisi di un blocco del CIF (per impostazione il primo che contiene atomi)."""
    limits = limits or load_limits()
    return analyse_cif(parse_cif(text, block=block, disorder=limits.disorder), name, limits, text, block)


def reanalyse(analysis: Analysis, limits: Limits) -> Analysis:
    """La stessa struttura con altre soglie, rifacendo solo gli stadi che le soglie cambiate toccano: con le soglie
    del criterio di stack (angolo, distanza interplanare, sovrapposizione, slittamento) nessuna geometria si
    rimisura. L'analisi di partenza resta com'è; se le soglie sono le stesse si riottiene lei. Le incertezze del
    Monte Carlo e i parametri del registro non passano alla nuova analisi (`pistak.api.reanalyse` li ricalcola)."""
    stage = first_stage(analysis.limits, limits)
    if stage is None:
        return analysis
    cache = analysis.cache
    if cache is None and stage == "pairs":
        stage = "contacts"  # nessuna misura da riusare
    cache = cache or StageCache()
    start = STAGES.index(stage)

    def redo(s: str) -> bool:
        return start <= STAGES.index(s)

    if redo("cif"):
        if cache.text is None:
            raise InputError("noSource")
        cif = parse_cif(cache.text, block=cache.block, disorder=limits.disorder)
    else:
        cif = analysis.cif
    structure = build_structure(cif, limits) if redo("structure") else analysis.structure
    fragments = _fragments(structure, limits) if redo("fragments") else (analysis.fragments or as_fragments(analysis.cores))
    new_cache = _contacts(structure, fragments, limits, cache.text, cache.block) if redo("contacts") else cache.derived()
    return _assemble(analysis.name, cif, structure, fragments, limits, new_cache)


def read_cif_file(path: str | Path) -> str:
    p = Path(path)
    if p.stat().st_size > MAX_BYTES:
        raise InputError("fileTooLarge", f"(> {MAX_BYTES // 2**20} MB): {p.name}")
    return p.read_text(encoding="utf-8", errors="replace")


def analyse_file(path: str | Path, limits: Limits | None = None, block: str | int | None = None) -> Analysis:
    p = Path(path)
    return analyse_text(read_cif_file(p), p.name, limits, block)
