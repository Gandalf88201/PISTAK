"""Catena completa: testo CIF → cella → molecole → nuclei π → coppie vicine."""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np

from pistak.core.cif import Cif, parse_cif
from pistak.core.config import Limits, load_limits
from pistak.core.cores import Core, pi_core
from pistak.core.pairs import Frame, Pair, collect_pairs, pair_frame
from pistak.core.stacks import StackColumn, find_stacks
from pistak.core.structure import Molecule, Structure, build_structure
from pistak.i18n import tf

MAX_BYTES = 50 * 1024 * 1024


@dataclass
class MoleculeGroup:
    formula: list[tuple[str, int]]
    count: int
    with_core: int

    @property
    def text(self) -> str:
        return "".join(f"{e}{n if n > 1 else ''}" for e, n in self.formula)


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
    #: Colonne di π-stack (molecole impilate lungo una direzione del cristallo).
    stacks: list[StackColumn] = field(default_factory=list)

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
        self.core_a: Core = analysis.cores[pair.a]  # type: ignore[assignment]
        self.core_b: Core = analysis.cores[pair.b].moved(self.shift_vector)  # type: ignore[union-attr]

    @cached_property
    def pos_a(self) -> np.ndarray:
        return self.A.pos

    @cached_property
    def pos_b(self) -> np.ndarray:
        """Posizioni degli atomi della seconda molecola, già traslata nella cella giusta."""
        return self.B.pos + self.shift_vector

    @cached_property
    def frame(self) -> Frame:
        return pair_frame(self.core_a, self.core_b)

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


def analyse_cif(cif: Cif, name: str = "", limits: Limits | None = None) -> Analysis:
    limits = limits or load_limits()
    structure = build_structure(cif, limits)
    cores = [None if m.extended else pi_core(m, limits) for m in structure.molecules]
    pairs, stack_edges = collect_pairs(structure, cores, limits)
    return Analysis(name or cif.name, cif, structure, cores, pairs, limits, stacks=find_stacks(structure, stack_edges))


def analyse_text(text: str, name: str = "", limits: Limits | None = None) -> Analysis:
    return analyse_cif(parse_cif(text), name, limits)


def read_cif_file(path: str | Path) -> str:
    p = Path(path)
    if p.stat().st_size > MAX_BYTES:
        raise ValueError(tf("file troppo grande (> {mb} MB): {name}", mb=MAX_BYTES // 2**20, name=p.name))
    return p.read_text(encoding="utf-8", errors="replace")


def analyse_file(path: str | Path, limits: Limits | None = None) -> Analysis:
    p = Path(path)
    return analyse_text(read_cif_file(p), p.name, limits)
