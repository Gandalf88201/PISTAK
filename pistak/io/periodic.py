"""La cella intera come input per calcoli periodici: VASP (POSCAR, INCAR, KPOINTS), Quantum ESPRESSO (pw.x), CP2K e
CRYSTAL.

Gli atomi sono quelli della cella ricostruita dall'analisi (stessa regola di disordine: una sola componente per sito,
le occupanze parziali non si possono dare a un codice periodico), in coordinate frazionarie riportate nella cella.
Gli idrogeni, come per il dimero (`pistak.io.qc`), si portano per impostazione alle distanze X–H da neutroni lungo il
legame; gli atomi pesanti restano dove li mette il cristallo. CRYSTAL riceve la cella in P1 (tutti gli atomi), così
non dipende dalle convenzioni di origine dei gruppi spaziali.

Punti k: una griglia centrata in Γ con n_i = ⌈L / a_i⌉ (L = `kpoint_length`, 20 Å per impostazione, adatto agli
isolanti molecolari); CP2K usa solo il punto Γ, e conviene una supercella se un lato è più corto di L/2.

Le parole chiave seguono le versioni correnti dei programmi (VASP 6, QE 7, CP2K 2024, CRYSTAL23): controllale con la
tua. Gli pseudopotenziali (VASP, QE) e le basi (CP2K, CRYSTAL) vanno scelti e copiati a parte: i file li nominano
con valori segnaposto commentati. Come gli altri export, i file sono ASCII con commenti in inglese.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pistak.core.analysis import Analysis
from pistak.core.elements import atomic_mass, atomic_number
from pistak.core.structure import to_frac
from pistak.io.qc import QCError, ascii_text, normalize_hydrogens

PROGRAMS = ("vasp", "qe", "cp2k", "crystal")
TASKS = ("scf", "relax", "vc-relax")
DISPERSIONS = ("d3bj", "d3zero", "d4", "none")
#: Estensione del file per i programmi con un solo file (VASP ne scrive tre in una cartella).
SUFFIX = {"qe": "pw.in", "cp2k": "inp", "crystal": "d12"}

_SUPPORTED_DISP = {"vasp": DISPERSIONS, "qe": ("d3bj", "d3zero", "none"), "cp2k": DISPERSIONS,
                   "crystal": ("d3bj", "d3zero", "none")}
_VASP_GGA = {"PBE": "PE", "PBESOL": "PS", "RPBE": "RP", "REVPBE": "RE"}
_VASP_IVDW = {"d3zero": 11, "d3bj": 12, "d4": 13}
#: Come CRYSTAL chiama i funzionali GGA (gli ibridi, B3LYP, PBE0, HSE06…, hanno già il nome giusto).
_CRYSTAL_DFT = {"PBE": "PBEXC", "PBESOL": "PBESOLXC", "BLYP": "BLYP"}


@dataclass
class PeriodicOptions:
    program: str = "vasp"  # vasp | qe | cp2k | crystal
    task: str = "scf"  # scf | relax (solo atomi) | vc-relax (atomi e cella)
    #: Funzionale: PBE per i codici a onde piane; per CRYSTAL anche ibridi (B3LYP, PBE0…).
    functional: str = "PBE"
    dispersion: str = "d3bj"
    #: Lunghezza (Å) per la griglia di punti k: n_i = ⌈L / a_i⌉.
    kpoint_length: float = 20.0
    normalize_h: bool = True
    #: Energia di taglio delle onde piane: eV per VASP, Ry per QE (0 = valore predefinito del programma qui sotto).
    cutoff: float = 0.0
    nprocs: int = 4

    def check(self) -> None:
        if self.program not in PROGRAMS:
            raise QCError("qcOptions", f"program={self.program}")
        if self.task not in TASKS:
            raise QCError("qcOptions", f"task={self.task}")
        if self.dispersion not in DISPERSIONS:
            raise QCError("qcOptions", f"dispersion={self.dispersion}")
        if self.dispersion not in _SUPPORTED_DISP[self.program]:
            raise QCError("qcUnsupported", f"{self.dispersion}: {self.program}")
        if self.program == "vasp" and self.functional.upper() not in _VASP_GGA:
            raise QCError("qcUnsupported", f"{self.functional}: vasp ({', '.join(_VASP_GGA)})")
        if self.kpoint_length <= 0:
            raise QCError("qcOptions", f"kpoint_length={self.kpoint_length}")


@dataclass
class PeriodicCell:
    vectors: np.ndarray  # righe a, b, c (Å)
    el: list[str]
    frac: np.ndarray  # N×3, nella cella [0, 1)
    notes: list[str]

    @property
    def lengths(self) -> np.ndarray:
        return np.linalg.norm(self.vectors, axis=1)

    def species(self) -> list[str]:
        """Gli elementi nell'ordine della prima comparsa, dopo l'ordinamento per elemento (C, H, poi alfabetico)."""
        return sorted(set(self.el), key=lambda e: (e != "C", e != "H", e))

    def kpoints(self, length: float) -> list[int]:
        return [max(1, math.ceil(length / x)) for x in self.lengths]


def periodic_cell(analysis: Analysis, options: PeriodicOptions | None = None) -> PeriodicCell:
    """Gli atomi della cella, ordinati per elemento, con gli idrogeni normalizzati se richiesto."""
    options = options or PeriodicOptions()
    s = analysis.structure
    cart = s.cell_frac @ s.vectors
    moved = 0
    if options.normalize_h:
        for m in s.molecules:  # la molecola ha posizioni contigue (non spezzate dalla cella): si normalizza lì
            pos, n = normalize_hydrogens(m.el, m.pos, m.bonds)
            moved += n
            delta = pos - m.pos
            cart[m.cell_index] += delta
    frac = to_frac(s.vectors, cart)
    frac -= np.floor(frac + 1e-9)
    order = sorted(range(len(s.cell_el)), key=lambda i: (s.cell_el[i] != "C", s.cell_el[i] != "H", s.cell_el[i], i))
    el = [s.cell_el[i] for i in order]
    c = analysis.cif.cell
    notes = [
        f"PISTAK: {analysis.name} [{analysis.cif.name}] space group {analysis.cif.space_group or '?'}, written in P1",
        f"cell a={c.a:.4f} b={c.b:.4f} c={c.c:.4f} alpha={c.alpha:.3f} beta={c.beta:.3f} gamma={c.gamma:.3f}, {len(el)} atoms",
        f"X-H distances set to neutron values ({moved} hydrogens moved)" if options.normalize_h
        else "X-H distances as in the CIF (X-ray: usually too short)",
        "experimental crystal geometry, not optimised",
    ]
    if analysis.cif.skipped_disorder:
        notes.append(f"disorder: {analysis.cif.skipped_disorder} atoms of minor components left out (policy {analysis.cif.disorder})")
    partial = sum(at.occ < 0.99 for at in analysis.cif.atoms)
    if partial:
        notes.append(f"WARNING: {partial} atoms with partial occupancy written as full atoms")
    if "H" not in el and "C" in el:
        notes.append("WARNING: no hydrogens in the CIF: add them before the calculation")
    return PeriodicCell(s.vectors.copy(), el, frac[order], notes)


# ---- scrittori ---------------------------------------------------------------------------------------------
def _vec(v: np.ndarray) -> str:
    return f"{v[0]:16.10f} {v[1]:16.10f} {v[2]:16.10f}"


def _counts(cell: PeriodicCell) -> list[tuple[str, int]]:
    return [(e, cell.el.count(e)) for e in cell.species()]


def _vasp(cell: PeriodicCell, o: PeriodicOptions, name: str) -> dict[str, str]:
    poscar = [cell.notes[0], "1.0", *(_vec(v) for v in cell.vectors)]
    counts = _counts(cell)
    poscar += ["  " + "  ".join(e for e, _ in counts), "  " + "  ".join(str(n) for _, n in counts), "Direct"]
    poscar += [f"{_vec(f)}  {e}" for e, f in zip(cell.el, cell.frac)]
    incar = [f"# {n}" for n in cell.notes]
    incar += [f"SYSTEM = {name}", "PREC = Accurate",
              f"ENCUT = {o.cutoff or 600:g}", f"GGA = {_VASP_GGA[o.functional.upper()]}", "EDIFF = 1E-6",
              "ISMEAR = 0", "SIGMA = 0.05", "LREAL = Auto", "LWAVE = .FALSE.", "LCHARG = .FALSE."]
    if o.dispersion != "none":
        incar.append(f"IVDW = {_VASP_IVDW[o.dispersion]}")
    if o.task != "scf":
        incar += ["IBRION = 2", f"ISIF = {3 if o.task == 'vc-relax' else 2}", "NSW = 200", "EDIFFG = -0.01"]
    incar.append("# POTCAR: concatenate the PAW potentials in the order of POSCAR line 6 (" + " ".join(e for e, _ in counts) + ")")
    k = cell.kpoints(o.kpoint_length)
    kpoints = [f"Gamma-centred grid from k-point length {o.kpoint_length:g} A", "0", "Gamma", " ".join(map(str, k)), "0 0 0"]
    return {"POSCAR": "\n".join(poscar) + "\n", "INCAR": "\n".join(incar) + "\n", "KPOINTS": "\n".join(kpoints) + "\n"}


def _qe(cell: PeriodicCell, o: PeriodicOptions, prefix: str) -> str:
    calc = {"scf": "scf", "relax": "relax", "vc-relax": "vc-relax"}[o.task]
    ecut = o.cutoff or 60
    lines = ["&CONTROL", f"  calculation = '{calc}'", f"  prefix = '{prefix}'", "  outdir = './tmp'", "  pseudo_dir = './pseudo'",
              "  tprnfor = .true.", "  tstress = .true.", "/", "&SYSTEM", "  ibrav = 0", f"  nat = {len(cell.el)}",
              f"  ntyp = {len(cell.species())}", f"  ecutwfc = {ecut:g}", f"  ecutrho = {8 * ecut:g}",
              f"  input_dft = '{o.functional.lower()}'"]
    if o.dispersion != "none":
        lines += ["  vdw_corr = 'grimme-d3'", f"  dftd3_version = {4 if o.dispersion == 'd3bj' else 3}"]
    lines += ["/", "&ELECTRONS", "  conv_thr = 1.0d-8", "/"]
    if o.task != "scf":
        lines += ["&IONS", "/"]
    if o.task == "vc-relax":
        lines += ["&CELL", "  cell_dofree = 'all'", "/"]
    # i commenti stanno tra le schede (pw.x salta le righe che iniziano con #), non prima delle namelist
    lines += [f"# {n}" for n in cell.notes]
    lines += ["# pseudopotentials: set the file names below (e.g. from the SSSP library) and pseudo_dir"]
    lines += ["ATOMIC_SPECIES"] + [f"  {e:<3} {atomic_mass(e):9.4f}  {e}.UPF" for e in cell.species()]
    lines += ["CELL_PARAMETERS angstrom"] + [_vec(v) for v in cell.vectors]
    lines += ["ATOMIC_POSITIONS crystal"] + [f"  {e:<3} {_vec(f)}" for e, f in zip(cell.el, cell.frac)]
    lines += ["K_POINTS automatic", " ".join(map(str, cell.kpoints(o.kpoint_length))) + " 0 0 0"]
    return "\n".join(lines) + "\n"


def _cp2k(cell: PeriodicCell, o: PeriodicOptions, prefix: str) -> str:
    run = {"scf": "ENERGY", "relax": "GEO_OPT", "vc-relax": "CELL_OPT"}[o.task]
    lines = [f"# {n}" for n in cell.notes]
    short = [n for n, x in zip("ABC", cell.lengths) if x < o.kpoint_length / 2]
    lines.append("# Gamma point only" + (f": cell edge(s) {', '.join(short)} shorter than {o.kpoint_length / 2:g} A, "
                                          "consider a supercell" if short else ""))
    lines += ["&GLOBAL", f"  PROJECT {prefix}", f"  RUN_TYPE {run}", "  PRINT_LEVEL LOW", "&END GLOBAL", "&FORCE_EVAL",
              "  METHOD Quickstep"]
    if o.task == "vc-relax":
        lines.append("  STRESS_TENSOR ANALYTICAL")
    lines += ["  &DFT", "    BASIS_SET_FILE_NAME BASIS_MOLOPT", "    POTENTIAL_FILE_NAME GTH_POTENTIALS",
              "    &MGRID", "      CUTOFF 600", "      REL_CUTOFF 60", "    &END MGRID",
              "    &SCF", "      EPS_SCF 1.0E-6", "      MAX_SCF 50", "      &OT", "        PRECONDITIONER FULL_SINGLE_INVERSE",
              "        MINIMIZER DIIS", "      &END OT", "      &OUTER_SCF", "        MAX_SCF 20", "        EPS_SCF 1.0E-6",
              "      &END OUTER_SCF", "    &END SCF", "    &XC", f"      &XC_FUNCTIONAL {o.functional.upper()}",
              "      &END XC_FUNCTIONAL"]
    if o.dispersion != "none":
        kind = {"d3bj": "DFTD3(BJ)", "d3zero": "DFTD3", "d4": "DFTD4"}[o.dispersion]
        lines += ["      &VDW_POTENTIAL", "        POTENTIAL_TYPE PAIR_POTENTIAL", "        &PAIR_POTENTIAL",
                  f"          TYPE {kind}", f"          REFERENCE_FUNCTIONAL {o.functional.upper()}"]
        if o.dispersion != "d4":
            lines.append("          PARAMETER_FILE_NAME dftd3.dat")
        lines += ["        &END PAIR_POTENTIAL", "      &END VDW_POTENTIAL"]
    lines += ["    &END XC", "  &END DFT", "  &SUBSYS", "    &CELL"]
    lines += [f"      {n} {_vec(v)}" for n, v in zip("ABC", cell.vectors)]
    lines += ["      PERIODIC XYZ", "    &END CELL", "    &COORD", "      SCALED"]
    lines += [f"      {e:<3} {_vec(f)}" for e, f in zip(cell.el, cell.frac)]
    lines += ["    &END COORD"]
    for e in cell.species():
        lines += [f"    &KIND {e}", "      BASIS_SET DZVP-MOLOPT-SR-GTH", "      POTENTIAL GTH-PBE", "    &END KIND"]
    lines += ["  &END SUBSYS", "&END FORCE_EVAL"]
    return "\n".join(lines) + "\n"


def _crystal(cell: PeriodicCell, o: PeriodicOptions, analysis: Analysis) -> tuple[str, str]:
    c = analysis.cif.cell
    lines = [cell.notes[0].replace("PISTAK: ", "")[:80], "CRYSTAL", "0 0 0", "1",
             f"{c.a:.6f} {c.b:.6f} {c.c:.6f} {c.alpha:.4f} {c.beta:.4f} {c.gamma:.4f}", str(len(cell.el))]
    # in P1 CRYSTAL costruisce la cella dai parametri con a lungo x: le coordinate frazionarie restano valide
    lines += [f"{atomic_number(e):3d} {f[0]:14.10f} {f[1]:14.10f} {f[2]:14.10f}" for e, f in zip(cell.el, cell.frac)]
    if o.task != "scf":
        lines += ["OPTGEOM"] + (["FULLOPTG"] if o.task == "vc-relax" else []) + ["END"]
    lines += ["END", "BASISSET", "POB-TZVP-REV2", "DFT", _CRYSTAL_DFT.get(o.functional.upper(), o.functional.upper()), "END"]
    if o.dispersion != "none":
        lines += ["DFTD3", "VERSION", "4" if o.dispersion == "d3bj" else "3", "END"]
    k = max(cell.kpoints(o.kpoint_length))
    lines += ["SHRINK", f"{k} {2 * k}", "TOLDEE", "8", "END"]
    head = [f"# {n}" for n in cell.notes]  # CRYSTAL non legge commenti: le note vanno in un file a parte
    return "\n".join(lines) + "\n", "\n".join(head) + "\n"


def stem(analysis: Analysis) -> str:
    return re.sub(r"[^\w.-]+", "_", Path(analysis.name).stem or "pistak")


def periodic_files(analysis: Analysis, options: PeriodicOptions | None = None) -> dict[str, str]:
    """Nome del file → contenuto (ASCII). VASP: POSCAR, INCAR, KPOINTS; gli altri: un file d'ingresso (CRYSTAL anche
    un .notes.txt con la provenienza, perché il formato non ammette commenti)."""
    options = options or PeriodicOptions()
    options.check()
    cell = periodic_cell(analysis, options)
    name = stem(analysis)
    if options.program == "vasp":
        files = _vasp(cell, options, name)
    elif options.program == "qe":
        files = {f"{name}.{SUFFIX['qe']}": _qe(cell, options, name)}
    elif options.program == "cp2k":
        files = {f"{name}.{SUFFIX['cp2k']}": _cp2k(cell, options, name)}
    else:
        body, notes = _crystal(cell, options, analysis)
        files = {f"{name}.{SUFFIX['crystal']}": body, f"{name}.notes.txt": notes}
    return {k: ascii_text(v) for k, v in files.items()}


def write_periodic(analysis: Analysis, folder: str | Path, options: PeriodicOptions | None = None) -> list[Path]:
    """Scrive i file in `folder` (creata se manca). Ritorna i percorsi scritti."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    out = []
    for name, text in periodic_files(analysis, options).items():
        p = folder / name
        p.write_text(text, encoding="utf-8")
        out.append(p)
    return out


__all__ = ["PeriodicOptions", "PeriodicCell", "periodic_cell", "periodic_files", "write_periodic", "PROGRAMS", "TASKS",
           "DISPERSIONS", "stem"]
