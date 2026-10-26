"""Il dimero di una coppia come input per calcoli quantistici: energia d'interazione DFT-D con correzione counterpoise
(ORCA, Psi4, Gaussian) o decomposizione SAPT (Psi4).

I monomeri sono le due molecole intere della coppia (per le coppie di frammenti π secondari, le molecole che li
contengono: tagliare un frammento richiederebbe di saturare le valenze), con B già nella posizione della coppia.
Le coppie intramolecolari non si esportano: i due frammenti sono legati covalentemente.

Idrogeni: per impostazione le distanze X–H si portano ai valori medi da diffrazione di neutroni (F. H. Allen,
Acta Cryst. B42, 515, 1986; Allen & Bruno, Acta Cryst. B66, 380, 2010), perché i raggi X le accorciano di ~0,1 Å
e l'energia d'interazione ne risente. Gli atomi pesanti restano dove li mette il cristallo (geometria sperimentale,
nessuna ottimizzazione).

Le parole chiave seguono le versioni correnti dei programmi (ORCA 5–6, Psi4 1.7–1.9, Gaussian 16): conviene
controllarle con quella installata.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from pistak.core.analysis import Analysis
from pistak.core.errors import PistakError
from pistak.core.pairs import Pair

#: Distanze X–H medie da diffrazione di neutroni (Å).
XH_NEUTRON = {"C": 1.083, "N": 1.009, "O": 0.983, "B": 1.19, "S": 1.338, "P": 1.42, "Si": 1.48}

PROGRAMS = ("orca", "psi4", "gaussian")
METHODS = ("dft-d", "sapt")
DISPERSIONS = ("d3bj", "d3zero", "d4", "none")
SUFFIX = {"orca": "inp", "psi4": "in", "gaussian": "gjf"}

_DISP = {  # come si scrive la dispersione in ciascun programma
    "orca": {"d3bj": "D3BJ", "d3zero": "D3ZERO", "d4": "D4", "none": ""},
    "psi4": {"d3bj": "-d3bj", "d3zero": "-d3zero", "d4": "-d4", "none": ""},
    "gaussian": {"d3bj": "EmpiricalDispersion=GD3BJ", "d3zero": "EmpiricalDispersion=GD3", "none": ""},
}


class QCError(PistakError):
    """Il dimero non si può esportare così. `code`: qcIntramolecular | qcUnsupported | qcOptions."""


@dataclass
class QCOptions:
    """Che cosa calcolare e con quale programma."""

    program: str = "orca"  # orca | psi4 | gaussian
    method: str = "dft-d"  # dft-d | sapt (solo Psi4)
    functional: str = "B3LYP"
    dispersion: str = "d3bj"  # d3bj | d3zero | d4 | none
    #: Base: vuota = def2-TZVP per DFT-D, jun-cc-pVDZ per SAPT (le scelte abituali).
    basis: str = ""
    sapt_level: str = "sapt0"  # sapt0 | ssapt0 | sapt2 | sapt2+ | sapt2+(3) | sapt2+3 …
    counterpoise: bool = True
    normalize_h: bool = True
    charges: tuple[int, int] = (0, 0)
    multiplicities: tuple[int, int] = (1, 1)
    nprocs: int = 4
    memory_gb: int = 4

    def resolved_basis(self) -> str:
        return self.basis or ("jun-cc-pVDZ" if self.method == "sapt" else "def2-TZVP")

    def check(self) -> None:
        if self.program not in PROGRAMS:
            raise QCError("qcOptions", f"program={self.program}")
        if self.method not in METHODS:
            raise QCError("qcOptions", f"method={self.method}")
        if self.dispersion not in DISPERSIONS:
            raise QCError("qcOptions", f"dispersion={self.dispersion}")
        if self.method == "sapt" and self.program != "psi4":
            raise QCError("qcUnsupported", f"SAPT: psi4 ({self.program})")
        if self.method == "dft-d" and self.dispersion not in _DISP[self.program]:
            raise QCError("qcUnsupported", f"{self.dispersion}: {self.program}")
        if any(m < 1 for m in self.multiplicities):
            raise QCError("qcOptions", f"multiplicities={self.multiplicities}")


@dataclass
class Monomer:
    el: list[str]
    pos: np.ndarray  # N×3, Å
    labels: list[str]
    charge: int = 0
    multiplicity: int = 1


@dataclass
class Dimer:
    a: Monomer
    b: Monomer
    #: Righe di provenienza e avvertenze, da mettere in testa all'input come commenti.
    notes: list[str] = field(default_factory=list)


def normalize_hydrogens(el: list[str], pos: np.ndarray, bonds: np.ndarray) -> tuple[np.ndarray, int]:
    """Sposta ogni H legato a un solo atomo X lungo la direzione X–H fino alla distanza da neutroni. Ritorna le nuove
    posizioni e quanti idrogeni sono stati spostati."""
    out = pos.copy()
    moved = 0
    neighbours: dict[int, list[int]] = {}
    for i, j in bonds:
        neighbours.setdefault(int(i), []).append(int(j))
        neighbours.setdefault(int(j), []).append(int(i))
    for h, e in enumerate(el):
        if e != "H" or len(neighbours.get(h, [])) != 1:
            continue
        x = neighbours[h][0]
        target = XH_NEUTRON.get(el[x])
        d = pos[h] - pos[x]
        n = np.linalg.norm(d)
        if target is None or n == 0:
            continue
        out[h] = pos[x] + d / n * target
        moved += 1
    return out, moved


def dimer(analysis: Analysis, pair: Pair, options: QCOptions | None = None) -> Dimer:
    """Le due molecole della coppia, pronte per l'input (B traslata, idrogeni normalizzati se richiesto)."""
    options = options or QCOptions()
    if pair.intramolecular:
        raise QCError("qcIntramolecular", f"#{pair.index + 1}")
    s = analysis.structure
    A, B = s.molecules[pair.a], s.molecules[pair.b]
    pos_a, pos_b = A.pos.copy(), B.pos + s.shift_vector(pair.shift)
    notes = [
        f"PISTAK: {analysis.name} [{analysis.cif.name}] coppia {pair.index + 1} ({pair.kind})",
        f"molecola {pair.a + 1} + molecola {pair.b + 1} traslata di [{' '.join(map(str, pair.shift))}]"
        + (f", frammenti π{pair.fa + 1}/π{pair.fb + 1}" if pair.fragmented else ""),
        f"interplanare {pair.interplanar:.3f} Å, slittamento L {pair.long_slip:.3f} Å, W {pair.lat_slip:.3f} Å, "
        f"angolo {pair.dihedral:.1f}°, contatto {pair.contact:.3f} Å",
    ]
    if options.normalize_h:
        pos_a, na = normalize_hydrogens(A.el, pos_a, A.bonds)
        pos_b, nb = normalize_hydrogens(B.el, pos_b, B.bonds)
        notes.append(f"distanze X-H portate ai valori da neutroni ({na + nb} idrogeni spostati)")
    else:
        notes.append("distanze X-H come nel CIF (raggi X: in genere troppo corte)")
    for name, m in (("A", A), ("B", B)):
        if "H" not in m.el:
            notes.append(f"ATTENZIONE: la molecola {name} non ha idrogeni nel CIF: vanno aggiunti prima del calcolo")
    notes.append("geometria sperimentale dal cristallo, nessuna ottimizzazione")
    qa, qb = options.charges
    ma, mb = options.multiplicities
    return Dimer(Monomer(list(A.el), pos_a, list(A.labels), qa, ma), Monomer(list(B.el), pos_b, list(B.labels), qb, mb), notes)


# ---- scrittori ---------------------------------------------------------------------------------------------
def _xyz(el: str, p: np.ndarray, prefix: str = "", suffix: str = "") -> str:
    return f"{prefix}{el}{suffix:<6} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}"


def _total(d: Dimer) -> tuple[int, int]:
    """Carica e molteplicità del dimero (accoppiamento ad alto spin degli spin dei monomeri)."""
    return d.a.charge + d.b.charge, d.a.multiplicity + d.b.multiplicity - 1


def _orca(d: Dimer, o: QCOptions) -> str:
    disp = _DISP["orca"][o.dispersion]
    keys = " ".join(x for x in (o.functional, disp, o.resolved_basis(), "TightSCF", "DefGrid2") if x)
    cp = "# DFT-D con correzione counterpoise (Boys-Bernardi): E_int = E(AB) - E(A in base AB) - E(B in base AB)"
    head = [f"# {n}" for n in d.notes] + [cp if o.counterpoise else "# DFT-D senza correzione counterpoise"]
    block = [f"%pal nprocs {o.nprocs} end", f"%maxcore {o.memory_gb * 1000 // max(1, o.nprocs)}"]

    def job(title: str, charge: int, mult: int, real_a: bool, real_b: bool, ghosts: bool) -> list[str]:
        lines = [f"# --- {title}", f"! {keys}", *block, f"* xyz {charge} {mult}"]
        for m, real in ((d.a, real_a), (d.b, real_b)):
            for el, p in zip(m.el, m.pos):
                if real:
                    lines.append(_xyz(el, p))
                elif ghosts:
                    lines.append(_xyz(el, p, suffix=" :"))  # atomo fantasma: solo funzioni di base
        lines.append("*")
        return lines

    q, mult = _total(d)
    out = head + job("dimero AB", q, mult, True, True, False)
    if o.counterpoise:
        out += ["", "$new_job"] + job("monomero A nella base del dimero (B fantasma)", d.a.charge, d.a.multiplicity, True, False, True)
        out += ["", "$new_job"] + job("monomero B nella base del dimero (A fantasma)", d.b.charge, d.b.multiplicity, False, True, True)
    return "\n".join(out) + "\n"


def _psi4(d: Dimer, o: QCOptions) -> str:
    lines = [f"# {n}" for n in d.notes]
    lines += [f"memory {o.memory_gb} GB", f"set_num_threads({o.nprocs})", "", "molecule dimer {"]
    for k, m in enumerate((d.a, d.b)):
        if k:
            lines.append("--")
        lines.append(f"{m.charge} {m.multiplicity}")
        lines += [_xyz(el, p, prefix="  ") for el, p in zip(m.el, m.pos)]
    lines += ["units angstrom", "no_reorient", "no_com", "symmetry c1", "}", "", "set {"]
    lines += [f"  basis {o.resolved_basis().lower()}", "  scf_type df"]
    if o.method == "sapt":
        lines += ["  freeze_core true", "}", "", "# SAPT: componenti elettrostatica, scambio, induzione, dispersione",
                  f"energy('{o.sapt_level.lower()}')"]
    else:
        name = f"{o.functional.lower()}{_DISP['psi4'][o.dispersion]}"
        lines += ["  dft_radial_points 99", "  dft_spherical_points 590", "}", ""]
        lines.append(f"energy('{name}', bsse_type='cp')" if o.counterpoise else f"energy('{name}')")
        lines.append("# con bsse_type='cp' Psi4 calcola dimero e monomeri nella base del dimero e stampa E_int (CP)")
    return "\n".join(lines) + "\n"


def _gaussian_basis(basis: str) -> str:
    """'def2-TZVP' → 'Def2TZVP' (così la scrive Gaussian); le altre basi restano come sono."""
    m = re.fullmatch(r"(?i)def2-?(\w+)", basis)
    return f"Def2{m.group(1).upper()}" if m else basis


def _gaussian(d: Dimer, o: QCOptions, name: str) -> str:
    disp = _DISP["gaussian"][o.dispersion]
    route = f"#P {o.functional}/{_gaussian_basis(o.resolved_basis())} {disp} SCF=Tight Integral=UltraFine".replace("  ", " ")
    if o.counterpoise:
        route += " Counterpoise=2"
    q, mult = _total(d)
    charge_line = f"{q} {mult}" + (f" {d.a.charge} {d.a.multiplicity} {d.b.charge} {d.b.multiplicity}" if o.counterpoise else "")
    lines = [f"! {n}" for n in d.notes]
    lines += [f"%nprocshared={o.nprocs}", f"%mem={o.memory_gb}GB", f"%chk={name}.chk", route.strip(), "",
              d.notes[0].replace("PISTAK: ", ""), "", charge_line]
    for k, m in enumerate((d.a, d.b), start=1):
        tag = f"(Fragment={k})" if o.counterpoise else ""
        lines += [_xyz(el, p, suffix=tag) for el, p in zip(m.el, m.pos)]
    lines.append("")
    return "\n".join(lines) + "\n"


_ASCII = {"Å": "A", "°": " deg", "π": "pi", "–": "-", "—": "-", "→": "->", "à": "a'", "è": "e'", "é": "e'",
          "ì": "i'", "ò": "o'", "ù": "u'", "·": "*"}


def ascii_text(text: str) -> str:
    """Solo caratteri ASCII: i programmi di chimica quantistica (Gaussian in particolare) non accettano altro."""
    for k, v in _ASCII.items():
        text = text.replace(k, v)
    return text.encode("ascii", "replace").decode("ascii")


def dimer_input(analysis: Analysis, pair: Pair, options: QCOptions | None = None) -> str:
    """Il testo dell'input per il programma scelto (solo ASCII)."""
    options = options or QCOptions()
    options.check()
    d = dimer(analysis, pair, options)
    if options.program == "orca":
        text = _orca(d, options)
    elif options.program == "psi4":
        text = _psi4(d, options)
    else:
        text = _gaussian(d, options, suggested_name(analysis, pair, options).rsplit(".", 1)[0])
    return ascii_text(text)


def suggested_name(analysis: Analysis, pair: Pair, options: QCOptions) -> str:
    stem = re.sub(r"[^\w.-]+", "_", Path(analysis.name).stem or "pistak")
    method = options.sapt_level if options.method == "sapt" else "dftd"
    return f"{stem}_coppia{pair.index + 1}_{method}.{SUFFIX[options.program]}"


def write_dimer_input(analysis: Analysis, pair: Pair, path: str | Path, options: QCOptions | None = None) -> Path:
    path = Path(path)
    options = options or QCOptions()
    if path.is_dir():  # una cartella: dentro, con il nome suggerito
        path = path / suggested_name(analysis, pair, options)
    path.write_text(dimer_input(analysis, pair, options), encoding="utf-8")
    return path


def exportable(pairs: list[Pair], kind: str | None = "stack") -> list[Pair]:
    """Le coppie che si possono esportare (non intramolecolari), di un tipo o tutte."""
    return [p for p in pairs if not p.intramolecular and (kind is None or p.kind == kind)]


__all__ = ["QCOptions", "QCError", "Dimer", "Monomer", "dimer", "dimer_input", "write_dimer_input", "suggested_name",
           "normalize_hydrogens", "exportable", "XH_NEUTRON", "replace"]
