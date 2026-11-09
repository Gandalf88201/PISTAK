"""Il dimero di una coppia come input per calcoli quantistici: energia d'interazione DFT-D con correzione counterpoise
(ORCA, Psi4, Gaussian) o decomposizione SAPT (Psi4).

I monomeri sono le due molecole intere della coppia (per le coppie di frammenti π secondari, le molecole che li
contengono: tagliare un frammento richiederebbe di saturare le valenze), con B già nella posizione della coppia.
Le coppie intramolecolari non si esportano: i due frammenti sono legati covalentemente.

Idrogeni: per impostazione le distanze X–H si portano ai valori medi da diffrazione di neutroni (Allen & Bruno,
Acta Cryst. B66, 380, 2010; il C–H secondo l'ibridazione del carbonio, vedi `XH_NEUTRON`), perché i raggi X le accorciano di ~0,1 Å
e l'energia d'interazione ne risente. Gli atomi pesanti restano dove li mette il cristallo (geometria sperimentale,
nessuna ottimizzazione).

Le parole chiave seguono le versioni correnti dei programmi (ORCA 5–6, Psi4 1.7–1.9, Gaussian 16): conviene
controllarle con quella installata.

Come gli altri export (`pistak.io.export`), i file hanno commenti e nomi fissi in inglese, uguali in ogni lingua
dell'interfaccia.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from pistak.core.analysis import Analysis
from pistak.core.errors import PistakError
from pistak.core.pairs import Pair

#: Distanze X–H medie da diffrazione di neutroni (Å), da Allen & Bruno, Acta Cryst. B66, 380 (2010). C, N, O: le
#: medie a T ≤ 60 K indicate nelle conclusioni (tab. 4), le meno accorciate dalla librazione; il C secondo
#: l'ibridazione (`xh_kind`), N⁺ = azoto con quattro vicini. O–H è quello degli ossidrili (gli acidi carbossilici,
#: 1,018 Å, non si distinguono). Csp–H 1,055 Å e B–H terminale 1,185 Å dalla tab. 2; S–H 1,338 e Si–H 1,506 Å da una
#: sola struttura ciascuno (tab. 2, NALCYS02 e COQYUC01). P–H non c'è nell'articolo: quegli H restano dove li mette il CIF.
XH_NEUTRON = {"Csp3": 1.095, "Csp2": 1.085, "Csp": 1.055, "N": 1.020, "N+": 1.040, "O": 0.989, "B": 1.185,
              "S": 1.338, "Si": 1.506}


def xh_kind(element: str, neighbours: int) -> str:
    """La voce di `XH_NEUTRON` per un atomo X con `neighbours` vicini (H compresi). Per il carbonio: 4 vicini = sp³,
    3 = sp² o aromatico, 2 = sp; per l'azoto 4 vicini = N⁺."""
    if element == "C":
        return {4: "Csp3", 3: "Csp2", 2: "Csp"}.get(neighbours, "")
    if element == "N":
        return "N+" if neighbours == 4 else "N"
    return element

PROGRAMS = ("orca", "psi4", "gaussian")
METHODS = ("dft-d", "sapt", "r2scan-3c")
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
    method: str = "dft-d"  # dft-d | sapt (solo Psi4) | r2scan-3c (solo ORCA: composito con base, D4 e gCP propri)
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
    #: Molteplicità del dimero: 0 = accoppiamento ad alto spin dei monomeri (m_A + m_B − 1). Per una coppia di radicali
    #: va scelta (singoletto o tripletto per due doppietti): l'energia d'interazione dipende da questa scelta.
    dimer_multiplicity: int = 0

    def total_multiplicity(self) -> int:
        ma, mb = self.multiplicities
        return self.dimer_multiplicity or ma + mb - 1

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
        if self.method == "r2scan-3c" and self.program != "orca":
            raise QCError("qcUnsupported", f"r2SCAN-3c: orca ({self.program})")
        if self.method == "dft-d" and self.dispersion not in _DISP[self.program]:
            raise QCError("qcUnsupported", f"{self.dispersion}: {self.program}")
        if any(m < 1 for m in self.multiplicities):
            raise QCError("qcOptions", f"multiplicities={self.multiplicities}")
        ma, mb = self.multiplicities
        m = self.dimer_multiplicity
        # spin totale da |S_A − S_B| a S_A + S_B a passi di 1: molteplicità da |m_A − m_B| + 1 a m_A + m_B − 1, di 2 in 2
        if m and not (abs(ma - mb) + 1 <= m <= ma + mb - 1 and (m - abs(ma - mb) - 1) % 2 == 0):
            raise QCError("qcOptions", f"dimer_multiplicity={m} (monomers {ma},{mb})")


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
    #: Molteplicità del dimero (0 = alto spin).
    multiplicity: int = 0


def normalize_hydrogens(el: list[str], pos: np.ndarray, bonds: np.ndarray) -> tuple[np.ndarray, int]:
    """Sposta ogni H legato a un solo atomo X lungo la direzione X–H fino alla distanza da neutroni. Ritorna le nuove
    posizioni e quanti idrogeni sono stati spostati. Restano dove sono gli H a ponte o senza vicini e quelli legati a un
    atomo senza distanza da neutroni (P, metalli): `unmoved_hydrogens` li conta."""
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
        target = XH_NEUTRON.get(xh_kind(el[x], len(neighbours[x])))
        d = pos[h] - pos[x]
        n = np.linalg.norm(d)
        if target is None or n == 0:
            continue
        out[h] = pos[x] + d / n * target
        moved += 1
    return out, moved


def unmoved_hydrogens(el: list[str], moved: int) -> int:
    return sum(e == "H" for e in el) - moved


def dimer(analysis: Analysis, pair: Pair, options: QCOptions | None = None) -> Dimer:
    """Le due molecole della coppia, pronte per l'input (B traslata, idrogeni normalizzati se richiesto)."""
    options = options or QCOptions()
    if pair.intramolecular:
        raise QCError("qcIntramolecular", f"#{pair.index + 1}")
    s = analysis.structure
    A, B = s.molecules[pair.a], s.molecules[pair.b]
    pos_a, pos_b = A.pos.copy(), B.pos + s.shift_vector(pair.shift)
    notes = [
        f"PISTAK: {analysis.name} [{analysis.cif.name}] pair {pair.index + 1} ({pair.kind})",
        f"molecule {pair.a + 1} + molecule {pair.b + 1} translated by [{' '.join(map(str, pair.shift))}]"
        + (f", fragments π{pair.fa + 1}/π{pair.fb + 1}" if pair.fragmented else ""),
        f"interplanar {pair.interplanar:.3f} Å, slip L {pair.long_slip:.3f} Å, W {pair.lat_slip:.3f} Å, "
        f"angle {pair.dihedral:.1f}°, contact {pair.contact:.3f} Å",
    ]
    if options.normalize_h:
        pos_a, na = normalize_hydrogens(A.el, pos_a, A.bonds)
        pos_b, nb = normalize_hydrogens(B.el, pos_b, B.bonds)
        left = unmoved_hydrogens(A.el, na) + unmoved_hydrogens(B.el, nb)
        notes.append(f"X-H distances set to neutron values ({na + nb} hydrogens moved"
                     + (f"; {left} left as in the CIF: bridging, unbonded or bonded to P/metals)" if left else ")"))
    else:
        notes.append("X-H distances as in the CIF (X-ray: usually too short)")
    for name, m in (("A", A), ("B", B)):
        if "H" not in m.el:
            notes.append(f"WARNING: molecule {name} has no hydrogens in the CIF: add them before the calculation")
    notes.append("experimental crystal geometry, not optimised")
    # revisione, voce 53: avvertenze sulla qualità del modello, che l'energia d'interazione eredita
    cif = analysis.cif
    partial = [at.label for at in cif.atoms if at.occ < 0.99]
    if partial or cif.skipped_disorder:
        notes.append(f"WARNING: disordered model ({len(partial)} atoms with partial occupancy, disorder rule "
                     f"'{cif.disorder}'): the dimer is one ordered component chosen by PISTAK, check it")
    if any(k.startswith(("_platon_squeeze", "_smtbx_masks")) for k in cif.items):
        notes.append("WARNING: the CIF reports solvent masking (SQUEEZE/masks): solvent molecules near the pair "
                     "are missing from the model")
    notes.append("input generated by PISTAK and checked for syntax only: inspect it before running")
    notes.append("monomers are the two whole molecules of the pair" + (" (the pair is between pi fragments of them)"
                                                                        if pair.fragmented else ""))
    qa, qb = options.charges
    ma, mb = options.multiplicities
    m = options.total_multiplicity()
    notes.append(f"charge {qa + qb}, multiplicity {m}" + (" (high-spin coupling of the monomers: set it explicitly for "
                                                        "radical pairs)" if not options.dimer_multiplicity and m > 1 else ""))
    return Dimer(Monomer(list(A.el), pos_a, list(A.labels), qa, ma), Monomer(list(B.el), pos_b, list(B.labels), qb, mb),
                 notes, m)


# ---- scrittori ---------------------------------------------------------------------------------------------
def _xyz(el: str, p: np.ndarray, prefix: str = "", suffix: str = "") -> str:
    return f"{prefix}{el}{suffix:<6} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}"


def _total(d: Dimer) -> tuple[int, int]:
    """Carica e molteplicità del dimero (per impostazione accoppiamento ad alto spin degli spin dei monomeri)."""
    return d.a.charge + d.b.charge, d.multiplicity or d.a.multiplicity + d.b.multiplicity - 1


def _orca(d: Dimer, o: QCOptions) -> str:
    if o.method == "r2scan-3c":
        # composito: base def2-mTZVPP, D4 e gCP inclusi; il gCP corregge già il BSSE, quindi niente counterpoise:
        # E_int = E(AB) - E(A) - E(B), ciascuno nella propria base
        keys = "r2SCAN-3c TightSCF DefGrid2"
        head = [f"# {n}" for n in d.notes] + ["# r2SCAN-3c (D4 and gCP included): E_int = E(AB) - E(A) - E(B), "
                                              "monomers in their own basis (gCP replaces the counterpoise correction)"]
    else:
        disp = _DISP["orca"][o.dispersion]
        keys = " ".join(x for x in (o.functional, disp, o.resolved_basis(), "TightSCF", "DefGrid2") if x)
        cp = "# DFT-D with counterpoise correction (Boys-Bernardi): E_int = E(AB) - E(A in AB basis) - E(B in AB basis)"
        head = [f"# {n}" for n in d.notes] + [cp if o.counterpoise else "# DFT-D without counterpoise correction"]
    # %maxcore è in MB per processo e ORCA a volte lo supera: il 75 % della memoria per processo, come consigliano i
    # centri di calcolo
    block = [f"%pal nprocs {o.nprocs} end", f"%maxcore {o.memory_gb * 750 // max(1, o.nprocs)}"]

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
    out = head + job("dimer AB", q, mult, True, True, False)
    if o.method == "r2scan-3c":
        out += ["", "$new_job"] + job("monomer A", d.a.charge, d.a.multiplicity, True, False, False)
        out += ["", "$new_job"] + job("monomer B", d.b.charge, d.b.multiplicity, False, True, False)
    elif o.counterpoise:
        out += ["", "$new_job"] + job("monomer A in the dimer basis (B ghost)", d.a.charge, d.a.multiplicity, True, False, True)
        out += ["", "$new_job"] + job("monomer B in the dimer basis (A ghost)", d.b.charge, d.b.multiplicity, False, True, True)
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
        lines += ["  freeze_core true", "}", "", "# SAPT: electrostatics, exchange, induction, dispersion",
                  f"energy('{o.sapt_level.lower()}')"]
    else:
        name = f"{o.functional.lower()}{_DISP['psi4'][o.dispersion]}"
        lines += ["  dft_radial_points 99", "  dft_spherical_points 590", "}", ""]
        lines.append(f"energy('{name}', bsse_type='cp')" if o.counterpoise else f"energy('{name}')")
        lines.append("# with bsse_type='cp' Psi4 computes dimer and monomers in the dimer basis and prints E_int (CP)")
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
    return f"{stem}_pair{pair.index + 1}_{method}.{SUFFIX[options.program]}"


def write_dimer_input(analysis: Analysis, pair: Pair, path: str | Path, options: QCOptions | None = None) -> Path:
    path = Path(path)
    options = options or QCOptions()
    if path.is_dir():  # una cartella: dentro, con il nome suggerito
        path = path / suggested_name(analysis, pair, options)
    path.write_text(dimer_input(analysis, pair, options), encoding="utf-8")
    return path


def exportable(pairs: list[Pair], kind: str | None = "stack") -> list[Pair]:
    """Le coppie che si possono esportare (non intramolecolari, non di polimeri: un'unità di una rete estesa non è
    una molecola da calcolare da sola), di un tipo o tutte."""
    return [p for p in pairs if not p.intramolecular and not p.polymer and (kind is None or p.kind == kind)]


__all__ = ["QCOptions", "QCError", "Dimer", "Monomer", "dimer", "dimer_input", "write_dimer_input", "suggested_name",
           "normalize_hydrogens", "exportable", "XH_NEUTRON", "xh_kind", "replace"]
