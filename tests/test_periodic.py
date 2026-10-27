"""Fase 4a: la cella intera come input per calcoli periodici (VASP, Quantum ESPRESSO, CP2K, CRYSTAL)."""
import re
from pathlib import Path

import numpy as np
import pytest

from pistak import api
from pistak.core.elements import atomic_number
from pistak.io import periodic
from pistak.io.qc import XH_NEUTRON, QCError

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
NUM = r"(-?\d+\.\d+)"


@pytest.fixture(scope="module")
def herring():
    return api.analyse_file(EX / "herringbone_anthracene.cif")


def test_cella_con_idrogeni_normalizzati(herring):
    s = herring.structure
    cell = periodic.periodic_cell(herring)
    assert len(cell.el) == s.n_cell_atoms and sorted(cell.el) == sorted(s.cell_el)
    assert cell.el == sorted(cell.el, key=lambda e: (e != "C", e != "H", e))  # raggruppati per elemento
    assert ((cell.frac >= 0) & (cell.frac < 1)).all()
    cart = cell.frac @ cell.vectors
    shifts = np.array([(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)]) @ cell.vectors
    h = [i for i, e in enumerate(cell.el) if e == "H"]
    c = [i for i, e in enumerate(cell.el) if e == "C"]
    for i in h:  # ogni H a 1,083 Å dal carbonio più vicino (immagini periodiche comprese)
        d = np.linalg.norm(cart[c][:, None, :] + shifts[None] - cart[i], axis=2).min()
        assert d == pytest.approx(XH_NEUTRON["C"], abs=1e-6)
    raw = periodic.periodic_cell(herring, periodic.PeriodicOptions(normalize_h=False))
    assert any("as in the CIF" in n for n in raw.notes)


def test_vasp(herring):
    files = periodic.periodic_files(herring)
    assert set(files) == {"POSCAR", "INCAR", "KPOINTS"}
    poscar = files["POSCAR"].splitlines()
    species, counts = poscar[5].split(), [int(x) for x in poscar[6].split()]
    assert species == ["C", "H"] and sum(counts) == herring.structure.n_cell_atoms and poscar[7] == "Direct"
    assert "IVDW = 12" in files["INCAR"] and "ISIF" not in files["INCAR"] and "GGA = PE" in files["INCAR"]
    relax = periodic.periodic_files(herring, periodic.PeriodicOptions(task="vc-relax", dispersion="none"))
    assert "ISIF = 3" in relax["INCAR"] and "IVDW" not in relax["INCAR"]
    k = [int(x) for x in files["KPOINTS"].splitlines()[3].split()]
    assert k == [max(1, int(np.ceil(20 / x))) for x in np.linalg.norm(herring.structure.vectors, axis=1)]


def test_quantum_espresso(herring):
    text = next(iter(periodic.periodic_files(herring, periodic.PeriodicOptions("qe", "relax")).values()))
    assert f"nat = {herring.structure.n_cell_atoms}" in text and "calculation = 'relax'" in text and "&IONS" in text
    assert "dftd3_version = 4" in text
    assert text.startswith("&CONTROL")  # niente commenti prima delle namelist
    assert len(re.findall(rf"^\s+[A-Z][a-z]?\s+{NUM}\s+{NUM}\s+{NUM}\s*$", text.split("ATOMIC_POSITIONS")[1], re.M)) == \
        herring.structure.n_cell_atoms


def test_cp2k_e_crystal(herring):
    cp2k = next(iter(periodic.periodic_files(herring, periodic.PeriodicOptions("cp2k", "vc-relax", dispersion="d4")).values()))
    assert "RUN_TYPE CELL_OPT" in cp2k and "TYPE DFTD4" in cp2k and "&KIND H" in cp2k and "SCALED" in cp2k
    files = periodic.periodic_files(herring, periodic.PeriodicOptions("crystal", "relax", functional="B3LYP"))
    d12 = next(v for k, v in files.items() if k.endswith(".d12")).splitlines()
    assert d12[1:4] == ["CRYSTAL", "0 0 0", "1"] and int(d12[5]) == herring.structure.n_cell_atoms
    assert {int(x.split()[0]) for x in d12[6:6 + int(d12[5])]} == {atomic_number("C"), atomic_number("H")}
    assert "OPTGEOM" in d12 and "B3LYP" in d12 and "PBEXC" not in d12
    pbe = periodic.periodic_files(herring, periodic.PeriodicOptions("crystal"))
    assert "PBEXC" in next(v for k, v in pbe.items() if k.endswith(".d12"))
    assert all(t.isascii() for f in (files, pbe) for t in f.values())


@pytest.mark.parametrize("options", [
    periodic.PeriodicOptions("qe", dispersion="d4"), periodic.PeriodicOptions("crystal", dispersion="d4"),
    periodic.PeriodicOptions("vasp", functional="B3LYP"), periodic.PeriodicOptions("abinit"),
    periodic.PeriodicOptions(task="md"), periodic.PeriodicOptions(kpoint_length=0),
])
def test_combinazioni_non_valide(herring, options):
    with pytest.raises(QCError):
        periodic.periodic_files(herring, options)


def test_cli(tmp_path, capsys):
    from pistak.cli import main

    assert main(["periodic", str(EX / "herringbone_anthracene.cif"), "--program", "qe"]) == 0
    assert "ATOMIC_POSITIONS crystal" in capsys.readouterr().out
    assert main(["periodic", str(EX / "herringbone_anthracene.cif"), "-o", str(tmp_path / "v")]) == 0
    assert {p.name for p in (tmp_path / "v").iterdir()} == {"POSCAR", "INCAR", "KPOINTS"}
