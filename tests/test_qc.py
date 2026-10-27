"""Export del dimero per calcoli quantistici: geometria, idrogeni normalizzati, input ORCA / Psi4 / Gaussian."""
import re
from pathlib import Path

import numpy as np
import pytest

from pistak import api
from pistak.core.config import Limits
from pistak.io import qc

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
COORD = re.compile(r"^\s*([A-Z][a-z]?)\S*\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", re.M)


@pytest.fixture(scope="module")
def herring():
    a = api.analyse_file(EX / "herringbone_anthracene.cif", Limits())
    return a, a.pairs[0]


def _coords(text: str) -> np.ndarray:
    return np.array([[float(x) for x in m.groups()[1:]] for m in COORD.finditer(text)])


def test_idrogeni_alle_distanze_da_neutroni(herring):
    a, _ = herring
    m = a.structure.molecules[0]
    pos, moved = qc.normalize_hydrogens(m.el, m.pos, m.bonds)
    assert moved == m.el.count("H")
    for i, j in m.bonds:
        i, j = (i, j) if m.el[j] == "H" else (j, i)
        if m.el[j] == "H":
            assert np.linalg.norm(pos[j] - pos[i]) == pytest.approx(qc.XH_NEUTRON[m.el[i]], abs=1e-9)
    heavy = [k for k, e in enumerate(m.el) if e != "H"]
    assert np.allclose(pos[heavy], m.pos[heavy])


def test_dimero(herring):
    a, p = herring
    d = qc.dimer(a, p)
    A, B = a.structure.molecules[p.a], a.structure.molecules[p.b]
    assert len(d.a.el) == len(A) and len(d.b.el) == len(B)
    heavy_b = [k for k, e in enumerate(B.el) if e != "H"]
    assert np.allclose(d.b.pos[heavy_b], (B.pos + a.structure.shift_vector(p.shift))[heavy_b])
    assert any("pair 1" in n for n in d.notes) and any("neutron" in n for n in d.notes)


def test_orca_con_counterpoise(herring):
    a, p = herring
    text = qc.dimer_input(a, p, qc.QCOptions("orca"))
    n = len(a.structure.molecules[p.a])
    assert text.count("$new_job") == 2 and text.count("! B3LYP D3BJ def2-TZVP") == 3
    jobs = text.split("$new_job")
    assert jobs[0].count(" :") == 0 and jobs[1].count(" :") == n and jobs[2].count(" :") == n  # atomi fantasma
    assert np.allclose(_coords(jobs[0]), np.vstack(_coords_of(qc.dimer(a, p))), atol=1e-7)
    single = qc.dimer_input(a, p, qc.QCOptions("orca", counterpoise=False))
    assert "$new_job" not in single and single.count("* xyz 0 1") == 1


def _coords_of(d):
    return [d.a.pos, d.b.pos]


def test_psi4_sapt_e_dft(herring):
    a, p = herring
    sapt = qc.dimer_input(a, p, qc.QCOptions("psi4", "sapt"))
    assert sapt.count("\n--\n") == 1 and "energy('sapt0')" in sapt and "basis jun-cc-pvdz" in sapt
    assert "no_com" in sapt and "symmetry c1" in sapt
    dft = qc.dimer_input(a, p, qc.QCOptions("psi4", "dft-d", dispersion="d3zero", basis="aug-cc-pVDZ"))
    assert "energy('b3lyp-d3zero', bsse_type='cp')" in dft and "basis aug-cc-pvdz" in dft


def test_gaussian(herring):
    a, p = herring
    text = qc.dimer_input(a, p, qc.QCOptions("gaussian"))
    n = len(a.structure.molecules[p.a])
    route = next(line for line in text.splitlines() if line.startswith("#P"))
    assert "B3LYP/Def2TZVP" in route and "EmpiricalDispersion=GD3BJ" in route and "Counterpoise=2" in route
    assert "\n0 1 0 1 0 1\n" in text
    assert text.count("(Fragment=1)") == n and text.count("(Fragment=2)") == n
    assert text.endswith("\n\n")  # Gaussian vuole una riga vuota in fondo


def test_cariche_e_molteplicita(herring):
    a, p = herring
    o = qc.QCOptions("orca", charges=(1, -1), multiplicities=(2, 2))
    text = qc.dimer_input(a, p, o)
    assert "* xyz 0 3" in text and "* xyz 1 2" in text and "* xyz -1 2" in text
    g = qc.dimer_input(a, p, qc.QCOptions("gaussian", charges=(1, -1), multiplicities=(2, 2)))
    assert "\n0 3 1 2 -1 2\n" in g


def test_solo_ascii(herring):
    a, p = herring
    for program in qc.PROGRAMS:
        text = qc.dimer_input(a, p, qc.QCOptions(program, "sapt" if program == "psi4" else "dft-d"))
        assert text.isascii()


def test_combinazioni_non_valide(herring):
    a, p = herring
    with pytest.raises(qc.QCError) as e:
        qc.dimer_input(a, p, qc.QCOptions("orca", "sapt"))
    assert e.value.code == "qcUnsupported"
    with pytest.raises(qc.QCError):
        qc.dimer_input(a, p, qc.QCOptions("gaussian", dispersion="d4"))
    with pytest.raises(qc.QCError):
        qc.dimer_input(a, p, qc.QCOptions("molpro"))


def test_coppia_intramolecolare_non_esportabile():
    from test_fragments import stacked_crystal

    a = api.analyse_text(stacked_crystal(), "fn", Limits())
    intra = next(p for p in a.pairs if p.intramolecular)
    with pytest.raises(qc.QCError) as e:
        qc.dimer_input(a, intra)
    assert e.value.code == "qcIntramolecular"
    assert intra not in qc.exportable(a.pairs, None)


def test_terminale(tmp_path, capsys):
    from pistak.cli import main

    out = tmp_path / "inputs"
    assert main(["qc", str(EX / "herringbone_anthracene.cif"), "--kind", "all", "--program", "psi4", "--method", "sapt",
                 "-o", str(out)]) == 0
    files = sorted(p.name for p in out.iterdir())
    assert files == ["herringbone_anthracene_pair1_sapt0.in", "herringbone_anthracene_pair2_sapt0.in"]
    assert main(["qc", str(EX / "demo_naphthalene.cif"), "--pair", "1", "--program", "gaussian"]) == 0
    assert "Counterpoise=2" in capsys.readouterr().out
    assert main(["qc", str(EX / "demo_naphthalene.cif"), "--program", "orca", "--method", "sapt"]) == 1


def test_commenti_in_inglese_in_ogni_lingua(herring):
    """Come gli altri export, l'input ha commenti fissi in inglese, qualunque sia la lingua dell'interfaccia."""
    from pistak import i18n

    a, p = herring
    texts = []
    for lang in ("it", "en"):
        i18n.set_language(lang)
        texts.append(qc.dimer_input(a, p, qc.QCOptions("orca")))
    i18n.set_language("it")
    assert texts[0] == texts[1] and "# PISTAK: herringbone_anthracene.cif [herringbone] pair 1 (stack)" in texts[0]
