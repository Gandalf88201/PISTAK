"""Strutture senza coppie: il motivo (`pistak.core.diagnosis`) al posto di viste vuote."""
import math
from pathlib import Path

import numpy as np
import pytest

from cifgen import make_cif, ring_system
from pistak import api, i18n
from pistak.core.cif import Cell
from pistak.core.diagnosis import ring_planarity, why_no_pairs
from pistak.labels import NO_PAIRS_TEXT, no_pairs_text

EX = Path(__file__).resolve().parents[1] / "examples"
BOX = Cell(9.0, 9.0, 9.0, 90.0, 90.0, 90.0)


def cyclohexane() -> list:
    """Cicloesano a sedia con gli idrogeni: anello saturo, scarto dal piano ~0,25 Å."""
    atoms = []
    for k in range(6):
        a = math.radians(60 * k)
        c = np.array([1.45 * math.cos(a), 1.45 * math.sin(a), 0.25 * (-1) ** k])
        out = np.array([math.cos(a), math.sin(a), 0.0])
        atoms += [("C", c), ("H", c + 1.09 * np.array([0, 0, (-1) ** k])), ("H", c + 1.03 * out - 0.35 * np.array([0, 0, (-1) ** k]))]
    return atoms


def analyse(atoms, cell=BOX, symops=("x, y, z",)):
    return api.analyse_text(make_cif(atoms, cell, list(symops), origin_frac=(0.5, 0.5, 0.5)), "t.cif")


def test_con_coppie_nessun_motivo():
    assert why_no_pairs(api.analyse_file(EX / "herringbone_anthracene.cif")) is None


def test_struttura_reale_con_anelli_saturi():
    # COD 1550133: cinque anelli, tutti saturi (scarto dal piano 0,17–0,21 Å): nessun nucleo π
    a = api.analyse_file(Path(__file__).parent / "data" / "cod" / "1550133.cif")
    assert not a.pairs
    r = why_no_pairs(a)
    assert r.code == "nonPlanarRings" and r.values["n"] == 5
    assert all(x.saturated for x in ring_planarity(a))


def test_anelli_saturi():
    a = analyse(cyclohexane())
    assert not a.pairs
    rings = ring_planarity(a)
    assert len(rings) == 1 and rings[0].saturated is True and rings[0].rms > 0.2
    r = why_no_pairs(a)
    assert r.code == "nonPlanarRings" and r.values["n"] == 1 and r.values["limit"] == "0.10"


def test_aromatico_distorto():
    atoms, k = [], 0  # benzene a sedia appena accennata: carboni alternati ±0,17 Å (oltre le soglie di anello e nucleo)
    for e, p in ring_system(1):
        if e == "C":
            k += 1
        atoms.append((e, p + np.array([0, 0, 0.17 * (-1) ** k])))
    a = analyse(atoms)
    r = why_no_pairs(a)
    assert r.code == "nearlyPlanarRings" and ring_planarity(a)[0].saturated is False
    assert "alza la soglia" in no_pairs_text(r)


def test_nessun_anello():
    c = np.zeros(3)  # neopentano: atomi pesanti a tetraedro, niente di planare
    atoms = [("C", c)] + [("C", c + 1.53 * np.array(v) / math.sqrt(3)) for v in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
    assert why_no_pairs(analyse(atoms)).code == "noRings"


def test_nuclei_pi_troppo_lontani():
    a = analyse(ring_system(2), Cell(30.0, 30.0, 30.0, 90.0, 90.0, 90.0))  # un naftalene solo, immagini a 30 Å
    r = why_no_pairs(a)
    assert r.code == "noContacts" and r.values["n"] == 1


def test_testi_in_due_lingue():
    r = why_no_pairs(analyse(cyclohexane()))
    it = no_pairs_text(r)
    try:
        i18n.set_language("en")
        en = no_pairs_text(r)
    finally:
        i18n.set_language("it")
    assert "saturi" in it and "saturated" in en and "{" not in it + en
    assert set(NO_PAIRS_TEXT) >= {"extendedOnly", "noContacts", "noRings", "nonPlanarRings", "nearlyPlanarRings", "coreNotPlanar"}


def test_cli(tmp_path, capsys):
    from pistak.cli import main

    f = tmp_path / "c.cif"
    f.write_text(make_cif(cyclohexane(), BOX, ["x, y, z"], origin_frac=(0.5, 0.5, 0.5)), encoding="utf-8")
    assert main(["analyze", str(f)]) == 0
    assert "anelli saturi" in capsys.readouterr().out


@pytest.fixture()
def app():
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


def test_viste_vuote_spiegano_perche(app):
    from pistak.gui.session import Session
    from pistak.gui.views.detail_view import DetailView
    from pistak.gui.views.table_view import TableView

    s = Session()
    table, detail = TableView(s, "A"), DetailView(s, "A")
    s._apply(analyse(cyclohexane()))
    assert not table.empty_notice.isHidden() and table.table.isHidden()
    assert "saturi" in table.empty_notice.text() and "Geometria" in table.empty_notice.text()
    assert detail.empty_notice.isHidden() and not detail.tabs.isHidden()  # il Dettaglio resta utilizzabile
    assert "Nessuna coppia di nuclei π" in detail.summary.toPlainText()
    s._apply(api.analyse_file(EX / "herringbone_anthracene.cif"))  # con coppie torna tutto come prima
    assert table.empty_notice.isHidden() and not table.table.isHidden()
    assert "Nessuna coppia" not in detail.summary.toPlainText()
