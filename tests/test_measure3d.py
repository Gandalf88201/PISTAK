"""Misure nella vista 3D: la scena di una misura della scheda Geometria e gli atomi scelti con il mouse."""
from pathlib import Path

import numpy as np
import pytest

from pistak import api
from pistak.core.intramolecular import molecule_geometry
from pistak.gui.scene import build_measure_scene, build_scene, nearest_atom, pick_measure, scene_frame

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
COD = ROOT / "tests" / "data" / "cod"


@pytest.fixture(scope="module")
def hydrate():
    return api.analyse_file(COD / "4500022.cif")


def test_scena_di_un_contatto(hydrate):
    c = next(g for g in molecule_geometry(hydrate.structure, hydrate.cif) if g.kind == "contact" and any(g.codes))
    sc = build_measure_scene(hydrate.structure, c)
    a, b = sc.highlight
    assert np.linalg.norm(sc.pos[a] - sc.pos[b]) == pytest.approx(c.value, abs=1e-9)  # proprio gli atomi misurati
    assert sc.names[a] == c.atoms[0] and sc.names[b].startswith(c.atoms[1] + "#")  # con il codice di simmetria
    assert sc.measure_paths[0].dashed and sc.measure_labels[0].text == f"{c.value:.3f} Å"
    assert set(sc.molecule) == {0, 1}  # la seconda molecola ha il carbonio ciano, come nelle coppie
    f = scene_frame(sc)
    assert abs(f.L @ (sc.pos[b] - sc.pos[a])) == pytest.approx(c.value)  # l'asse L lungo il contatto


def test_scena_di_un_angolo_senza_idrogeni(hydrate):
    g = next(g for g in molecule_geometry(hydrate.structure, hydrate.cif, kinds=("angle",)) if "H" in g.atoms[0])
    sc = build_measure_scene(hydrate.structure, g, show_h=False)
    assert len(sc.highlight) == 3 and sc.el[sc.highlight[0]] == "H"  # l'idrogeno misurato resta
    assert sum(e == "H" for e in sc.el) == 1
    assert sc.measure_labels[0].text.endswith("°") and not sc.measure_paths[0].dashed


def test_atomo_sotto_il_clic_e_misura():
    a = api.analyse_file(EX / "herringbone_anthracene.cif")
    sc = build_scene(a.pair_context(a.pairs[0]))
    assert len(sc.names) == len(sc.pos) and any(n.endswith("(B)") for n in sc.names)
    assert nearest_atom(sc, sc.pos[2] + np.array([0.0, 0.0, 0.5])) == 2
    assert nearest_atom(sc, sc.pos.max(axis=0) + 10) is None  # lontano da tutto: nessun atomo
    i, j, k = 0, 1, len(sc.pos) - 1
    paths, labels, text = pick_measure(sc, [i, j])
    d = np.linalg.norm(sc.pos[i] - sc.pos[j])
    assert len(paths) == 1 and labels[0].text == f"{d:.3f} Å" and f"{d:.3f} Å" in text
    paths, labels, text = pick_measure(sc, [i, j, k])
    assert len(paths) == 2 and labels[-1].text.endswith("°") and text.count("·") >= 6
    assert pick_measure(sc, []) == ([], [], "")


@pytest.fixture()
def gui(monkeypatch):
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    pytest.importorskip("PySide6")
    pytest.importorskip("pyvistaqt")
    monkeypatch.setenv("PISTAK_NO_3D", "1")  # la logica della vista senza OpenGL
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


def test_dalla_scheda_geometria_alla_vista_3d(gui, hydrate):
    from pistak.gui.session import Session
    from pistak.gui.views.detail_view import DetailView
    from pistak.gui.views.view3d import View3D

    s = Session()
    v, d = View3D(s, "A"), DetailView(s, "A")
    s._apply(hydrate)
    d.tabs.setCurrentIndex(6)
    d.geom_kind.setCurrentIndex(d.geom_kind.findData("contact"))
    d._geometry_chosen(0, 0)
    assert v.measure is d._geometry_rows[0] and v.ctx is None and v.scene.highlight
    assert "Geometria" in v.pick_info.text()
    # con Misura: clic vicino a due atomi
    a, b = v.scene.highlight
    v.pick_point(v.scene.pos[a] + np.array([0.3, 0.0, 0.0]))
    v.pick_point(v.scene.pos[b])
    assert v.picks == [a, b] and f"{np.linalg.norm(v.scene.pos[a] - v.scene.pos[b]):.3f} Å" in v.pick_info.text()
    v.pick_point(v.scene.pos.max(axis=0) + 20)  # nel vuoto: niente cambia
    assert v.picks == [a, b]
    # scegliendo una coppia si torna alla coppia, senza misure
    s.select("A", 0)
    assert v.measure is None and v.ctx is not None and v.picks == []


def test_misura_anche_senza_coppie(gui):
    from pistak.gui.session import Session
    from pistak.gui.views.view3d import View3D

    s = Session()
    v = View3D(s, "A")
    s._apply(api.analyse_file(COD / "1550133.cif"))  # anelli saturi: nessuna coppia, la vista mostra il perché
    assert not v.empty_notice.isHidden()
    bond = molecule_geometry(s.analysis.structure, s.analysis.cif, kinds=("bond",))[0]
    s.show_measure("A", bond)
    assert v.empty_notice.isHidden() and v.scene is not None and len(v.scene.highlight) == 2
    s.show_measure("B", bond)  # un altro gruppo: questa vista non cambia
    assert v.measure is bond


def test_la_misura_resta_se_la_vista_si_aggiorna(gui, hydrate):
    from pistak.gui.session import Session
    from pistak.gui.views.view3d import View3D

    s = Session()
    v = View3D(s, "A")
    s._apply(hydrate)
    c = next(g for g in molecule_geometry(hydrate.structure, hydrate.cif) if g.kind == "contact")
    s.show_measure("A", c)
    v.refresh()  # per esempio quando il pannello torna visibile
    assert v.measure is c and v.scene.highlight
    s._apply(api.analyse_file(EX / "herringbone_anthracene.cif"))  # un'altra analisi: si torna alla coppia
    assert v.measure is None and v.ctx is not None
