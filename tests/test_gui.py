"""Prove della GUI senza schermo (piattaforma Qt "offscreen"): sincronizzazione delle viste e scena 3D."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_API", "pyside6")

from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("PySide6")
from PySide6.QtWidgets import QApplication

from pistak import api
from pistak.core.config import Limits
from pistak.gui.scene import build_scene
from pistak.gui.session import Session
from pistak.gui.views.detail_view import DetailView
from pistak.gui.views.hist_view import HistogramView
from pistak.gui.views.projection_view import ProjectionView
from pistak.gui.views.scatter_view import ScatterView
from pistak.gui.views.table_view import TableView

EX = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def world(app):
    s = Session()
    s.limits = Limits.site()  # 3 coppie (a→b e b→a separate): servono a provare selezione e filtri
    a = api.analyse_file(EX / "herringbone_anthracene.cif", s.limits)
    views = {
        "tableA": TableView(s, "A"), "scatterA": ScatterView(s, "A"), "projA": ProjectionView(s, "A"),
        "histA": HistogramView(s, "A"), "detailA": DetailView(s, "A"),
        "tableB": TableView(s, "B"), "free": ScatterView(s, None),
    }
    s._apply(a)
    return s, a, views


def test_tutte_le_viste_partono_dalla_prima_coppia(world):
    s, a, v = world
    assert len(a.pairs) == 3
    assert all(view.index == 0 for view in v.values() if view.group)


def test_la_selezione_segue_il_gruppo(world):
    s, a, v = world
    v["scatterA"].choose(2)
    assert v["tableA"].index == v["projA"].index == v["detailA"].index == 2
    assert v["tableA"].pair is a.pairs[2]
    assert v["tableB"].index == 0  # un altro gruppo non si muove
    assert v["free"].index == 0


def test_vista_non_sincronizzata_ha_selezione_propria(world):
    s, a, v = world
    v["free"].choose(1)
    assert v["free"].index == 1
    assert v["tableA"].index == 0  # e non sposta le altre
    s.select("A", 2)
    assert v["free"].index == 1  # né le segue


def test_cambio_di_gruppo_di_una_vista(world):
    s, a, v = world
    s.select("B", 1)
    v["tableA"].sync_combo.setCurrentIndex(v["tableA"].sync_combo.findData("B"))
    assert v["tableA"].group == "B" and v["tableA"].index == 1


def test_filtro_condiviso_dal_gruppo(world):
    s, a, v = world
    v["scatterA"].filter_to({1, 2})
    assert v["tableA"].proxy.rowCount() == 2
    assert v["tableB"].proxy.rowCount() == 3
    v["scatterA"].filter_to(None)
    assert v["tableA"].proxy.rowCount() == 3


def test_scena_della_coppia_stack(world):
    s, a, v = world
    pair = a.pairs[0]
    assert pair.kind == "stack"
    scene = build_scene(a.pair_context(pair), show_h=True, show_column=True)
    n = len(a.structure.molecules[0])
    assert len(scene.pos) == 2 * n and len(scene.polygons) == 2
    assert len(scene.ghost_pos) > 0  # le altre molecole della colonna
    heavy = build_scene(a.pair_context(pair), show_h=False)
    assert len(heavy.pos) < len(scene.pos) and all(e != "H" for e in heavy.el)
    # il vettore tra i centroidi è la somma di L, W e n
    steps = [p.points[1] - p.points[0] for p in scene.paths[:3]]
    d = scene.paths[3].points[1] - scene.paths[3].points[0]
    assert np.allclose(sum(steps), d)


def test_scena_fuori_schermo_con_pyvista(world):
    pv = pytest.importorskip("pyvista")
    from pistak.gui.views.view3d import Options, populate

    s, a, _ = world
    scene = build_scene(a.pair_context(a.pairs[0]), show_column=True)
    try:
        pl = pv.Plotter(off_screen=True, window_size=(300, 200))
        for style in (0, 1, 2):
            populate(pl, scene, Options(style=style))
        pl.close()
    except Exception as e:  # nessun contesto OpenGL disponibile in questo ambiente
        pytest.skip(f"rendering non disponibile: {e}")
