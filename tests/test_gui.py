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


def test_scena_misure_sugli_assi(world):
    s, a, _ = world
    pair = a.pairs[0]
    scene = build_scene(a.pair_context(pair))
    assert [p.color for p in scene.axes] == [p.color for p in scene.paths[:3]]  # L, W, n con i colori dei gradini
    texts = [lb.text for lb in scene.labels]
    assert {"L", "W", "n"} <= set(texts)
    for name, value in (("L", pair.long_slip), ("W", pair.lat_slip), ("n", pair.interplanar)):
        assert f"{name} {value:.2f} Å" in texts  # gli stessi valori della tabella


def test_vista_dall_alto_come_la_proiezione(world):
    """La camera «dall'alto» (quella di partenza) mette a schermo L a destra e W in alto, con la stessa scala sui due
    assi: ogni atomo cade dove lo disegna la Proiezione, a meno di traslazione e ingrandimento."""
    pv = pytest.importorskip("pyvista")
    from pistak.gui.views.view3d import view_directions

    s, a, _ = world
    ctx = a.pair_context(a.pairs[0])
    scene = build_scene(ctx, show_h=False)
    direction, up = view_directions(ctx.frame)["top"]
    try:
        pl = pv.Plotter(off_screen=True, window_size=(400, 300))
    except Exception as e:  # noqa: BLE001
        pytest.skip(f"rendering non disponibile: {e}")
    pl.enable_parallel_projection()
    c = scene.center
    pl.camera_position = [tuple(c - direction * 30), tuple(c), tuple(up)]
    pl.camera.parallel_scale = 8
    ren = pl.renderer
    screen = []
    for p in scene.pos:
        ren.SetWorldPoint(*p, 1.0)
        ren.WorldToDisplay()
        screen.append(ren.GetDisplayPoint()[:2])
    pl.close()
    screen = np.array(screen)
    q = ctx.frame.project(scene.pos)[:, :2]  # coordinate (L, W) della Proiezione
    # schermo = k · (L, W) + t, con k > 0 uguale sui due assi (niente rotazioni né specchiature)
    A = np.column_stack([q, np.ones(len(q))])
    coef, *_ = np.linalg.lstsq(A, screen, rcond=None)
    assert np.allclose(A @ coef, screen, atol=0.5)
    k = coef[0, 0]
    assert k > 0 and np.isclose(coef[1, 1], k, rtol=1e-3) and abs(coef[1, 0]) < 1e-3 * k and abs(coef[0, 1]) < 1e-3 * k


def test_incertezze_nel_dettaglio(world):
    s, a, v = world
    from test_uncertainty import with_su

    from pistak.core.analysis import read_cif_file

    b = api.analyse_text(with_su(read_cif_file(EX / "herringbone_anthracene.cif"), 20), "h", s.limits)
    s._apply(b)
    api.uncertainty(b, samples=10, keys=["interplanar"])
    s.uncertaintyChanged.emit()
    tree = v["detailA"].params
    texts = [tree.topLevelItem(i).text(1) for i in range(tree.topLevelItemCount())]
    assert any("10" in t for t in texts)  # la riga «Incertezze: Monte Carlo, 10 campioni»
    values = [tree.topLevelItem(i).child(j).text(1) for i in range(tree.topLevelItemCount())
              for j in range(tree.topLevelItem(i).childCount())]
    assert any("±" in t for t in values)


def test_finestra_calcoli_quantistici(app):
    from pistak.gui.qc_dialog import QCDialog
    from pistak.io.qc import QCOptions

    d = QCDialog(QCOptions("gaussian", dispersion="d3zero", charges=(1, -1), multiplicities=(2, 2)))
    o = d.options()
    assert (o.program, o.dispersion, o.charges, o.multiplicities) == ("gaussian", "d3zero", (1, -1), (2, 2))
    d.method.setCurrentIndex(d.method.findData("sapt"))  # SAPT esiste solo in Psi4
    assert d.options().program == "psi4" and not d.functional.isEnabled() and d.sapt_level.isEnabled()
    d.method.setCurrentIndex(d.method.findData("dft-d"))
    d.dispersion.setCurrentIndex(d.dispersion.findData("d4"))
    d.program.setCurrentIndex(d.program.findData("gaussian"))  # Gaussian non ha D4
    assert d.options().dispersion == "d3bj"
    d.options().check()


def test_dettaglio_geometria_e_controlli(world):
    s, a, v = world
    d = v["detailA"]
    assert d.geometry.rowCount() == 0  # la geometria si calcola solo quando la scheda è aperta
    d.tabs.setCurrentIndex(6)
    n_bonds = d.geometry.rowCount()
    assert n_bonds > 0 and d.geometry.columnCount() == 5
    d.geom_kind.setCurrentIndex(d.geom_kind.findData("torsion"))
    assert d.geometry.columnCount() == 7 and d.geometry.rowCount() > 0
    d.geom_kind.setCurrentIndex(d.geom_kind.findData("bond"))
    d.geom_h.setChecked(False)
    assert 0 < d.geometry.rowCount() < n_bonds
    assert d.checks.rowCount() == len(a.alerts) and d.checks.item(0, 1).text()
    assert "Volume della cella" in d.summary.toPlainText()


def test_finestra_calcoli_periodici(app):
    from pistak.gui.periodic_dialog import PeriodicDialog
    from pistak.io.periodic import PeriodicOptions

    d = PeriodicDialog(PeriodicOptions("cp2k", "vc-relax", dispersion="d4", kpoint_length=25))
    o = d.options()
    assert (o.program, o.task, o.dispersion, o.kpoint_length) == ("cp2k", "vc-relax", "d4", 25)
    d.program.setCurrentIndex(d.program.findData("qe"))  # Quantum ESPRESSO non ha D4
    assert d.options().dispersion == "d3bj"
    d.options().check()


@pytest.fixture()
def view3d(world):
    pytest.importorskip("pyvistaqt")
    from pistak.gui.views.view3d import View3D, _NullPlotter

    s, a, _ = world
    v = View3D(s, "A")
    if isinstance(v.plotter, _NullPlotter):
        pytest.skip("rendering non disponibile")
    v.resize(600, 400)
    v.refresh()
    yield s, a, v
    v.close()


def test_png_della_3d_in_alta_risoluzione(view3d, tmp_path):
    from PySide6.QtGui import QImage

    _, _, v = view3d
    out = tmp_path / "v.png"
    v.screenshot(str(out), scale=3)
    im = QImage(str(out))
    w, h = max(v.plotter.interactor.width(), 100), max(v.plotter.interactor.height(), 100)  # senza schermo il pannello è minuscolo
    assert (im.width(), im.height()) == (3 * w, 3 * h)


def test_mantieni_la_vista_rispetto_alla_coppia(view3d):
    s, a, v = view3d
    v.lock_cam.setChecked(True)
    v.look("side")
    v.plotter.camera.Azimuth(25)  # l'utente ruota: la vista non è più quella automatica
    before = v._relative_camera()
    s.select("A", 2)
    after = v._relative_camera()
    assert v.ctx.pair is a.pairs[2]
    for k in ("pos", "focal", "up"):
        assert np.allclose(before[k], after[k], atol=1e-6)
    v.lock_cam.setChecked(False)
    s.select("A", 0)
    assert v._still_auto() and v._auto_view[0] == "top"  # senza: si torna alla vista dall'alto
