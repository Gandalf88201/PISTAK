"""Session (save/restore) and comparison window, with Qt without a screen and without OpenGL."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_API", "pyside6")
os.environ["PISTAK_NO_3D"] = "1"

import json
import time
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from pistak.core.config import Limits
from pistak.gui import session_file
from pistak.gui.collection_window import CollectionWindow
from pistak.gui.main_window import MainWindow
from pistak.gui.session import Session
from pistak.gui.views.hist_view import HistogramView
from pistak.gui.views.scatter_view import ScatterView
from pistak.gui.views.table_view import TableView

EX = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def wait(app, cond, timeout=40):
    t = time.time()
    while not cond() and time.time() - t < timeout:
        app.processEvents()
        time.sleep(0.02)
    assert cond(), "timed out"


def new_window(app) -> MainWindow:
    win = MainWindow(Session())
    win.session.limits = Limits.site()  # 3 pairs in the herringbone: needed to test selection and filters
    win.autosave = False
    return win


def view(win, cls, n=0):
    return [d.widget() for d in win.docks if isinstance(d.widget(), cls)][n]


def test_save_and_restore(app, tmp_path):
    win = new_window(app)
    win.session.limits = Limits.site(max_interplanar=4.5)
    win.session.load_file(EX / "herringbone_anthracene.cif")
    wait(app, lambda: win.session.analysis is not None)
    s = win.session
    s.select("A", 2)
    s.select("B", 1)
    s.set_subset("A", {1, 2})
    sc = view(win, ScatterView)
    sc.x_combo.setCurrentIndex(sc.x_combo.findData("dihedral"))
    sc.y_combo.setCurrentIndex(sc.y_combo.findData("contact"))
    table = view(win, TableView)
    table.kind_boxes["parallel"].setChecked(False)
    dock = win.add_view("histogram", "B")
    win.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
    h = dock.widget()
    h.apply_state({"key": "lat_slip", "bins": 9, "weighted": False})
    n_views = len(win.docks)

    out = win.save_session(tmp_path / "test.pistak")
    data = json.loads(out.read_text())
    assert data["version"] == 1 and data["app"] == "PISTAK"
    assert data["source"]["name"] == "herringbone_anthracene.cif" and "data_" in data["source"]["text"]
    assert data["selection"]["A"] == 2 and data["subset"]["A"] == [1, 2] and data["subset"]["B"] is None
    assert len(data["views"]) == n_views

    win2 = new_window(app)
    assert win2.load_session(out)
    wait(app, lambda: win2.session.analysis is not None and win2.session.selection("A") == 2)
    s2 = win2.session
    assert s2.limits.max_interplanar == 4.5
    assert s2.source_name == "herringbone_anthracene.cif"
    assert (s2.selection("A"), s2.selection("B")) == (2, 1)
    assert s2.subset("A") == {1, 2} and s2.subset("B") is None
    assert len(win2.docks) == n_views
    sc2 = view(win2, ScatterView)
    assert (sc2.x_combo.currentData(), sc2.y_combo.currentData()) == ("dihedral", "contact")
    assert view(win2, TableView).kind_boxes["parallel"].isChecked() is False
    h2 = view(win2, HistogramView, 1)
    assert h2.group == "B" and h2.bins.value() == 9 and h2.combo.currentData() == "lat_slip" and not h2.weighted.isChecked()
    # the panel layout: in the default layout scatter and histogram are tabs of the same group
    dock_of = lambda w: next(d for d in win2.docks if d.widget() is w)  # noqa: E731
    assert dock_of(h2 if False else view(win2, HistogramView, 0)) in win2.tabifiedDockWidgets(dock_of(sc2))
    assert view(win2, TableView).index == 2  # the view follows the restored selection of its group
    assert h2.index == 1


def test_session_without_embedded_cif_and_missing_file(app, tmp_path):
    win = new_window(app)
    errors = []
    win._show_error = errors.append
    p = tmp_path / "x.pistak"
    session_file.write(p, {"source": {"name": "missing.cif", "path": str(tmp_path / "missing.cif")}, "views": []})
    assert win.load_session(p) is False
    assert errors and "non c'è più" in errors[0]


def test_invalid_session_file(app, tmp_path):
    win = new_window(app)
    errors = []
    win._show_error = errors.append
    bad = tmp_path / "bad.pistak"
    bad.write_text("{non json")
    assert win.load_session(bad) is False
    other = tmp_path / "other.pistak"
    other.write_text(json.dumps({"app": "other"}))
    assert win.load_session(other) is False
    newer = tmp_path / "newer.pistak"
    newer.write_text(json.dumps({"app": "PISTAK", "version": 99}))
    assert win.load_session(newer) is False
    assert len(errors) == 3 and "più nuova" in errors[2]


def test_unknown_keys_and_invalid_thresholds_are_ignored(app, tmp_path):
    win = new_window(app)
    p = tmp_path / "f.pistak"
    session_file.write(p, {"limits": {"max_dihedral": 12.0, "future_threshold": 1}, "views": [{"kind": "FutureView", "state": {}}],
                           "source": {"name": "d.cif", "text": (EX / "demo_naphthalene.cif").read_text()}, "news": {"x": 1}})
    assert win.load_session(p)
    wait(app, lambda: win.session.analysis is not None)
    assert win.session.limits.max_dihedral == 12.0 and len(win.docks) > 0  # default layout instead of the unknown views


def test_comparison_window(app):
    win = CollectionWindow(Limits())
    opened = []
    win.openRequested.connect(lambda path, i, block: opened.append((Path(path).name, i)))
    win.add_paths([EX])
    wait(app, lambda: len(win.entries) == 6)
    assert win.struct.rowCount() == 6 and len(win.table) > 0
    assert win.describe_rows() and win.describe_rows()[0]["n"] > 0
    assert win.cloud.data["x"].size > 0  # the scatter has the points
    names = [win.struct.item(r, 1).text() for r in range(6)]
    assert win.struct.item(0, 0).text() == "★"  # the first one is the reference
    win._open_entry(names.index("demo_naphthalene.cif"))
    assert opened == [("demo_naphthalene.cif", 0)]
    win._set_reference(names.index("demo_naphthalene.cif"))
    assert win.reference == names.index("demo_naphthalene.cif")
    win.kind_combo.setCurrentIndex(win.kind_combo.findData("tilted"))
    assert win.describe_rows()[0]["n"] > 0
    win.add_paths([EX])  # files already present: no duplicates
    app.processEvents()
    assert len(win.entries) == 6
    win.clear()
    assert win.struct.rowCount() == 0


def _loaded(app, limits=None):
    win = new_window(app)
    if limits is not None:
        win.session.limits = limits
    win.session.load_file(EX / "herringbone_anthracene.cif")
    wait(app, lambda: win.session.analysis is not None)
    return win


def test_threshold_sliders_in_real_time(app):
    win = _loaded(app)
    s, panel = win.session, win.limits_dock.widget()
    a0 = s.analysis
    assert [p.kind for p in a0.pairs].count("stack") == 1
    s.select("A", 1)
    key = a0.pairs[1].key
    sc = view(win, ScatterView)
    sc.y_combo.setCurrentIndex(sc.y_combo.findData("interplanar"))
    row = panel.rows["max_interplanar"]
    row.slider.setValue(round((3.0 - row.lo) / row.step))  # like dragging the slider
    assert row.spin.value() == pytest.approx(3.0)
    panel._commit()  # without waiting for the timer
    a1 = s.analysis
    assert a1 is not a0 and s.limits.max_interplanar == pytest.approx(3.0)
    assert "stack" not in [p.kind for p in a1.pairs]  # reclassified at once, in the interface thread
    assert a1.structure is a0.structure and a1.cache.contacts is a0.cache.contacts  # no geometry redone
    assert a1.pairs[s.selection("A")].key == key  # the selection stays on the same pair
    assert "0" in panel.summary.text()
    lines = [ln.value() for ln in sc.lines if ln.angle == 0]
    assert lines == [pytest.approx(3.0)]  # the threshold line follows the slider
    win.close()


def test_remeasuring_thresholds_in_a_thread(app):
    win = _loaded(app)
    s, panel = win.session, win.limits_dock.widget()
    a0 = s.analysis
    assert panel.rows["contact_cutoff"].heavy and not panel.rows["max_slip"].heavy
    for v in (3.5, 3.6, 3.7):  # close requests: they merge, the last one wins
        s.set_limits(s.limits.with_values(contact_cutoff=v))
    wait(app, lambda: s.analysis is not a0 and not s._limits_task and s.analysis.limits.contact_cutoff == pytest.approx(3.7))
    assert s.analysis.structure is a0.structure and s.analysis.cache.contacts is not a0.cache.contacts
    assert panel.rows["contact_cutoff"].value() == pytest.approx(3.7)  # the panel follows the session thresholds
    win.close()


def test_thresholds_from_elsewhere_update_the_panel(app):
    win = _loaded(app)
    s, panel = win.session, win.limits_dock.widget()
    s.set_limits(s.limits.with_values(max_dihedral=12.0, overlap_measure="area"))
    assert panel.rows["max_dihedral"].value() == pytest.approx(12.0)
    assert panel.measure.currentData() == "area"
    assert s.analysis.limits.max_dihedral == pytest.approx(12.0)
    panel._site()
    assert s.limits.max_dihedral == Limits.site().max_dihedral and s.limits.pair_mode == "site"
    win.close()


def test_manual_from_help_menu(app):
    win = new_window(app)
    win.show()
    m = win.show_manual("8-soglie-in-tempo-reale")
    app.processEvents()
    assert m.lang == "it" and m.toc.topLevelItemCount() == 21  # index + 20 chapters
    assert m.text.verticalScrollBar().value() > 0
    for i in range(m.toc.topLevelItemCount()):  # every index entry leads to its heading
        assert m.go_to(m.toc.topLevelItem(i).data(0, Qt.ItemDataRole.UserRole))
    widths = []  # wide figures fit the text (long code blocks scroll instead)
    block = m.text.document().begin()
    while block.isValid():
        it = block.begin()
        while not it.atEnd():
            fmt = it.fragment().charFormat()
            if fmt.isImageFormat():
                widths.append(fmt.toImageFormat().width())
            it += 1
        block = block.next()
    assert len(widths) >= 12 and 0 < max(widths) <= m.text.viewport().width()
    m.lang_combo.setCurrentIndex(m.lang_combo.findData("en"))
    assert m.lang == "en" and m.toc.topLevelItem(1).text(0) == "1. Basic concepts"
    assert win.show_manual() is not m or m.lang == "it"  # with the interface in Italian it reopens in Italian
    m.close()
    win.close()
