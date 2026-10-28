"""Rigenera le figure del manuale (`pistak/manual/img/<lingua>/*.png`) dalla finestra vera, senza schermo.

    .venv/bin/python tools/manual_figures.py            # italiano e inglese
    .venv/bin/python tools/manual_figures.py --lang en

Le figure usano i CIF sintetici di `examples/` e, per la geometria, una struttura reale della COD in
`tests/data/cod/` (dominio pubblico). Qt gira «offscreen»; la vista 3D si disegna fuori schermo
con PyVista e compare come immagine al posto della finestra OpenGL, con la stessa scena e la stessa camera.
Dopo un cambiamento dell'interfaccia basta rilanciare lo script e controllare le immagini.
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ["QT_QPA_PLATFORM"] = "offscreen"
os.environ["PISTAK_NO_3D"] = "1"
os.environ["PISTAK_CONFIG_DIR"] = tempfile.mkdtemp(prefix="pistak-figures-")  # nessuna impostazione personale

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
EX = ROOT / "examples"
COD = ROOT / "tests" / "data" / "cod"
OUT = ROOT / "pistak" / "manual" / "img"

import pyvista as pv  # noqa: E402
from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtGui import QImage, QPixmap  # noqa: E402
from PySide6.QtWidgets import QApplication, QLabel, QWidget  # noqa: E402

from pistak import i18n  # noqa: E402
from pistak.core.config import Limits  # noqa: E402


class ShotPlotter:
    """Un plotter PyVista fuori schermo che mostra ogni render come immagine in un QLabel."""

    def __init__(self, parent: QWidget):
        object.__setattr__(self, "pl", pv.Plotter(off_screen=True, window_size=(900, 640)))
        label = QLabel(parent)
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setMinimumSize(200, 150)
        object.__setattr__(self, "interactor", label)

    def __getattr__(self, name):
        return getattr(self.pl, name)

    def __setattr__(self, name, value):
        setattr(self.pl, name, value)

    def render(self) -> None:
        size = self.interactor.size()
        self.pl.window_size = [max(300, size.width()), max(200, size.height())]
        img = self.pl.screenshot(return_img=True)
        h, w, _ = img.shape
        qimg = QImage(img.tobytes(), w, h, 3 * w, QImage.Format.Format_RGB888).copy()
        self.interactor.setPixmap(QPixmap.fromImage(qimg))


def wait(app: QApplication, cond, timeout: float = 60) -> None:
    t = time.time()
    while not cond():
        if time.time() - t > timeout:
            raise TimeoutError("l'analisi non è arrivata")
        app.processEvents()
        time.sleep(0.02)
    for _ in range(5):
        app.processEvents()


def save(widget: QWidget, path: Path) -> None:
    """PNG a 256 colori (le schermate dell'interfaccia restano nitide e pesano un quarto)."""
    from PIL import Image

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.png")
    widget.grab().save(str(tmp))
    with Image.open(tmp) as im:
        im.convert("RGB").quantize(256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(path, optimize=True)
    tmp.unlink()
    print("  ", path.relative_to(ROOT), f"{path.stat().st_size // 1024} kB")


def figures(app: QApplication, lang: str) -> None:
    from pistak.gui.collection_window import CollectionWindow
    from pistak.gui.main_window import MainWindow
    from pistak.gui.periodic_dialog import PeriodicDialog
    from pistak.gui.qc_dialog import QCDialog
    from pistak.gui.session import Session
    from pistak.gui.settings_dialog import SettingsDialog
    from pistak.gui.views.detail_view import DetailView
    from pistak.gui.views.hist_view import HistogramView
    from pistak.gui.views.projection_view import ProjectionView
    from pistak.gui.views.scatter_view import ScatterView
    from pistak.gui.views.table_view import TableView
    from pistak.gui.views.view3d import View3D

    i18n.set_language(lang)
    out = OUT / lang
    print(f"[{lang}]")
    win = MainWindow(Session())
    win.autosave = False
    for d in win.docks:
        v = d.widget()
        if isinstance(v, View3D):
            old = v.plotter.interactor
            v.plotter = ShotPlotter(v)
            v.body.replaceWidget(old, v.plotter.interactor)
            old.hide()
            old.deleteLater()
    win.resize(1360, 840)
    win.show()

    def load(name: str | Path, limits: Limits | None = None) -> None:
        if limits is not None:
            win.session.limits = limits
        before = win.session.analysis
        win.session.load_file(EX / name)  # un percorso assoluto resta quello
        wait(app, lambda: win.session.analysis is not before and win.session.analysis is not None)

    def view(cls):
        return next(d.widget() for d in win.docks if isinstance(d.widget(), cls))

    def dock_of(widget):
        return next(d for d in win.docks if d.widget() is widget)

    def alone(widget, w: int, h: int):
        """La vista staccata in una finestra propria di dimensione fissa, per una figura leggibile."""
        d = dock_of(widget)
        d.setFloating(True)
        d.resize(w, h)
        d.show()
        wait(app, lambda: True)  # prima la nuova dimensione, poi il disegno (la 3D si rende alla misura del riquadro)
        widget.refresh()
        if hasattr(widget, "_rebuild"):  # la 3D non ricostruisce una scena invariata: qui serve, cambia la misura
            widget._rebuild(keep_camera=False)
        wait(app, lambda: True)
        return d

    # finestra principale e viste: l'herringbone ha uno stack (colonna) e una coppia inclinata
    load("herringbone_anthracene.cif", Limits())
    win.session.select("A", 0)
    wait(app, lambda: True)
    save(win, out / "main_window.png")
    v3 = view(View3D)
    alone(v3, 720, 540)
    save(v3, out / "view3d.png")
    proj = view(ProjectionView)
    alone(proj, 720, 520)
    save(proj, out / "projection.png")
    det = view(DetailView)
    alone(det, 680, 560)
    for i, name in ((0, "detail_structure"), (1, "detail_pair"), (3, "detail_columns"), (7, "detail_checks")):
        det.tabs.setCurrentIndex(i)
        wait(app, lambda: True)
        save(det, out / f"{name}.png")

    # geometria: i contatti di van der Waals di una struttura reale con legami a idrogeno (COD 4500022)
    load(COD / "4500022.cif")
    det.tabs.setCurrentIndex(6)
    det.geom_kind.setCurrentIndex(det.geom_kind.findData("contact"))
    wait(app, lambda: det.geometry.rowCount() > 0)
    save(det, out / "detail_geometry.png")
    # lo stesso legame a idrogeno nella 3D (doppio clic sulla riga), con la molecola d'acqua
    det._geometry_chosen(0, 0)
    v3.picks = []
    alone(v3, 720, 540)
    v3.look("top")
    wait(app, lambda: True)
    save(v3, out / "view3d_measure.png")
    det.tabs.setCurrentIndex(0)

    # tabella, scatter e istogramma: i dimeri del naftalene hanno uno stack e tre coppie parallele
    load("dimers_naphthalene.cif")
    table = view(TableView)
    alone(table, 900, 260)
    save(table, out / "table.png")
    sc = view(ScatterView)
    sc.x_combo.setCurrentIndex(sc.x_combo.findData("long_slip"))
    sc.y_combo.setCurrentIndex(sc.y_combo.findData("interplanar"))
    alone(sc, 720, 500)
    save(sc, out / "scatter.png")
    hist = view(HistogramView)
    alone(hist, 720, 440)
    save(hist, out / "histogram.png")

    # soglie: pannello in tempo reale, finestra delle soglie
    win.show_limits_panel()
    wait(app, lambda: True)
    save(win.limits_dock, out / "limits_panel.png")
    dlg = SettingsDialog(win.session.limits, win)
    dlg.resize(560, 640)
    dlg.show()
    wait(app, lambda: True)
    save(dlg, out / "settings.png")
    dlg.close()
    qc = QCDialog(parent=win)
    qc.show()
    wait(app, lambda: True)
    save(qc, out / "qc_dialog.png")
    qc.close()
    per = PeriodicDialog(parent=win)
    per.show()
    wait(app, lambda: True)
    save(per, out / "periodic_dialog.png")
    per.close()

    # confronto e statistiche su tutti gli esempi
    col: CollectionWindow = win.show_collection()
    col.resize(1150, 720)
    col.add_paths([EX])
    wait(app, lambda: len(col.entries) > 0 and col._job is None)
    for i, name in ((0, "collection_structures"), (2, "collection_scatter"), (4, "collection_summary")):
        col.tabs.setCurrentIndex(i)
        wait(app, lambda: True)
        save(col, out / f"{name}.png")
    col.close()
    win.limits_dock.hide()
    win.close()


def main() -> int:
    p = argparse.ArgumentParser(description="Figure del manuale")
    p.add_argument("--lang", choices=["it", "en"], action="append")
    args = p.parse_args()
    app = QApplication.instance() or QApplication([])
    for lang in args.lang or ["it", "en"]:
        figures(app, lang)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
