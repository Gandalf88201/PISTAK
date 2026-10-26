"""Proiezione sul piano medio della coppia: atomi pesanti e legami delle due molecole, contorni dei nuclei,
area di sovrapposizione, baricentri e scomposizione dello slittamento in L e W (come la figura del sito)."""
from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QBrush, QColor, QPen, QPolygonF
from PySide6.QtWidgets import QCheckBox, QGraphicsPolygonItem, QLabel

from pistak.core.elements import CYAN_CARBON, color
from pistak.core.geometry import clip_convex, convex_hull
from pistak.gui import theme
from pistak.gui.scene import COLOR_L, COLOR_PLANE_A, COLOR_PLANE_B, COLOR_W
from pistak.gui.session import Session
from pistak.gui.views.base import PairView
from pistak.i18n import tf, tr


def _pen(rgb, w=2.0, style=Qt.PenStyle.SolidLine, alpha=255) -> QPen:
    pen = QPen(theme.qcolor(rgb, alpha))
    pen.setWidthF(w)
    pen.setStyle(style)
    pen.setCosmetic(True)
    return pen


class ProjectionView(PairView):
    TITLE = "Proiezione"

    def __init__(self, session: Session, group: str | None = "A", parent=None):
        super().__init__(session, group, parent)
        self.show_h = QCheckBox("H")
        self.show_h.toggled.connect(lambda _: self._redraw())
        self.show_overlap = QCheckBox(tr("Sovrapposizione"))
        self.show_overlap.setChecked(True)
        self.show_overlap.toggled.connect(lambda _: self._redraw())
        self.show_atoms = QCheckBox(tr("Atomi"))
        self.show_atoms.setChecked(True)
        self.show_atoms.toggled.connect(lambda _: self._redraw())
        for w in (self.show_atoms, self.show_h, self.show_overlap):
            self.toolbar.addWidget(w)
        self.plot = pg.PlotWidget(background=theme.background())
        self.plot.setAspectLocked(True)
        self.plot.showGrid(x=True, y=True, alpha=0.15)
        self.plot.setLabel("bottom", tr("L (Å) — asse lungo della prima molecola"))
        self.plot.setLabel("left", "W (Å)")
        self.plot.getPlotItem().getViewBox().invertY(False)
        self.caption = QLabel("")
        self.caption.setWordWrap(True)
        self.body.addWidget(self.plot, 1)
        self.body.addWidget(self.caption)

    def on_analysis(self) -> None:
        if not self.refreshing:  # durante refresh() ridisegna on_pair
            self._redraw()

    def on_pair(self, index: int) -> None:
        self._redraw()

    def _redraw(self) -> None:
        self.plot.clear()
        self.plot.setBackground(theme.background())
        pair = self.pair
        if not (self.analysis and pair):
            self.caption.setText("")
            return
        ctx = self.analysis.pair_context(pair)
        qa, qb = ctx.proj_a, ctx.proj_b
        for mol, q, other in ((ctx.A, qa, 0), (ctx.B, qb, 1)):
            base = COLOR_PLANE_A if other == 0 else COLOR_PLANE_B
            keep = np.ones(len(mol), bool) if self.show_h.isChecked() else mol.heavy
            segs = [(i, j) for i, j in mol.bonds if keep[i] and keep[j]]
            if segs:
                xs = np.array([[q[i, 0], q[j, 0]] for i, j in segs]).ravel()
                ys = np.array([[q[i, 1], q[j, 1]] for i, j in segs]).ravel()
                self.plot.addItem(pg.PlotCurveItem(xs, ys, connect="pairs", pen=_pen(base, 2.2)))
            if self.show_atoms.isChecked():
                idx = np.flatnonzero(keep)
                brushes = [QBrush(theme.qcolor(CYAN_CARBON if (other == 1 and mol.el[i] == "C") else color(mol.el[i]))) for i in idx]
                self.plot.addItem(pg.ScatterPlotItem(q[idx, 0], q[idx, 1], size=7 if self.show_h.isChecked() else 8, brush=brushes,
                                                     pen=_pen(base, 1.0), pxMode=True))
        ha = convex_hull(qa[ctx.core_a.atoms, :2])
        hb = convex_hull(qb[ctx.core_b.atoms, :2])
        for h, rgb in ((ha, COLOR_PLANE_A), (hb, COLOR_PLANE_B)):
            if len(h) >= 3:
                closed = np.vstack([h, h[:1]])
                self.plot.addItem(pg.PlotCurveItem(closed[:, 0], closed[:, 1], pen=_pen(rgb, 1.6, Qt.PenStyle.DashLine)))
        if self.show_overlap.isChecked() and len(ha) >= 3 and len(hb) >= 3:
            inter = clip_convex(ha, hb)
            if len(inter) >= 3:
                item = QGraphicsPolygonItem(QPolygonF([QPointF(float(x), float(y)) for x, y in inter]))
                item.setBrush(QBrush(QColor(46, 160, 67, 70)))
                item.setPen(_pen((46, 160, 67), 1.0))
                self.plot.addItem(item)
        cb = qb[:, :2][ctx.core_b.atoms].mean(axis=0)
        # slittamento: prima lungo L (rosso), poi lungo W (verde)
        self._arrow((0, 0), (cb[0], 0), COLOR_L)
        self._arrow((cb[0], 0), (cb[0], cb[1]), COLOR_W)
        self.plot.addItem(pg.ScatterPlotItem([0], [0], size=11, brush=theme.qcolor(COLOR_PLANE_A), pen=_pen((0, 0, 0), 1.0)))
        self.plot.addItem(pg.ScatterPlotItem([cb[0]], [cb[1]], size=11, brush=theme.qcolor(COLOR_PLANE_B), pen=_pen((0, 0, 0), 1.0)))
        self.plot.autoRange()
        self.caption.setText(
            f"<span style='color:rgb{COLOR_L}'>■</span> L {pair.long_slip:.2f} Å &nbsp; "
            f"<span style='color:rgb{COLOR_W}'>■</span> W {pair.lat_slip:.2f} Å &nbsp; · "
            + tf("interplanare {ip:.2f} Å · angolo {ang:.1f}° · sovrapposizione {ov:.1f} Å²", ip=pair.interplanar, ang=pair.dihedral,
                 ov=pair.get("overlap_area"))
        )

    def _arrow(self, a, b, rgb) -> None:
        if abs(a[0] - b[0]) + abs(a[1] - b[1]) < 1e-6:
            return
        self.plot.addItem(pg.PlotCurveItem([a[0], b[0]], [a[1], b[1]], pen=_pen(rgb, 3.0)))
        angle = float(np.degrees(np.arctan2(b[1] - a[1], b[0] - a[0])))
        head = pg.ArrowItem(pos=(b[0], b[1]), angle=180 - angle, headLen=16, tipAngle=30, brush=theme.qcolor(rgb), pen=_pen(rgb, 1.0))
        self.plot.addItem(head)

    def extra_state(self) -> dict:
        return {"atoms": self.show_atoms.isChecked(), "h": self.show_h.isChecked(), "overlap": self.show_overlap.isChecked()}

    def apply_state(self, st: dict) -> None:
        for box, key, default in ((self.show_atoms, "atoms", True), (self.show_h, "h", False), (self.show_overlap, "overlap", True)):
            box.blockSignals(True)
            box.setChecked(bool(st.get(key, default)))
            box.blockSignals(False)
