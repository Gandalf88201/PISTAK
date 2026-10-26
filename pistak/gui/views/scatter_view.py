"""Proiezione dei parametri: ogni coppia è un punto in un grafico a due assi scelti tra tutti i parametri
del registro. Un clic seleziona la coppia (e muove le viste sincronizzate); con "Selezione a riquadro" si
trascina un rettangolo e le coppie dentro diventano il filtro del gruppo (la tabella le mostra da sole).
Se un asse mostra un parametro confrontato da una soglia del criterio di stack, la soglia è una linea tratteggiata
che segue il pannello «Soglie in tempo reale»."""
from __future__ import annotations

import math

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QPen
from PySide6.QtWidgets import QCheckBox, QComboBox, QLabel, QPushButton

from pistak.gui import theme
from pistak.gui.session import Session
from pistak.gui.views.base import PairView
from pistak.i18n import tf, tr
from pistak.params import get_spec, specs

KIND_COLOR = "kind"  # valore interno per «colora per tipo di coppia» (non tradotto)
#: Parametro → la soglia che lo confronta (None se con queste soglie non si applica).
THRESHOLDS = {
    "interplanar": lambda lim: lim.max_interplanar,
    "dihedral": lambda lim: lim.max_dihedral,
    "slip": lambda lim: lim.max_slip if lim.max_slip > 0 else None,
    "overlap": lambda lim: lim.min_overlap if lim.overlap_measure != "area" else None,
    "overlap_ratio": lambda lim: lim.min_overlap if lim.overlap_measure == "area" else None,
    "contact": lambda lim: lim.contact_cutoff,
}


class ScatterView(PairView):
    TITLE = "Scatter parametri"

    def __init__(self, session: Session, group: str | None = "A", parent=None,
                 x: str = "long_slip", y: str = "interplanar"):
        super().__init__(session, group, parent)
        self.x_combo, self.y_combo, self.c_combo = QComboBox(), QComboBox(), QComboBox()
        for lab, combo in (("X", self.x_combo), ("Y", self.y_combo), (tr("Colore"), self.c_combo)):
            self.toolbar.addWidget(QLabel(lab))
            self.toolbar.addWidget(combo)
        self._want = {"x": x, "y": y, "c": KIND_COLOR}
        self.box_mode = QCheckBox(tr("Selezione a riquadro"))
        self.box_mode.setToolTip(tr("Trascina un rettangolo: le coppie dentro diventano il filtro del gruppo"))
        self.box_mode.toggled.connect(self._box_toggled)
        self.toolbar.addWidget(self.box_mode)
        self.clear_btn = QPushButton(tr("Azzera filtro"))
        self.clear_btn.clicked.connect(lambda: self.filter_to(None))
        self.toolbar.addWidget(self.clear_btn)

        self.plot = pg.PlotWidget(background=theme.background())
        self.plot.showGrid(x=True, y=True, alpha=0.2)
        self.body.addWidget(self.plot, 1)
        self.info = QLabel("")
        self.body.addWidget(self.info)

        self.cloud = pg.ScatterPlotItem(size=10, pxMode=True, hoverable=True, hoverSize=14, tip=self._tip)
        self.cloud.sigClicked.connect(self._clicked)
        self.marker = pg.ScatterPlotItem(size=20, pxMode=True, brush=QBrush(Qt.BrushStyle.NoBrush), pen=pg.mkPen(theme.foreground(), width=2.5))
        self.plot.addItem(self.cloud)
        self.plot.addItem(self.marker)
        self.roi: pg.RectROI | None = None
        self.lines: list[pg.InfiniteLine] = []
        self._fill_combos()
        for c in (self.x_combo, self.y_combo, self.c_combo):
            c.currentIndexChanged.connect(self._axes_changed)

    # ---- assi ----------------------------------------------------------------------
    def _fill_combos(self) -> None:
        for combo, want in ((self.x_combo, self._want["x"]), (self.y_combo, self._want["y"])):
            combo.blockSignals(True)
            combo.clear()
            for s in specs():
                combo.addItem(s.title, s.key)
            combo.setCurrentIndex(max(0, combo.findData(want)))
            combo.blockSignals(False)
        self.c_combo.blockSignals(True)
        self.c_combo.clear()
        self.c_combo.addItem(tr("tipo di coppia"), None)
        for s in specs():
            self.c_combo.addItem(s.title, s.key)
        self.c_combo.setCurrentIndex(max(0, self.c_combo.findData(None if self._want["c"] == KIND_COLOR else self._want["c"])))
        self.c_combo.blockSignals(False)

    def _axes_changed(self) -> None:
        self._want = {"x": self.x_combo.currentData(), "y": self.y_combo.currentData(), "c": self.c_combo.currentData() or KIND_COLOR}
        self._redraw()

    def on_params(self) -> None:
        self._fill_combos()

    def on_analysis(self) -> None:
        if not self.refreshing:  # durante refresh() ridisegna on_subset
            self._redraw()

    def on_pair(self, index: int) -> None:
        self._update_marker()

    def on_subset(self, indices) -> None:
        self._redraw()

    # ---- disegno -------------------------------------------------------------------
    def _values(self, key: str) -> np.ndarray:
        a = self.analysis
        return np.array([p.get(key) for p in a.pairs]) if a else np.zeros(0)

    def _redraw(self) -> None:
        a = self.analysis
        kx, ky = self.x_combo.currentData(), self.y_combo.currentData()
        sx, sy = get_spec(kx), get_spec(ky)
        self.plot.setLabel("bottom", sx.title if sx else "")
        self.plot.setLabel("left", sy.title if sy else "")
        self.plot.setBackground(theme.background())
        if not a or not a.pairs or not sx or not sy:
            self.cloud.setData([], [])
            self.marker.setData([], [])
            self.info.setText("")
            return
        x, y = self._values(kx), self._values(ky)
        sub = self.subset
        brushes, pens = [], []
        ckey = self.c_combo.currentData()
        cvals = self._values(ckey) if ckey else None
        if cvals is not None and np.isfinite(cvals).any():
            lo, hi = np.nanmin(cvals), np.nanmax(cvals)
            cmap = pg.colormap.get("viridis")
        for p in a.pairs:
            if cvals is not None and np.isfinite(cvals).any():
                t = 0.0 if hi == lo else (cvals[p.index] - lo) / (hi - lo)
                rgb = cmap.map(0.0 if not math.isfinite(t) else float(t), mode="byte")[:3]
            else:
                rgb = theme.KIND_COLORS[p.kind]
            alpha = 230 if (sub is None or p.index in sub) else 40
            brushes.append(QBrush(theme.qcolor(tuple(int(v) for v in rgb), alpha)))
            pens.append(QPen(theme.qcolor((0, 0, 0), alpha // 2)))
        size = [8 + min(10, 2.5 * math.log2(p.count)) if p.count > 0 else 8 for p in a.pairs]
        ok = np.isfinite(x) & np.isfinite(y)
        self.cloud.setData(x[ok], y[ok], brush=[b for b, k in zip(brushes, ok) if k], pen=[q for q, k in zip(pens, ok) if k],
                           size=[s for s, k in zip(size, ok) if k], data=[p.index for p, k in zip(a.pairs, ok) if k])
        self._update_marker()
        self._threshold_lines(kx, ky)
        self.plot.autoRange()
        n_sub = len(sub) if sub is not None else len(a.pairs)
        self.info.setText(tf("{n} / {tot} coppie in evidenza · dimensione del punto ∝ coppie equivalenti", n=n_sub, tot=len(a.pairs)))

    def _threshold_lines(self, kx: str, ky: str) -> None:
        for line in self.lines:
            self.plot.removeItem(line)
        self.lines = []
        a = self.analysis
        if a is None:
            return
        pen = pg.mkPen(theme.foreground(), width=1, style=Qt.PenStyle.DashLine)
        for key, angle in ((kx, 90), (ky, 0)):
            v = THRESHOLDS[key](a.limits) if key in THRESHOLDS else None
            if v is not None:
                line = pg.InfiniteLine(pos=v, angle=angle, pen=pen, movable=False)
                line.setToolTip(tr("Soglia del criterio di stack (Analisi ▸ Soglie in tempo reale)"))
                self.plot.addItem(line)
                self.lines.append(line)

    def _update_marker(self) -> None:
        pair = self.pair
        kx, ky = self.x_combo.currentData(), self.y_combo.currentData()
        if pair and kx and ky and math.isfinite(pair.get(kx)) and math.isfinite(pair.get(ky)):
            self.marker.setData([pair.get(kx)], [pair.get(ky)])
        else:
            self.marker.setData([], [])

    def _tip(self, x, y, data) -> str:
        a = self.analysis
        p = a.pairs[data] if a and data is not None and data < len(a.pairs) else None
        return f"#{p.index + 1} {theme.kind_label(p.kind)}  ×{p.count}\nX = {x:.3g}   Y = {y:.3g}" if p else ""

    # ---- interazione ---------------------------------------------------------------
    def _clicked(self, _item, points, *_) -> None:
        if len(points):
            self.choose(int(points[0].data()))

    def _box_toggled(self, on: bool) -> None:
        if self.roi is not None:
            self.plot.removeItem(self.roi)
            self.roi = None
        if on:
            (x0, x1), (y0, y1) = self.plot.viewRange()
            w, h = (x1 - x0) / 3, (y1 - y0) / 3
            self.roi = pg.RectROI((x0 + w, y0 + h), (w, h), pen=pg.mkPen((220, 120, 0), width=2))
            self.roi.addScaleHandle((1, 1), (0, 0))
            self.roi.addScaleHandle((0, 0), (1, 1))
            self.roi.sigRegionChangeFinished.connect(self._roi_changed)
            self.plot.addItem(self.roi)
            self._roi_changed()
        else:
            self.filter_to(None)

    def _roi_changed(self) -> None:
        a = self.analysis
        if not (self.roi and a):
            return
        pos, size = self.roi.pos(), self.roi.size()
        x, y = self._values(self.x_combo.currentData()), self._values(self.y_combo.currentData())
        x0, x1 = sorted((pos.x(), pos.x() + size.x()))
        y0, y1 = sorted((pos.y(), pos.y() + size.y()))
        inside = {p.index for p in a.pairs if x0 <= x[p.index] <= x1 and y0 <= y[p.index] <= y1}
        self.filter_to(inside)

    def extra_state(self) -> dict:
        return {"x": self.x_combo.currentData(), "y": self.y_combo.currentData(), "color": self.c_combo.currentData()}

    def apply_state(self, st: dict) -> None:
        self._want = {"x": st.get("x", self._want["x"]), "y": st.get("y", self._want["y"]), "c": st.get("color") or KIND_COLOR}
        for c in (self.x_combo, self.y_combo, self.c_combo):
            c.blockSignals(True)
        self._fill_combos()
        for c in (self.x_combo, self.y_combo, self.c_combo):
            c.blockSignals(False)
