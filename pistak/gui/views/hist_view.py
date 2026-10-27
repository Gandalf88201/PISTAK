"""Istogramma di un parametro, con le coppie impilate per tipo. Un clic su una barra filtra le coppie di quel
intervallo (nel gruppo della vista); la riga tratteggiata è la coppia selezionata."""
from __future__ import annotations

import math

import numpy as np
import pyqtgraph as pg
from PySide6.QtWidgets import QCheckBox, QComboBox, QLabel, QSpinBox

from pistak.gui import theme
from pistak.gui.session import Session
from pistak.gui.views.base import PairView
from pistak.i18n import tr
from pistak.params import get_spec, specs


class HistogramView(PairView):
    TITLE = "Istogramma"

    def __init__(self, session: Session, group: str | None = "A", parent=None, key: str = "interplanar"):
        super().__init__(session, group, parent)
        self._want = key
        self.combo = QComboBox()
        self.bins = QSpinBox()
        self.bins.setRange(3, 80)
        self.bins.setValue(15)
        self.weighted = QCheckBox(tr("Pesa per coppie equivalenti"))
        self.weighted.setChecked(True)
        self.toolbar.addWidget(self.combo)
        self.toolbar.addWidget(QLabel(tr("Intervalli")))
        self.toolbar.addWidget(self.bins)
        self.toolbar.addWidget(self.weighted)
        self.plot = theme.plot_widget()
        self.plot.showGrid(x=False, y=True, alpha=0.2)
        self.body.addWidget(self.plot, 1)
        self.line = pg.InfiniteLine(angle=90, pen=pg.mkPen(theme.foreground(), width=2, style=pg.QtCore.Qt.PenStyle.DashLine))
        self._edges = np.zeros(0)
        self._fill()
        self.combo.currentIndexChanged.connect(self._changed)
        self.bins.valueChanged.connect(self._redraw)
        self.weighted.toggled.connect(self._redraw)
        self.plot.scene().sigMouseClicked.connect(self._click)

    def _fill(self) -> None:
        self.combo.blockSignals(True)
        self.combo.clear()
        for s in specs():
            self.combo.addItem(s.title, s.key)
        self.combo.setCurrentIndex(max(0, self.combo.findData(self._want)))
        self.combo.blockSignals(False)

    def _changed(self) -> None:
        self._want = self.combo.currentData()
        self._redraw()

    def on_params(self) -> None:
        self._fill()

    def on_analysis(self) -> None:
        self._redraw()

    def on_subset(self, indices) -> None:
        self._redraw()

    def on_pair(self, index: int) -> None:
        self._update_line()

    def _redraw(self) -> None:
        self.plot.clear()
        self.plot.setBackground(theme.background())
        a = self.analysis
        key = self.combo.currentData()
        spec = get_spec(key)
        if not (a and a.pairs and spec):
            return
        vals = np.array([p.get(key) for p in a.pairs])
        w = np.array([p.count if self.weighted.isChecked() else 1 for p in a.pairs], dtype=float)
        ok = np.isfinite(vals)
        if not ok.any():
            return
        edges = np.histogram_bin_edges(vals[ok], bins=self.bins.value())
        self._edges = edges
        width = (edges[1] - edges[0]) * 0.92 if len(edges) > 1 else 1.0
        centers = (edges[:-1] + edges[1:]) / 2
        bottom = np.zeros(len(centers))
        for kind in ("stack", "parallel", "tilted"):
            sel = ok & np.array([p.kind == kind for p in a.pairs])
            h, _ = np.histogram(vals[sel], bins=edges, weights=w[sel])
            if h.sum():
                self.plot.addItem(pg.BarGraphItem(x=centers, height=h, width=width, y0=bottom, brush=theme.qcolor(theme.KIND_COLORS[kind]),
                                                  pen=pg.mkPen(None), name=theme.kind_label(kind)))
            bottom = bottom + h
        self.plot.setLabel("bottom", spec.title)
        self.plot.setLabel("left", tr("coppie") + (tr(" (con equivalenti)") if self.weighted.isChecked() else ""))
        self.plot.addItem(self.line)
        self._update_line()

    def _update_line(self) -> None:
        pair, key = self.pair, self.combo.currentData()
        v = pair.get(key) if pair and key else math.nan
        self.line.setVisible(math.isfinite(v))
        if math.isfinite(v):
            self.line.setValue(v)

    def _click(self, ev) -> None:
        a = self.analysis
        if not (a and len(self._edges) > 1) or not self.plot.sceneBoundingRect().contains(ev.scenePos()):
            return
        x = self.plot.getPlotItem().vb.mapSceneToView(ev.scenePos()).x()
        i = int(np.searchsorted(self._edges, x, side="right")) - 1
        if not 0 <= i < len(self._edges) - 1:
            self.filter_to(None)
            return
        lo, hi = self._edges[i], self._edges[i + 1]
        key = self.combo.currentData()
        self.filter_to({p.index for p in a.pairs if lo <= p.get(key) <= hi})

    def extra_state(self) -> dict:
        return {"key": self.combo.currentData(), "bins": self.bins.value(), "weighted": self.weighted.isChecked()}

    def apply_state(self, st: dict) -> None:
        self._want = st.get("key", self._want)
        for w in (self.combo, self.bins, self.weighted):
            w.blockSignals(True)
        self._fill()
        self.bins.setValue(int(st.get("bins", self.bins.value())))
        self.weighted.setChecked(bool(st.get("weighted", True)))
        for w in (self.combo, self.bins, self.weighted):
            w.blockSignals(False)
