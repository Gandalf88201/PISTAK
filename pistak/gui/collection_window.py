"""Finestra «Confronto e statistiche»: un insieme di strutture analizzate insieme.

Schede: Strutture (il motivo di stack migliore di ciascuna e le differenze dal riferimento), Distribuzioni
(istogrammi su tutte le coppie), Dispersione (due parametri, ogni punto una coppia di un file), Correlazioni
(matrice di Pearson) e Riepilogo (descrizione statistica). Un clic su un punto o un doppio clic su una struttura
chiedono alla finestra principale di aprire quel file e quella coppia.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import QObject, QRunnable, Qt, QThreadPool, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMenu,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pistak.core.config import Limits
from pistak.gui import theme
from pistak.gui.flow_layout import FlowLayout
from pistak.gui.screen import fit_to_screen
from pistak.i18n import tf, tr
from pistak.io import compare as cmp
from pistak.io.batch import StructureEntry, find_cifs, run_batch, write_summary_csv
from pistak.io.stats import PairTable, correlation, describe, summary_rows, write_describe_csv
from pistak.params import get_spec

KIND_FILTER = [("Solo π-stack", "stack"), ("Parallele", "parallel"), ("Inclinate", "tilted"), ("Tutti i tipi", None)]
HEATMAP_DEFAULT = ["interplanar", "long_slip", "lat_slip", "dihedral", "centroid_distance", "contact", "overlap_area", "beta_a"]
BY_KIND, BY_FILE = "kind", "file"  # valori interni (non tradotti) di «colora/impila per…»
DIVERGING = pg.ColorMap([0.0, 0.5, 1.0], [(46, 111, 214), (245, 245, 245), (214, 69, 65)])


class _Signals(QObject):
    progress = Signal(int, int, str)
    done = Signal(object)


class _Job(QRunnable):
    def __init__(self, files: list[Path], limits: Limits):
        super().__init__()
        self.files, self.limits = files, limits
        self.signals = _Signals()

    def run(self) -> None:
        try:
            res = run_batch(self.files, self.limits, 0, progress=lambda i, n, name: self.signals.progress.emit(i, n, name))
        except Exception as e:  # noqa: BLE001  (la finestra non deve restare ad aspettare per sempre)
            res = [StructureEntry(str(f), f.name, error=f"{type(e).__name__}: {e}") for f in self.files]
        self.signals.done.emit(res)


def _title(key: str) -> str:
    s = get_spec(key)
    return s.title if s else key


def _num(v: float, digits: int = 2) -> str:
    return "–" if v is None or not math.isfinite(v) else f"{v:.{digits}f}"


class CollectionWindow(QMainWindow):
    #: percorso del file e indice della coppia che l'utente vuole vedere nella finestra principale
    openRequested = Signal(str, int, str)  # file, indice della coppia, blocco data_ ("" = il primo)
    changed = Signal()  # l'elenco dei file è cambiato (per il salvataggio della sessione)

    def __init__(self, limits: Limits, parent=None):
        super().__init__(parent)
        self.setWindowTitle(tr("Confronto e statistiche"))
        fit_to_screen(self, 1150, 760, min_width=560, min_height=380)
        self.limits = limits
        self.entries: list[StructureEntry] = []
        self.reference = 0
        self.table = PairTable.from_entries([])
        self._job: _Job | None = None
        self._picked: tuple[str, int] | None = None

        # barra superiore ------------------------------------------------------------
        top = FlowLayout()  # a capo automatico se la finestra è stretta
        add = QPushButton(tr("Aggiungi CIF…"))
        add.clicked.connect(self.add_files_dialog)
        add_dir = QPushButton(tr("Aggiungi cartella…"))
        add_dir.clicked.connect(self.add_folder_dialog)
        rm = QPushButton(tr("Rimuovi selezionati"))
        rm.clicked.connect(self.remove_selected)
        clear = QPushButton(tr("Svuota"))
        clear.clicked.connect(self.clear)
        redo = QPushButton(tr("Rianalizza"))
        redo.setToolTip(tr("Rifà l'analisi di tutti i file con le soglie correnti"))
        redo.clicked.connect(self.reanalyse)
        self.kind_combo = QComboBox()
        for label, k in KIND_FILTER:
            self.kind_combo.addItem(tr(label), k)
        self.kind_combo.currentIndexChanged.connect(self.refresh_plots)
        self.weighted = QCheckBox(tr("Pesa per coppie equivalenti"))
        self.weighted.setChecked(True)
        self.weighted.setToolTip(tr("Ogni coppia vale quante coppie equivalenti per simmetria ha raccolto"))
        self.weighted.toggled.connect(self.refresh_plots)
        export = QToolButton()
        export.setText(tr("Esporta ▾"))
        export.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QMenu(export)
        menu.addAction(tr("Confronto delle strutture (CSV)…"), lambda: self._export("compare"))
        menu.addAction(tr("Tutte le coppie (CSV)…"), lambda: self._export("pairs"))
        menu.addAction(tr("Descrizione statistica (CSV)…"), lambda: self._export("describe"))
        menu.addAction(tr("Riepilogo per file (CSV)…"), lambda: self._export("summary"))
        export.setMenu(menu)
        for w in (add, add_dir, rm, clear, redo):
            top.addWidget(w)
        top.addWidget(QLabel(tr("Coppie")))
        top.addWidget(self.kind_combo)
        top.addWidget(self.weighted)
        top.addWidget(export)
        self.bar = QProgressBar()
        self.bar.setTextVisible(True)
        self.bar.hide()
        self.info = QLabel(tr("Aggiungi file o cartelle di CIF: le coppie di tutti i file finiscono nella stessa analisi."))

        # schede ----------------------------------------------------------------------
        self.tabs = QTabWidget()
        self._build_structures_tab()
        self._build_hist_tab()
        self._build_scatter_tab()
        self._build_corr_tab()
        self._build_summary_tab()
        root = QWidget()
        lay = QVBoxLayout(root)
        lay.addLayout(top)
        lay.addWidget(self.bar)
        lay.addWidget(self.tabs, 1)
        lay.addWidget(self.info)
        self.setCentralWidget(root)
        self.setAcceptDrops(True)

    # ---- dati ------------------------------------------------------------------------
    @property
    def paths(self) -> list[str]:
        return [e.path for e in self.entries]

    def add_paths(self, items: list[str | Path]) -> None:
        """Aggiunge file o cartelle (analisi in parallelo, senza bloccare la finestra)."""
        known = set(self.paths)
        files: list[Path] = []
        for it in items:
            files += [f for f in find_cifs(it) if str(f) not in known and f not in files]
        if not files or self._job is not None:
            return
        self.bar.setRange(0, len(files))
        self.bar.setValue(0)
        self.bar.show()
        self._job = _Job(files, self.limits)
        self._job.signals.progress.connect(lambda i, n, name: (self.bar.setValue(i), self.bar.setFormat(f"{i}/{n}  {name}")))
        self._job.signals.done.connect(self._finished)
        QThreadPool.globalInstance().start(self._job)

    def _finished(self, results: list[StructureEntry]) -> None:
        self._job = None
        self.bar.hide()
        self.entries += results
        bad = [r for r in results if r.error]
        self._rebuild()
        if bad:
            self.info.setText(tf("{n} file non letti: ", n=len(bad)) + "; ".join(f"{r.name} ({r.error})" for r in bad[:3]))

    def clear(self) -> None:
        self.entries, self.reference = [], 0
        self._rebuild()

    def remove_selected(self) -> None:
        rows = {i.row() for i in self.struct.selectedIndexes()}
        if rows:
            self.entries = [e for i, e in enumerate(self.entries) if i not in rows]
            self.reference = min(self.reference, max(0, len(self.entries) - 1))
            self._rebuild()

    def reanalyse(self) -> None:
        files = [Path(p) for p in self.paths]
        self.entries = []
        self._rebuild()
        self.add_paths(files)

    def set_limits(self, limits: Limits) -> None:
        self.limits = limits

    def add_files_dialog(self) -> None:
        names, _ = QFileDialog.getOpenFileNames(self, tr("Aggiungi file CIF"), "", tr("CIF (*.cif *.CIF);;Tutti i file (*)"))
        if names:
            self.add_paths(names)

    def add_folder_dialog(self) -> None:
        d = QFileDialog.getExistingDirectory(self, tr("Cartella con file CIF"))
        if d:
            self.add_paths([d])

    def dragEnterEvent(self, e) -> None:  # noqa: N802
        if e.mimeData().hasUrls():
            e.acceptProposedAction()

    def dropEvent(self, e) -> None:  # noqa: N802
        self.add_paths([u.toLocalFile() for u in e.mimeData().urls() if u.isLocalFile()])

    def _rebuild(self) -> None:
        ok = [e for e in self.entries if not e.error]
        self.table = PairTable.from_entries(ok)
        keys = self.table.keys
        for combo, want in ((self.hist_key, "interplanar"), (self.sx, "long_slip"), (self.sy, "interplanar")):
            cur = combo.currentData() or want
            combo.blockSignals(True)
            combo.clear()
            for k in keys:
                combo.addItem(_title(k), k)
            combo.setCurrentIndex(max(0, combo.findData(cur)))
            combo.blockSignals(False)
        cur_c = self.sc.currentData()
        self.sc.blockSignals(True)
        self.sc.clear()
        self.sc.addItem(tr("tipo di coppia"), BY_KIND)
        self.sc.addItem(tr("file"), BY_FILE)
        for k in keys:
            self.sc.addItem(_title(k), k)
        self.sc.setCurrentIndex(max(0, self.sc.findData(cur_c)))
        self.sc.blockSignals(False)
        checked = {self.corr_list.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.corr_list.count()) if self.corr_list.item(i).checkState() == Qt.CheckState.Checked}
        self.corr_list.blockSignals(True)
        self.corr_list.clear()
        for k in keys:
            it = QListWidgetItem(_title(k))
            it.setData(Qt.ItemDataRole.UserRole, k)
            it.setCheckState(Qt.CheckState.Checked if (k in checked if checked else k in HEATMAP_DEFAULT) else Qt.CheckState.Unchecked)
            self.corr_list.addItem(it)
        self.corr_list.blockSignals(False)
        self._fill_structures()
        self.refresh_plots()
        self.changed.emit()

    def refresh_plots(self) -> None:
        self._draw_hist()
        self._draw_scatter()
        self._draw_corr()
        self._fill_summary()
        n = len(self.table)
        self.info.setText(tf("{s} strutture · {n} coppie", s=len(self.entries), n=n)
                          + (tf(" (riferimento: {name})", name=self.entries[self.reference].name) if self.entries else ""))

    # ---- scheda Strutture --------------------------------------------------------------
    def _build_structures_tab(self) -> None:
        self.struct = QTableWidget()
        self.struct.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.struct.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.struct.setSortingEnabled(False)
        self.struct.cellDoubleClicked.connect(lambda r, _c: self._open_entry(r))
        self.struct.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.struct.customContextMenuRequested.connect(self._struct_menu)
        w = QWidget()
        lay = QVBoxLayout(w)
        hint = QLabel(tr("Per ogni struttura, la coppia π-stack più vicina. Doppio clic: apre il file nella finestra principale. "
                      "Tasto destro: imposta come riferimento. «Δ» è la differenza dal riferimento; «distanza motivo» è la distanza (Å) "
                      "in (interplanare, slittamento L, slittamento W) tra i due motivi di stack più simili."))
        hint.setWordWrap(True)
        lay.addWidget(hint)
        lay.addWidget(self.struct)
        self.tabs.addTab(w, tr("Strutture"))

    def _struct_menu(self, pos) -> None:
        row = self.struct.rowAt(pos.y())
        if row < 0:
            return
        m = QMenu(self)
        m.addAction(tr("Imposta come riferimento"), lambda: self._set_reference(row))
        m.addAction(tr("Apri nella finestra principale"), lambda: self._open_entry(row))
        m.exec(self.struct.viewport().mapToGlobal(pos))

    def _set_reference(self, row: int) -> None:
        self.reference = row
        self._fill_structures()
        self.changed.emit()

    def _open_entry(self, row: int) -> None:
        e = self.entries[row]
        if not e.error:
            best = next((p.index for p in e.pairs if p.kind == "stack"), 0)
            self.openRequested.emit(e.path, best, e.block or "")

    def _fill_structures(self) -> None:
        rows = cmp.compare(self.entries, self.reference) if self.entries else []
        keys = ["interplanar", "long_slip", "lat_slip", "dihedral"]
        heads = (["", "File", tr("Gruppo"), tr("Molecole"), "π-stack", tr("Colonne"), tr("Interpl. (Å)"), tr("Slitt. L (Å)"), tr("Slitt. W (Å)"), tr("Angolo (°)"),
                  tr("Δ interpl."), tr("Δ slitt. L"), tr("Δ slitt. W"), tr("Distanza motivo (Å)")])
        self.struct.clear()
        self.struct.setColumnCount(len(heads))
        self.struct.setHorizontalHeaderLabels(heads)
        self.struct.setRowCount(len(rows))
        for r, row in enumerate(rows):
            cells = ["★" if r == self.reference else "", row.name, row.space_group, row.formulas, str(row.n_stack), str(row.n_columns)]
            if row.error:
                cells += [row.error] + [""] * (len(heads) - len(cells) - 1)
            else:
                cells += [_num(row.best.get(k, math.nan)) for k in keys]
                cells += [("–" if k not in row.delta else f"{row.delta[k]:+.2f}") for k in keys[:3]]
                cells += [_num(row.motif_distance)]
            for c, v in enumerate(cells):
                it = QTableWidgetItem(v)
                if c >= 4:
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                if row.error:
                    it.setForeground(QBrush(QColor(200, 60, 60)))
                self.struct.setItem(r, c, it)
        self.struct.resizeColumnsToContents()

    # ---- raggruppamenti per i grafici ----------------------------------------------------
    def _mask(self) -> np.ndarray:
        return self.table.mask(self.kind_combo.currentData())

    def _groups(self, by: str, mask: np.ndarray) -> list[tuple[str, tuple, np.ndarray]]:
        """(etichetta, colore RGB, maschera) per ogni gruppo di coppie."""
        t = self.table
        if by == BY_FILE:
            names = list(dict.fromkeys(t.files))
            return [(n, pg.intColor(i, hues=max(len(names), 6)).getRgb()[:3], mask & (np.array(t.files, dtype=str) == n)) for i, n in enumerate(names)]
        return [(theme.kind_label(k), theme.KIND_COLORS[k], mask & (t.kinds == k)) for k in ("stack", "parallel", "tilted")]

    # ---- scheda Distribuzioni ---------------------------------------------------------
    def _build_hist_tab(self) -> None:
        self.hist_key = QComboBox()
        self.hist_bins = QSpinBox()
        self.hist_bins.setRange(3, 100)
        self.hist_bins.setValue(20)
        self.hist_by = QComboBox()
        self.hist_by.addItem(tr("tipo di coppia"), BY_KIND)
        self.hist_by.addItem(tr("file"), BY_FILE)
        for w in (self.hist_key, self.hist_by):
            w.currentIndexChanged.connect(self._draw_hist)
        self.hist_bins.valueChanged.connect(self._draw_hist)
        self.hist_plot = pg.PlotWidget(background=theme.background())
        self.hist_plot.showGrid(x=False, y=True, alpha=0.2)
        w = QWidget()
        lay = QVBoxLayout(w)
        bar = QHBoxLayout()
        for lab, wid in ((tr("Parametro"), self.hist_key), (tr("Intervalli"), self.hist_bins), (tr("Impila per"), self.hist_by)):
            bar.addWidget(QLabel(lab))
            bar.addWidget(wid)
        bar.addStretch(1)
        lay.addLayout(bar)
        lay.addWidget(self.hist_plot, 1)
        self.tabs.addTab(w, tr("Distribuzioni"))

    def _draw_hist(self) -> None:
        pl = self.hist_plot
        pl.clear()
        pl.setBackground(theme.background())
        key = self.hist_key.currentData()
        if not key or not len(self.table):
            return
        mask = self._mask()
        x = self.table.column(key)
        ok = mask & np.isfinite(x)
        if not ok.any():
            return
        edges = np.histogram_bin_edges(x[ok], bins=self.hist_bins.value())
        centers, width = (edges[:-1] + edges[1:]) / 2, (edges[1] - edges[0]) * 0.92
        w = self.table.weights(self.weighted.isChecked())
        bottom = np.zeros(len(centers))
        pl.addLegend(offset=(10, 10))
        for label, rgb, gm in self._groups(self.hist_by.currentData(), mask):
            sel = gm & np.isfinite(x)
            h, _ = np.histogram(x[sel], bins=edges, weights=w[sel])
            if h.sum():
                item = pg.BarGraphItem(x=centers, height=h, width=width, y0=bottom, brush=theme.qcolor(rgb), pen=pg.mkPen(None), name=label)
                pl.addItem(item)  # con `name` entra da sola nella legenda
            bottom = bottom + h
        pl.setLabel("bottom", _title(key))
        pl.setLabel("left", tr("coppie") + (tr(" (con equivalenti)") if self.weighted.isChecked() else ""))

    # ---- scheda Dispersione -----------------------------------------------------------
    def _build_scatter_tab(self) -> None:
        self.sx, self.sy, self.sc = QComboBox(), QComboBox(), QComboBox()
        for c in (self.sx, self.sy, self.sc):
            c.currentIndexChanged.connect(self._draw_scatter)
        self.scatter_plot = pg.PlotWidget(background=theme.background())
        self.scatter_plot.showGrid(x=True, y=True, alpha=0.2)
        self.cloud = pg.ScatterPlotItem(size=10, pxMode=True, hoverable=True, hoverSize=14, tip=self._tip)
        self.cloud.sigClicked.connect(self._picked_point)
        self.pick_marker = pg.ScatterPlotItem(size=20, pxMode=True, brush=QBrush(Qt.BrushStyle.NoBrush), pen=pg.mkPen(theme.foreground(), width=2.5))
        self.scatter_plot.addItem(self.cloud)
        self.scatter_plot.addItem(self.pick_marker)
        self.pick_label = QLabel(tr("Clic su un punto per sceglierlo."))
        self.open_btn = QPushButton(tr("Apri nella finestra principale"))
        self.open_btn.setEnabled(False)
        self.open_btn.clicked.connect(self._open_picked)
        w = QWidget()
        lay = QVBoxLayout(w)
        bar = QHBoxLayout()
        for lab, wid in (("X", self.sx), ("Y", self.sy), (tr("Colore"), self.sc)):
            bar.addWidget(QLabel(lab))
            bar.addWidget(wid)
        bar.addStretch(1)
        lay.addLayout(bar)
        lay.addWidget(self.scatter_plot, 1)
        foot = QHBoxLayout()
        foot.addWidget(self.pick_label, 1)
        foot.addWidget(self.open_btn)
        lay.addLayout(foot)
        self.tabs.addTab(w, tr("Dispersione"))

    def _tip(self, x, y, data) -> str:
        if data is None or data >= len(self.table):
            return ""
        t = self.table
        return (f"{t.files[data]} · " + tf("coppia {n}", n=t.index[data] + 1) + f" · {theme.kind_label(str(t.kinds[data]))}\n"
                f"X = {x:.3g}   Y = {y:.3g}")

    def _draw_scatter(self) -> None:
        self.pick_marker.setData([], [])
        self.scatter_plot.setBackground(theme.background())
        kx, ky, by = self.sx.currentData(), self.sy.currentData(), self.sc.currentData()
        if not (kx and ky and len(self.table)):
            self.cloud.setData([], [])
            return
        t = self.table
        mask = self._mask()
        x, y = t.column(kx), t.column(ky)
        ok = mask & np.isfinite(x) & np.isfinite(y)
        idx = np.flatnonzero(ok)
        if by in (BY_KIND, BY_FILE):
            colors = np.zeros((len(t), 3), int)
            for _label, rgb, gm in self._groups(by, mask):
                colors[gm] = rgb
        else:
            c = t.column(by)
            lo, hi = np.nanmin(c[ok]) if ok.any() else 0, np.nanmax(c[ok]) if ok.any() else 1
            cmap = pg.colormap.get("viridis")
            tt = np.where(np.isfinite(c), (c - lo) / ((hi - lo) or 1), 0)
            colors = np.array([cmap.map(float(v), mode="byte")[:3] for v in tt]).astype(int)
        brushes = [QBrush(theme.qcolor(tuple(int(v) for v in colors[i]), 220)) for i in idx]
        sizes = [8 + min(10, 2.5 * math.log2(max(1, int(t.counts[i])))) for i in idx]
        self.cloud.setData(x[idx], y[idx], brush=brushes, pen=[theme.point_pen(tuple(int(v) for v in colors[i]), 220) for i in idx], size=sizes, data=[int(i) for i in idx])
        self.scatter_plot.setLabel("bottom", _title(kx))
        self.scatter_plot.setLabel("left", _title(ky))
        self.scatter_plot.autoRange()

    def _picked_point(self, _item, points, *_) -> None:
        if not len(points):
            return
        i = int(points[0].data())
        t = self.table
        self._picked = (t.paths[i], int(t.index[i]), t.blocks[i] if t.blocks else "")
        self.pick_marker.setData([points[0].pos().x()], [points[0].pos().y()])
        self.pick_label.setText(f"<b>{t.files[i]}</b> · " + tf("coppia {n}", n=t.index[i] + 1)
                                + f" · {theme.kind_label(str(t.kinds[i]))} · ×{t.counts[i]}")
        self.open_btn.setEnabled(True)

    def _open_picked(self) -> None:
        if self._picked:
            self.openRequested.emit(*self._picked)

    # ---- scheda Correlazioni ----------------------------------------------------------
    def _build_corr_tab(self) -> None:
        self.corr_list = QListWidget()
        self.corr_list.setMaximumWidth(260)
        self.corr_list.itemChanged.connect(lambda _i: self._draw_corr())
        self.corr_plot = pg.PlotWidget(background=theme.background())
        self.corr_plot.setAspectLocked(True)
        self.corr_plot.invertY(True)
        w = QWidget()
        lay = QHBoxLayout(w)
        side = QVBoxLayout()
        side.addWidget(QLabel(tr("Parametri da correlare")))
        side.addWidget(self.corr_list, 1)
        lay.addLayout(side)
        lay.addWidget(self.corr_plot, 1)
        self.tabs.addTab(w, tr("Correlazioni"))

    def _draw_corr(self) -> None:
        pl = self.corr_plot
        pl.clear()
        pl.setBackground(theme.background())
        keys = [self.corr_list.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.corr_list.count())
                if self.corr_list.item(i).checkState() == Qt.CheckState.Checked]
        if len(keys) < 2 or not len(self.table):
            return
        c = correlation(self.table, keys, self.kind_combo.currentData())
        n = len(keys)
        img = pg.ImageItem(np.nan_to_num(c, nan=0.0), levels=(-1, 1))
        img.setLookupTable(DIVERGING.getLookupTable(0.0, 1.0, 256))
        pl.addItem(img)
        for i in range(n):
            for j in range(n):
                v = c[i, j]
                txt = pg.TextItem("–" if not math.isfinite(v) else f"{v:.2f}", color=(20, 20, 20) if not math.isfinite(v) or abs(v) < 0.75 else (255, 255, 255), anchor=(0.5, 0.5))
                txt.setPos(i + 0.5, j + 0.5)
                pl.addItem(txt)
        # sotto i numeri (le etichette lunghe si sovrapporrebbero), a sinistra «numero · nome»
        pl.getAxis("bottom").setTicks([[(i + 0.5, str(i + 1)) for i in range(n)]])
        pl.getAxis("left").setTicks([[(i + 0.5, f"{i + 1} · {_title(k)[:26]}") for i, k in enumerate(keys)]])
        pl.getAxis("left").setWidth(190)
        pl.setXRange(0, n)
        pl.setYRange(0, n)

    # ---- scheda Riepilogo --------------------------------------------------------------
    def _build_summary_tab(self) -> None:
        self.sum_table = QTableWidget()
        self.sum_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.file_table = QTableWidget()
        self.file_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.addWidget(QLabel(tr("<b>Descrizione statistica</b> delle coppie del tipo scelto (tutte le strutture insieme)")))
        lay.addWidget(self.sum_table, 3)
        lay.addWidget(QLabel(tr("<b>Per file</b>")))
        lay.addWidget(self.file_table, 2)
        self.tabs.addTab(w, tr("Riepilogo"))

    def describe_rows(self) -> list[dict]:
        return describe(self.table, None, self.kind_combo.currentData(), weighted=self.weighted.isChecked()) if len(self.table) else []

    def _fill_summary(self) -> None:
        rows = self.describe_rows()
        heads = [tr("Parametro"), "n", tr("Media"), tr("Dev. st."), tr("Min"), "25%", tr("Mediana"), "75%", tr("Max")]
        self.sum_table.clear()
        self.sum_table.setColumnCount(len(heads))
        self.sum_table.setHorizontalHeaderLabels(heads)
        self.sum_table.setRowCount(len(rows))
        for r, d in enumerate(rows):
            vals = [_title(d["key"]), str(d["n"])] + [_num(d.get(k, math.nan), 3) for k in ("mean", "std", "min", "q25", "median", "q75", "max")]
            for c, v in enumerate(vals):
                it = QTableWidgetItem(v)
                if c:
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.sum_table.setItem(r, c, it)
        self.sum_table.resizeColumnsToContents()
        srows = summary_rows([e for e in self.entries])
        heads2 = ["File", tr("Coppie"), "π-stack", "parallele", "inclinate", tr("Colonne"), tr("Interpl. medio"), tr("Slitt. L medio"), tr("Slitt. W medio"), tr("Angolo medio")]
        self.file_table.clear()
        self.file_table.setColumnCount(len(heads2))
        self.file_table.setHorizontalHeaderLabels(heads2)
        self.file_table.setRowCount(len(srows))
        for r, d in enumerate(srows):
            vals = [d["file"] if not d["error"] else f"{d['file']} — {d['error']}", str(d["pairs"]), str(d["n_stack"]), str(d["n_parallel"]),
                    str(d["n_tilted"]), str(d["columns"])] + [_num(d[f"mean_stack_{k}"]) for k in ("interplanar", "long_slip", "lat_slip", "dihedral")]
            for c, v in enumerate(vals):
                it = QTableWidgetItem(v)
                if c:
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.file_table.setItem(r, c, it)
        self.file_table.resizeColumnsToContents()

    # ---- esportazione -----------------------------------------------------------------
    def _export(self, what: str) -> None:
        if not self.entries:
            QMessageBox.information(self, "PISTAK", tr("Aggiungi prima qualche file."))
            return
        name = {"compare": tr("confronto.csv"), "pairs": tr("coppie.csv"), "describe": tr("descrizione.csv"), "summary": tr("riepilogo.csv")}[what]
        path, _ = QFileDialog.getSaveFileName(self, tr("Esporta"), name, "CSV (*.csv)")
        if not path:
            return
        if what == "compare":
            cmp.write_csv(cmp.compare(self.entries, self.reference), path)
        elif what == "pairs":
            self.table.write_csv(path)
        elif what == "summary":
            write_summary_csv(self.entries, path)
        else:
            write_describe_csv(self.describe_rows(), path)
        self.info.setText(tf("Scritto {path}", path=path))
