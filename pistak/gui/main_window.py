"""Finestra principale: un insieme di viste agganciabili (dock) che condividono la `Session`.

Ogni vista può essere trascinata fuori dalla finestra, affiancata, impilata a schede o duplicata
(Vista ▸ Nuova vista) e ha il suo selettore di sincronizzazione. La finestra si adatta allo schermo e si può
ridimensionare liberamente; la lingua (italiano/inglese) si cambia dal menu e la finestra si ricostruisce da sola.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import QByteArray, QSettings, Qt, QTimer
from PySide6.QtGui import QAction, QActionGroup, QDragEnterEvent, QDropEvent, QKeySequence, QResizeEvent, QShowEvent
from PySide6.QtWidgets import (
    QApplication,
    QDialog,
    QDockWidget,
    QFileDialog,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QTabBar,
    QTextBrowser,
    QVBoxLayout,
)

from pistak import i18n
from pistak.core.analysis import read_cif_file
from pistak.core.cif import block_names
from pistak.core.config import ROOT, Limits, user_config_dir
from pistak.gui import session_file, theme
from pistak.gui.screen import clamp_into_screen, fit_to_screen
from pistak.gui.session import GROUPS, Session
from pistak.gui.settings_dialog import SettingsDialog
from pistak.gui.views import LEGACY_KINDS, VIEW_TYPES, view_title
from pistak.gui.views.base import PairView
from pistak.i18n import tf, tr
from pistak.io import export
from pistak.params import load_errors, plugin_dirs, specs

APP_NAME = "PISTAK"
DEFAULT_SIZE = (1500, 920)


def last_session_path() -> Path:
    """Dove si salva da sola l'ultima sessione, alla chiusura (cartella di configurazione dell'utente)."""
    return user_config_dir() / "last.pistak"


class MainWindow(QMainWindow):
    def __init__(self, session: Session):
        super().__init__()
        self.session = session
        self.setWindowTitle(tf("{app} — analisi dei π-stack", app=APP_NAME))
        self.setDockNestingEnabled(True)
        self.setAcceptDrops(True)
        self.setCorner(Qt.Corner.BottomLeftCorner, Qt.DockWidgetArea.LeftDockWidgetArea)
        self.docks: list[QDockWidget] = []
        self._counter: dict[str, int] = {}
        self.settings = QSettings("PISTAK", "PISTAK")
        self.collection = None  # finestra «Confronto e statistiche», creata alla prima richiesta
        self.session_path: str | None = None
        self.autosave = True
        app = QApplication.instance()
        if app is not None:
            app.main_window = self  # type: ignore[attr-defined]  # tiene viva la finestra dopo un cambio di lingua

        self.status_label = QLabel("")
        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setMaximumWidth(140)
        self.progress.hide()
        self.statusBar().addWidget(self.status_label, 1)
        self.statusBar().addPermanentWidget(self.progress)
        session.message.connect(self.status_label.setText)
        session.error.connect(self._show_error)
        session.busyChanged.connect(self._busy)
        session.analysisChanged.connect(self._analysis_changed)
        self._tabbar_timer = QTimer(self)
        self._tabbar_timer.setSingleShot(True)
        self._tabbar_timer.timeout.connect(self._hide_stray_tab_bars)

        self._build_menus()
        self.default_layout()
        # dopo la disposizione: aggiunto prima, Qt lasciava visibile una barra di schede vuota sopra lo scatter
        self.limits_dock = self._limits_dock()
        self.fit_to_screen()  # dopo aver creato i pannelli: aggiungerli può far crescere la finestra oltre lo schermo

    # ---- viste e disposizione ------------------------------------------------------
    def add_view(self, kind: str, group: str | None = "A", object_name: str | None = None, **kw) -> QDockWidget:
        kind = LEGACY_KINDS.get(kind, kind)
        view: PairView = VIEW_TYPES[kind](self.session, group, **kw)
        n = self._counter[kind] = self._counter.get(kind, 0) + 1
        dock = QDockWidget(view_title(kind) if n == 1 else f"{view_title(kind)} {n}", self)
        dock.setObjectName(object_name or f"{kind}-{n}")
        dock.setProperty("kind", kind)
        dock.setWidget(view)
        dock.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
        dock.destroyed.connect(lambda _=None, d=dock: self.docks.remove(d) if d in self.docks else None)
        dock.dockLocationChanged.connect(self._check_tab_bars)
        dock.topLevelChanged.connect(self._check_tab_bars)
        self.docks.append(dock)
        return dock

    def _limits_dock(self) -> QDockWidget:
        """Il pannello dei cursori delle soglie: fuori da `self.docks` (non è una vista delle coppie), nascosto
        all'inizio, si apre da Analisi ▸ Soglie in tempo reale."""
        from pistak.gui.limits_panel import LimitsPanel

        dock = QDockWidget(tr("Soglie in tempo reale"), self)
        dock.setObjectName("limits-panel")
        dock.setWidget(LimitsPanel(self.session))
        dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)
        dock.setFloating(True)
        dock.hide()
        return dock

    def show_limits_panel(self) -> None:
        d = self.limits_dock
        if d.isFloating() and not d.isVisible():
            d.resize(420, 360)
            clamp_into_screen(d)
        d.show()
        d.raise_()

    # ---- barre di schede orfane (difetto di Qt 6) ---------------------------------------
    def _check_tab_bars(self, *_) -> None:
        self._tabbar_timer.start(0)

    def _hide_stray_tab_bars(self) -> None:
        """Con due gruppi di pannelli a schede Qt 6 crea a volte una barra di schede doppia che, dopo un
        ridimensionamento, resta ferma dove stava, sopra un altro pannello. Si nasconde ogni barra che non è attaccata
        (sopra o sotto, allineata a sinistra) a uno dei pannelli che elenca: quelle vere lo sono sempre."""
        docks = [d for d in self.findChildren(QDockWidget) if d.isVisible() and not d.isFloating()]
        for bar in self.findChildren(QTabBar, options=Qt.FindChildOption.FindDirectChildrenOnly):
            if not bar.isVisible() or bar.count() == 0:
                continue
            titles = {bar.tabText(i) for i in range(bar.count())}
            g = bar.geometry()
            attached = any(abs(g.left() - d.geometry().left()) <= 6
                           and (abs(g.top() - d.geometry().bottom()) <= 6 or abs(g.bottom() - d.geometry().top()) <= 6)
                           for d in docks if d.windowTitle() in titles)
            if not attached:
                bar.hide()

    def showEvent(self, e: QShowEvent) -> None:  # noqa: N802
        super().showEvent(e)
        self._check_tab_bars()

    def resizeEvent(self, e: QResizeEvent) -> None:  # noqa: N802
        super().resizeEvent(e)
        self._check_tab_bars()

    def clear_docks(self) -> None:
        for d in list(self.docks):
            self.removeDockWidget(d)
            # deleteLater è differito: finché esiste, un vecchio pannello con lo stesso objectName intercetta
            # restoreState() e il layout salvato finisce su un pannello che sta per sparire
            d.setObjectName("")
            d.setParent(None)
            d.deleteLater()
        self.docks.clear()
        self._counter.clear()

    def fit_to_screen(self) -> None:
        """Porta la finestra alla dimensione preferita, ma mai oltre lo schermo; poi ridistribuisce i pannelli in proporzione."""
        self._size = fit_to_screen(self, *DEFAULT_SIZE)
        self._proportional_docks()

    def default_layout(self) -> None:
        self.clear_docks()
        table, scatter, hist = self.add_view("table"), self.add_view("scatter"), self.add_view("histogram")
        view3d, proj, detail = self.add_view("3d"), self.add_view("projection"), self.add_view("detail")
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, table)
        self.splitDockWidget(table, scatter, Qt.Orientation.Vertical)
        self.tabifyDockWidget(scatter, hist)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, view3d)
        self.splitDockWidget(view3d, proj, Qt.Orientation.Vertical)
        self.tabifyDockWidget(proj, detail)
        scatter.raise_()
        proj.raise_()
        self._layout_docks = (table, scatter, view3d, proj)
        self._proportional_docks()
        self._check_tab_bars()

    def _proportional_docks(self) -> None:
        """Dimensioni dei pannelli in proporzione alla finestra, non in pixel assoluti: vale per qualunque schermo."""
        if not getattr(self, "_layout_docks", None) or any(d not in self.docks for d in self._layout_docks):
            return
        table, scatter, view3d, proj = self._layout_docks
        w, h = getattr(self, "_size", (self.width(), self.height()))
        self.resizeDocks([table, view3d], [int(w * 0.42), int(w * 0.58)], Qt.Orientation.Horizontal)
        self.resizeDocks([table, scatter], [int(h * 0.36), int(h * 0.52)], Qt.Orientation.Vertical)
        self.resizeDocks([view3d, proj], [int(h * 0.54), int(h * 0.34)], Qt.Orientation.Vertical)

    def refresh_views(self) -> None:
        """Dà alle viste i dati già presenti (dopo averle create a analisi già fatta)."""
        for d in self.docks:
            v = d.widget()
            if isinstance(v, PairView):
                v.refresh()

    def _new_view(self, kind: str) -> None:
        dock = self.add_view(kind, "A")
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        dock.setFloating(True)
        dock.resize(min(640, self.width() - 40), min(480, self.height() - 40))
        dock.show()
        if self.session.analysis:
            dock.widget().refresh()

    # ---- menu ----------------------------------------------------------------------
    def _act(self, menu, text: str, slot, shortcut: str | QKeySequence.StandardKey | None = None) -> QAction:
        a = QAction(text, self)
        if shortcut:
            a.setShortcut(shortcut)
        a.triggered.connect(lambda _=False: slot())
        menu.addAction(a)
        return a

    def _build_menus(self) -> None:
        mb = self.menuBar()
        f = mb.addMenu(tr("&File"))
        self._act(f, tr("Apri CIF…"), self.open_dialog, QKeySequence.StandardKey.Open)
        self._act(f, tr("Apri il cristallo dimostrativo"), self.session.load_demo)
        self.recent_menu = f.addMenu(tr("Recenti"))
        self._refresh_recent()
        f.addSeparator()
        self._act(f, tr("Salva la sessione…"), self.save_session_dialog, QKeySequence.StandardKey.Save)
        self._act(f, tr("Apri una sessione…"), self.open_session_dialog, "Ctrl+Shift+O")
        self._act(f, tr("Riapri l'ultima sessione"), lambda: self.load_session(last_session_path()))
        f.addSeparator()
        ex = f.addMenu(tr("Esporta"))
        self._act(ex, tr("Tabella delle coppie (CSV)…"), lambda: self._export("csv"))
        self._act(ex, tr("Risultati completi (JSON)…"), lambda: self._export("json"))
        self._act(ex, tr("Coppia selezionata (XYZ)…"), lambda: self._export("xyz"))
        self._act(ex, tr("Coppia selezionata per calcoli quantistici (ORCA, Psi4, Gaussian)…"), self._export_qc)
        self._act(ex, tr("Cella per calcoli periodici (VASP, Quantum ESPRESSO, CP2K, CRYSTAL)…"), self._export_periodic)
        self._act(ex, tr("Geometria delle molecole (CSV)…"), lambda: self._export("geom"))
        self._act(ex, tr("Avvisi dei controlli rapidi (CSV)…"), lambda: self._export("alerts"))
        self._act(ex, tr("Vista 3D (PNG)…"), lambda: self._export("png"))
        f.addSeparator()
        self._act(f, tr("Esci"), self.close, QKeySequence.StandardKey.Quit)

        v = mb.addMenu(tr("&Vista"))
        new = v.addMenu(tr("Nuova vista (finestra staccata)"))
        for kind in VIEW_TYPES:
            self._act(new, view_title(kind), lambda k=kind: self._new_view(k))
        self._act(v, tr("Reimposta la disposizione"), lambda: (self.default_layout(), self.fit_to_screen()))
        theme_menu = v.addMenu(tr("Tema"))
        group = QActionGroup(self)
        group.setExclusive(True)
        for mode, name in theme.MODES.items():
            act = QAction(tr(name), self, checkable=True)
            act.setChecked(mode == theme.mode())
            act.triggered.connect(lambda _=False, m=mode: QTimer.singleShot(0, lambda: self.change_theme(m)))
            group.addAction(act)
            theme_menu.addAction(act)

        a = mb.addMenu(tr("&Analisi"))
        self._act(a, tr("Soglie…"), self.edit_limits, "Ctrl+,")
        self._act(a, tr("Soglie in tempo reale"), self.show_limits_panel, "Ctrl+Shift+L")
        self._act(a, tr("Rianalizza"), self.session.reanalyse, "Ctrl+R")
        self._act(a, tr("Incertezze (Monte Carlo)…"), self._uncertainty)
        a.addSeparator()
        self._act(a, tr("Ricarica i plugin"), self._reload_plugins)
        self._act(a, tr("Apri la cartella dei plugin"), self._open_plugins)
        self._act(a, tr("Elenco dei parametri…"), self._show_params)
        a.addSeparator()
        self._act(a, tr("Confronto e statistiche…"), self.show_collection, "Ctrl+Shift+C")
        self._act(a, tr("Elabora una cartella di CIF…"), lambda: self.show_collection(ask_folder=True))

        lang_menu = mb.addMenu(tr("&Lingua"))
        group = QActionGroup(self)
        group.setExclusive(True)
        for code, name in i18n.LANGS.items():
            act = QAction(name, self, checkable=True)
            act.setChecked(code == i18n.get_language())
            # differito: la finestra sta per essere sostituita e non può farlo dall'interno del proprio menu
            act.triggered.connect(lambda _=False, c=code: QTimer.singleShot(0, lambda: self.change_language(c)))
            group.addAction(act)
            lang_menu.addAction(act)

        h = mb.addMenu(tr("&Aiuto"))
        self._act(h, tr("Manuale"), self.show_manual, QKeySequence.StandardKey.HelpContents)
        self._act(h, tr("Informazioni"), self._about)

    # ---- file ----------------------------------------------------------------------
    def open_dialog(self) -> None:
        start = str(self.settings.value("last_dir", str(ROOT / "examples")))
        path, _ = QFileDialog.getOpenFileName(self, tr("Apri un file CIF"), start, tr("CIF (*.cif *.CIF);;Tutti i file (*)"))
        if path:
            self.open_file(path)

    def open_file(self, path: str, block: str | None = None) -> None:
        self.settings.setValue("last_dir", str(Path(path).parent))
        recent = [p for p in self.settings.value("recent", [], type=list) if p != path][:7]
        self.settings.setValue("recent", [path, *recent])
        self._refresh_recent()
        if block is None:
            block = self._choose_block(path)
        self.session.load_file(path, block)

    def _choose_block(self, path: str) -> str | None:
        """Se il file ha più blocchi data_ con atomi, chiede quale analizzare (None = il primo)."""
        try:
            names = block_names(read_cif_file(path))
        except Exception:  # noqa: BLE001  (l'errore lo mostrerà l'analisi)
            return None
        if len(names) <= 1 or not self.isVisible():
            return None
        name, ok = QInputDialog.getItem(self, APP_NAME, tf("{file} contiene {n} strutture: quale analizzare?",
                                                            file=Path(path).name, n=len(names)), names, 0, False)
        return name if ok else names[0]

    def _refresh_recent(self) -> None:
        self.recent_menu.clear()
        for p in self.settings.value("recent", [], type=list):
            if Path(p).exists():
                self._act(self.recent_menu, Path(p).name, lambda q=p: self.open_file(q))

    def dragEnterEvent(self, e: QDragEnterEvent) -> None:  # noqa: N802
        if e.mimeData().hasUrls():
            e.acceptProposedAction()

    def dropEvent(self, e: QDropEvent) -> None:  # noqa: N802
        for u in e.mimeData().urls():
            if u.isLocalFile():
                f = u.toLocalFile()
                self.load_session(f) if f.endswith(session_file.SUFFIX) else self.open_file(f)
                break

    def _export(self, kind: str) -> None:
        a = self.session.analysis
        if not a:
            return
        base = Path(a.name).stem or "pistak"
        if kind == "png":
            for d in self.docks:
                w = d.widget()
                if hasattr(w, "screenshot"):
                    path, _ = QFileDialog.getSaveFileName(self, tr("Salva la vista 3D"), f"{base}_3d.png", "PNG (*.png)")
                    if path:
                        w.screenshot(path)
                    return
            return
        suffix = {"csv": "csv", "json": "json", "xyz": "xyz", "geom": "csv", "alerts": "csv"}[kind]
        name = {"geom": f"{base}_geometry", "alerts": f"{base}_alerts"}.get(kind, base)
        path, _ = QFileDialog.getSaveFileName(self, tr("Esporta"), f"{name}.{suffix}", f"*.{suffix}")
        if not path:
            return
        if kind == "csv":
            export.write_csv(a, path)
        elif kind == "json":
            export.write_json(a, path)
        elif kind == "geom":
            export.write_geometry_csv(a, path)
        elif kind == "alerts":
            export.write_alerts_csv(a, path)
        else:
            idx = self.session.selection("A")
            if not 0 <= idx < len(a.pairs):
                QMessageBox.information(self, APP_NAME, tr("Seleziona prima una coppia."))
                return
            export.write_pair_xyz(a, a.pairs[idx], path)
        self.status_label.setText(tf("Scritto {path}", path=path))

    def _export_qc(self) -> None:
        from pistak.errors import describe_error
        from pistak.gui.qc_dialog import QCDialog
        from pistak.io import qc

        a = self.session.analysis
        idx = self.session.selection("A")
        if not a or not 0 <= idx < len(a.pairs):
            QMessageBox.information(self, APP_NAME, tr("Seleziona prima una coppia."))
            return
        dlg = QCDialog(getattr(self, "_qc_options", None), self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        options = self._qc_options = dlg.options()
        pair = a.pairs[idx]
        try:
            text = qc.dimer_input(a, pair, options)
        except Exception as e:  # noqa: BLE001
            QMessageBox.warning(self, APP_NAME, describe_error(e))
            return
        suffix = qc.SUFFIX[options.program]
        path, _ = QFileDialog.getSaveFileName(self, tr("Esporta"), qc.suggested_name(a, pair, options), f"*.{suffix}")
        if path:
            Path(path).write_text(text, encoding="utf-8")
            self.status_label.setText(tf("Scritto {path}", path=path))

    def _export_periodic(self) -> None:
        from pistak.errors import describe_error
        from pistak.gui.periodic_dialog import PeriodicDialog
        from pistak.io import periodic

        a = self.session.analysis
        if not a:
            return
        dlg = PeriodicDialog(getattr(self, "_periodic_options", None), self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        options = self._periodic_options = dlg.options()
        try:
            files = periodic.periodic_files(a, options)
        except Exception as e:  # noqa: BLE001
            QMessageBox.warning(self, APP_NAME, describe_error(e))
            return
        start = str(self.settings.value("last_dir", str(ROOT)))
        folder = QFileDialog.getExistingDirectory(self, tr("Cartella in cui scrivere i file"), start)
        if not folder:
            return
        target = Path(folder) / f"{periodic.stem(a)}_{options.program}" if len(files) > 1 else Path(folder)
        target.mkdir(parents=True, exist_ok=True)
        for name, text in files.items():
            (target / name).write_text(text, encoding="utf-8")
        self.status_label.setText(tf("Scritto {path}", path=target if len(files) > 1 else target / next(iter(files))))

    # ---- analisi -------------------------------------------------------------------
    def edit_limits(self) -> None:
        dlg = SettingsDialog(self.session.limits, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self.session.set_limits(dlg.limits(), save=dlg.save)
            if self.collection is not None:
                self.collection.set_limits(dlg.limits())

    def _uncertainty(self) -> None:
        if self.session.analysis is None or self.session.busy:
            return
        n, ok = QInputDialog.getInt(self, APP_NAME, tr("Campioni del Monte Carlo (le σ vengono dalle s.u. del CIF):"),
                                    200, 10, 5000, 50)
        if ok:
            self.session.compute_uncertainty(n)

    def _reload_plugins(self) -> None:
        loaded = self.session.reload_plugins()
        if load_errors:
            QMessageBox.warning(self, APP_NAME, tr("Alcuni plugin non si caricano:") + "\n\n" + "\n".join(load_errors))
        else:
            self.status_label.setText(tf("Plugin ricaricati: {names}", names=", ".join(loaded) or tr("nessuno")))

    def _open_plugins(self) -> None:
        folder = plugin_dirs()[0]  # plugins/ del progetto se si lavora dai sorgenti, altrimenti quella dell'utente
        folder.mkdir(parents=True, exist_ok=True)
        if sys.platform == "darwin":
            subprocess.Popen(["open", str(folder)])
        elif os.name == "nt":
            os.startfile(folder)  # type: ignore[attr-defined]
        else:
            subprocess.Popen(["xdg-open", str(folder)])

    def _show_params(self) -> None:
        dlg = QDialog(self)
        dlg.setWindowTitle(tr("Parametri disponibili"))
        fit_to_screen(dlg, 760, 560, min_width=360, min_height=240, avail=None)
        box = QTextBrowser()
        rows = "".join(f"<tr><td><b>{s.key}</b></td><td>{s.title}</td><td>{s.group_name}</td><td>{s.source}</td><td>{s.info}</td></tr>"
                       for s in specs())
        head = "".join(f"<th align=left>{h}</th>" for h in (tr("chiave"), tr("nome"), tr("gruppo"), tr("origine"), tr("descrizione")))
        box.setHtml(f"<table cellpadding=4 border=0><tr>{head}</tr>{rows}</table>")
        lay = QVBoxLayout(dlg)
        lay.addWidget(box)
        dlg.exec()

    # ---- confronto e statistiche -----------------------------------------------------
    def show_collection(self, ask_folder: bool = False):
        from pistak.gui.collection_window import CollectionWindow

        if self.collection is None:
            self.collection = CollectionWindow(self.session.limits, self)
            self.collection.openRequested.connect(self.open_pair)
        self.collection.set_limits(self.session.limits)
        self.collection.show()
        self.collection.raise_()
        self.collection.activateWindow()
        if ask_folder:
            self.collection.add_folder_dialog()
        return self.collection

    def open_pair(self, path: str, index: int, block: str = "") -> None:
        """Apre un file (o un suo blocco data_) e seleziona una coppia (richiesta della finestra di confronto)."""
        same = (self.session.source_path and Path(self.session.source_path) == Path(path)
                and (self.session.source_block or "") == block and self.session.analysis)
        if same:
            self.session.select("A", min(index, len(self.session.analysis.pairs) - 1))
        else:
            self.open_file(path, block or None)
            self.after_analysis(lambda: self.session.select("A", min(index, max(0, len(self.session.analysis.pairs) - 1))))
        self.raise_()
        self.activateWindow()

    def after_analysis(self, fn) -> None:
        """Esegue `fn` una sola volta, dopo che le viste hanno ricevuto la prossima analisi."""

        def once() -> None:
            self.session.analysisChanged.disconnect(once)
            fn()

        self.session.analysisChanged.connect(once)

    # ---- sessione ---------------------------------------------------------------------
    def state(self) -> dict:
        s = self.session
        views = [{"kind": d.property("kind"), "object": d.objectName(), "title": d.windowTitle(), "state": d.widget().state()}
                 for d in self.docks if isinstance(d.widget(), PairView)]
        st: dict = {
            "source": ({**session_file.embed_source(s.source_name, s.source_path, s.source_text), "block": s.source_block}
                       if s.source_text else None),
            "limits": s.limits.to_dict(),
            "selection": {g: s.selection(g) for g in GROUPS},
            "subset": {g: (sorted(s.subset(g)) if s.subset(g) is not None else None) for g in GROUPS},
            "views": views,
            "window": {"geometry": session_file.b64(bytes(self.saveGeometry())), "dock_state": session_file.b64(bytes(self.saveState()))},
            "language": i18n.get_language(),
        }
        if self.collection is not None and self.collection.entries:
            st["collection"] = {"paths": self.collection.paths, "reference": self.collection.reference}
        return st

    def save_session(self, path: str | Path) -> Path:
        out = session_file.write(path, self.state())
        self.session_path = str(out)
        return out

    def save_session_dialog(self) -> None:
        stem = Path(self.session.source_name or tr("sessione")).stem
        start = self.session_path or str(Path(self.settings.value("last_dir", str(ROOT))) / f"{stem}{session_file.SUFFIX}")
        path, _ = QFileDialog.getSaveFileName(self, tr("Salva la sessione"), start, tf("Sessione PISTAK (*{ext})", ext=session_file.SUFFIX))
        if path:
            if not path.endswith(session_file.SUFFIX):
                path += session_file.SUFFIX
            self.save_session(path)
            self.status_label.setText(tf("Sessione salvata in {path}", path=path))

    def open_session_dialog(self) -> None:
        start = str(Path(self.session_path).parent) if self.session_path else str(ROOT)
        path, _ = QFileDialog.getOpenFileName(self, tr("Apri una sessione"), start, tf("Sessione PISTAK (*{ext})", ext=session_file.SUFFIX))
        if path:
            self.load_session(path)

    def _build_docks_from(self, st: dict) -> None:
        """Ricrea i pannelli salvati (nell'ordine e con i nomi originali), poi disposizione e impostazioni di ciascuno."""
        self.clear_docks()
        for v in st.get("views", []):
            kind = LEGACY_KINDS.get(v.get("kind"), v.get("kind"))
            if kind in VIEW_TYPES:
                v["kind"] = kind
                dock = self.add_view(kind, v.get("state", {}).get("group"), object_name=v.get("object"))
                self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)
        if not self.docks:
            self.default_layout()
            return
        win = st.get("window", {})
        if win.get("geometry"):
            self.restoreGeometry(QByteArray(session_file.unb64(win["geometry"])))
        if win.get("dock_state"):
            self.restoreState(QByteArray(session_file.unb64(win["dock_state"])))
            self._check_tab_bars()
        clamp_into_screen(self)  # la sessione può venire da uno schermo più grande
        for v in st.get("views", []):
            dock = next((d for d in self.docks if d.objectName() == v.get("object")), None)
            if dock is not None:
                dock.widget().restore(v.get("state", {}))

    def _apply_selection(self, st: dict) -> None:
        """Riapplica selezione e filtri di ogni gruppo (a analisi pronta) e le impostazioni che dipendono dai dati."""
        selection, subset = st.get("selection", {}), st.get("subset", {})
        n = len(self.session.analysis.pairs) if self.session.analysis else 0
        for g in GROUPS:
            i = selection.get(g, 0)
            self.session.select(g, i if 0 <= i < n else (0 if n else -1))
            sub = subset.get(g)
            self.session.set_subset(g, {i for i in sub if i < n} if sub is not None else None)
        for d in self.docks:
            st_v = next((v for v in st.get("views", []) if v.get("object") == d.objectName()), None)
            if st_v:
                d.widget().apply_state(st_v.get("state", {}).get("view", {}))
                d.widget().on_pair(d.widget().index)

    def load_session(self, path: str | Path) -> bool:
        """Ripristina viste, disposizione, soglie, file, selezione e filtri salvati. Ritorna False se il file non si legge."""
        try:
            st = session_file.read(path)
        except session_file.SessionFileError as e:
            self._show_error(str(e))
            return False
        self.session_path = str(path)
        self._build_docks_from(st)
        if st.get("limits"):
            self.session.limits = Limits.from_dict(st["limits"], self.session.limits)
            self.session.limitsChanged.emit()
        src = st.get("source") or {}
        loaded = True
        if src.get("text"):
            self.after_analysis(lambda: self._apply_selection(st))
            self.session.load_text(src["text"], src.get("name", "sessione.cif"), src.get("path"), src.get("block"))
        elif src.get("path") and Path(src["path"]).exists():
            self.after_analysis(lambda: self._apply_selection(st))
            self.session.load_file(src["path"], src.get("block"))
        else:
            loaded = False
            if src:
                self._show_error(tf("Il file della sessione non c'è più: {f}", f=src.get("path") or src.get("name")))
        self._restore_collection(st)
        self.status_label.setText(tf("Sessione ripristinata da {path}", path=path))
        return loaded

    def _restore_collection(self, st: dict) -> None:
        col = st.get("collection")
        if col and col.get("paths"):
            win_c = self.show_collection()
            win_c.reference = int(col.get("reference", 0))
            win_c.add_paths([p for p in col["paths"] if Path(p).exists()])

    # ---- lingua -----------------------------------------------------------------------
    def change_language(self, lang: str) -> None:
        """Passa a un'altra lingua: salva la scelta, ricostruisce la finestra con gli stessi pannelli e la stessa
        selezione, e chiude quella vecchia. I dati (analisi, soglie, file) restano nella `Session`."""
        lang = i18n.resolve(lang)
        if lang == i18n.get_language():
            return
        st = self.state()
        i18n.set_language(lang)
        try:
            i18n.save_language(lang)
        except OSError:
            pass  # la scelta vale comunque per questa esecuzione
        self._replace_window(st, tf("Lingua: {name}", name=i18n.LANGS[lang]))

    def change_theme(self, mode: str) -> None:
        """Passa al tema chiaro, scuro o di sistema: salva la scelta e ricostruisce la finestra come per la lingua
        (i grafici prendono i colori quando nascono)."""
        if mode == theme.mode():
            return
        st = self.state()
        theme.set_mode(mode)
        try:
            theme.save_mode(mode)
        except OSError:
            pass
        self._replace_window(st, tf("Tema: {name}", name=tr(theme.MODES[mode])))

    def _replace_window(self, st: dict, message: str) -> None:
        """Una finestra nuova con i pannelli, la disposizione e la selezione di `st` al posto di questa."""
        new = MainWindow(self.session)
        new.autosave, new.session_path = self.autosave, self.session_path
        new._build_docks_from(st)
        if self.session.analysis:
            new.refresh_views()
            new._apply_selection(st)
        new.show()
        new._restore_collection(st)
        new.status_label.setText(message)
        self.autosave = False
        if self.collection is not None:
            self.collection.close()
        self.close()

    def closeEvent(self, event) -> None:  # noqa: N802
        if self.autosave and self.session.source_text:
            try:
                path = last_session_path()
                path.parent.mkdir(parents=True, exist_ok=True)
                session_file.write(path, self.state())  # per «Riapri l'ultima sessione»
            except OSError:
                pass
        if self.collection is not None:
            self.collection.close()
        super().closeEvent(event)

    def show_manual(self, anchor: str | None = None):
        """Il manuale nella lingua dell'interfaccia (una finestra sola, riusata)."""
        from pistak.gui.manual_window import ManualWindow

        win = getattr(self, "_manual", None)
        if win is None or win.lang != i18n.get_language():
            win = self._manual = ManualWindow(self)
        win.show()
        win.raise_()
        win.activateWindow()
        if anchor:
            win.go_to(anchor)
        return win

    def _about(self) -> None:
        QMessageBox.about(self, APP_NAME, tf(
            "<b>{app}</b> — analisi dei π-stack da file CIF.<br>Basato sull'algoritmo di moleculesfriends.eu.<br>"
            "Codice in <code>pistak/</code>; parametri e viste si estendono con i plugin.", app=APP_NAME))

    # ---- stato ---------------------------------------------------------------------
    def _busy(self, busy: bool, label: str) -> None:
        self.progress.setVisible(busy)
        if busy:
            self.status_label.setText(label)
        QApplication.setOverrideCursor(Qt.CursorShape.BusyCursor) if busy else QApplication.restoreOverrideCursor()

    def _analysis_changed(self) -> None:
        a = self.session.analysis
        if a:
            self.setWindowTitle(f"{APP_NAME} — {a.name}")

    def _show_error(self, msg: str) -> None:
        self.status_label.setText(msg)
        QMessageBox.warning(self, APP_NAME, msg)
