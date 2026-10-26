"""Avvio della finestra grafica."""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("QT_API", "pyside6")  # per pyvistaqt/QtPy


def run(file: str | None = None, screenshot: str | None = None) -> int:
    import pyqtgraph as pg
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication

    from pistak.api import ensure_plugins
    from pistak.gui import theme
    from pistak.gui.main_window import MainWindow
    from pistak.gui.session import Session

    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("PISTAK")
    pg.setConfigOptions(antialias=True, background=theme.background(), foreground=theme.foreground())
    ensure_plugins()
    session = Session()
    win = MainWindow(session)
    win.show()
    if file and file.endswith(".pistak"):
        if not win.load_session(file):
            session.load_demo()
    elif file:
        win.open_file(file)
    else:
        session.load_demo()

    if screenshot:
        # modalità di prova: attende la prima analisi, cattura la finestra e le viste 3D ed esce
        def snap() -> None:
            out = Path(screenshot)
            win.grab().save(str(out))
            for i, dock in enumerate(win.docks):
                if hasattr(dock.widget(), "screenshot"):
                    dock.widget().screenshot(str(out.with_name(f"{out.stem}_3d{i}{out.suffix}")))
            app.quit()

        session.analysisChanged.connect(lambda: QTimer.singleShot(2500, snap))
    return app.exec()
