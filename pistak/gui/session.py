"""Stato condiviso dell'applicazione: analisi corrente, soglie e selezione, con i segnali che tengono
allineate le viste.

Sincronizzazione: le viste appartengono a un *gruppo* ("A", "B", "C") oppure a nessuno. La selezione e il
filtro (sottoinsieme di coppie) sono per gruppo: scegliere una coppia in una vista del gruppo A sposta tutte
le altre viste del gruppo A, e solo quelle. Una vista senza gruppo tiene una selezione propria.
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal

from pistak import api
from pistak.core.analysis import Analysis, read_cif_file
from pistak.core.cif import parse_cif
from pistak.core.config import Limits, load_limits, save_user_config
from pistak.core.demo import demo_cif
from pistak.errors import describe_error
from pistak.i18n import tf

GROUPS = ("A", "B", "C")


class _Signals(QObject):
    done = Signal(object)
    failed = Signal(str)


class _Task(QRunnable):
    def __init__(self, fn: Callable[[], object]):
        super().__init__()
        self.fn = fn
        self.signals = _Signals()

    def run(self) -> None:
        try:
            res = self.fn()
        except Exception as e:  # noqa: BLE001
            self.signals.failed.emit(describe_error(e))
        else:
            self.signals.done.emit(res)


class Session(QObject):
    analysisChanged = Signal()
    selectionChanged = Signal(str, int)  # gruppo, indice della coppia (-1 = nessuna)
    subsetChanged = Signal(str, object)  # gruppo, insieme di indici o None
    paramsChanged = Signal()  # i parametri del registro sono cambiati (plugin ricaricati)
    busyChanged = Signal(bool, str)
    message = Signal(str)
    error = Signal(str)

    def __init__(self) -> None:
        super().__init__()
        self.analysis: Analysis | None = None
        self.limits: Limits = load_limits()
        self.source_text = ""
        self.source_name = ""
        self.source_path: str | None = None
        self._selection: dict[str, int] = {g: 0 for g in GROUPS}
        self._subset: dict[str, set[int] | None] = {g: None for g in GROUPS}
        self._pool = QThreadPool.globalInstance()
        self._busy = False
        self._tasks: list[_Task] = []

    # ---- selezione e filtro per gruppo -------------------------------------------------
    def selection(self, group: str) -> int:
        return self._selection.get(group, -1)

    def select(self, group: str, index: int) -> None:
        if group in self._selection:
            self._selection[group] = index
            self.selectionChanged.emit(group, index)

    def subset(self, group: str) -> set[int] | None:
        return self._subset.get(group)

    def set_subset(self, group: str, indices: set[int] | None) -> None:
        if group in self._subset:
            self._subset[group] = indices
            self.subsetChanged.emit(group, indices)

    # ---- caricamento e analisi ---------------------------------------------------------
    @property
    def busy(self) -> bool:
        return self._busy

    def load_file(self, path: str | Path) -> None:
        p = Path(path)
        try:
            text = read_cif_file(p)
        except Exception as e:  # noqa: BLE001
            self.error.emit(f"{p.name}: {describe_error(e)}")
            return
        self.load_text(text, p.name, str(p.resolve()))

    def load_demo(self) -> None:
        self.load_text(demo_cif(), "demo.cif")

    def load_text(self, text: str, name: str, path: str | None = None) -> None:
        self.source_text, self.source_name, self.source_path = text, name, path
        self._run(name)

    def reanalyse(self) -> None:
        if self.source_text:
            self._run(self.source_name)

    def _run(self, name: str) -> None:
        text, limits = self.source_text, self.limits

        def work() -> Analysis:
            cif = parse_cif(text)
            from pistak.core.analysis import analyse_cif

            return api.finish(analyse_cif(cif, name, limits))

        self._submit(work, tf("Analisi di {name}…", name=name), self._apply)

    def _apply(self, analysis: Analysis) -> None:
        self.analysis = analysis
        for g in GROUPS:
            self._selection[g] = 0 if analysis.pairs else -1
            self._subset[g] = None
        self.analysisChanged.emit()
        n = len(analysis.pairs)
        self.message.emit(tf("{name}: {n} coppie di molecole vicine", name=analysis.name, n=n) if n
                          else tf("{name}: nessuna coppia con nuclei π", name=analysis.name))

    def _submit(self, fn: Callable[[], object], label: str, on_done: Callable[[object], None]) -> None:
        self._set_busy(True, label)
        task = _Task(fn)
        self._tasks.append(task)

        def finish(res: object) -> None:
            self._tasks.remove(task)
            self._set_busy(False, "")
            on_done(res)

        def fail(msg: str) -> None:
            self._tasks.remove(task)
            self._set_busy(False, "")
            self.error.emit(msg)

        task.signals.done.connect(finish)
        task.signals.failed.connect(fail)
        self._pool.start(task)

    def _set_busy(self, busy: bool, label: str) -> None:
        self._busy = busy
        self.busyChanged.emit(busy, label)

    # ---- soglie e plugin ---------------------------------------------------------------
    def set_limits(self, limits: Limits, save: bool = False) -> None:
        self.limits = limits
        if save:
            save_user_config(limits)
        self.reanalyse()

    def reload_plugins(self) -> list[str]:
        loaded = api.ensure_plugins(force=True)
        self.paramsChanged.emit()
        self.reanalyse()
        return loaded
