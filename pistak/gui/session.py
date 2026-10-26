"""Stato condiviso dell'applicazione: analisi corrente, soglie e selezione, con i segnali che tengono
allineate le viste.

Sincronizzazione: le viste appartengono a un *gruppo* ("A", "B", "C") oppure a nessuno. La selezione e il
filtro (sottoinsieme di coppie) sono per gruppo: scegliere una coppia in una vista del gruppo A sposta tutte
le altre viste del gruppo A, e solo quelle. Una vista senza gruppo tiene una selezione propria.

Soglie: `set_limits` rifà solo gli stadi dell'analisi che le soglie cambiate toccano (`pistak.core.stages`). Il
criterio di stack si applica subito, nel thread dell'interfaccia (pochi millisecondi: i cursori del pannello
«Soglie in tempo reale»); gli altri stadi in un thread a parte, uno alla volta, e le richieste arrivate nel frattempo
si fondono nell'ultima. Dopo una rianalisi per nuove soglie selezione e filtri restano sulle stesse coppie.
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal

from pistak import api
from pistak.core.analysis import Analysis, read_cif_file
from pistak.core.config import Limits, load_limits, save_user_config
from pistak.core.demo import demo_cif
from pistak.core.stages import first_stage
from pistak.errors import describe_error
from pistak.i18n import tf, tr

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
    uncertaintyChanged = Signal()  # sono state calcolate le incertezze dell'analisi corrente
    limitsChanged = Signal()  # le soglie della sessione sono cambiate (l'analisi le segue subito o poco dopo)
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
        #: Blocco data_ analizzato (None = il primo con atomi).
        self.source_block: str | None = None
        self._selection: dict[str, int] = {g: 0 for g in GROUPS}
        self._subset: dict[str, set[int] | None] = {g: None for g in GROUPS}
        self._pool = QThreadPool.globalInstance()
        self._busy = False
        self._tasks: list[_Task] = []
        #: Una rianalisi per nuove soglie è in corso in un thread; `_limits_pending`: ne è arrivata un'altra intanto.
        self._limits_task = False
        self._limits_pending = False

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

    def load_file(self, path: str | Path, block: str | None = None) -> None:
        p = Path(path)
        try:
            text = read_cif_file(p)
        except Exception as e:  # noqa: BLE001
            self.error.emit(f"{p.name}: {describe_error(e)}")
            return
        self.load_text(text, p.name, str(p.resolve()), block)

    def load_demo(self) -> None:
        self.load_text(demo_cif(), "demo.cif")

    def load_text(self, text: str, name: str, path: str | None = None, block: str | None = None) -> None:
        self.source_text, self.source_name, self.source_path, self.source_block = text, name, path, block
        self._run(name)

    def reanalyse(self) -> None:
        if self.source_text:
            self._run(self.source_name)

    def _run(self, name: str) -> None:
        text, limits, block = self.source_text, self.limits, self.source_block

        def work() -> Analysis:
            return api.analyse_text(text, name, limits, block)

        self._submit(work, tf("Analisi di {name}…", name=name), self._apply)

    def _apply(self, analysis: Analysis, keep: bool = False) -> None:
        """Rende corrente `analysis`. Con `keep` (stessa struttura, altre soglie) selezione e filtri restano sulle
        stesse coppie, riconosciute dai contatti che contengono; altrimenti si riparte dalla prima coppia."""
        old = self.analysis
        where = {m: p.index for p in analysis.pairs for m in p.members} if keep and old is not None else None

        def moved(i: int) -> int:
            return where.get(old.pairs[i].key, -1) if where is not None and old is not None and 0 <= i < len(old.pairs) else -1

        self.analysis = analysis
        first = 0 if analysis.pairs else -1
        for g in GROUPS:
            if where is None:
                self._selection[g], self._subset[g] = first, None
                continue
            i = moved(self._selection[g])
            self._selection[g] = i if i >= 0 else first
            sub = self._subset[g]
            self._subset[g] = None if sub is None else {j for j in map(moved, sub) if j >= 0}
        self.analysisChanged.emit()
        if not keep:
            n = len(analysis.pairs)
            self.message.emit(tf("{name}: {n} coppie di molecole vicine", name=analysis.name, n=n) if n
                              else tf("{name}: nessuna coppia con nuclei π", name=analysis.name))

    def _report(self, analysis: Analysis) -> None:
        kinds = [p.kind for p in analysis.pairs]
        self.message.emit(tf("{name}: {stack} stack, {parallel} parallele, {tilted} inclinate · {cols} colonne",
                             name=analysis.name, stack=kinds.count("stack"), parallel=kinds.count("parallel"),
                             tilted=kinds.count("tilted"), cols=len(analysis.stacks)))

    def _submit(self, fn: Callable[[], object], label: str, on_done: Callable[[object], None],
                on_fail: Callable[[], None] | None = None) -> None:
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
            if on_fail is not None:
                on_fail()
            self.error.emit(msg)

        task.signals.done.connect(finish)
        task.signals.failed.connect(fail)
        self._pool.start(task)

    def _set_busy(self, busy: bool, label: str) -> None:
        self._busy = busy
        self.busyChanged.emit(busy, label)

    # ---- soglie e plugin ---------------------------------------------------------------
    def set_limits(self, limits: Limits, save: bool = False) -> None:
        """Nuove soglie: l'analisi corrente le segue rifacendo solo gli stadi che toccano."""
        self.limits = limits
        if save:
            save_user_config(limits)
        self.limitsChanged.emit()
        self.update_analysis()

    def update_analysis(self) -> None:
        """Porta l'analisi corrente alle soglie della sessione (`self.limits`), dallo stadio che serve."""
        base = self.analysis
        if base is None:
            self.reanalyse()
            return
        stage = first_stage(base.limits, self.limits)
        if stage is None:
            return
        if stage == "pairs":  # solo il criterio di stack: pochi millisecondi, subito
            try:
                new = api.reanalyse(base, self.limits)
            except Exception as e:  # noqa: BLE001
                self.error.emit(describe_error(e))
                return
            self._apply(new, keep=True)
            self._report(new)
            return
        if self._limits_task:
            self._limits_pending = True
            return
        self._limits_task = True
        limits = self.limits

        def done(new: Analysis) -> None:
            self._limits_task = False
            if self.analysis is base:  # nel frattempo non è stato aperto un altro file
                self._apply(new, keep=True)
                self._report(new)
            if self._limits_pending:
                self._limits_pending = False
                self.update_analysis()

        def failed() -> None:
            self._limits_task = self._limits_pending = False

        self._submit(lambda: api.reanalyse(base, limits), tf("Rianalisi di {name} con le nuove soglie…", name=base.name),
                     done, failed)

    def compute_uncertainty(self, samples: int, include_cell: bool = True) -> None:
        """σ dei parametri dell'analisi corrente con il Monte Carlo, in un thread a parte."""
        analysis = self.analysis
        if analysis is None:
            return

        def work():
            return api.uncertainty(analysis, samples, include_cell=include_cell)

        def done(u) -> None:
            if self.analysis is analysis:
                self.uncertaintyChanged.emit()
                self.message.emit(tf("Incertezze calcolate su {n} campioni", n=samples) if u.has_su
                                  else tr("Il CIF non ha incertezze standard (s.u.): tutte le σ sono 0"))

        self._submit(work, tf("Monte Carlo delle incertezze ({n} campioni)…", n=samples), done)

    def reload_plugins(self) -> list[str]:
        loaded = api.ensure_plugins(force=True)
        self.paramsChanged.emit()
        self.reanalyse()
        return loaded
