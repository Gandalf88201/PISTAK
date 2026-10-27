"""Base comune delle viste: collegamento alla sessione e sincronizzazione per gruppi.

Per scrivere una vista nuova basta ereditare da `PairView`, riempire `self.body` con i widget e
implementare gli hook (`on_analysis`, `on_pair`, `on_subset`, `on_params`). La sincronizzazione è già
gestita: la vista chiama `self.choose(indice)` quando l'utente sceglie una coppia e riceve `on_pair`
quando cambia la selezione del suo gruppo. Poi si registra in `views/__init__.py` (VIEW_TYPES).

Una nuova analisi arriva con `refresh()`: on_analysis, on_subset, on_pair in quest'ordine. Durante `refresh()`
`self.refreshing` è vero, così una vista che ridisegna tutto in on_subset/on_pair può saltare on_analysis. Una
vista nascosta (dietro un'altra scheda) rimanda il lavoro a quando torna visibile: con i cursori delle soglie in
tempo reale le analisi si susseguono e ridisegnare ciò che non si vede rallenterebbe tutto.
"""
from __future__ import annotations

import html

from PySide6.QtCore import Qt
from PySide6.QtGui import QShowEvent
from PySide6.QtWidgets import QComboBox, QLabel, QSizePolicy, QVBoxLayout, QWidget

from pistak.core.analysis import Analysis
from pistak.core.diagnosis import why_no_pairs
from pistak.core.pairs import Pair
from pistak.gui.flow_layout import FlowLayout
from pistak.gui.session import GROUPS, Session
from pistak.i18n import tf, tr
from pistak.labels import NO_PAIRS_HINT, no_pairs_text


def sync_choices() -> list[tuple[str, str | None]]:
    return [(tf("Gruppo {g}", g=g), g) for g in GROUPS] + [(tr("Non sincronizzata"), None)]


class PairView(QWidget):
    TITLE = "Vista"
    #: Senza coppie mostra il motivo al posto del contenuto (il Dettaglio no: ha altro da mostrare).
    EMPTY_NOTICE = True

    def __init__(self, session: Session, group: str | None = "A", parent: QWidget | None = None):
        super().__init__(parent)
        self.session = session
        self.group: str | None = group
        self._local = -1  # selezione propria, usata solo se group è None
        self._local_subset: set[int] | None = None
        #: Vero durante `refresh()` (on_pair e on_subset stanno per arrivare).
        self.refreshing = False
        #: Qualcosa è cambiato mentre la vista era nascosta: si rifà tutto quando torna visibile.
        self._stale = False

        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 2, 4, 4)
        outer.setSpacing(3)
        bar = FlowLayout()  # a capo automatico: il pannello può restringersi senza forzare la larghezza della finestra
        self.toolbar = bar
        self.sync_combo = QComboBox()
        for label, g in sync_choices():
            self.sync_combo.addItem(label, g)
        self.sync_combo.setCurrentIndex(max(0, self.sync_combo.findData(group)))
        self.sync_combo.setToolTip(tr("Le viste dello stesso gruppo seguono la stessa coppia selezionata.\n"
                                   "'Non sincronizzata': questa vista tiene una selezione propria."))
        self.sync_combo.currentIndexChanged.connect(self._sync_changed)
        bar.addWidget(QLabel("🔗"))
        bar.addWidget(self.sync_combo)
        outer.addLayout(bar)
        # al posto del contenuto quando l'analisi non ha coppie: il perché, invece di una vista vuota
        self.empty_notice = QLabel()
        self.empty_notice.setWordWrap(True)
        self.empty_notice.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.empty_notice.setTextFormat(Qt.TextFormat.RichText)
        self.empty_notice.setMargin(16)
        self.empty_notice.hide()
        outer.addWidget(self.empty_notice, 1)
        self.body = QVBoxLayout()
        self.body.setContentsMargins(0, 0, 0, 0)
        outer.addLayout(self.body, 1)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMinimumSize(160, 110)  # i pannelli si possono rimpicciolire molto: la finestra deve stare in schermi piccoli

        session.analysisChanged.connect(self._analysis_changed)
        session.selectionChanged.connect(self._selection_changed)
        session.subsetChanged.connect(self._subset_changed)
        session.paramsChanged.connect(self.on_params)
        session.uncertaintyChanged.connect(lambda: self.on_pair(self.index))

    # ---- stato ---------------------------------------------------------------------
    @property
    def analysis(self) -> Analysis | None:
        return self.session.analysis

    @property
    def index(self) -> int:
        return self.session.selection(self.group) if self.group else self._local

    @property
    def pair(self) -> Pair | None:
        a = self.analysis
        i = self.index
        return a.pairs[i] if a and 0 <= i < len(a.pairs) else None

    @property
    def subset(self) -> set[int] | None:
        return self.session.subset(self.group) if self.group else self._local_subset

    # ---- azioni dell'utente --------------------------------------------------------
    def choose(self, index: int) -> None:
        """Da chiamare quando l'utente sceglie una coppia in questa vista."""
        if self.group:
            self.session.select(self.group, index)
        elif index != self._local:
            self._local = index
            self.on_pair(index)

    def filter_to(self, indices: set[int] | None) -> None:
        """Da chiamare quando l'utente restringe l'insieme di coppie (None = tutte)."""
        if self.group:
            self.session.set_subset(self.group, indices)
        else:
            self._local_subset = indices
            self.on_subset(indices)

    # ---- eventi della sessione -----------------------------------------------------
    def _sync_changed(self) -> None:
        self.group = self.sync_combo.currentData()
        if self.group is None:  # parte dalla selezione che aveva
            self._local = self.session.selection(next(iter(GROUPS)))
        self.on_pair(self.index)
        self.on_subset(self.subset)

    def refresh(self) -> None:
        """Tutto da capo: dati, filtro e coppia selezionata."""
        if self._hidden():
            self._stale = True
            return
        self._stale = False
        self._update_empty_notice()
        self.refreshing = True
        try:
            self.on_analysis()
            self.on_subset(self.subset)
        finally:
            self.refreshing = False
        self.on_pair(self.index)

    def _update_empty_notice(self) -> None:
        """Senza coppie il contenuto della vista lascia il posto al motivo (`EMPTY_NOTICE = False` per non farlo)."""
        reason = why_no_pairs(self.analysis) if self.analysis is not None and self.EMPTY_NOTICE else None
        if reason is not None:
            self.empty_notice.setText(f"<p><b>{tr('Nessuna coppia di nuclei π da mostrare')}</b></p>"
                                      f"<p>{html.escape(no_pairs_text(reason))}</p><p><i>{tr(NO_PAIRS_HINT)}</i></p>")
        self.empty_notice.setVisible(reason is not None)
        for i in range(self.body.count()):
            w = self.body.itemAt(i).widget()
            if w is not None:
                w.setVisible(reason is None)

    def _hidden(self) -> bool:
        """Nascosta in una finestra visibile (scheda non in primo piano, pannello chiuso). Una vista mai mostrata
        (prove, finestre in costruzione) non conta come nascosta."""
        w = self.window()
        return w is not self and w.isVisible() and not self.isVisible()

    def showEvent(self, e: QShowEvent) -> None:  # noqa: N802
        super().showEvent(e)
        if self._stale:
            self.refresh()

    def _analysis_changed(self) -> None:
        self._local = 0 if self.analysis and self.analysis.pairs else -1
        self._local_subset = None
        self.refresh()

    def _selection_changed(self, group: str, index: int) -> None:
        if group == self.group:
            if self._hidden():
                self._stale = True
            else:
                self.on_pair(index)

    def _subset_changed(self, group: str, indices) -> None:
        if group == self.group:
            if self._hidden():
                self._stale = True
            else:
                self.on_subset(indices)

    # ---- salvataggio della sessione ------------------------------------------------
    def state(self) -> dict:
        """Impostazioni della vista da salvare con la sessione (tutto JSON)."""
        return {"group": self.group, "local": self._local, "view": self.extra_state()}

    def restore(self, st: dict) -> None:
        """Ripristina `state()`. La coppia selezionata del gruppo viene riapplicata dalla sessione, non qui."""
        self.group = st.get("group", self.group)
        self._local = int(st.get("local", -1))
        self.sync_combo.blockSignals(True)
        self.sync_combo.setCurrentIndex(max(0, self.sync_combo.findData(self.group)))
        self.sync_combo.blockSignals(False)
        self.apply_state(st.get("view", {}))

    def extra_state(self) -> dict:
        """Le sottoclassi aggiungono qui le proprie impostazioni (assi scelti, opzioni…)."""
        return {}

    def apply_state(self, st: dict) -> None:
        """Le sottoclassi rileggono qui ciò che hanno scritto in `extra_state`; ignorare le chiavi sconosciute."""

    # ---- hook per le sottoclassi ---------------------------------------------------
    def on_analysis(self) -> None:
        """È arrivata una nuova analisi: aggiorna liste e dati."""

    def on_pair(self, index: int) -> None:
        """La coppia selezionata (per il gruppo di questa vista) è cambiata; -1 = nessuna."""

    def on_subset(self, indices: set[int] | None) -> None:
        """L'insieme di coppie in evidenza è cambiato (None = tutte)."""

    def on_params(self) -> None:
        """L'elenco dei parametri del registro è cambiato (plugin ricaricati)."""
