"""Pannello «Soglie in tempo reale»: cursori per il criterio di stack, con le coppie riclassificate mentre si muovono.

Le soglie del criterio di stack toccano solo lo stadio `pairs` dell'analisi (`pistak.core.stages`): le coppie già
misurate si riclassificano in pochi millisecondi, senza rifare la geometria, e tabella, scatter e colonne seguono
il cursore. Il contatto massimo invece decide quali coppie misurare: si applica al rilascio del cursore, in un
thread a parte. I campi vengono da `Limits` (etichette, intervalli, passi), come nella finestra «Soglie».
"""
from __future__ import annotations

from dataclasses import fields

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from pistak.core.config import Limits, load_limits, save_user_config, user_config_file
from pistak.core.stages import stage_of
from pistak.gui.session import Session
from pistak.i18n import tf, tr

#: Le soglie con il cursore, nell'ordine del pannello.
LIVE_FIELDS = ("max_dihedral", "max_interplanar", "min_overlap", "max_slip", "contact_cutoff")
#: Passo del cursore quando serve più fine di quello della finestra «Soglie».
FINE_STEP = {"min_overlap": 0.01}
#: Attesa dopo l'ultimo movimento prima di rianalizzare (ms): fonde gli scatti ravvicinati del cursore.
DEBOUNCE_MS = 30

_META = {f.name: f.metadata for f in fields(Limits)}


class _Row:
    """Un cursore e la sua casella numerica, sempre allineati."""

    def __init__(self, name: str, on_change):
        m = _META[name]
        self.name = name
        self.lo, self.hi = float(m["min"]), float(m["max"])
        self.step = FINE_STEP.get(name, float(m["step"]))
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(0, round((self.hi - self.lo) / self.step))
        self.slider.setMinimumWidth(120)
        self.spin = QDoubleSpinBox()
        self.spin.setRange(self.lo, self.hi)
        self.spin.setSingleStep(self.step)
        self.spin.setDecimals(0 if self.step >= 1 else 2)
        self.spin.setKeyboardTracking(False)
        self.spin.setFixedWidth(self.spin.fontMetrics().horizontalAdvance("000.00") + 36)  # cursori allineati
        # le soglie che rimisurano le coppie si applicano al rilascio, le altre mentre il cursore si muove
        self.heavy = stage_of(name) != "pairs"
        self.slider.setTracking(not self.heavy)
        self.slider.sliderMoved.connect(self._moving)
        self.slider.valueChanged.connect(self._from_slider)
        self.spin.valueChanged.connect(self._from_spin)
        self._on_change = on_change
        self.widget = QWidget()
        lay = QHBoxLayout(self.widget)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.slider, 1)
        lay.addWidget(self.spin)

    def value(self) -> float:
        return float(self.spin.value())

    def set_value(self, v: float) -> None:
        """Dall'esterno (soglie cambiate altrove): nessuna notifica."""
        for w in (self.slider, self.spin):
            w.blockSignals(True)
        self.spin.setValue(v)
        self.slider.setValue(round((v - self.lo) / self.step))
        for w in (self.slider, self.spin):
            w.blockSignals(False)

    def _moving(self, i: int) -> None:
        if self.heavy:  # senza tracking: si mostra il valore, si applica al rilascio
            self.spin.blockSignals(True)
            self.spin.setValue(self.lo + i * self.step)
            self.spin.blockSignals(False)

    def _from_slider(self, i: int) -> None:
        self.spin.blockSignals(True)
        self.spin.setValue(min(self.hi, self.lo + i * self.step))
        self.spin.blockSignals(False)
        self._on_change()

    def _from_spin(self, v: float) -> None:
        self.slider.blockSignals(True)
        self.slider.setValue(round((v - self.lo) / self.step))
        self.slider.blockSignals(False)
        self._on_change()


class LimitsPanel(QWidget):
    def __init__(self, session: Session, parent: QWidget | None = None):
        super().__init__(parent)
        self.session = session
        self.rows = {name: _Row(name, self._schedule) for name in LIVE_FIELDS}
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        for name, row in self.rows.items():
            label = QLabel(tr(_META[name]["label"]))
            label.setWordWrap(True)
            if row.heavy:
                tip = tr("Rimisura le coppie: si applica quando lasci il cursore.")
                label.setToolTip(tip)
                row.slider.setToolTip(tip)
            form.addRow(label, row.widget)
        self.measure = QComboBox()
        for key, text in _META["overlap_measure"]["choices"].items():
            self.measure.addItem(tr(text), key)
        self.measure.currentIndexChanged.connect(self._schedule)
        form.addRow(tr(_META["overlap_measure"]["label"]), self.measure)

        self.summary = QLabel("")
        self.summary.setWordWrap(True)
        note = QLabel(tr("Le coppie già misurate si riclassificano mentre muovi i cursori: la geometria non si rifà. "
                         "Le incertezze del Monte Carlo vanno ricalcolate."))
        note.setWordWrap(True)
        note.setEnabled(False)
        save = QPushButton(tr("Salva come predefinite"))
        save.setToolTip(tf("Scrive le soglie in {path}", path=user_config_file()))
        save.clicked.connect(self._save)
        restore = QPushButton(tr("Soglie salvate"))
        restore.setToolTip(tr("Torna alle soglie di user.toml (o a quelle predefinite)"))
        restore.clicked.connect(lambda: self.session.set_limits(load_limits()))
        site = QPushButton(tr("Valori del sito"))
        site.clicked.connect(self._site)
        buttons = QHBoxLayout()
        for b in (save, restore, site):
            buttons.addWidget(b)

        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addWidget(self.summary)
        lay.addWidget(note)
        lay.addLayout(buttons)
        lay.addStretch(1)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(DEBOUNCE_MS)
        self._timer.timeout.connect(self._commit)
        session.limitsChanged.connect(self._sync)
        session.analysisChanged.connect(self._sync)
        self._sync()

    # ---- dalla sessione ai cursori -----------------------------------------------------
    def _sync(self) -> None:
        if self._timer.isActive():  # un movimento non ancora applicato: i cursori hanno ragione loro
            self._update_summary()
            return
        lim = self.session.limits
        for name, row in self.rows.items():
            if row.slider.isSliderDown():
                continue
            v = float(getattr(lim, name))
            if abs(row.value() - v) > 1e-9 or row.slider.value() != round((v - row.lo) / row.step):
                row.set_value(v)
        self.measure.blockSignals(True)
        self.measure.setCurrentIndex(max(0, self.measure.findData(lim.overlap_measure)))
        self.measure.blockSignals(False)
        self._update_summary()

    def _update_summary(self) -> None:
        a = self.session.analysis
        if a is None:
            self.summary.setText("")
            return
        kinds = [p.kind for p in a.pairs]
        self.summary.setText(tf("<b>{stack}</b> stack · {parallel} parallele · {tilted} inclinate · {cols} colonne",
                                stack=kinds.count("stack"), parallel=kinds.count("parallel"),
                                tilted=kinds.count("tilted"), cols=len(a.stacks)))

    # ---- dai cursori alla sessione -----------------------------------------------------
    def values(self) -> dict:
        out: dict = {name: row.value() for name, row in self.rows.items()}
        out["overlap_measure"] = self.measure.currentData()
        return out

    def _schedule(self, *_) -> None:
        self._timer.start()

    def _commit(self) -> None:
        limits = self.session.limits.with_values(**self.values())
        if limits != self.session.limits:
            self.session.set_limits(limits)

    def _save(self) -> None:
        self._timer.stop()
        self._commit()
        path = save_user_config(self.session.limits)
        self.session.message.emit(tf("Soglie salvate in {path}", path=path))

    def _site(self) -> None:
        """Le soglie del criterio di stack del sito, lasciando com'erano le altre (modo delle coppie, frammenti…)."""
        site = Limits.site()
        self.session.set_limits(self.session.limits.with_values(
            **{name: getattr(site, name) for name in (*LIVE_FIELDS, "overlap_measure")}))
