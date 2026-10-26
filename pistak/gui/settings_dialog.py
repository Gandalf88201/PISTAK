"""Finestra delle soglie: i campi sono generati da `Limits`, quindi una soglia nuova aggiunta lì compare qui."""
from __future__ import annotations

from dataclasses import fields

from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from pistak.core.config import Limits
from pistak.gui.screen import fit_to_screen
from pistak.i18n import tr


class SettingsDialog(QDialog):
    def __init__(self, limits: Limits, parent=None):
        super().__init__(parent)
        self.setWindowTitle(tr("Soglie dell'analisi"))
        self.editors: dict[str, QSpinBox | QDoubleSpinBox] = {}
        form = QFormLayout()
        for f in fields(Limits):
            m = f.metadata
            if f.type == "int":
                w: QSpinBox | QDoubleSpinBox = QSpinBox()
                w.setRange(int(m["min"]), int(m["max"]))
                w.setSingleStep(int(m["step"]))
            else:
                w = QDoubleSpinBox()
                w.setRange(float(m["min"]), float(m["max"]))
                w.setSingleStep(float(m["step"]))
                w.setDecimals(2)
            w.setValue(getattr(limits, f.name))
            self.editors[f.name] = w
            form.addRow(tr(m["label"]), w)
        note = QLabel(tr("Cambiare una soglia rifà l'analisi. «Salva» le scrive in config/user.toml."))
        note.setWordWrap(True)
        reset = QPushButton(tr("Valori del sito"))
        reset.clicked.connect(self._reset)
        self.buttons = QDialogButtonBox()
        self.apply_btn = self.buttons.addButton(tr("Applica"), QDialogButtonBox.ButtonRole.AcceptRole)
        self.save_btn = self.buttons.addButton(tr("Applica e salva"), QDialogButtonBox.ButtonRole.AcceptRole)
        self.buttons.addButton(QDialogButtonBox.StandardButton.Cancel)
        self.apply_btn.clicked.connect(lambda: self._done(False))
        self.save_btn.clicked.connect(lambda: self._done(True))
        self.buttons.rejected.connect(self.reject)
        # il modulo sta in un'area scorrevole: con schermi bassi la finestra non esce dallo schermo
        body = QWidget()
        body.setLayout(form)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setWidget(body)
        lay = QVBoxLayout(self)
        lay.addWidget(scroll, 1)
        lay.addWidget(note)
        lay.addWidget(reset)
        lay.addWidget(self.buttons)
        self.save = False
        fit_to_screen(self, 560, 620, min_width=320, min_height=260)

    def _reset(self) -> None:
        default = Limits()
        for name, w in self.editors.items():
            w.setValue(getattr(default, name))

    def _done(self, save: bool) -> None:
        self.save = save
        self.accept()

    def limits(self) -> Limits:
        return Limits(**{n: (int(w.value()) if isinstance(w, QSpinBox) else float(w.value())) for n, w in self.editors.items()})
