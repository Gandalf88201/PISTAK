"""Finestra delle soglie: i campi sono generati da `Limits`, quindi una soglia nuova aggiunta lì compare qui."""
from __future__ import annotations

from dataclasses import fields

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from pistak.core.config import Limits, user_config_file
from pistak.gui.screen import fit_to_screen
from pistak.i18n import tf, tr


class SettingsDialog(QDialog):
    def __init__(self, limits: Limits, parent=None):
        super().__init__(parent)
        self.setWindowTitle(tr("Soglie dell'analisi"))
        self.editors: dict[str, QSpinBox | QDoubleSpinBox | QComboBox | QLineEdit] = {}
        form = QFormLayout()
        for f in fields(Limits):
            m = f.metadata
            if f.type == "str" and m.get("free"):
                w: QSpinBox | QDoubleSpinBox | QComboBox | QLineEdit = QLineEdit()
            elif f.type == "str":
                w = QComboBox()
                for key, text in m["choices"].items():
                    w.addItem(tr(text), key)
            elif f.type == "int":
                w = QSpinBox()
                w.setRange(int(m["min"]), int(m["max"]))
                w.setSingleStep(int(m["step"]))
            else:
                w = QDoubleSpinBox()
                w.setRange(float(m["min"]), float(m["max"]))
                w.setSingleStep(float(m["step"]))
                w.setDecimals(2)
            self.editors[f.name] = w
            form.addRow(tr(m["label"]), w)
        self._set_values(limits)
        note = QLabel(tf("Cambiare una soglia rifà solo la parte dell'analisi che ne dipende. «Salva» le scrive in {path}.",
                      path=user_config_file()))
        note.setWordWrap(True)
        reset = QPushButton(tr("Valori predefiniti"))
        reset.clicked.connect(lambda: self._set_values(Limits()))
        site = QPushButton(tr("Valori del sito"))
        site.clicked.connect(lambda: self._set_values(Limits.site()))
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
        resets = QHBoxLayout()
        resets.addWidget(reset)
        resets.addWidget(site)
        lay.addLayout(resets)
        lay.addWidget(self.buttons)
        self.save = False
        fit_to_screen(self, 560, 620, min_width=320, min_height=260)

    def _set_values(self, limits: Limits) -> None:
        for name, w in self.editors.items():
            v = getattr(limits, name)
            if isinstance(w, QLineEdit):
                w.setText(v)
            elif isinstance(w, QComboBox):
                if w.findData(v) < 0:  # un valore valido ma senza voce nel menu (es. disorder = "2" da user.toml)
                    w.addItem(v, v)
                w.setCurrentIndex(w.findData(v))
            else:
                w.setValue(v)

    def _done(self, save: bool) -> None:
        self.save = save
        self.accept()

    def limits(self) -> Limits:
        def value(w):
            if isinstance(w, QLineEdit):
                return w.text().strip()
            if isinstance(w, QComboBox):
                return w.currentData()
            return int(w.value()) if isinstance(w, QSpinBox) else float(w.value())

        return Limits(**{n: value(w) for n, w in self.editors.items()})
