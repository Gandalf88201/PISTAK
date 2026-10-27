"""Finestra delle opzioni per esportare il dimero di una coppia verso ORCA, Psi4 o Gaussian (`pistak.io.qc`)."""
from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from pistak.gui.screen import fit_to_screen
from pistak.i18n import tr
from pistak.io.qc import QCOptions

METHOD_TEXT = {"dft-d": "DFT-D con correzione counterpoise", "sapt": "SAPT (solo Psi4)"}
DISP_TEXT = {"d3bj": "D3(BJ)", "d3zero": "D3 (zero)", "d4": "D4", "none": "nessuna"}


def _pair(minimum: int, value: tuple[int, int]) -> tuple[QWidget, QSpinBox, QSpinBox]:
    """Due caselle numeriche affiancate, una per monomero (A e B)."""
    box = QWidget()
    lay = QHBoxLayout(box)
    lay.setContentsMargins(0, 0, 0, 0)
    spins = []
    for name, v in zip("AB", value):
        s = QSpinBox()
        s.setRange(minimum, 9)
        s.setValue(v)
        lay.addWidget(QLabel(name))
        lay.addWidget(s)
        spins.append(s)
    return box, spins[0], spins[1]


class QCDialog(QDialog):
    def __init__(self, options: QCOptions | None = None, parent=None):
        super().__init__(parent)
        o = options or QCOptions()
        self.setWindowTitle(tr("Dimero per calcoli quantistici"))
        form = QFormLayout()
        self.program = QComboBox()
        for p in ("orca", "psi4", "gaussian"):
            self.program.addItem({"orca": "ORCA", "psi4": "Psi4", "gaussian": "Gaussian"}[p], p)
        self.method = QComboBox()
        for k, text in METHOD_TEXT.items():
            self.method.addItem(tr(text), k)
        self.functional = QLineEdit(o.functional)
        self.dispersion = QComboBox()
        for k, text in DISP_TEXT.items():
            self.dispersion.addItem(tr(text), k)
        self.basis = QLineEdit(o.basis)
        self.basis.setPlaceholderText(tr("def2-TZVP (DFT-D) · jun-cc-pVDZ (SAPT)"))
        self.sapt_level = QComboBox()
        for lev in ("sapt0", "ssapt0", "sapt2", "sapt2+", "sapt2+(3)", "sapt2+3"):
            self.sapt_level.addItem(lev, lev)
        self.counterpoise = QCheckBox(tr("Correzione counterpoise (BSSE)"))
        self.counterpoise.setChecked(o.counterpoise)
        self.normalize_h = QCheckBox(tr("Distanze X–H da neutroni (C–H 1,083 Å, N–H 1,009 Å, O–H 0,983 Å)"))
        self.normalize_h.setChecked(o.normalize_h)
        charges, self.qa, self.qb = _pair(-9, o.charges)
        mults, self.ma, self.mb = _pair(1, o.multiplicities)
        self.nprocs = QSpinBox()
        self.nprocs.setRange(1, 256)
        self.nprocs.setValue(o.nprocs)
        self.memory = QSpinBox()
        self.memory.setRange(1, 2048)
        self.memory.setSuffix(" GB")
        self.memory.setValue(o.memory_gb)
        for combo, value in ((self.program, o.program), (self.method, o.method), (self.dispersion, o.dispersion),
                             (self.sapt_level, o.sapt_level)):
            combo.setCurrentIndex(max(0, combo.findData(value)))
        form.addRow(tr("Programma"), self.program)
        form.addRow(tr("Metodo"), self.method)
        form.addRow(tr("Funzionale"), self.functional)
        form.addRow(tr("Correzione di dispersione"), self.dispersion)
        form.addRow(tr("Base"), self.basis)
        form.addRow(tr("Livello SAPT"), self.sapt_level)
        form.addRow("", self.counterpoise)
        form.addRow("", self.normalize_h)
        form.addRow(tr("Cariche"), charges)
        form.addRow(tr("Molteplicità"), mults)
        form.addRow(tr("Processori"), self.nprocs)
        form.addRow(tr("Memoria"), self.memory)
        note = QLabel(tr("I monomeri sono le due molecole intere della coppia, nella geometria del cristallo. "
                         "Controlla le parole chiave con la versione del programma che usi."))
        note.setWordWrap(True)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addWidget(note)
        lay.addWidget(buttons)
        self.method.currentIndexChanged.connect(self._sync)
        self.program.currentIndexChanged.connect(self._sync)
        self._sync()
        fit_to_screen(self, 520, 560, min_width=320, min_height=300)

    def _sync(self) -> None:
        """SAPT esiste solo in Psi4; con SAPT funzionale, dispersione e counterpoise non servono."""
        sapt = self.method.currentData() == "sapt"
        if sapt and self.program.currentData() != "psi4":
            self.program.setCurrentIndex(self.program.findData("psi4"))
        for w in (self.functional, self.dispersion, self.counterpoise):
            w.setEnabled(not sapt)
        self.sapt_level.setEnabled(sapt)
        d4 = self.dispersion.findData("d4")
        item = self.dispersion.model().item(d4)
        if item is not None:  # Gaussian non ha D4
            item.setEnabled(self.program.currentData() != "gaussian")
        if self.program.currentData() == "gaussian" and self.dispersion.currentData() == "d4":
            self.dispersion.setCurrentIndex(self.dispersion.findData("d3bj"))

    def options(self) -> QCOptions:
        return QCOptions(
            program=self.program.currentData(), method=self.method.currentData(), functional=self.functional.text().strip() or "B3LYP",
            dispersion=self.dispersion.currentData(), basis=self.basis.text().strip(), sapt_level=self.sapt_level.currentData(),
            counterpoise=self.counterpoise.isChecked(), normalize_h=self.normalize_h.isChecked(),
            charges=(self.qa.value(), self.qb.value()), multiplicities=(self.ma.value(), self.mb.value()),
            nprocs=self.nprocs.value(), memory_gb=self.memory.value())
