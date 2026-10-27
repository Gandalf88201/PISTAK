"""Finestra delle opzioni per esportare la cella intera verso VASP, Quantum ESPRESSO, CP2K o CRYSTAL
(`pistak.io.periodic`)."""
from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QVBoxLayout,
)

from pistak.gui.screen import fit_to_screen
from pistak.i18n import tr
from pistak.io.periodic import PeriodicOptions

PROGRAM_TEXT = {"vasp": "VASP (POSCAR, INCAR, KPOINTS)", "qe": "Quantum ESPRESSO (pw.x)", "cp2k": "CP2K", "crystal": "CRYSTAL"}
TASK_TEXT = {"scf": "energia (singolo punto)", "relax": "ottimizza gli atomi", "vc-relax": "ottimizza atomi e cella"}
DISP_TEXT = {"d3bj": "D3(BJ)", "d3zero": "D3 (zero)", "d4": "D4", "none": "nessuna"}
#: Dispersioni che ciascun programma accetta da queste opzioni.
SUPPORTED = {"vasp": ("d3bj", "d3zero", "d4", "none"), "qe": ("d3bj", "d3zero", "none"),
             "cp2k": ("d3bj", "d3zero", "d4", "none"), "crystal": ("d3bj", "d3zero", "none")}


class PeriodicDialog(QDialog):
    def __init__(self, options: PeriodicOptions | None = None, parent=None):
        super().__init__(parent)
        o = options or PeriodicOptions()
        self.setWindowTitle(tr("Cella per calcoli periodici"))
        form = QFormLayout()
        self.program = QComboBox()
        for k, text in PROGRAM_TEXT.items():
            self.program.addItem(text, k)
        self.task = QComboBox()
        for k, text in TASK_TEXT.items():
            self.task.addItem(tr(text), k)
        self.functional = QLineEdit(o.functional)
        self.functional.setToolTip(tr("VASP e Quantum ESPRESSO: PBE o PBEsol. CP2K e CRYSTAL accettano anche ibridi (B3LYP, PBE0)."))
        self.dispersion = QComboBox()
        for k, text in DISP_TEXT.items():
            self.dispersion.addItem(tr(text), k)
        self.kpoints = QDoubleSpinBox()
        self.kpoints.setRange(1, 200)
        self.kpoints.setDecimals(0)
        self.kpoints.setSuffix(" Å")
        self.kpoints.setValue(o.kpoint_length)
        self.kpoints.setToolTip(tr("Punti k lungo ciascun asse: n = ⌈L / a⌉. 20 Å va bene per gli isolanti molecolari."))
        self.normalize_h = QCheckBox(tr("Distanze X–H da neutroni (C–H 1,083 Å, N–H 1,009 Å, O–H 0,983 Å)"))
        self.normalize_h.setChecked(o.normalize_h)
        for combo, value in ((self.program, o.program), (self.task, o.task), (self.dispersion, o.dispersion)):
            combo.setCurrentIndex(max(0, combo.findData(value)))
        form.addRow(tr("Programma"), self.program)
        form.addRow(tr("Calcolo"), self.task)
        form.addRow(tr("Funzionale"), self.functional)
        form.addRow(tr("Correzione di dispersione"), self.dispersion)
        form.addRow(tr("Lunghezza per i punti k"), self.kpoints)
        form.addRow("", self.normalize_h)
        note = QLabel(tr("Tutti gli atomi della cella in P1, con la regola di disordine dell'analisi. Pseudopotenziali e basi "
                         "vanno scelti a parte; controlla le parole chiave con la versione del programma che usi."))
        note.setWordWrap(True)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay = QVBoxLayout(self)
        lay.addLayout(form)
        lay.addWidget(note)
        lay.addWidget(buttons)
        self.program.currentIndexChanged.connect(self._sync)
        self._sync()
        fit_to_screen(self, 640, 400, min_width=320, min_height=280)

    def _sync(self) -> None:
        """Le dispersioni che il programma scelto non ha sono disattivate."""
        ok = SUPPORTED[self.program.currentData()]
        for i in range(self.dispersion.count()):
            item = self.dispersion.model().item(i)
            if item is not None:
                item.setEnabled(self.dispersion.itemData(i) in ok)
        if self.dispersion.currentData() not in ok:
            self.dispersion.setCurrentIndex(self.dispersion.findData("d3bj"))

    def options(self) -> PeriodicOptions:
        return PeriodicOptions(program=self.program.currentData(), task=self.task.currentData(),
                               functional=self.functional.text().strip() or "PBE", dispersion=self.dispersion.currentData(),
                               kpoint_length=self.kpoints.value(), normalize_h=self.normalize_h.isChecked())
