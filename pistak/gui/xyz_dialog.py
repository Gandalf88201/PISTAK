"""Scelta dell'esportazione XYZ di tutta la struttura: un file per dimero o un file per coppia della tabella, con gli
avvisi su che cosa contiene la seconda (`pistak.io.export.xyz_counts`)."""
from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QLabel, QRadioButton, QVBoxLayout

from pistak.i18n import tf, tr
from pistak.io.export import XyzCounts
from pistak.labels import xyz_dimers_text, xyz_pairs_warnings

DIMERS, PAIRS = "dimers", "pairs"


def _note(text: str, warn: bool = False) -> QLabel:
    label = QLabel(text)
    label.setTextFormat(Qt.TextFormat.PlainText)  # nomi come «<struttura>» non sono tag HTML
    label.setWordWrap(True)
    label.setContentsMargins(24, 0, 0, 6)  # rientrata sotto il suo pulsante
    if warn:
        label.setStyleSheet("color: #c26a00;")
    return label


class XyzExportDialog(QDialog):
    def __init__(self, name: str, counts: XyzCounts, choice: str = DIMERS, parent=None, stem: str = "pistak"):
        super().__init__(parent)
        self.setWindowTitle(tf("Esporta XYZ — {name}", name=name))
        self.setMinimumWidth(520)
        lay = QVBoxLayout(self)
        self.dimers = QRadioButton(tf("Un file per dimero — {n} file", n=counts.dimers))
        self.pairs = QRadioButton(tf("Un file per coppia della tabella — {n} file", n=counts.pairs))
        lay.addWidget(self.dimers)
        lay.addWidget(_note(xyz_dimers_text(counts)))
        lay.addWidget(self.pairs)
        warnings = xyz_pairs_warnings(counts)
        for w in warnings:
            lay.addWidget(_note("⚠ " + w, warn=True))
        if not warnings:
            lay.addWidget(_note(tr("Ogni coppia è un dimero diverso: i due insiemi di file hanno le stesse molecole.")))
        lay.addWidget(_note(tf("I file per coppia vanno nella cartella {pairs}, quelli per dimero in {dimers}: i due "
                               "insiemi non si mescolano.", pairs=f"{stem}_pairs_xyz", dimers=f"{stem}_xyz")))
        self.dimers.setEnabled(counts.dimers > 0)
        (self.pairs if choice == PAIRS or counts.dimers == 0 else self.dimers).setChecked(True)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel)
        buttons.addButton(tr("Esporta…"), QDialogButtonBox.ButtonRole.AcceptRole)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)

    def choice(self) -> str:
        return PAIRS if self.pairs.isChecked() else DIMERS
