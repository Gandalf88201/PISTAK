"""Errori del calcolo. Ognuno ha un `code` stabile (non tradotto) e un dettaglio libero: il testo per l'utente,
nella sua lingua, lo costruisce `pistak.errors.describe_error` fuori dal nucleo."""
from __future__ import annotations


class PistakError(Exception):
    """Errore previsto (file illeggibile, struttura non trattabile…). `code` identifica il caso."""

    def __init__(self, code: str, detail: str = ""):
        super().__init__(f"{code}: {detail}" if detail else code)
        self.code = code
        self.detail = detail


class InputError(PistakError):
    """Il file in ingresso non si può leggere. `code`: fileTooLarge."""


class CifError(PistakError):
    """Errore di lettura del CIF. `code` è uno tra: empty, noBlock, noCell, noAtoms, noSymops, badSymop."""


class StructureError(PistakError):
    """La cella non si può ricostruire. `code`: tooLarge."""
