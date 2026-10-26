"""Testi degli errori per l'utente, nella lingua corrente. Le classi stanno in `pistak.core.errors` (il nucleo
non dipende dalle traduzioni); qui si trasformano in messaggi per terminale, finestra e riepiloghi."""
from __future__ import annotations

from pistak.core.errors import CifError, InputError, PistakError, StructureError  # noqa: F401  (comodità per chi importa da qui)
from pistak.i18n import tr

ERROR_TEXT = {  # testi italiani: si traducono al momento dell'uso, con tr()
    "empty": "il file è vuoto", "noCell": "mancano i parametri di cella", "noAtoms": "non ci sono atomi con coordinate frazionarie",
    "noSymops": "mancano le operazioni di simmetria (e il gruppo non è P1)", "badSymop": "operazione di simmetria illeggibile",
    "tooLarge": "troppi atomi nella cella", "fileTooLarge": "file troppo grande", "noBlock": "blocco data_ non trovato:",
}


def error_code(e: Exception) -> str:
    """Codice stabile dell'errore: quello di `PistakError`, altrimenti il nome della classe."""
    return e.code if isinstance(e, PistakError) else type(e).__name__


def describe_error(e: Exception) -> str:
    """Messaggio leggibile per un errore qualunque."""
    if isinstance(e, PistakError):
        return f"{tr(ERROR_TEXT.get(e.code, e.code))} {e.detail}".strip()
    return f"{type(e).__name__}: {e}"
