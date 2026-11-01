"""PISTAK: analisi dei π-stack nei cristalli molecolari a partire da file CIF.

Il pacchetto è diviso in strati che dipendono solo da quelli sotto di loro:

    pistak.core    calcolo puro (numpy/scipy), nessuna interfaccia: usabile anche come libreria
    pistak.params  registro dei parametri misurabili su una coppia di molecole (estendibile con plugin)
    pistak.io      esportazione dei risultati ed elaborazione di cartelle di file
    pistak.gui     finestra grafica (PySide6, PyVista, pyqtgraph)
    pistak.cli     comandi da terminale

Esempio d'uso come libreria::

    from pistak import analyse_file
    result = analyse_file("struttura.cif")
    for pair in result.pairs:
        print(pair.kind, pair.interplanar, pair.long_slip, pair.lat_slip, pair.get("beta_a"))
"""
__version__ = "0.1.0"

#: Caricati alla prima richiesta (PEP 562): il nucleo porta con sé scipy, e la finestra grafica deve potersi
#: mostrare prima di averlo importato.
_LAZY = {"analyse_file": "pistak.api", "analyse_text": "pistak.api", "Analysis": "pistak.core.analysis",
         "Limits": "pistak.core.config"}
__all__ = ["Analysis", "Limits", "analyse_file", "analyse_text", "__version__"]


def __getattr__(name: str):
    if name in _LAZY:
        import importlib

        value = getattr(importlib.import_module(_LAZY[name]), name)
        globals()[name] = value
        return value
    raise AttributeError(f"module 'pistak' has no attribute {name!r}")
