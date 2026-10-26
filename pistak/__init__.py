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
from pistak.api import analyse_file, analyse_text  # noqa: F401
from pistak.core.analysis import Analysis  # noqa: F401
from pistak.core.config import Limits  # noqa: F401

__version__ = "0.1.0"
