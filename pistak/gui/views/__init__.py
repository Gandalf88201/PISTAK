"""Registro delle viste disponibili. Per aggiungerne una: crea la classe (vedi `base.PairView`, con un attributo
TITLE in italiano da passare a tr()) e inseriscila qui: comparirà da sola nel menu Vista ▸ Nuova vista.

Le chiavi sono identificatori stabili (non tradotti), scritti anche nei file di sessione."""
from pistak.gui.views.base import PairView
from pistak.gui.views.detail_view import DetailView
from pistak.gui.views.hist_view import HistogramView
from pistak.gui.views.projection_view import ProjectionView
from pistak.gui.views.scatter_view import ScatterView
from pistak.gui.views.table_view import TableView
from pistak.gui.views.view3d import View3D
from pistak.i18n import tr

VIEW_TYPES: dict[str, type[PairView]] = {
    "3d": View3D,
    "projection": ProjectionView,
    "scatter": ScatterView,
    "histogram": HistogramView,
    "table": TableView,
    "detail": DetailView,
}

#: Nomi usati dalle sessioni salvate con la prima versione (in italiano) → identificatori attuali.
LEGACY_KINDS = {"3D": "3d", "Proiezione": "projection", "Scatter parametri": "scatter", "Istogramma": "histogram",
                "Tabella": "table", "Dettaglio": "detail"}


def view_title(kind: str) -> str:
    """Il nome della vista nella lingua corrente."""
    return tr(VIEW_TYPES[kind].TITLE)
