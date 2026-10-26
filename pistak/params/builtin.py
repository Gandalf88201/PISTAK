"""I parametri dell'analizzatore del sito, già calcolati da `measure_pair`: qui sono solo registrati."""
from __future__ import annotations

from pistak.params.registry import parameter


def _attr(key: str, label: str, unit: str, attr: str, description: str, digits: int = 2) -> None:
    parameter(key, label, unit, description, digits, group="sito", source="builtin")(lambda ctx, a=attr: getattr(ctx.pair, a))


_attr("interplanar", "Distanza interplanare", "Å", "interplanar", "|d·n|: distanza tra i piani medi misurata lungo la normale.")
_attr("long_slip", "Slittamento longitudinale", "Å", "long_slip", "Componente dello slittamento lungo l'asse lungo della prima molecola.")
_attr("lat_slip", "Slittamento laterale", "Å", "lat_slip", "Componente dello slittamento lungo l'asse corto.")
_attr("slip", "Slittamento totale", "Å", "slip", "Modulo della componente di d nel piano.")
_attr("dihedral", "Angolo tra i piani", "°", "dihedral", "Angolo tra le normali dei due nuclei π (0–90°).", 1)
_attr("centroid_distance", "Distanza tra i centroidi", "Å", "centroid_distance", "|d|: distanza tra i baricentri dei due nuclei π.")
_attr("contact", "Contatto più corto", "Å", "contact", "Distanza minima tra atomi pesanti delle due molecole.")
_attr("overlap", "Sovrapposizione", "", "overlap", "Frazione degli atomi di un nucleo che cade, in proiezione, dentro il contorno dell'altro.")
_attr("count", "Coppie equivalenti", "", "count", "Quante coppie con gli stessi valori, legate dalla simmetria, sono state raccolte.", 0)
