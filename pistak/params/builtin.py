"""I parametri dell'analizzatore del sito, già calcolati da `measure_pair`: qui sono solo registrati."""
from __future__ import annotations

from pistak.params.registry import parameter


def _attr(key: str, label: str, unit: str, attr: str, description: str, digits: int = 2) -> None:
    parameter(key, label, unit, description, digits, group="sito", source="builtin")(lambda ctx, a=attr: getattr(ctx.pair, a))


_attr("interplanar", "Distanza interplanare", "Å", "interplanar", "|d·n|: distanza tra i piani medi misurata lungo la normale.")
_attr("long_slip", "Slittamento longitudinale", "Å", "long_slip", "Componente dello slittamento lungo L, la bisettrice degli assi lunghi dei due frammenti proiettati sul piano medio (nel modo del sito: l'asse lungo della prima molecola). Non definito se un asse lungo non lo è o se gli assi sono ruotati più di max_axis_twist.")
_attr("lat_slip", "Slittamento laterale", "Å", "lat_slip", "Componente dello slittamento lungo W = n × L.")
_attr("long_slip_a", "Slittamento longitudinale visto da A", "Å", "long_slip_a", "Componente lungo l'asse lungo di A della posizione del centroide di B proiettata sul piano di A (distanza di «pitch» di Curtis et al., 2004). Non dipende dalla rotazione tra i due assi.")
_attr("lat_slip_a", "Slittamento laterale visto da A", "Å", "lat_slip_a", "Componente lungo l'asse corto di A (n_A × l_A) della posizione del centroide di B proiettata sul piano di A (distanza di «roll» di Curtis et al., 2004).")
_attr("long_slip_b", "Slittamento longitudinale visto da B", "Å", "long_slip_b", "Come «visto da A», nel sistema del frammento B.")
_attr("lat_slip_b", "Slittamento laterale visto da B", "Å", "lat_slip_b", "Come «visto da A», nel sistema del frammento B.")
_attr("slip", "Slittamento totale", "Å", "slip", "Modulo della componente di d nel piano.")
_attr("dihedral", "Angolo tra i piani", "°", "dihedral", "Angolo tra le normali dei due nuclei π (0–90°).", 1)
_attr("centroid_distance", "Distanza tra i centroidi", "Å", "centroid_distance", "|d|: distanza tra i baricentri dei due nuclei π.")
_attr("contact", "Contatto più corto", "Å", "contact", "Distanza minima tra atomi pesanti delle due molecole.")
_attr("overlap", "Sovrapposizione", "", "overlap", "Frazione degli atomi di un nucleo che cade, in proiezione, dentro il contorno dell'altro.")
_attr("count", "Coppie equivalenti", "", "count", "Quante coppie con gli stessi valori, legate dalla simmetria, sono state raccolte.", 0)
_attr("per_molecule", "Contatti per molecola", "", "per_molecule", "Quanti contatti di questa classe ha una molecola (il frammento A): a differenza delle coppie equivalenti per cella non dipende da Z. Una colonna di molecole traslate: 2.", 0)
