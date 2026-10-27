"""Testi per l'utente costruiti dai risultati del calcolo, nella lingua corrente.

Il nucleo (`pistak.core`) restituisce solo dati e identificatori stabili; terminale e finestra li trasformano qui
in parole. Gli export su file non usano queste funzioni: le loro colonne hanno chiavi fisse.
"""
from __future__ import annotations

from pistak.core.pairs import Pair
from pistak.core.stacks import RANK_NAMES
from pistak.i18n import tr


def relation(pair: Pair) -> str:
    """Come la seconda molecola è legata alla prima: «stessa [0 1 0]» o «molecola 2 [1 0 0]»; con più frammenti π
    anche quali («· π1–π2»), e per le coppie nella stessa molecola «intramolecolare · π1–π3»."""
    frag = f"π{pair.fa + 1}–π{pair.fb + 1}"
    if pair.intramolecular:
        return f"{tr('intramolecolare')} · {frag}"
    t = "[" + " ".join(map(str, pair.shift)) + "]"
    base = f"{tr('stessa')} {t}" if pair.a == pair.b else f"{tr('molecola')} {pair.b + 1} {t}"
    return f"{base} · {frag}" if pair.fragmented else base


def rank_name(rank: int) -> str:
    """Nome del tipo di colonna di stack (0 = aggregato finito … 3 = rete 3D)."""
    return tr(RANK_NAMES[rank])


#: Testi degli avvisi dei controlli rapidi (`pistak.core.validation`), in italiano: si traducono al momento dell'uso.
ALERT_TEXT = {
    "overlap": "{n} coppie di atomi quasi sovrapposti (sotto 0,7 Å): {atoms}. Disordine non risolto o atomi doppi?",
    "formulaHeavy": "La formula del modello per unità formula ({model}) non coincide con _chemical_formula_sum ({cif}): "
                    "{diff} (modello/CIF). Solvente tolto con SQUEEZE, atomi mancanti o Z sbagliato?",
    "formulaH": "Gli idrogeni del modello non coincidono con la formula dichiarata: {diff} (modello/CIF).",
    "formulaUnchecked": "Formula non controllata: nel CIF mancano _chemical_formula_sum o Z.",
    "density": "Densità calcolata dal modello {calc} g/cm³, dichiarata nel CIF {cif} g/cm³.",
    "noHydrogens": "Nessun idrogeno nel modello: contatti e dimeri per i calcoli non sono realistici.",
    "missingH": "{n} atomi di carbonio terminali senza idrogeni (idrogeni mancanti?): {atoms}.",
    "bondLength": "{n} legami di lunghezza insolita: {bonds}.",
    "shortContact": "{n} contatti intermolecolari {pair} più corti di {limit} Å: {contacts}.",
    "highR1": "R1 alto: {value}.",
    "highWR2": "wR2 alto: {value}.",
    "goof": "GoF lontano da 1: {value}.",
    "residualDensity": "Densità elettronica residua alta per una struttura di atomi leggeri: {value} e/Å³.",
    "noRefinement": "Il CIF non riporta R1, wR2 né GoF.",
    "disorderSkipped": "{n} atomi di componenti disordinate esclusi dall'analisi (regola «{policy}»).",
    "partialOccupancy": "{n} atomi con occupanza parziale nel modello analizzato: {atoms}.",
    "symopsDerived": "Operazioni di simmetria ricavate dal simbolo del gruppo ({source}), non elencate nel CIF.",
    "missedSymmetry": "Simmetria probabilmente mancata: entro {tolerance} Å il modello ha la simmetria di {found}, "
                      "il CIF dichiara {declared}.",
    "pseudoSymmetry": "Pseudosimmetria: entro {tolerance} Å il modello ha la simmetria di {found} (dichiarato {declared}).",
    "addsymUnavailable": "Ricerca della simmetria mancata non eseguita: serve il pacchetto spglib.",
    "checkFailed": "Il controllo «{check}» non è riuscito: {error}",
}
#: Che cosa vuol dire ciascun livello degli avvisi.
ALERT_LEVELS = {"A": "problema serio", "B": "possibile problema", "C": "da controllare", "G": "informazione"}


def alert_text(alert) -> str:
    """Il testo di un avviso (`pistak.core.validation.Alert`) nella lingua corrente."""
    template = ALERT_TEXT.get(alert.code)
    if template is None:
        return f"{alert.code}: {alert.values}"
    try:
        return tr(template).format(**alert.values)
    except (KeyError, IndexError, ValueError):
        return f"{alert.code}: {alert.values}"
