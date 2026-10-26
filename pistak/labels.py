"""Testi per l'utente costruiti dai risultati del calcolo, nella lingua corrente.

Il nucleo (`pistak.core`) restituisce solo dati e identificatori stabili; terminale e finestra li trasformano qui
in parole. Gli export su file non usano queste funzioni: le loro colonne hanno chiavi fisse.
"""
from __future__ import annotations

from pistak.core.pairs import Pair
from pistak.core.stacks import RANK_NAMES
from pistak.i18n import tr


def relation(pair: Pair) -> str:
    """Come la seconda molecola è legata alla prima: «stessa [0 1 0]» o «molecola 2 [1 0 0]»."""
    t = "[" + " ".join(map(str, pair.shift)) + "]"
    return f"{tr('stessa')} {t}" if pair.a == pair.b else f"{tr('molecola')} {pair.b + 1} {t}"


def rank_name(rank: int) -> str:
    """Nome del tipo di colonna di stack (0 = aggregato finito … 3 = rete 3D)."""
    return tr(RANK_NAMES[rank])
