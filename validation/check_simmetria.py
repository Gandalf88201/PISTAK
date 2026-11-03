"""Checklist, sezione 1: operazioni di simmetria di PISTAK contro gemmi (oracolo esterno, `pip install gemmi`).

1. Per tutte le impostazioni della tabella di gemmi, `hall_operations(simbolo di Hall)` deve dare lo stesso insieme di
   operazioni (rotazione e traslazione modulo il reticolo) di gemmi.
2. Ogni voce di `HM_TO_HALL` deve dare le operazioni del gruppo che gemmi trova con lo stesso simbolo H-M.
3. `operations()` con un simbolo H-M in forme alternative (completa, senza spazi, con :1/:2/:H/:R). Senza
   suffisso, per i gruppi con due origini, il riferimento è l'origine 2 (la convenzione dichiarata da PISTAK).

    .venv/bin/python validation/check_simmetria.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pistak.core.spacegroups import HM_TO_HALL, SpaceGroupError, hall_operations, operations  # noqa: E402


def gemmi_ops(ops: gemmi.GroupOps) -> set:
    """Operazioni di gemmi come (R, t in dodicesimi modulo 12): gemmi tiene R×24 e t in ventiquattresimi."""
    out = set()
    for op in ops:
        r = tuple(tuple(x // 24 for x in row) for row in op.rot)
        assert all(x % 2 == 0 for x in op.tran), op.triplet()
        t = tuple((x // 2) % 12 for x in op.tran)
        out.add((r, t))
    return out


def ours(hall: str) -> set:
    return {(r, tuple(x % 12 for x in t)) for r, t in hall_operations(hall)}


def parse_triplets(triplets: list[str]) -> set:
    return gemmi_ops(gemmi.GroupOps([gemmi.Op(t) for t in triplets]))


def check_hall() -> None:
    print("== hall_operations contro gemmi, tutte le impostazioni ==")
    n = bad = err = 0
    for sg in gemmi.spacegroup_table():
        n += 1
        try:
            mine = ours(sg.hall)
        except SpaceGroupError as e:
            err += 1
            print(f"  NON LETTO  {sg.number:3d} {sg.xhm():22s} Hall «{sg.hall}»: {e}")
            continue
        ref = gemmi_ops(sg.operations())
        if mine != ref:
            bad += 1
            print(f"  DIVERSO    {sg.number:3d} {sg.xhm():22s} Hall «{sg.hall}»: {len(mine)} operazioni contro {len(ref)}")
    print(f"  {n} impostazioni: {n - bad - err} uguali, {bad} diverse, {err} non lette")


def reference(hm: str):
    """Il gruppo di gemmi per un simbolo H-M, con la convenzione di PISTAK: senza suffisso, origine 2."""
    sg = gemmi.find_spacegroup_by_name(hm)
    if sg is not None and sg.ext == "1" and ":" not in hm:
        sg = gemmi.find_spacegroup_by_name(f"{hm}:2")
    return sg


def check_table() -> None:
    print("\n== HM_TO_HALL contro gemmi (stesso simbolo H-M) ==")
    bad = 0
    for hm, hall in HM_TO_HALL.items():
        sg = reference(hm)
        if sg is None:
            print(f"  gemmi non conosce «{hm}»")
            continue
        if ours(hall) != gemmi_ops(sg.operations()):
            bad += 1
            print(f"  DIVERSO  «{hm}» → «{hall}»; gemmi: {sg.xhm()} = «{sg.hall}»")
    print(f"  {len(HM_TO_HALL)} voci, {bad} diverse")


def check_lookup() -> None:
    """Il percorso usato davvero da un CIF: `operations(hm=...)` con le scritture che si trovano nei file."""
    print("\n== operations(hm=...) con scritture alternative ==")
    cases = ["P 1 21/c 1", "P21/c", "P 1 21/n 1", "P21/n", "P 1 21/a 1", "C 1 2/c 1", "C2/c", "I 1 2/a 1", "P b c a",
             "Pbca", "P n m a", "Pna21", "P 21 21 21", "P212121", "C m c a", "C m c e", "Cmce", "R -3 :H", "R -3",
             "R -3 :R", "R 3 c :H", "P 63/m m c", "F d d d :1", "F d d d :2", "F d d d", "I 41/a :1", "I 41/a :2",
             "I 41/a", "P 4/n", "P 4/n :1", "P -1", "P1", "P 2/c", "P 1 2/c 1", "P c c n", "R 3 :R", "R -3 c :R", "R 3 2 :R", "F d -3 m :1", "P n n n",
             "P n n n :1", "P 42/n", "P 42/n :1"]
    bad = 0
    for hm in cases:
        res = operations(hm=hm)
        sg = reference(hm)
        if res is None or sg is None:
            print(f"  {hm:14s} PISTAK {'nessuna' if res is None else res[1]}, gemmi {'nessuna' if sg is None else sg.xhm()}")
            continue
        triplets, source = res
        same = parse_triplets(triplets) == gemmi_ops(sg.operations())
        bad += not same
        print(f"  {hm:14s} {source:6s} → {'uguale' if same else 'DIVERSO'} a gemmi {sg.xhm()}")
    print(f"  {len(cases)} simboli, {bad} diversi")


if __name__ == "__main__":
    check_hall()
    check_table()
    check_lookup()
