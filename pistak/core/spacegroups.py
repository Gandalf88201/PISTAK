"""Operazioni di simmetria dal simbolo del gruppo spaziale, per i CIF che non le elencano.

Si usa, in ordine: il simbolo di Hall (`_space_group_name_Hall`), che descrive i generatori del gruppo e si
interpreta senza tabelle (S. R. Hall, Acta Cryst. A37, 517, 1981; International Tables B, §1.4); poi il simbolo
Hermann–Mauguin, tradotto in Hall con la tabella `HM_TO_HALL` (i gruppi e le impostazioni più comuni nei cristalli
molecolari); infine, se è installato, `gemmi`, che conosce tutte le impostazioni.

Le traslazioni si tengono in dodicesimi interi: tutte quelle dei gruppi spaziali lo sono, e così il gruppo si chiude
senza errori di arrotondamento.
"""
from __future__ import annotations

import re
from fractions import Fraction

import numpy as np

Op = tuple[tuple[tuple[int, int, int], ...], tuple[int, int, int]]  # (R 3×3 interi, t in dodicesimi)

_I = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

# Rotazioni proprie lungo gli assi principali (righe: x', y', z' in funzione di x, y, z).
_ROT = {
    ("1", "z"): _I, ("1", "x"): _I, ("1", "y"): _I,
    ("2", "z"): ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),
    ("3", "z"): ((0, -1, 0), (1, -1, 0), (0, 0, 1)),
    ("4", "z"): ((0, -1, 0), (1, 0, 0), (0, 0, 1)),
    ("6", "z"): ((1, -1, 0), (1, 0, 0), (0, 0, 1)),
    ("2", "x"): ((1, 0, 0), (0, -1, 0), (0, 0, -1)),
    ("3", "x"): ((1, 0, 0), (0, 0, -1), (0, 1, -1)),
    ("4", "x"): ((1, 0, 0), (0, 0, -1), (0, 1, 0)),
    ("6", "x"): ((1, 0, 0), (0, 1, -1), (0, 1, 0)),
    ("2", "y"): ((-1, 0, 0), (0, 1, 0), (0, 0, -1)),
    ("3", "y"): ((-1, 0, 1), (0, 1, 0), (-1, 0, 0)),
    ("4", "y"): ((0, 0, 1), (0, 1, 0), (-1, 0, 0)),
    ("6", "y"): ((0, 0, 1), (0, 1, 0), (-1, 0, 1)),
    ("3", "*"): ((0, 0, 1), (1, 0, 0), (0, 1, 0)),
}
# Assi binari lungo le diagonali di faccia, rispetto all'asse della rotazione precedente: ' = a−b (per z), " = a+b.
_DIAG2 = {
    ("z", "'"): ((0, -1, 0), (-1, 0, 0), (0, 0, -1)), ("z", '"'): ((0, 1, 0), (1, 0, 0), (0, 0, -1)),
    ("x", "'"): ((-1, 0, 0), (0, 0, -1), (0, -1, 0)), ("x", '"'): ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
    ("y", "'"): ((0, 0, -1), (0, -1, 0), (-1, 0, 0)), ("y", '"'): ((0, 0, 1), (0, -1, 0), (1, 0, 0)),
}
_TRANS = {"a": (6, 0, 0), "b": (0, 6, 0), "c": (0, 0, 6), "n": (6, 6, 6), "u": (3, 0, 0), "v": (0, 3, 0),
          "w": (0, 0, 3), "d": (3, 3, 3)}
_CENTRING = {
    "P": [], "A": [(0, 6, 6)], "B": [(6, 0, 6)], "C": [(6, 6, 0)], "I": [(6, 6, 6)],
    "R": [(8, 4, 4), (4, 8, 8)], "S": [(4, 4, 8), (8, 8, 4)], "T": [(4, 8, 4), (8, 4, 8)],
    "F": [(0, 6, 6), (6, 0, 6), (6, 6, 0)],
}
_MATRIX = re.compile(r"^(-?)([12346])([1-5]?)([xyz'\"*]?)([1-5]?)([abcnuvwd]*)$")


class SpaceGroupError(ValueError):
    pass


def _mul(a: Op, b: Op) -> Op:
    ra, ta = np.array(a[0]), np.array(a[1])
    rb, tb = np.array(b[0]), np.array(b[1])
    r = ra @ rb
    t = (ra @ tb + ta) % 12
    return tuple(map(tuple, r.tolist())), tuple(int(x) for x in t)


def _neg(r) -> tuple:
    return tuple(tuple(-x for x in row) for row in r)


def hall_operations(symbol: str) -> list[Op]:
    """Tutte le operazioni (modulo le traslazioni del reticolo) del gruppo descritto dal simbolo di Hall."""
    s = symbol.strip().replace("_", " ")
    shift = (0, 0, 0)
    m = re.search(r"\(([^)]*)\)\s*$", s)
    if m:
        try:
            shift = tuple(int(x) for x in m.group(1).split())
        except ValueError:
            raise SpaceGroupError(symbol) from None
        if len(shift) != 3:
            raise SpaceGroupError(symbol)
        s = s[:m.start()].strip()
    tokens = s.split()
    if not tokens:
        raise SpaceGroupError(symbol)
    lat = tokens[0]
    centro = lat.startswith("-")
    lat = lat.lstrip("-").upper()
    if lat not in _CENTRING:
        raise SpaceGroupError(symbol)
    gens: list[Op] = [(_I, c) for c in _CENTRING[lat]]
    if centro:
        gens.append((_neg(_I), (0, 0, 0)))
    prev_n, prev_axis = None, None
    for k, tok in enumerate(tokens[1:]):
        mm = _MATRIX.match(tok)
        if not mm:
            raise SpaceGroupError(f"{symbol}: {tok}")
        improper, n, screw1, axis, screw2, trans = mm.groups()
        screw = screw1 or screw2
        if not axis:
            if k == 0:
                axis = "z"
            elif k == 1 and n == "2":
                axis = "x" if prev_n in ("2", "4") else "'"
            elif k == 2 and n == "3":
                axis = "*"
            else:
                axis = "z" if n == "1" else None
            if axis is None:
                raise SpaceGroupError(f"{symbol}: {tok}")
        if axis in ("'", '"'):
            if n != "2":
                raise SpaceGroupError(f"{symbol}: {tok}")
            r = _DIAG2[(prev_axis or "z", axis)]
        else:
            if (n, axis) not in _ROT:
                raise SpaceGroupError(f"{symbol}: {tok}")
            r = _ROT[(n, axis)]
        t = np.zeros(3, dtype=int)
        for ch in trans:
            t += _TRANS[ch]
        if screw:
            if axis not in "xyz":
                raise SpaceGroupError(f"{symbol}: {tok}")
            t["xyz".index(axis)] += 12 * int(screw) // int(n)
        if improper:
            r = _neg(r)
        gens.append((r, tuple(int(x) % 12 for x in t)))
        prev_n = n
        if axis in "xyz":
            prev_axis = axis
    # origine spostata di v (in dodicesimi): W' = V W V⁻¹, cioè t' = t + (I − R) v
    if any(shift):
        v = np.array(shift)
        gens = [(r, tuple(int(x) for x in (np.array(t) + (np.eye(3, dtype=int) - np.array(r)) @ v) % 12)) for r, t in gens]
    return _close(gens)


def _close(gens: list[Op]) -> list[Op]:
    ops: list[Op] = [(_I, (0, 0, 0))]
    seen = set(ops)
    frontier = list(ops)
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                c = _mul(a, g)
                if c not in seen:
                    seen.add(c)
                    new.append(c)
        ops += new
        frontier = new
        if len(ops) > 192:
            raise SpaceGroupError("il gruppo non si chiude (più di 192 operazioni)")
    return ops


def triplet(op: Op) -> str:
    """Operazione nella forma dei CIF: '-x+1/2,y+1/2,-z+1/2'."""
    out = []
    for row, t in zip(op[0], op[1]):
        s = ""
        for c, v in zip(row, "xyz"):
            if c:
                s += ("+" if c > 0 and s else "") + ("-" if c < 0 else "") + (str(abs(c)) if abs(c) != 1 else "") + v
        if t:
            f = Fraction(t, 12)
            s += f"+{f.numerator}/{f.denominator}" if f.denominator != 1 else f"+{f.numerator}"
        out.append(s or "0")
    return ",".join(out)


def _hm_key(symbol: str) -> str:
    """Forma confrontabile di un simbolo H-M: 'P 1 21/c 1' → 'p21/c', 'P21/n' → 'p21/n', 'R -3 :H' → 'r-3'."""
    s = re.sub(r":\s*[hHrR12]\s*$", "", symbol.strip())
    s = s.replace("_", "").strip()
    tokens = s.split()
    if len(tokens) == 4 and tokens[0][:1].isalpha() and sum(t == "1" for t in tokens[1:]) == 2:
        tokens = [tokens[0]] + [t for t in tokens[1:] if t != "1"]  # simbolo monoclino completo
    return "".join(tokens).lower()


#: Simboli Hermann–Mauguin (impostazione standard, o quella indicata) → simbolo di Hall.
HM_TO_HALL = {
    "P 1": "P 1", "P -1": "-P 1",
    "P 2": "P 2y", "P 21": "P 2yb", "C 2": "C 2y", "P m": "P -2y", "P c": "P -2yc", "P n": "P -2yac", "P a": "P -2ya",
    "C m": "C -2y", "C c": "C -2yc", "P 2/m": "-P 2y", "P 21/m": "-P 2yb", "C 2/m": "-C 2y", "P 2/c": "-P 2yc",
    "P 21/c": "-P 2ybc", "P 21/n": "-P 2yn", "P 21/a": "-P 2yab", "C 2/c": "-C 2yc", "I 2/a": "-I 2ya",
    "P 2 2 2": "P 2 2", "P 2 2 21": "P 2c 2", "P 21 21 2": "P 2 2ab", "P 21 21 21": "P 2ac 2ab", "C 2 2 21": "C 2c 2",
    "C 2 2 2": "C 2 2", "F 2 2 2": "F 2 2", "I 2 2 2": "I 2 2", "I 21 21 21": "I 2b 2c",
    "P m m 2": "P 2 -2", "P c a 21": "P 2c -2ac", "P n a 21": "P 2c -2n", "P b a 2": "P 2 -2ab", "P n n 2": "P 2 -2n",
    "C m c 21": "C 2c -2", "F d d 2": "F 2 -2d", "I b a 2": "I 2 -2c",
    "P m m m": "-P 2 2", "P b c a": "-P 2ac 2ab", "P n m a": "-P 2ac 2n", "P b c n": "-P 2n 2ab", "P c c n": "-P 2ab 2ac",
    "C m c m": "-C 2c 2", "C m c a": "-C 2bc 2", "C m c e": "-C 2bc 2", "I b a m": "-I 2 2c", "I b c a": "-I 2b 2c",
    "F d d d": "-F 2uv 2vw",
    "P 4": "P 4", "P 41": "P 4w", "P 42": "P 4c", "P 43": "P 4cw", "I 4": "I 4", "I 41": "I 4bw", "P -4": "P -4",
    "I -4": "I -4", "P 4/m": "-P 4", "P 42/m": "-P 4c", "I 4/m": "-I 4", "I 41/a": "-I 4ad",
    "P 41 21 2": "P 4abw 2nw", "P 43 21 2": "P 4nw 2abw",
    "P 3": "P 3", "P 31": "P 31", "P 32": "P 32", "R 3": "R 3", "P -3": "-P 3", "R -3": "-R 3",
    "P 3 2 1": 'P 3 2"', "P 31 2 1": 'P 31 2"', "P 32 2 1": 'P 32 2"', "R 3 2": 'R 3 2"', "R 3 c": 'R 3 -2"c',
    "R -3 c": '-R 3 2"c',
    "P 6": "P 6", "P 61": "P 61", "P 65": "P 65", "P 63": "P 6c", "P 6/m": "-P 6", "P 63/m": "-P 6c",
    "P 61 2 2": "P 61 2 (0 0 -1)", "P 63/m m c": "-P 6c 2c",
    "P 21 3": "P 2ac 2ab 3", "P a -3": "-P 2ac 2ab 3", "I 21 3": "I 2b 2c 3", "P m -3 m": "-P 4 2 3",
    "F m -3 m": "-F 4 2 3", "F d -3 m": "-F 4vw 2vw 3", "I a -3 d": "-I 4bd 2c 3",
}
_HM_INDEX = {_hm_key(k): v for k, v in HM_TO_HALL.items()}


def operations(hall: str | None = None, hm: str | None = None) -> tuple[list[str], str] | None:
    """Le operazioni del gruppo (come triplette CIF) e da dove vengono ("hall", "table", "gemmi"); None se il
    simbolo è sconosciuto."""
    if hall and hall.strip() not in ("?", "."):
        try:
            return [triplet(o) for o in hall_operations(hall)], "hall"
        except SpaceGroupError:
            pass
    if hm and hm.strip() not in ("?", "."):
        h = _HM_INDEX.get(_hm_key(hm))
        if h:
            return [triplet(o) for o in hall_operations(h)], "table"
    try:  # facoltativo: gemmi conosce tutte le impostazioni
        import gemmi  # type: ignore[import-not-found]
    except ImportError:
        return None
    try:
        if hall and hall.strip() not in ("?", "."):
            return [op.triplet() for op in gemmi.symops_from_hall(hall)], "gemmi"
        sg = gemmi.find_spacegroup_by_name(hm or "")
        if sg is not None:
            return [op.triplet() for op in sg.operations()], "gemmi"
    except Exception:  # noqa: BLE001  (simbolo che gemmi non accetta)
        return None
    return None
