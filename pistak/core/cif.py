"""Lettura dei file CIF (Crystallographic Information File) quanto basta per ricostruire il cristallo:
cella, operazioni di simmetria e posizioni frazionarie degli atomi del primo blocco di dati.
Il resto del file (fattori di spostamento, riflessi, geometria) viene ignorato.

Porting fedele di `src/lib/crystal/cif.ts` del sito.
"""
from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from pistak.core.elements import is_known


@dataclass(frozen=True)
class Cell:
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float


@dataclass(frozen=True)
class SymOp:
    """Operazione di simmetria x' = R·x + t, in coordinate frazionarie."""

    r: np.ndarray  # 3×3
    t: np.ndarray  # 3
    text: str

    def apply(self, frac: np.ndarray) -> np.ndarray:
        return self.r @ frac + self.t


@dataclass(frozen=True)
class CifAtom:
    label: str
    el: str
    frac: np.ndarray
    occ: float = 1.0


@dataclass
class Cif:
    name: str
    space_group: str
    cell: Cell
    symops: list[SymOp]
    atoms: list[CifAtom]
    skipped_disorder: int = 0
    items: dict[str, str] = field(default_factory=dict)  # tutti i valori semplici del blocco, per i plugin


class CifError(Exception):
    """Errore di lettura. `code` è uno tra: empty, noCell, noAtoms, noSymops, badSymop."""

    def __init__(self, code: str, detail: str = ""):
        super().__init__(f"{code}: {detail}" if detail else code)
        self.code = code
        self.detail = detail


@dataclass
class _Loop:
    tags: list[str]
    rows: list[list[str]]


@dataclass
class _Block:
    name: str = ""
    items: dict[str, str] = field(default_factory=dict)
    loops: list[_Loop] = field(default_factory=list)


def tokenize(text: str) -> list[str]:
    """Divide il testo in valori CIF: parole, stringhe tra apici e campi di testo tra righe che iniziano con ';'."""
    out: list[str] = []
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            body = [line[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                body.append(lines[i])
                i += 1
            out.append("\n".join(body).strip())
            i += 1
            continue
        j, n = 0, len(line)
        while j < n:
            ch = line[j]
            if ch in " \t":
                j += 1
                continue
            if ch == "#":
                break
            if ch in "'\"":
                # l'apice chiude solo se seguito da uno spazio o dalla fine della riga ("O'Brien" resta intero)
                k = j + 1
                while k < n and not (line[k] == ch and (k + 1 == n or line[k + 1].isspace())):
                    k += 1
                out.append(line[j + 1:k])
                j = k + 1
                continue
            k = j
            while k < n and not line[k].isspace():
                k += 1
            out.append(line[j:k])
            j = k
        i += 1
    return out


_KEYWORD = re.compile(r"^(loop_|data_|save_|global_|stop_)", re.I)


def _is_tag(t: str) -> bool:
    return t.startswith("_")


def _is_keyword(t: str) -> bool:
    return bool(_KEYWORD.match(t))


def _norm_tag(t: str) -> str:
    """_symmetry_equiv_pos_as_xyz e _space_group_symop.operation_xyz (stile CIF2) diventano confrontabili."""
    return t.lower().replace(".", "_")


def _first_block(tokens: list[str]) -> _Block:
    block = _Block()
    i, n = 0, len(tokens)
    while i < n and not re.match(r"^data_", tokens[i], re.I):
        i += 1
    if i < n:
        block.name = tokens[i][5:]
        i += 1
    while i < n:
        t = tokens[i]
        if re.match(r"^data_", t, re.I):
            break  # solo il primo blocco
        if t.lower() == "loop_":
            i += 1
            tags: list[str] = []
            while i < n and _is_tag(tokens[i]):
                tags.append(_norm_tag(tokens[i]))
                i += 1
            values: list[str] = []
            while i < n and not _is_tag(tokens[i]) and not _is_keyword(tokens[i]):
                values.append(tokens[i])
                i += 1
            if tags:
                k = len(tags)
                rows = [values[s:s + k] for s in range(0, len(values) - k + 1, k)]
                block.loops.append(_Loop(tags, rows))
            continue
        if _is_tag(t):
            v = tokens[i + 1] if i + 1 < n else None
            if v is not None and not _is_tag(v) and not _is_keyword(v):
                block.items[_norm_tag(t)] = v
                i += 2
            else:
                i += 1
            continue
        i += 1
    return block


_UNCERTAINTY = re.compile(r"\(\d+\)$")


def cif_number(v: str | None) -> float:
    """'10.3528(16)' → 10.3528; '?' e '.' → NaN."""
    if v is None or v in ("?", "."):
        return float("nan")
    try:
        return float(_UNCERTAINTY.sub("", v))
    except ValueError:
        return float("nan")


def _loop_with(b: _Block, *tags: str) -> _Loop | None:
    return next((lp for lp in b.loops if any(t in lp.tags for t in tags)), None)


_TERM = re.compile(r"^(-?)(\d*\.?\d*(?:/\d+)?)\*?([xyz]?)$")


def _parse_component(raw: str) -> tuple[np.ndarray, float]:
    """Una componente di un'operazione, '-x+1/2' o '1/2+y' o 'x-y': coefficienti di x, y, z e costante."""
    s = re.sub(r"\s+", "", raw.lower())
    if not s:
        raise CifError("badSymop", raw)
    r = np.zeros(3)
    t = 0.0
    for term in filter(None, s.replace("-", "+-").split("+")):
        m = _TERM.match(term)
        if not m or (not m.group(2) and not m.group(3)):
            raise CifError("badSymop", raw)
        sign = -1.0 if m.group(1) else 1.0
        mag = 1.0
        if m.group(2):
            p, _, q = m.group(2).partition("/")
            try:
                mag = float(p) / float(q) if q else float(p)
            except (ValueError, ZeroDivisionError):
                raise CifError("badSymop", raw) from None
        if m.group(3):
            r["xyz".index(m.group(3))] += sign * mag
        else:
            t += sign * mag
    return r, t


def parse_symop(text: str) -> SymOp:
    parts = text.split(",")
    if len(parts) != 3:
        raise CifError("badSymop", text)
    comps = [_parse_component(p) for p in parts]
    return SymOp(np.array([c[0] for c in comps]), np.array([c[1] for c in comps]), re.sub(r"\s+", "", text))


def element_of(symbol: str | None, label: str, known: Callable[[str], bool] = is_known) -> str:
    """Simbolo dell'elemento da type_symbol ('O2-', 'C') o dall'etichetta ('C12A', 'Cl1')."""
    src = symbol if symbol and symbol not in ("?", ".") else label
    m = re.match(r"^([A-Za-z])([A-Za-z]?)", src)
    if not m:
        return ""
    one = m.group(1).upper()
    two = one + m.group(2).lower()
    if m.group(2) and known(two) and (symbol or m.group(2) == m.group(2).lower()):
        return two
    if one == "D":
        return "H"
    return one


_SYM_TAGS = ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz")


def parse_cif(text: str, known: Callable[[str], bool] = is_known) -> Cif:
    if not text.strip():
        raise CifError("empty")
    b = _first_block(tokenize(text))
    cell = Cell(
        cif_number(b.items.get("_cell_length_a")), cif_number(b.items.get("_cell_length_b")),
        cif_number(b.items.get("_cell_length_c")), cif_number(b.items.get("_cell_angle_alpha")),
        cif_number(b.items.get("_cell_angle_beta")), cif_number(b.items.get("_cell_angle_gamma")),
    )
    if any(not np.isfinite(x) or x <= 0 for x in vars(cell).values()):
        raise CifError("noCell")

    space_group = (b.items.get("_symmetry_space_group_name_h-m") or b.items.get("_space_group_name_h-m_alt") or "").strip()
    sym_loop = _loop_with(b, *_SYM_TAGS)
    if sym_loop:
        col = next(i for i, t in enumerate(sym_loop.tags) if t in _SYM_TAGS)
        symops = [parse_symop(row[col]) for row in sym_loop.rows]
    else:
        single = next((b.items[t] for t in _SYM_TAGS if b.items.get(t)), None)
        if single:
            symops = [parse_symop(single)]
        # senza operazioni il file descrive già tutta la cella solo se il gruppo è P1
        elif not space_group or re.fullmatch(r"p\s*1", space_group, re.I):
            symops = [parse_symop("x,y,z")]
        else:
            raise CifError("noSymops", space_group)

    atom_loop = _loop_with(b, "_atom_site_fract_x")
    if not atom_loop:
        raise CifError("noAtoms")

    def col_of(tag: str) -> int:
        return atom_loop.tags.index(tag) if tag in atom_loop.tags else -1

    i_label, i_type, i_x, i_y, i_z, i_occ, i_dis = (col_of(t) for t in (
        "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
        "_atom_site_fract_z", "_atom_site_occupancy", "_atom_site_disorder_group"))
    atoms: list[CifAtom] = []
    skipped = 0
    for row in atom_loop.rows:
        label = row[i_label] if i_label >= 0 else f"X{len(atoms) + 1}"
        # dei gruppi di disordine si tiene solo il primo (componente principale)
        dis = row[i_dis] if i_dis >= 0 else "."
        if dis not in (".", "?", "1", "-1"):
            skipped += 1
            continue
        frac = np.array([cif_number(row[i_x]), cif_number(row[i_y]), cif_number(row[i_z])])
        if not np.all(np.isfinite(frac)):
            continue
        el = element_of(row[i_type] if i_type >= 0 else None, label, known)
        if not el:
            continue
        occ = cif_number(row[i_occ]) if i_occ >= 0 else 1.0
        atoms.append(CifAtom(label, el, frac, occ if np.isfinite(occ) else 1.0))
    if not atoms:
        raise CifError("noAtoms")
    return Cif(b.name, space_group, cell, symops, atoms, skipped, dict(b.items))
