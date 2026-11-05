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
from pistak.core.errors import CifError


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
    #: Incertezze standard (s.u.) delle coordinate frazionarie, dalle cifre tra parentesi (0 se mancano).
    su: np.ndarray = field(default_factory=lambda: np.zeros(3))
    #: Gruppo di disordine come nel CIF ('.' se l'atomo non è disordinato). Un gruppo negativo (convenzione SHELXL
    #: PART −n) è disordinato su un elemento di simmetria: le sue immagini per simmetria sono l'altra orientazione.
    group: str = "."
    #: Assemblaggio di disordine ('.' se manca): con il gruppo identifica la parte disordinata.
    assembly: str = "."


@dataclass
class Cif:
    name: str
    space_group: str
    cell: Cell
    symops: list[SymOp]
    atoms: list[CifAtom]
    skipped_disorder: int = 0
    items: dict[str, str] = field(default_factory=dict)  # tutti i valori semplici del blocco, per i plugin
    #: Da dove vengono le operazioni: cif (elencate nel file) | hall | table | gemmi (dal simbolo) | p1 (nessuna).
    symops_source: str = "cif"
    #: Regola usata per gli atomi disordinati (vedi `parse_cif`).
    disorder: str = "first"
    #: s.u. dei parametri di cella (a, b, c in Å; α, β, γ in gradi), 0 se mancano.
    cell_su: tuple[float, ...] = (0.0,) * 6
    #: Gli atomi delle componenti disordinate scartate dalla regola `disorder` (servono al contenuto della cella).
    omitted: list[CifAtom] = field(default_factory=list)


@dataclass
class _Loop:
    tags: list[str]
    rows: list[list[str]]


@dataclass
class _Block:
    name: str = ""
    items: dict[str, str] = field(default_factory=dict)
    loops: list[_Loop] = field(default_factory=list)


class _Quoted(str):
    """Valore scritto tra apici o in un campo di testo ';': è sempre un valore, mai un tag o una parola chiave
    (un FCF incollato in `_iucr_refine_fcf_details` comincia con 'data_I' e non deve aprire un blocco)."""


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
            out.append(_Quoted("\n".join(body).strip()))
            i += 1
            continue
        j, n = 0, len(line)
        while j < n:
            ch = line[j]
            if ch.isspace():  # non solo spazio e tab: con NBSP o form feed la parola sotto sarebbe vuota e j fermo
                j += 1
                continue
            if ch == "#":
                break
            if ch in "'\"":
                # l'apice chiude solo se seguito da uno spazio o dalla fine della riga ("O'Brien" resta intero)
                k = j + 1
                while k < n and not (line[k] == ch and (k + 1 == n or line[k + 1].isspace())):
                    k += 1
                out.append(_Quoted(line[j + 1:k]))
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
    return not isinstance(t, _Quoted) and t.startswith("_")


def _is_keyword(t: str) -> bool:
    return not isinstance(t, _Quoted) and bool(_KEYWORD.match(t))


def _is_data(t: str) -> bool:
    return not isinstance(t, _Quoted) and bool(re.match(r"^data_", t, re.I))


def _norm_tag(t: str) -> str:
    """_symmetry_equiv_pos_as_xyz e _space_group_symop.operation_xyz (stile CIF2) diventano confrontabili."""
    return t.lower().replace(".", "_")


def _blocks(tokens: list[str]) -> list[_Block]:
    """Tutti i blocchi `data_` del file, nell'ordine (i valori prima del primo `data_` sono ignorati)."""
    blocks: list[_Block] = []
    i, n = 0, len(tokens)
    while i < n and not _is_data(tokens[i]):
        i += 1
    block: _Block | None = None
    while i < n:
        t = tokens[i]
        if _is_data(t):
            block = _Block(t[5:])
            blocks.append(block)
            i += 1
            continue
        assert block is not None
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
    return blocks


def _first_block(tokens: list[str]) -> _Block:
    blocks = _blocks(tokens)
    return blocks[0] if blocks else _Block()


_UNCERTAINTY = re.compile(r"\(\d+\)$")


def cif_su(v: str | None) -> float:
    """Incertezza standard scritta tra parentesi, nelle unità dell'ultima cifra: '10.3528(16)' → 0.0016,
    '1.5(12)' → 1.2, '123(4)' → 4; 0 se manca (o il valore non è un numero)."""
    if not v:
        return 0.0
    m = re.fullmatch(r"\s*[-+]?(\d*)(?:\.(\d*))?(?:[eE][-+]?\d+)?\((\d+)\)\s*", v)
    if not m:
        return 0.0
    decimals = len(m.group(2) or "")
    exp = re.search(r"[eE]([-+]?\d+)", v)
    return int(m.group(3)) * 10.0 ** (-decimals + (int(exp.group(1)) if exp else 0))


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


#: Senza type_symbol, un'etichetta tutta maiuscola si legge a una lettera ('CA1' è un C-alfa, 'CD1' un C-delta, non
#: Ca o Cd), tranne questi alogeni, che come etichette di carbonio o boro non si usano.
_UPPER_TWO_LETTER = frozenset({"Cl", "Br"})
_CENTROID = re.compile(r"^Cg\d", re.I)


def element_of(symbol: str | None, label: str, known: Callable[[str], bool] = is_known) -> str:
    """Simbolo dell'elemento da type_symbol ('O2-', 'C') o dall'etichetta ('C12A', 'Cl1'). Vuoto per i siti che non
    sono atomi: centroidi ('Cg1', convenzione di PLATON) ed elementi sconosciuti (picchi 'Q1', 'X')."""
    has_symbol = bool(symbol and symbol not in ("?", "."))
    src = symbol if has_symbol else label
    if not has_symbol and _CENTROID.match(label):
        return ""
    m = re.match(r"^([A-Za-z])([A-Za-z]?)", src)
    if not m:
        return ""
    one = m.group(1).upper()
    two = one + m.group(2).lower()
    if m.group(2) and known(two) and (has_symbol or m.group(2).islower() or two in _UPPER_TWO_LETTER):
        return two
    if one == "D":
        return "H"
    return one if known(one) else ""


_SYM_TAGS = ("_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz")
_HALL_TAGS = ("_symmetry_space_group_name_hall", "_space_group_name_hall")
_ATOM_TAG = "_atom_site_fract_x"


def _has_atoms(b: _Block) -> bool:
    return _loop_with(b, _ATOM_TAG) is not None


def block_names(text: str) -> list[str]:
    """I nomi dei blocchi `data_` che contengono atomi (quelli analizzabili), nell'ordine del file."""
    return [b.name for b in _blocks(tokenize(text)) if _has_atoms(b)]


def _symops(b: _Block, space_group: str) -> tuple[list[SymOp], str]:
    sym_loop = _loop_with(b, *_SYM_TAGS)
    if sym_loop:
        col = next(i for i, t in enumerate(sym_loop.tags) if t in _SYM_TAGS)
        return [parse_symop(row[col]) for row in sym_loop.rows], "cif"
    single = next((b.items[t] for t in _SYM_TAGS if b.items.get(t)), None)
    if single:
        return [parse_symop(single)], "cif"
    hall = next((b.items[t] for t in _HALL_TAGS if b.items.get(t)), None)
    from pistak.core.spacegroups import operations

    found = operations(hall=hall, hm=space_group)
    if found:
        ops, source = found
        return [parse_symop(o) for o in ops], source
    # senza operazioni il file descrive già tutta la cella solo se il gruppo è P1
    if not space_group or re.fullmatch(r"p\s*1", space_group, re.I):
        return [parse_symop("x,y,z")], "p1"
    raise CifError("noSymops", space_group)


#: Distanza sotto la quale due atomi di gruppi di disordine diversi non possono esistere insieme (Å): nessun legame è
#: così corto tra atomi pesanti, e due alternative dello stesso sito distano tipicamente 0,3–1 Å.
ALTERNATIVE_DISTANCE = 1.0


def _overlap(a: list[np.ndarray], b: list[np.ndarray], vectors: np.ndarray) -> bool:
    """Vero se un atomo di `a` dista meno di `ALTERNATIVE_DISTANCE` da uno di `b` (immagine minima)."""
    d = np.array(a)[:, None, :] - np.array(b)[None, :, :]
    d -= np.round(d)
    return bool((np.linalg.norm(d @ vectors, axis=2) < ALTERNATIVE_DISTANCE).any())


def _merged_assemblies(groups: list[str], assemblies: list[str], fracs: list[np.ndarray], cell: Cell) -> list[str]:
    """Gli assemblaggi di disordine, uniti quando i gruppi di due assemblaggi diversi si sovrappongono.

    Gruppi che si sovrappongono non possono esistere insieme: sono alternative dello stesso sito anche se il CIF li
    mette in assemblaggi diversi (COD 7208500: acetonitrile con le due orientazioni in A e in B, a 0,31 Å). Uniti gli
    assemblaggi, la regola `major` ne tiene una sola."""
    from pistak.core.structure import cell_vectors  # import qui: structure dipende da cif

    vectors = cell_vectors(cell)
    pos: dict[tuple[str, str], list[np.ndarray]] = {}
    for g, a, f in zip(groups, assemblies, fracs):
        if g not in (".", "?") and np.all(np.isfinite(f)):
            pos.setdefault((a, g), []).append(f)
    parent = {a: a for a, _g in pos}

    def root(a: str) -> str:
        while parent[a] != a:
            a = parent[a]
        return a

    keys = sorted(pos)
    for i, (a1, g1) in enumerate(keys):
        for a2, g2 in keys[i + 1:]:
            if a1 != a2 and root(a1) != root(a2) and _overlap(pos[a1, g1], pos[a2, g2], vectors):
                r1, r2 = sorted((root(a1), root(a2)))
                parent[r2] = r1
    return [root(a) if a in parent else a for a in assemblies]


def _signed_alternatives(groups: list[str], assemblies: list[str], fracs: list[np.ndarray], cell: Cell) -> set[tuple[str, str]]:
    """(assemblaggio, gruppo) dei gruppi negativi −n che nello stesso assemblaggio sono l'alternativa del gruppo n.

    Di norma −n è lo stesso gruppo n su un elemento di simmetria (SHELXL), e si tratta come n. Alcuni CIF però usano
    n e −n come due alternative distinte (COD 4113834: 1 con occupanza 0,61, −1 con 0,39, atomi a 0,4 Å): lo si
    riconosce perché i loro atomi si sovrappongono, e allora −n è un gruppo a sé."""
    from pistak.core.structure import cell_vectors  # import qui: structure dipende da cif

    vectors = cell_vectors(cell)
    out = set()
    pos: dict[tuple[str, str], list[np.ndarray]] = {}
    for g, a, f in zip(groups, assemblies, fracs):
        if g not in (".", "?") and np.all(np.isfinite(f)):
            pos.setdefault((a, g), []).append(f)
    for (a, g), fs in pos.items():
        other = pos.get((a, g[1:])) if g.startswith("-") else None
        if other and _overlap(fs, other, vectors):
            out.add((a, g))
    return out


def _disorder_keep(groups: list[str], assemblies: list[str], occs: list[float], policy: str,
                   split: set[tuple[str, str]] = frozenset()) -> list[bool]:
    """Quali atomi tenere secondo la regola di disordine. `split`: gruppi negativi da non unire al gruppo positivo
    (vedi `_signed_alternatives`).

    first   gruppi '.', '?', '1', '-1' (il comportamento del sito)
    major   per ogni assemblaggio, il gruppo con l'occupanza media più alta (a parità, il numero più basso)
    all     tutti
    <n>     il gruppo n (e −n), più gli atomi non disordinati
    """
    none = [g in (".", "?") for g in groups]
    if policy == "all":
        return [True] * len(groups)
    if policy == "first":
        return [g in (".", "?", "1", "-1") for g in groups]
    keys = [g if (a, g) in split else g.lstrip("-") for g, a in zip(groups, assemblies)]
    if policy.isdigit():
        return [no or k == policy for no, k in zip(none, keys)]
    # major
    stats: dict[tuple[str, str], list[float]] = {}
    for no, k, a, o in zip(none, keys, assemblies, occs):
        if not no:
            stats.setdefault((a, k), []).append(o)

    def order(item):
        (a, g), vals = item
        num = int(g) if g.lstrip("-").isdigit() else 10**6
        return (-sum(vals) / len(vals), abs(num), g)

    best: dict[str, str] = {}
    for (a, g), _vals in sorted(stats.items(), key=order):
        best.setdefault(a, g)
    return [no or best.get(a) == k for no, k, a in zip(none, keys, assemblies)]


def _cif_from_block(b: _Block, known: Callable[[str], bool], disorder: str) -> Cif:
    cell = Cell(
        cif_number(b.items.get("_cell_length_a")), cif_number(b.items.get("_cell_length_b")),
        cif_number(b.items.get("_cell_length_c")), cif_number(b.items.get("_cell_angle_alpha")),
        cif_number(b.items.get("_cell_angle_beta")), cif_number(b.items.get("_cell_angle_gamma")),
    )
    if any(not np.isfinite(x) or x <= 0 for x in (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)):
        raise CifError("noCell")
    cell_su = tuple(cif_su(b.items.get(t)) for t in ("_cell_length_a", "_cell_length_b", "_cell_length_c",
                                                         "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"))

    space_group = (b.items.get("_symmetry_space_group_name_h-m") or b.items.get("_space_group_name_h-m_alt") or "").strip()
    symops, source = _symops(b, space_group)

    atom_loop = _loop_with(b, _ATOM_TAG)
    if not atom_loop:
        raise CifError("noAtoms")

    def col_of(tag: str) -> int:
        return atom_loop.tags.index(tag) if tag in atom_loop.tags else -1

    i_label, i_type, i_x, i_y, i_z, i_occ, i_dis, i_asm, i_flag = (col_of(t) for t in (
        "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
        "_atom_site_fract_z", "_atom_site_occupancy", "_atom_site_disorder_group", "_atom_site_disorder_assembly",
        "_atom_site_calc_flag"))
    rows = atom_loop.rows
    if i_flag >= 0:  # siti fittizi (dummy): non sono atomi
        rows = [r for r in rows if r[i_flag].lower() != "dum"]
    occs = [cif_number(r[i_occ]) if i_occ >= 0 else 1.0 for r in rows]
    occs = [o if np.isfinite(o) else 1.0 for o in occs]
    groups = [r[i_dis] if i_dis >= 0 else "." for r in rows]
    assemblies = [r[i_asm] if i_asm >= 0 else "." for r in rows]
    fracs = [np.array([cif_number(r[i_x]), cif_number(r[i_y]), cif_number(r[i_z])]) for r in rows]
    merged = _merged_assemblies(groups, assemblies, fracs, cell)
    keep = _disorder_keep(groups, merged, occs, disorder, _signed_alternatives(groups, merged, fracs, cell))
    atoms: list[CifAtom] = []
    omitted: list[CifAtom] = []
    skipped = 0
    for row, ok, occ, group, assembly in zip(rows, keep, occs, groups, assemblies):
        label = row[i_label] if i_label >= 0 else f"X{len(atoms) + 1}"
        if not ok:
            skipped += 1
        frac = np.array([cif_number(row[i_x]), cif_number(row[i_y]), cif_number(row[i_z])])
        su = np.array([cif_su(row[i_x]), cif_su(row[i_y]), cif_su(row[i_z])])
        if not np.all(np.isfinite(frac)):
            continue
        el = element_of(row[i_type] if i_type >= 0 else None, label, known)
        if not el:
            continue
        (atoms if ok else omitted).append(CifAtom(label, el, frac, occ, su, group, assembly))
    if not atoms:
        raise CifError("noAtoms")
    return Cif(b.name, space_group, cell, symops, atoms, skipped, dict(b.items), source, disorder, cell_su, omitted)


def parse_cif(text: str, known: Callable[[str], bool] = is_known, block: str | int | None = None,
              disorder: str = "first") -> Cif:
    """Un blocco del CIF: per nome, per posizione (da 0, tra tutti i blocchi) o, con `block=None`, il primo che
    contiene atomi. `disorder` è la regola per gli atomi disordinati (vedi `_disorder_keep`)."""
    if not text.strip():
        raise CifError("empty")
    blocks = _blocks(tokenize(text))
    if isinstance(block, int):
        if not 0 <= block < len(blocks):
            raise CifError("noBlock", str(block))
        b = blocks[block]
    elif block is not None:
        b = next((x for x in blocks if x.name.lower() == block.lower()), None)
        if b is None:
            raise CifError("noBlock", block)
    else:
        # il primo blocco con atomi; se non ce n'è nessuno, il primo (che darà l'errore giusto)
        b = next((x for x in blocks if _has_atoms(x)), blocks[0] if blocks else _Block())
    return _cif_from_block(b, known, disorder)
