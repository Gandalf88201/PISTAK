"""Dal contenuto di un CIF alle molecole del cristallo: si applicano le operazioni di simmetria,
si eliminano i doppioni sulla cella, si riconoscono i legami (anche attraverso i bordi della cella)
e si ricompongono le molecole intere in coordinate cartesiane (Å).

Porting di `src/lib/crystal/structure.ts`; i cicli annidati del sito sono qui vettorizzati con numpy/scipy.
Convenzione: una matrice di cella ha i vettori a, b, c come RIGHE, quindi cart = frac @ M.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from pistak.core.cif import ALTERNATIVE_DISTANCE, Cell, Cif, SymOp
from pistak.core.config import Limits
from pistak.core.elements import NONMETALS, covalent_radius
from pistak.core.errors import StructureError


def cell_vectors(c: Cell) -> np.ndarray:
    """Vettori di cella in coordinate cartesiane, con a lungo x e b nel piano xy (convenzione abituale)."""
    ca, cb, cg = np.cos(np.radians([c.alpha, c.beta, c.gamma]))
    sg = np.sin(np.radians(c.gamma))
    cx = c.c * cb
    cy = c.c * (ca - cb * cg) / sg
    cz = np.sqrt(max(0.0, c.c ** 2 - cx ** 2 - cy ** 2))
    return np.array([[c.a, 0.0, 0.0], [c.b * cg, c.b * sg, 0.0], [cx, cy, cz]])


def to_cart(m: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Coordinate frazionarie → cartesiane (accetta un punto o un array N×3)."""
    return np.asarray(f) @ m


def to_frac(m: np.ndarray, p: np.ndarray) -> np.ndarray:
    return np.asarray(p) @ np.linalg.inv(m)


@dataclass
class Molecule:
    el: list[str]
    labels: list[str]
    pos: np.ndarray  # N×3, Å
    bonds: np.ndarray  # M×2, indici in `el`
    centroid: np.ndarray
    #: Vero se la rete di legami non si chiude: polimero o solido esteso, non una molecola.
    extended: bool
    #: Per ogni atomo: indice dell'atomo nella cella (`Structure.cell_el`) e cella di appartenenza.
    cell_index: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    cell_shift: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), dtype=int))
    #: Rete estesa: base (righe intere, forma a scala) del reticolo delle sue periodicità, le traslazioni t per cui un
    #: atomo dell'unità è legato, attraverso la rete, alla sua copia traslata di t. Una catena lungo a ha la sola riga
    #: (1, 0, 0); una molecola discreta nessuna riga. Due copie dell'unità traslate di t sono la stessa rete se e solo
    #: se t sta in questo reticolo (`same_network`), qualunque sia il taglio dell'unità ripetuta (revisione 5, N1).
    periods: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), dtype=int))
    #: Rete estesa: per ogni atomo dell'unità, i suoi vicini legati nella rete infinita come (atomo dell'unità,
    #: traslazione in celle della copia dell'unità in cui sta). Serve a riconoscere un atomo che lega due frammenti
    #: anche quando sta in un'altra copia dell'unità ripetuta (`pairs.fragments_link`, revisione 5, N1). Vuoto per le
    #: molecole discrete, i cui legami sono tutti in `bonds`.
    links: list = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.el)

    @property
    def heavy(self) -> np.ndarray:
        """Maschera booleana degli atomi che non sono idrogeno."""
        return np.array([e != "H" for e in self.el])

    @property
    def extent(self) -> float:
        """Distanza massima di un atomo dal baricentro: serve a scartare in fretta le molecole lontane."""
        return float(np.linalg.norm(self.pos - self.centroid, axis=1).max(initial=0.0))

    def formula(self) -> list[tuple[str, int]]:
        """Formula in ordine di Hill: C, H, poi gli altri in ordine alfabetico."""
        n: dict[str, int] = {}
        for e in self.el:
            n[e] = n.get(e, 0) + 1
        order = sorted(n, key=lambda e: (e != "C", e != "H", e))
        return [(e, n[e]) for e in order]


#: Dove finisce una molecola sotto un'operazione: (indice della molecola immagine, traslazione in celle).
Image = tuple[int, tuple[int, int, int]]


@dataclass
class Structure:
    vectors: np.ndarray
    cell_el: list[str]
    cell_labels: list[str]
    cell_frac: np.ndarray
    molecules: list[Molecule]
    #: Le operazioni di simmetria del CIF (servono a riconoscere le coppie equivalenti).
    symops: list[SymOp] = field(default_factory=list)
    #: Per ogni atomo della cella: (indice dell'atomo in `Cif.atoms`, indice dell'operazione che lo genera).
    cell_origin: np.ndarray = field(default_factory=lambda: np.zeros((0, 2), dtype=int))
    _images: dict[float, list[list[Image | None]]] = field(default_factory=dict, repr=False)
    _atom_maps: dict[float, list[list[np.ndarray | None]]] = field(default_factory=dict, repr=False)

    @property
    def n_cell_atoms(self) -> int:
        return len(self.cell_el)

    def shift_vector(self, shift) -> np.ndarray:
        return to_cart(self.vectors, np.asarray(shift, dtype=float))

    def molecule_images(self, tol: float = 0.3) -> list[list[Image | None]]:
        """Per ogni operazione di simmetria e ogni molecola della cella, la molecola in cui va a finire:
        g(m nella cella 0) = m' traslata di T celle. None se l'immagine non coincide con una molecola (simmetria
        non esatta, atomi disordinati tolti in modo non simmetrico…): quell'operazione non si usa per m."""
        if tol not in self._images:
            self._images[tol] = _molecule_images(self, tol)
        return self._images[tol]

    def atom_maps(self, tol: float = 0.3) -> list[list[np.ndarray | None]]:
        """Come `molecule_images`, atomo per atomo: per l'operazione k e la molecola m, l'array che dà per ogni atomo
        di m l'indice dell'atomo corrispondente nella molecola immagine (None se l'immagine non è esatta)."""
        if tol not in self._atom_maps:
            imgs = self.molecule_images(tol)
            out = []
            for k, row in enumerate(imgs):
                out.append([_atom_map(self, self.symops[k], m, img, tol) if img is not None else None
                            for m, img in zip(self.molecules, row)])
            self._atom_maps[tol] = out
        return self._atom_maps[tol]

    def rotation(self, k: int) -> np.ndarray:
        """Parte rotazionale (intera) dell'operazione k, in coordinate frazionarie."""
        return np.rint(self.symops[k].r).astype(int)


def _wrap(x: np.ndarray) -> np.ndarray:
    return x - np.floor(x + 1e-9)


def special_disorder_copies(cif: Cif, vectors: np.ndarray, limits: Limits) -> dict[tuple[str, str], set[int]]:
    """Per ogni gruppo di disordine negativo (disordinato su un elemento di simmetria), le operazioni di simmetria le
    cui copie del gruppo si tengono.

    Il gruppo nel CIF ha occupanza parziale e la sua immagine per l'elemento di simmetria è l'altra orientazione: in un
    sito del cristallo ce n'è una sola. Si scorrono le operazioni in ordine e si scarta intera ogni copia del gruppo
    con un atomo a meno di `ALTERNATIVE_DISTANCE` (ma oltre `same_site`, che è lo stesso atomo) da una copia già
    tenuta: così un sito ha un'orientazione sola e le orientazioni non si mescolano."""
    groups: dict[tuple[str, str], list[int]] = {}
    for ia, at in enumerate(cif.atoms):
        if at.group.startswith("-"):
            groups.setdefault((at.assembly, at.group), []).append(ia)
    out: dict[tuple[str, str], set[int]] = {}
    for key, idx in groups.items():
        frac = np.array([cif.atoms[i].frac for i in idx])
        kept: list[np.ndarray] = []
        out[key] = set()
        for k, op in enumerate(cif.symops):
            f = _wrap(frac @ op.r.T + op.t)
            if kept:
                d = f[:, None, :] - np.vstack(kept)[None, :, :]
                d -= np.round(d)
                dist = np.linalg.norm(d @ vectors, axis=2)
                if ((dist >= limits.same_site) & (dist < ALTERNATIVE_DISTANCE)).any():
                    continue
            kept.append(f)
            out[key].add(k)
    return out


DUPLICATE = 0.05  # Å: due atomi del CIF più vicini di così sono lo stesso atomo (nessuna posizione sdoppiata raffinabile)


def fill_cell(cif: Cif, vectors: np.ndarray, limits: Limits,
              origins: list[tuple[int, int]] | None = None) -> tuple[list[str], list[str], np.ndarray]:
    """Tutti gli atomi della cella: ogni operazione di simmetria applicata a ogni atomo, senza doppioni. Se si passa
    `origins`, vi si aggiunge per ogni atomo della cella (indice in `cif.atoms`, indice dell'operazione). Dei gruppi
    di disordine su un elemento di simmetria si tiene un'orientazione per sito (`special_disorder_copies`)."""
    special = special_disorder_copies(cif, vectors, limits)
    els: list[str] = []
    labels: list[str] = []
    buf = np.zeros((16, 3))  # posizioni frazionarie accettate finora (cresce al bisogno)
    by_el: dict[str, list[int]] = {}
    src: list[int] = []  # per ogni posizione accettata, l'atomo del CIF da cui viene

    def ordered(i: int) -> bool:  # atomo a occupanza piena e fuori da gruppi di disordine
        a = cif.atoms[i]
        return a.occ >= 0.99 and a.group in (".", "?", "", "0")

    for ia, at in enumerate(cif.atoms):
        allowed = special.get((at.assembly, at.group))
        for k, op in enumerate(cif.symops):
            if allowed is not None and k not in allowed:
                continue  # l'altra orientazione di un gruppo disordinato su un elemento di simmetria
            f = _wrap(op.apply(at.frac))
            # si fondono le immagini per simmetria dello stesso sito (atomo su una posizione speciale), gli atomi
            # ordinati ripetuti nel CIF (molecola intera su un centro di inversione) e le posizioni che coincidono entro
            # DUPLICATE (lo stesso atomo scritto due volte, anche a occupanza parziale: COD 2100449 elenca le immagini
            # per inversione con occupanza ½); mai due posizioni alternative distinte di un atomo sdoppiato, che
            # distano più di DUPLICATE: revisione 2, voci 60 e 65
            idx = by_el.get(at.el)
            if idx:
                d = f - buf[idx]
                d -= np.round(d)
                dist = np.linalg.norm(d @ vectors, axis=1)
                close = np.flatnonzero(dist < limits.same_site)
                if any(src[idx[c]] == ia or dist[c] < DUPLICATE or (ordered(ia) and ordered(src[idx[c]]))
                       for c in close):
                    continue
            if len(els) == len(buf):
                buf = np.vstack([buf, np.zeros_like(buf)])
            buf[len(els)] = f
            by_el.setdefault(at.el, []).append(len(els))
            src.append(ia)
            els.append(at.el)
            labels.append(at.label)
            if origins is not None:
                origins.append((ia, k))
            if len(els) > limits.max_cell_atoms:
                raise StructureError("tooLarge", str(limits.max_cell_atoms))
    return els, labels, buf[:len(els)].copy()


_SHIFTS = np.array([(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)], dtype=int)


def cell_bonds(els: list[str], pos: np.ndarray, vectors: np.ndarray, limits: Limits) -> list[list[tuple[int, tuple[int, int, int]]]]:
    """Legami tra atomi della cella, anche verso le celle vicine: per ogni atomo, lista di (vicino, traslazione in celle)."""
    n = len(els)
    r = np.array([covalent_radius(e) for e in els])
    is_h = np.array([e == "H" for e in els])
    max_cut = 2 * r.max() + limits.bond_tolerance
    tree = cKDTree(pos)
    rows: list[np.ndarray] = []
    for s in _SHIFTS:
        v = s @ vectors
        other = cKDTree(pos + v)
        # coppie (i, j) con |pos[j] + v − pos[i]| < max_cut
        sm = tree.sparse_distance_matrix(other, max_cut, output_type="ndarray")
        if not len(sm):
            continue
        i, j, d = sm["i"], sm["j"], sm["v"]
        keep = (d < r[i] + r[j] + limits.bond_tolerance) & (d * d >= 0.16)
        keep &= ~(is_h[i] & is_h[j])  # due idrogeni non sono mai legati tra loro in un cristallo molecolare
        if not s.any():
            keep &= i != j
        i, j = i[keep], j[keep]
        rows.append(np.column_stack([i, j, np.broadcast_to(s, (len(i), 3))]))
    edges: list[list[tuple[int, tuple[int, int, int]]]] = [[] for _ in range(n)]
    if rows:
        allrows = np.concatenate(rows)
        order = np.lexsort((allrows[:, 4], allrows[:, 3], allrows[:, 2], allrows[:, 1], allrows[:, 0]))
        for i, j, sx, sy, sz in allrows[order]:
            edges[i].append((int(j), (int(sx), int(sy), int(sz))))
    return edges


#: Gli anelli fino a questa taglia (di soli non metalli) restano interi nell'unità ripetuta di una rete estesa.
UNIT_RING_SIZE = 10

Shift = tuple[int, int, int]


def lattice_basis(vectors) -> np.ndarray:
    """Base intera in forma a scala (righe con il primo elemento non nullo positivo e colonne pivot crescenti) del
    reticolo generato dai vettori interi `vectors`: eliminazione di Gauss sugli interi (algoritmo di Euclide per
    colonna), quindi il reticolo è esattamente quello generato, non il sottospazio razionale."""
    rows = [[int(x) for x in t] for t in vectors if any(int(x) for x in t)]
    basis: list[list[int]] = []
    for col in range(3):
        while True:
            nz = [r for r in rows if r[col] != 0]
            if len(nz) <= 1:
                break
            piv = min(nz, key=lambda r: abs(r[col]))
            out = [piv]
            for r in rows:
                if r is piv:
                    continue
                if r[col] != 0:
                    q = r[col] // piv[col]
                    r = [r[k] - q * piv[k] for k in range(3)]
                if any(r):
                    out.append(r)
            rows = out
        nz = [r for r in rows if r[col] != 0]
        if nz:
            piv = nz[0]
            rows = [r for r in rows if r is not piv]
            basis.append(piv if piv[col] > 0 else [-x for x in piv])
    return np.array(basis, dtype=int).reshape(-1, 3)


def in_lattice(basis: np.ndarray, t) -> bool:
    """Vero se il vettore intero `t` sta nel reticolo della base a scala `basis` (`lattice_basis`)."""
    v = [int(x) for x in t]
    for row in basis:
        col = next(k for k in range(3) if row[k] != 0)
        if v[col] % int(row[col]):
            return False
        q = v[col] // int(row[col])
        v = [v[k] - q * int(row[k]) for k in range(3)]
    return not any(v)


def same_network(s: Structure, a: int, b: int, shift) -> bool:
    """Vero se la molecola `a` (cella 0) e la molecola `b` traslata di `shift` sono lo stesso insieme di atomi legati:
    la stessa molecola discreta, o la stessa rete estesa (catena, strato, rete 3D). Non dipende da come è stata
    tagliata l'unità ripetuta di una rete (revisione 5, N1): due unità traslate di una catena lungo a sono la stessa
    catena se la traslazione è un multiplo di a, due catene parallele altrimenti."""
    if a != b:
        return False
    if not any(int(x) for x in shift):
        return True
    return bool(len(s.molecules[a].periods)) and in_lattice(s.molecules[a].periods, shift)


def network_links(cell_index, cell_shift, edges) -> list[list[tuple[int, Shift]]]:
    """I vicini legati di ogni atomo di un'unità ripetuta nella rete infinita: (posizione nell'unità del vicino,
    traslazione della copia dell'unità che lo contiene), dai legami della cella `edges` (`cell_bonds`)."""
    where = {int(c): k for k, c in enumerate(cell_index)}
    out: list[list[tuple[int, Shift]]] = [[] for _ in range(len(cell_index))]
    for k, c in enumerate(cell_index):
        for j, es in edges[int(c)]:
            kk = where.get(int(j))
            if kk is None:
                continue
            out[k].append((kk, tuple(int(cell_shift[k][d]) + int(es[d]) - int(cell_shift[kk][d]) for d in range(3))))
    return out


def _add(a: Shift, b: Shift) -> Shift:
    return a[0] + b[0], a[1] + b[1], a[2] + b[2]


def _ring_bond(i: int, j: int, s: Shift, edges, organic: set[int]) -> bool:
    """Il legame i→j (j nella cella s) chiude un anello di non metalli di al massimo `UNIT_RING_SIZE` atomi: da j si
    torna a i nella stessa cella senza ripassare per questo legame. Un giro che torna a i in un'altra cella non è un
    anello ma la periodicità della rete."""
    frontier = [(j, s)]
    seen = {(j, s)}
    for step in range(UNIT_RING_SIZE - 1):  # un anello di n atomi: n − 1 passi da j a i
        nxt = []
        for a, sa in frontier:
            for b, es in edges[a]:
                sb = _add(sa, es)
                if b not in organic or (step == 0 and b == i and sb == (0, 0, 0)):  # non lo stesso legame all'indietro
                    continue
                if b == i and sb == (0, 0, 0):
                    return True
                if (b, sb) not in seen:
                    seen.add((b, sb))
                    nxt.append((b, sb))
        frontier = nxt
    return False


def _extended_unit(members: list[tuple[int, Shift]], edges, els: list[str]
                   ) -> tuple[list[tuple[int, Shift]], list[tuple[int, int]]]:
    """L'unità ripetuta di una rete estesa, con gli anelli e i leganti interi.

    La visita in ampiezza mette ogni atomo nella cella da cui lo raggiunge per prima: se arriva a un anello da due lati
    lungo cammini che differiscono di una traslazione, l'anello resta spezzato e i suoi frammenti π non si trovano.
    Qui le celle si riassegnano unendo i legami per livelli (union-find con le traslazioni): prima i legami degli
    anelli di non metalli, poi gli altri legami tra non metalli (i leganti), infine quelli con i metalli. Un legame che
    chiude un giro con una traslazione diversa da zero è la periodicità della rete e si scarta; quelli dei livelli
    precedenti sono già tutti coerenti, così anelli e leganti finiti restano interi. Il primo atomo resta nella cella
    (0, 0, 0)."""
    atoms = [i for i, _ in members]
    inside = set(atoms)
    organic = {i for i in atoms if els[i] in NONMETALS}
    ranked = []
    for i in atoms:
        for j, es in edges[i]:
            if j not in inside or j < i or (j == i and es < (0, 0, 0)):  # ogni legame una volta sola
                continue
            if i in organic and j in organic:
                level = 0 if _ring_bond(i, j, es, edges, organic) else 1
            else:
                level = 2
            ranked.append((level, i, j, es))
    ranked.sort()
    parent = {i: i for i in atoms}
    off: dict[int, Shift] = {i: (0, 0, 0) for i in atoms}  # cella di i rispetto al suo genitore

    def find(i: int) -> tuple[int, Shift]:
        path = []
        while parent[i] != i:
            path.append(i)
            i = parent[i]
        root, acc = i, (0, 0, 0)
        for k in reversed(path):  # compressione: ognuno punta alla radice con la traslazione accumulata
            acc = _add(off[k], acc)
            parent[k], off[k] = root, acc
        return root, acc

    for _level, i, j, es in ranked:
        ri, oi = find(i)
        rj, oj = find(j)
        if ri != rj:  # cella di j = cella di i + es
            parent[rj] = ri
            want = _add(oi, es)
            off[rj] = (want[0] - oj[0], want[1] - oj[1], want[2] - oj[2])
    _, o0 = find(atoms[0])
    cell = {}
    for i in atoms:
        _, oi = find(i)
        cell[i] = (oi[0] - o0[0], oi[1] - o0[1], oi[2] - o0[2])
    index = {i: k for k, i in enumerate(atoms)}
    bonds = [(k, index[j]) for k, i in enumerate(atoms) for j, es in edges[i]
             if index[j] > k and cell[j] == _add(cell[i], es)]
    return [(i, cell[i]) for i in atoms], bonds


def build_structure(cif: Cif, limits: Limits | None = None) -> Structure:
    """Ricompone le molecole intere partendo dagli atomi della cella e dai loro legami."""
    limits = limits or Limits()
    vectors = cell_vectors(cif.cell)
    origins: list[tuple[int, int]] = []
    els, labels, frac = fill_cell(cif, vectors, limits, origins)
    cart = to_cart(vectors, frac)
    edges = cell_bonds(els, cart, vectors, limits)
    seen: dict[int, tuple[int, int, int]] = {}  # atomo della cella → cella in cui è stato messo
    molecules: list[Molecule] = []
    inv = np.linalg.inv(vectors)
    for start in range(len(els)):
        if start in seen:
            continue
        members: list[tuple[int, tuple[int, int, int]]] = [(start, (0, 0, 0))]
        index = {start: 0}
        bonds: list[tuple[int, int]] = []
        extended = False
        periods: list[tuple[int, int, int]] = []
        seen[start] = (0, 0, 0)
        k = 0
        while k < len(members):
            i, shift = members[k]
            for j, es in edges[i]:
                s = (shift[0] + es[0], shift[1] + es[1], shift[2] + es[2])
                prev = seen.get(j)
                if prev is None:
                    seen[j] = s
                    index[j] = len(members)
                    members.append((j, s))
                    bonds.append((k, len(members) - 1))
                elif prev != s:
                    extended = True  # lo stesso atomo raggiunto in due celle diverse: la rete continua all'infinito
                    # la differenza è una periodicità della rete; quelle dei cicli chiusi dalla visita (gli archi fuori
                    # dall'albero) generano tutto il reticolo delle periodicità
                    periods.append((s[0] - prev[0], s[1] - prev[1], s[2] - prev[2]))
                elif index[j] > k:
                    bonds.append((k, index[j]))
            k += 1
        if extended:  # la visita in ampiezza può spezzare un anello tra due celle: si ricolloca l'unità ripetuta
            members, bonds = _extended_unit(members, edges, els)
        ci = np.array([m[0] for m in members])
        cs = np.array([m[1] for m in members], dtype=int)
        pos = to_cart(vectors, frac[ci] + cs)
        # la molecola si sposta di celle intere finché il suo baricentro non cade nella cella di partenza
        c = pos.mean(axis=0)
        fc = np.floor(c @ inv + 1e-6)
        back = -(fc @ vectors)
        pos = pos + back
        molecules.append(Molecule(
            el=[els[i] for i in ci], labels=[labels[i] for i in ci], pos=pos,
            bonds=_dedupe(bonds), centroid=c + back, extended=extended,
            cell_index=ci, cell_shift=cs - fc.astype(int), periods=lattice_basis(periods),
            links=network_links(ci, cs, edges) if extended else [],
        ))
    return Structure(vectors, els, labels, frac, molecules, list(cif.symops),
                     cell_origin=np.array(origins, dtype=int).reshape(-1, 2))


def _molecule_images(s: Structure, tol: float) -> list[list[Image | None]]:
    owner: dict[int, tuple[int, np.ndarray]] = {}  # atomo della cella → (molecola, cella in cui sta in quella molecola)
    for mi, m in enumerate(s.molecules):
        for ci, cs in zip(m.cell_index, m.cell_shift):
            owner[int(ci)] = (mi, np.asarray(cs))
    by_el: dict[str, np.ndarray] = {}
    for e in set(s.cell_el):
        by_el[e] = np.array([i for i, x in enumerate(s.cell_el) if x == e])
    centroids = [(s.cell_frac[m.cell_index] + m.cell_shift).mean(axis=0) for m in s.molecules]
    out: list[list[Image | None]] = []
    for op in s.symops:
        row: list[Image | None] = []
        for mi, m in enumerate(s.molecules):
            row.append(_image_of(s, op, mi, m, owner, by_el, centroids, tol))
        out.append(row)
    return out


def _image_of(s: Structure, op: SymOp, mi: int, m: Molecule, owner, by_el, centroids, tol: float) -> Image | None:
    if not len(m):
        return None
    f0 = s.cell_frac[m.cell_index[0]] + m.cell_shift[0]
    img = op.apply(f0)
    cand = by_el[m.el[0]]
    d = img - s.cell_frac[cand]
    d -= np.round(d)
    dist = np.linalg.norm(d @ s.vectors, axis=1)
    k = int(np.argmin(dist))
    if dist[k] > tol:
        return None
    j = int(cand[k])
    if j not in owner:
        return None
    mj, sj = owner[j]
    t = np.rint(img - s.cell_frac[j] - sj).astype(int)
    # controllo: anche il baricentro deve andare nel baricentro dell'immagine (stessa molecola intera)
    gc = op.apply(centroids[mi])
    if np.linalg.norm((gc - centroids[mj] - t) @ s.vectors) > tol or len(s.molecules[mj]) != len(m):
        return None
    return mj, (int(t[0]), int(t[1]), int(t[2]))


def _atom_map(s: Structure, op: SymOp, m: Molecule, img: Image, tol: float) -> np.ndarray | None:
    mj, t = img
    target = s.molecules[mj]
    f_src = s.cell_frac[m.cell_index] + m.cell_shift
    f_img = f_src @ op.r.T + op.t
    f_tgt = s.cell_frac[target.cell_index] + target.cell_shift + np.asarray(t)
    d = np.linalg.norm((f_img[:, None, :] - f_tgt[None, :, :]) @ s.vectors, axis=2)
    j = d.argmin(axis=1)
    if (d[np.arange(len(j)), j] > tol).any() or len(set(j.tolist())) != len(j):
        return None
    return j


def _dedupe(b: list[tuple[int, int]]) -> np.ndarray:
    if not b:
        return np.zeros((0, 2), dtype=int)
    arr = np.sort(np.array(b, dtype=int), axis=1)
    return np.unique(arr, axis=0)
