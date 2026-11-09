"""Anelli e nuclei π di una molecola.

Per ogni molecola si cerca il nucleo π: il sistema di anelli condensati planari più grande
(in mancanza, tutta la molecola senza idrogeni, se è planare). Del nucleo si calcolano il piano medio
(minimi quadrati), il baricentro e l'asse lungo (direzione di massima estensione nel piano).

Porting di `piCore`/`smallRings` di `src/lib/crystal/pistack.ts`.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, replace

import numpy as np

from pistak.core.config import Limits
from pistak.core.elements import NONMETALS
from pistak.core.geometry import Plane, fit_plane
from pistak.core.structure import Molecule

#: Sotto questo rapporto √(λ2/λ3) la nuvola di punti è quasi una retta: il piano medio (e la sua normale) non è definito.
MIN_PLANE_ASPECT = 0.05
#: Anelli più grandi cercati nei frammenti uniti da `fragment_merge = "conjugated"` (macrocicli porfirinici: 16).
MAX_MACROCYCLE = 20


@dataclass(frozen=True)
class Core:
    """Nucleo π: piano medio e atomi (indici nella molecola) che lo compongono."""

    plane: Plane
    atoms: np.ndarray
    kind: str  # "rings" | "molecule"
    #: Da dove viene l'asse lungo: "pca" (direzione di massima estensione) | "atoms" (definito dall'utente).
    axis_source: str = "pca"
    #: Gli anelli del frammento (indici nella molecola): i loro poligoni, proiettati, danno la sovrapposizione delle
    #: coppie (unione di poligoni convessi, non l'involucro convesso del frammento). Vuoto: si usa l'involucro.
    rings: tuple[tuple[int, ...], ...] = ()

    @property
    def centroid(self) -> np.ndarray:
        return self.plane.centroid

    @property
    def normal(self) -> np.ndarray:
        return self.plane.normal

    @property
    def long(self) -> np.ndarray:
        return self.plane.long

    @property
    def rms(self) -> float:
        return self.plane.rms

    def moved(self, v: np.ndarray) -> Core:
        """Lo stesso nucleo traslato di `v` (la molecola vicina in un'altra cella)."""
        return replace(self, plane=replace(self.plane, centroid=self.plane.centroid + v))

    def with_long(self, direction: np.ndarray) -> Core | None:
        """Lo stesso nucleo con l'asse lungo lungo `direction` (proiettata sul piano); None se è quasi normale al piano."""
        n = self.plane.normal
        v = np.asarray(direction, dtype=float) - n * np.dot(direction, n)
        norm = np.linalg.norm(v)
        if norm < 1e-3 * max(1.0, float(np.linalg.norm(direction))):
            return None
        return replace(self, plane=replace(self.plane, long=v / norm), axis_source="atoms")


def small_rings(n: int, bonds: np.ndarray, allowed, max_size: int = 8) -> list[list[int]]:
    """Anelli fino a `max_size` atomi: per ogni legame, tutti i cicli più corti che lo contengono.

    Tenerne uno solo per legame dipende dall'ordine degli atomi e perde gli anelli con tutti i legami condivisi con
    anelli della stessa taglia (l'anello centrale del coronene); con tutti i cicli più corti il risultato non dipende
    dall'ordine e per i policicli benzenoidi dà esattamente gli esagoni."""
    adj: list[list[int]] = [[] for _ in range(n)]
    ok = [(int(p), int(q)) for p, q in bonds if allowed(p) and allowed(q)]
    for p, q in ok:
        adj[p].append(q)
        adj[q].append(p)
    rings: dict[tuple[int, ...], list[int]] = {}
    for p, q in ok:
        # distanze da p senza usare il legame p–q, con tutti i predecessori sui cammini più corti
        dist = {p: 0}
        preds: dict[int, list[int]] = {p: []}
        queue = deque([p])
        while queue:
            u = queue.popleft()
            if q in dist and dist[u] >= dist[q] or dist[u] >= max_size - 1:
                break
            for w in adj[u]:
                if u == p and w == q:
                    continue
                if w not in dist:
                    dist[w] = dist[u] + 1
                    preds[w] = [u]
                    queue.append(w)
                elif dist[w] == dist[u] + 1:
                    preds[w].append(u)
        if q not in dist or not 3 <= dist[q] + 1 <= max_size:
            continue
        stack = [[q]]  # tutti i cammini più corti q → p, all'indietro sui predecessori
        while stack:
            path = stack.pop()
            if path[-1] == p:
                rings.setdefault(tuple(sorted(path)), path)
                continue
            stack.extend(path + [u] for u in sorted(preds[path[-1]], reverse=True))
    return list(rings.values())


def planar_rings(m: Molecule, limits: Limits) -> list[tuple[list[int], Plane]]:
    """Gli anelli planari della molecola (senza idrogeni), ciascuno con il suo piano medio.

    Due passate: prima gli anelli tra i soli atomi non metallici, poi quelli che contengono un metallo. In una sola
    passata i cicli più corti attraverso un legame C–C di un Cp o di un arene η⁵/η⁶ sono i triangoli M–C–C, e
    l'anello aromatico non verrebbe mai trovato (ferrocene, Cr(η⁶-C₆H₆)(CO)₃: revisione, voce 51)."""
    heavy = m.heavy
    metal = [_is_metal(e) for e in m.el]
    found: dict[tuple[int, ...], list[int]] = {}
    for r in small_rings(len(m), m.bonds, lambda i: heavy[i] and not metal[i], limits.max_ring_size):
        found.setdefault(tuple(sorted(r)), r)
    if any(metal):
        for r in small_rings(len(m), m.bonds, lambda i: heavy[i], limits.max_ring_size):
            if any(metal[i] for i in r):
                found.setdefault(tuple(sorted(r)), r)
    out = []
    for r in found.values():
        pl = fit_plane(m.pos[r])
        if pl.rms < limits.ring_rms:
            out.append((r, pl))
    return out


def pi_core(m: Molecule, limits: Limits) -> Core | None:
    """Nucleo π principale di una molecola, o None se non ha parti planari adatte."""
    frags = pi_fragments(m, limits.with_values(fragments="largest"))
    return frags[0] if frags else None


def _is_metal(el: str) -> bool:
    return el not in NONMETALS


def _fuse(rings: list[list[int]]) -> list[tuple[set[int], list[list[int]]]]:
    """Anelli condensati (che condividono almeno due atomi) riuniti in sistemi: (atomi, anelli) per sistema."""
    groups: list[tuple[set[int], list[list[int]]]] = []
    for r in rings:
        touching = [g for g in groups if sum(1 for i in r if i in g[0]) >= 2]
        atoms, members = set(r), [r]
        for g in touching:
            atoms |= g[0]
            members += g[1]
            groups.remove(g)
        groups.append((atoms, members))
    return groups


def _planar_systems(m: Molecule, limits: Limits) -> list[tuple[list[int], list[list[int]]]]:
    """I sistemi di anelli planari condensati che passano la soglia di planarità `molecule_rms`.

    Un sistema che non la passa e contiene anelli con un metallo si riprova senza quegli anelli: in un metallocene o in
    un complesso η⁶-arene i triangoli M–C–C (legami M–C entro la tolleranza) sono «anelli» planari condensati con il Cp
    o l'arene, e il sistema intero, con il metallo a ~1,6 Å dal piano, non è planare; senza di loro resta l'anello
    aromatico. Un chelato planare con il metallo (metallo-porfirina) passa la soglia così com'è e resta intero."""
    rings = [r for r, _ in planar_rings(m, limits)]
    out = []
    for atoms, members in _fuse(rings):
        g = sorted(atoms)
        if fit_plane(m.pos[g]).rms < limits.molecule_rms:
            out.append((g, members))
            continue
        organic = [r for r in members if not any(_is_metal(m.el[i]) for i in r)]
        if len(organic) < len(members):
            for atoms2, members2 in _fuse(organic):
                g2 = sorted(atoms2)
                if fit_plane(m.pos[g2]).rms < limits.molecule_rms:
                    out.append((g2, members2))
    return out


def _conjugated(m: Molecule, systems: list[tuple[list[int], list[list[int]]]], limits: Limits) -> list[tuple[list[int], list[list[int]]]]:
    """`fragment_merge = "conjugated"`: unisce i sistemi di anelli complanari collegati da un legame o da atomi sp/sp²
    (non metalli, al più tre vicini: carboni meso delle porfirine, ponti vinilici, alchini) e aggiunge gli atomi
    sp/sp² esociclici complanari. Si cresce un frammento alla volta, dal sistema più grande, accettando un atomo o un
    sistema adiacente solo se l'insieme resta entro `molecule_rms`: un biarile ruotato resta in due frammenti."""
    heavy = m.heavy
    nbr: list[list[int]] = [[] for _ in range(len(m))]
    for p, q in m.bonds:
        nbr[int(p)].append(int(q))
        nbr[int(q)].append(int(p))
    owner = {i: k for k, (g, _) in enumerate(systems) for i in g}
    linker = {i for i in range(len(m)) if heavy[i] and i not in owner and not _is_metal(m.el[i]) and len(nbr[i]) <= 3}
    used: set[int] = set()
    out = []
    for k in sorted(range(len(systems)), key=lambda k: -len(systems[k][0])):
        if k in used:
            continue
        used.add(k)
        atoms, rings = set(systems[k][0]), list(systems[k][1])
        grown = True
        while grown:
            grown = False
            frontier = sorted({j for i in atoms for j in nbr[i] if j not in atoms and (j in linker or (j in owner and owner[j] not in used))})
            for j in frontier:
                if j in atoms or (j in owner and owner[j] in used):
                    continue  # sistema già unito attraverso un altro atomo della frontiera (anelli doppi: revisione 2)
                add = set(systems[owner[j]][0]) if j in owner else {j}
                trial = sorted(atoms | add)
                if fit_plane(m.pos[trial]).rms < limits.molecule_rms:
                    atoms |= add
                    if j in owner:
                        rings += systems[owner[j]][1]
                        used.add(owner[j])
                    grown = True
        if len(rings) > len(systems[k][1]):
            # macrocicli (anello interno di 16 atomi di porfirine e ftalocianine): non sono anelli «piccoli», ma l'area
            # che racchiudono fa parte del frammento e conta nella sovrapposizione delle coppie
            known = {tuple(sorted(r)) for r in rings}
            rings += [r for r in small_rings(len(m), m.bonds, lambda i, s=frozenset(atoms): i in s, MAX_MACROCYCLE)
                      if len(r) > limits.max_ring_size and tuple(sorted(r)) not in known]
        out.append((sorted(atoms), rings))
    return out


def pi_fragments(m: Molecule, limits: Limits) -> list[Core]:
    """I frammenti π della molecola: i sistemi di anelli planari condensati (che condividono almeno un legame), dal
    più grande. Con `limits.fragments == "largest"` solo il primo, il nucleo π del sito; con "all" anche gli altri
    con almeno `min_fragment_atoms` atomi (la soglia non vale per il primo). Con `fragment_merge == "conjugated"` i
    sistemi complanari collegati da legami o atomi sp/sp² diventano un frammento solo (`_conjugated`). Senza anelli
    adatti: tutta la molecola senza idrogeni, se è planare e non lineare."""
    systems = _planar_systems(m, limits)
    if getattr(limits, "fragment_merge", "fused") == "conjugated" and systems:
        systems = _conjugated(m, systems, limits)
    systems.sort(key=lambda s: len(s[0]), reverse=True)  # ordinamento stabile come nel sito
    if systems:
        keep = systems[:1] if limits.fragments == "largest" else \
            systems[:1] + [s for s in systems[1:] if len(s[0]) >= limits.min_fragment_atoms]
        return [Core(fit_plane(m.pos[g]), np.array(g), "rings", rings=tuple(tuple(r) for r in rs)) for g, rs in keep]
    allh = np.flatnonzero(m.heavy)
    if len(allh) < 4:
        return []
    pl = fit_plane(m.pos[allh])
    if pl.rms >= limits.molecule_rms or pl.aspect < MIN_PLANE_ASPECT:  # lineare: normale non definita
        return []
    return [Core(pl, allh, "molecule")]


def axis_atoms(m: Molecule, spec: str) -> tuple[int, int] | None:
    """Gli atomi che definiscono l'asse lungo secondo `spec` ("C1,C9": etichette del CIF), in questa molecola: tra
    tutti gli atomi con quelle etichette (una molecola su un elemento di simmetria ne ha più copie, "C1,C1" vale) la
    coppia più lontana. None se la molecola non li contiene."""
    parts = [p.strip() for p in spec.replace(";", ",").split(",") if p.strip()]
    if len(parts) != 2:
        return None
    ia = [i for i, lab in enumerate(m.labels) if lab == parts[0]]
    ib = [i for i, lab in enumerate(m.labels) if lab == parts[1]]
    pairs = [(i, j) for i in ia for j in ib if i != j]
    if not pairs:
        return None
    return max(pairs, key=lambda ij: (float(np.linalg.norm(m.pos[ij[1]] - m.pos[ij[0]])), -ij[0], -ij[1]))


def apply_axis(m: Molecule, frags: list[Core], spec: str) -> list[Core]:
    """Asse lungo definito dall'utente: si applica ai frammenti che contengono entrambi gli atomi, o, se nessuno li
    contiene, al frammento principale. Il verso va dal primo atomo al secondo."""
    ij = axis_atoms(m, spec) if spec else None
    if ij is None or not frags:
        return frags
    i, j = ij
    direction = m.pos[j] - m.pos[i]
    targets = [k for k, f in enumerate(frags) if i in f.atoms and j in f.atoms] or [0]
    out = list(frags)
    for k in targets:
        out[k] = frags[k].with_long(direction) or frags[k]
    return out
