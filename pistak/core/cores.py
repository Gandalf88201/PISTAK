"""Anelli e nuclei π di una molecola.

Per ogni molecola si cerca il nucleo π: il sistema di anelli condensati planari più grande
(in mancanza, tutta la molecola senza idrogeni, se è planare e sp²/sp). Del nucleo si calcolano il piano medio
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
    #: Da dove viene l'asse lungo: "pca" (direzione di massima estensione degli atomi del frammento) | "atoms"
    #: (definito dall'utente con due atomi) | "molecule" (asse lungo di tutta la molecola, proiettato sul piano del
    #: frammento) | "substituents" (asse del frammento con i suoi sostituenti): gli ultimi due sono il ripiego per i
    #: frammenti quasi isotropi (`fallback_axis`, revisione 4, M1.1).
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

    def with_long(self, direction: np.ndarray, source: str = "atoms") -> Core | None:
        """Lo stesso nucleo con l'asse lungo lungo `direction` (proiettata sul piano); None se è quasi normale al piano.
        `source` dice da dove viene l'asse (`axis_source`)."""
        n = self.plane.normal
        v = np.asarray(direction, dtype=float) - n * np.dot(direction, n)
        norm = np.linalg.norm(v)
        if norm < 1e-3 * max(1.0, float(np.linalg.norm(direction))):
            return None
        return replace(self, plane=replace(self.plane, long=v / norm), axis_source=source)


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


#: Elementi di cui un anello deve contenerne almeno uno per contare in `pi_scope = "organic"`: tutti i carbocicli ed
#: eterocicli, le azine e le triazine (melammina, acido cianurico), tetrazoli e pentazoli, borazine e boroxine,
#: fosfazeni, anelli S–N. Restano fuori gli anelli di soli O, S, Se, Te, P, As, Si, Ge (silicati, fosfati, solfati
#: ciclici, anelli e gabbie di calcogeni e pnicogeni).
PI_RING_ELEMENTS = frozenset({"B", "C", "N"})
#: Non metalli che non chiudono anelli π: idrogeno, alogeni e gas nobili (monovalenti o senza legami; un «anello» con
#: loro viene da contatti entro la tolleranza dei legami, come i triangoli Br₃ di un polibromuro).
_NOT_IN_PI_RINGS = frozenset({"H", "D", "F", "Cl", "Br", "I", "At", "He", "Ne", "Ar", "Kr", "Xe", "Rn"})


def has_pi_ring(m: Molecule, limits: Limits) -> bool:
    """Vero se la molecola (o l'unità ripetuta) ha un anello planare (`ring_rms`, al più `max_ring_size` atomi) fatto
    di soli non metalli e con almeno un atomo di B, C o N: il criterio ad anelli di `pi_scope = "organic"` (revisione 3,
    N1). Gli anelli con un metallo non contano: i chelati M–O–C–O di un carbonato, i quadrati M–C≡N–M dei cianuri e gli
    anelli Ta/V/C di un carburo hanno due non metalli adiacenti ma non sono sistemi π."""
    allowed = [bool(h) and e in NONMETALS and e not in _NOT_IN_PI_RINGS for h, e in zip(m.heavy, m.el)]
    for r in small_rings(len(m), m.bonds, lambda i: allowed[i], limits.max_ring_size):
        if any(m.el[i] in PI_RING_ELEMENTS for i in r) and fit_plane(m.pos[r]).rms < limits.ring_rms:
            return True
    return False


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


#: Somma minima (°) dei tre angoli attorno a un atomo dell'anello con tre vicini non metallici perché sia sp² (360° se
#: planare; ~328° per un sp³ con un idrogeno mancante o un'ammina piramidale). Con `ring_atoms = "sp2"`.
MIN_SP2_ANGLE_SUM = 345.0


def _neighbours(m: Molecule) -> list[list[int]]:
    nbr: list[list[int]] = [[] for _ in range(len(m))]
    for p, q in m.bonds:
        nbr[int(p)].append(int(q))
        nbr[int(q)].append(int(p))
    return nbr


def sp2_ring(m: Molecule, ring: list[int], nbr: list[list[int]] | None = None) -> bool:
    """Vero se ogni atomo non metallico dell'anello può essere sp² o sp: al più tre vicini non metallici (idrogeni
    compresi; i legami ai metalli non contano, come in un Cp o in una piridina coordinati) e, con tre vicini, nel loro
    piano (somma degli angoli ≥ `MIN_SP2_ANGLE_SUM`). Esclude dai frammenti π gli anelli saturi che sono planari per
    simmetria o per caso (ciclobutano su un centro d'inversione, anelli di chelati sp³), i triangoli delle gabbie
    (carborani: boro con cinque o sei vicini) e i ciclopropani. Un atomo sp³ senza idrogeni nel CIF e con due soli vicini
    non si riconosce (revisione del paper, 13 novembre 2026, punto 1)."""
    if nbr is None:
        nbr = _neighbours(m)
    for i in ring:
        if _is_metal(m.el[i]):
            continue
        nm = [j for j in nbr[i] if not _is_metal(m.el[j])]
        if len(nm) >= 4:
            return False
        if len(nm) == 3:
            v = [m.pos[j] - m.pos[i] for j in nm]
            v = [x / np.linalg.norm(x) for x in v]
            total = sum(np.degrees(np.arccos(np.clip(np.dot(v[a], v[b]), -1.0, 1.0)))
                        for a, b in ((0, 1), (0, 2), (1, 2)))
            if total < MIN_SP2_ANGLE_SUM:
                return False
    return True


def _planar_systems(m: Molecule, limits: Limits) -> list[tuple[list[int], list[list[int]]]]:
    """I sistemi di anelli planari condensati che passano la soglia di planarità `molecule_rms`.

    Con `ring_atoms = "sp2"` (predefinito dalla v1.4.0) contano solo gli anelli di almeno quattro atomi i cui atomi non
    metallici sono sp² o sp (`sp2_ring`): i triangoli (M–C–C di un η²-legante, cluster, X₃) non sono sistemi π; resta
    fuori anche il catione ciclopropenio, l'unico triangolo aromatico. Un sistema che non passa `molecule_rms` e contiene anelli con un metallo si riprova senza quegli
    anelli: in un metallocene o in un complesso η⁶-arene i triangoli M–C–C (legami M–C entro la tolleranza) sono
    «anelli» planari condensati con il Cp o l'arene, e il sistema intero, con il metallo a ~1,6 Å dal piano, non è
    planare; senza di loro resta l'anello aromatico. Un chelato planare con il metallo (metallo-porfirina) passa la
    soglia così com'è e resta intero."""
    rings = [r for r, _ in planar_rings(m, limits)]
    if getattr(limits, "ring_atoms", "sp2") == "sp2":
        nbr = _neighbours(m)
        rings = [r for r in rings if len(r) >= 4 and sp2_ring(m, r, nbr)]
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
    adatti e solo con `ringless_fragments = "planar"`: tutta la molecola senza idrogeni, se è planare, non lineare e
    sp²/sp (nessun atomo con quattro vicini o più, `_sp2_skeleton`)."""
    systems = _planar_systems(m, limits)
    if getattr(limits, "fragment_merge", "fused") == "conjugated" and systems:
        systems = _conjugated(m, systems, limits)
    systems.sort(key=lambda s: len(s[0]), reverse=True)  # ordinamento stabile come nel sito
    if systems:
        keep = systems[:1] if limits.fragments == "largest" else \
            systems[:1] + [s for s in systems[1:] if len(s[0]) >= limits.min_fragment_atoms]
        return [Core(fit_plane(m.pos[g]), np.array(g), "rings", rings=tuple(tuple(r) for r in rs)) for g, rs in keep]
    if getattr(limits, "ringless_fragments", "none") != "planar":
        return []  # predefinito dalla revisione 4: senza anelli nessun frammento π (vedi `Limits.ringless_fragments`)
    allh = np.flatnonzero(m.heavy)
    if len(allh) < 4 or not _sp2_skeleton(m, allh):
        return []
    pl = fit_plane(m.pos[allh])
    if pl.rms >= limits.molecule_rms or pl.aspect < MIN_PLANE_ASPECT:  # lineare: normale non definita
        return []
    return [Core(pl, allh, "molecule")]


#: Nel ripiego dello scheletro planare un atomo con questi vicini o più (H compresi) è sp³ (o ipervalente): niente π.
MAX_FALLBACK_NEIGHBOURS = 4


def _sp2_skeleton(m: Molecule, atoms: np.ndarray) -> bool:
    """Vero se nessun atomo dello scheletro ha `MAX_FALLBACK_NEIGHBOURS` vicini legati o più, contando tutti i vicini
    (idrogeni e metalli compresi): il ripiego tiene solo scheletri sp²/sp (nitrato, carbonato, guanidinio, urea) e
    scarta le catene sp³ che per caso giacciono in un piano (1,2-dicloroetano, etilendiammina, n-pentano, acido
    glutarico, etere etilico, THF: revisione 3, N5). Scarta anche le molecole con un metile (acetone, DMF, acetato): il
    CH₃ ha quattro vicini. Limite: se il CIF non ha gli idrogeni, un CH₂ ha due vicini e passa come sp²; il criterio è
    affidabile solo con gli H nel modello (lo sono nella quasi totalità delle strutture molecolari depositate)."""
    count = np.zeros(len(m), dtype=int)
    for i, j in m.bonds:
        count[int(i)] += 1
        count[int(j)] += 1
    return bool(np.all(count[atoms] < MAX_FALLBACK_NEIGHBOURS))


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


#: Origini dell'asse lungo che non dipendono dall'isotropia del frammento (`pairs.axis_defined`).
FALLBACK_SOURCES = ("molecule", "substituents")


def fallback_atoms(m: Molecule, core: Core, mode: str) -> np.ndarray:
    """Gli atomi (non H) che definiscono l'asse di ripiego di un frammento: tutta la molecola (`mode = "molecule"`)
    oppure il frammento con gli atomi legati direttamente ai suoi (`"substituents"`)."""
    heavy = m.heavy
    if mode == "molecule":
        return np.flatnonzero(heavy)
    own = set(int(i) for i in core.atoms)
    out = set(own)
    for i, j in m.bonds:
        i, j = int(i), int(j)
        if i in own and heavy[j]:
            out.add(j)
        elif j in own and heavy[i]:
            out.add(i)
    return np.array(sorted(out))


def fallback_axis(m: Molecule, core: Core, limits: Limits) -> tuple[np.ndarray, float] | None:
    """Asse lungo di ripiego per un frammento quasi isotropo (revisione 4, M1.1): la direzione di massima estensione,
    nel piano del frammento, degli atomi di `fallback_atoms` proiettati su quel piano, con la sua isotropia
    (λ2/λ3)^½ in proiezione. None se anche questa supera `max_axis_isotropy` (benzene, mesitilene, C₆F₆: nessun asse
    lungo nemmeno con i sostituenti) o se la molecola è una rete estesa nel modo "molecule" (la sua unità ripetuta è un
    taglio arbitrario)."""
    mode = getattr(limits, "axis_fallback", "molecule")
    if mode not in FALLBACK_SOURCES or (mode == "molecule" and m.extended):
        return None
    idx = fallback_atoms(m, core, mode)
    if len(idx) <= len(core.atoms):
        return None
    n = core.normal
    d = m.pos[idx] - m.pos[idx].mean(axis=0)
    d = d - np.outer(d @ n, n)  # nel piano del frammento
    values, vectors = np.linalg.eigh(d.T @ d)  # crescente: ~0 (lungo n), λ2, λ3
    if values[2] <= 0:
        return None
    aspect = float(np.sqrt(max(0.0, values[1]) / values[2]))
    if aspect > limits.max_axis_isotropy:
        return None
    v = vectors[:, 2]
    return (v if v[np.argmax(np.abs(v))] >= 0 else -v), aspect


def apply_fallback_axis(m: Molecule, frags: list[Core], limits: Limits) -> list[Core]:
    """Dà l'asse di ripiego (`fallback_axis`) ai frammenti il cui asse automatico non è definito (isotropia sopra
    `max_axis_isotropy`); gli altri, e quelli con l'asse definito dall'utente, restano come sono. `axis_source` dice
    quale asse è stato usato, così valori basati sull'anello e sulla molecola non si mescolano in silenzio."""
    out = []
    for f in frags:
        if f.axis_source == "pca" and f.plane.aspect > limits.max_axis_isotropy:
            hit = fallback_axis(m, f, limits)
            if hit is not None:
                f = f.with_long(hit[0], source=getattr(limits, "axis_fallback", "molecule")) or f
        out.append(f)
    return out
