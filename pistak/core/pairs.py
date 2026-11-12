"""Coppie di frammenti π vicini e geometria del loro π-stacking.

Per ogni coppia, con d il vettore tra i baricentri dei due nuclei (frammenti π) e n la normale media:
- distanza interplanare = |d·n|;
- slittamento = componente di d nel piano, s = d − (d·n)n;
- slittamento longitudinale = |s·L| lungo l'asse lungo L, slittamento laterale = |s·W| lungo W = n × L.
Una coppia è un π-stack se i piani sono quasi paralleli, abbastanza vicini e i nuclei, proiettati sul piano,
si sovrappongono quanto chiede il criterio di stack (`Limits.overlap_measure`, `min_overlap`, `max_slip`).

Ogni molecola ha uno o più frammenti π (`cores.pi_fragments`, `Limits.fragments`): una coppia è tra il frammento
`fa` della molecola `a` e il frammento `fb` della molecola `b` (traslata di `shift`). Tra due molecole a contatto si
misurano il frammento principale di ciascuna (0–0) e le altre combinazioni di frammenti abbastanza vicini; con più
frammenti nella stessa molecola anche le coppie intramolecolari (a = b, shift = 0).

Due modi di raccogliere le coppie (`Limits.pair_mode`):

* ``site`` (porting di `measurePair`/`neighbourPairs` di `src/lib/crystal/pistack.ts`): ogni contatto compare due
  volte, a→b e b→a; l'asse L è quello della prima molecola, quindi i due versi possono dare slittamenti diversi;
  le coppie «equivalenti» si riconoscono dai valori arrotondati; si esplorano sempre ±`neighbour_shells` celle.
* ``canonical``: ogni contatto compare una volta (a < b, oppure a = b con la traslazione di segno negativo);
  L è la bisettrice degli assi lunghi dei due frammenti, quindi i valori non dipendono da quale si chiama A;
  due contatti sono equivalenti se un'operazione del gruppo spaziale porta l'uno nell'altro (molecole e frammenti);
  le celle esplorate bastano sempre a trovare tutti i contatti (almeno ±`neighbour_shells`). `count` = contatti
  per cella.
"""
from __future__ import annotations

import itertools
import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.spatial.distance import cdist

from pistak.core.config import Limits
from pistak.core.cores import FALLBACK_SOURCES, Core
from pistak.core.elements import covalent_radius
from pistak.core.geometry import clip_convex, convex_hull, points_inside, polygon_area, unit
from pistak.core.structure import Structure, same_network

KINDS = ("stack", "parallel", "tilted")  # in ordine di interesse
Shift = tuple[int, int, int]
#: Limite alle celle esplorate per asse in modalità canonical (celle molto piccole o molto oblique).
MAX_SHELLS = 8
#: Due frammenti secondari formano una coppia se i loro atomi sono entro `contact_cutoff` + questo margine (Å).
FRAGMENT_MARGIN = 1.0
#: (a, b, shift, fa, fb): un contatto tra il frammento fa della molecola a e il frammento fb di b traslata di shift.
Member = tuple[int, int, Shift, int, int]
Fragments = list[list[Core]]
#: Chiavi dei parametri rinominate: il nome vecchio vale ancora in `Pair.get` e nel registro.
RENAMED_KEYS = {"beta_a": "beta", "beta_b": "gamma"}


@dataclass
class Frame:
    """Sistema di riferimento della coppia: origine nel baricentro del primo nucleo,
    assi L, W nel piano medio, n normale."""

    origin: np.ndarray
    L: np.ndarray
    W: np.ndarray
    n: np.ndarray

    def project(self, p: np.ndarray) -> np.ndarray:
        """Coordinate (L, W, n) di uno o più punti 3D."""
        return (np.asarray(p) - self.origin) @ np.column_stack([self.L, self.W, self.n])


def pair_frame(ca: Core, cb: Core, symmetric: bool = False) -> Frame:
    """n = normale media; L = asse lungo di A proiettato sul piano medio, oppure (``symmetric``) la bisettrice degli
    assi lunghi di A e B, che non cambia scambiando le due molecole."""
    nb = -cb.normal if np.dot(ca.normal, cb.normal) < 0 else cb.normal
    n = unit(ca.normal + nb)
    la = ca.long - n * np.dot(ca.long, n)
    if symmetric:
        lb = cb.long - n * np.dot(cb.long, n)
        if np.dot(la, lb) < 0:
            lb = -lb
        if np.linalg.norm(la + lb) > 1e-6:
            la = la + lb
    L = unit(la)
    return Frame(ca.centroid, L, np.cross(n, L), n)


def frame_for(ca: Core, cb: Core, limits: Limits) -> Frame:
    return pair_frame(ca, cb, symmetric=limits.pair_mode != "site")


@dataclass
class Pair:
    """Una coppia di frammenti π vicini. `a`, `b` sono indici in `Structure.molecules`, `fa`, `fb` i frammenti."""

    a: int
    b: int
    shift: Shift  # traslazione (in celle) della seconda molecola
    kind: str  # stack | parallel | tilted
    interplanar: float
    long_slip: float
    lat_slip: float
    slip: float
    dihedral: float
    centroid_distance: float
    #: Contatto più corto tra atomi pesanti delle due molecole (intramolecolare: tra i due frammenti).
    contact: float
    #: Frazione degli atomi del nucleo di una molecola che cade, in proiezione, dentro il contorno dell'altra.
    overlap: float
    #: Quante coppie equivalenti sono state raccolte in questa (canonical: contatti per cella).
    count: int = 1
    #: Tutte le coppie raccolte (a, b, shift, fa, fb), la prima è questa.
    members: list[Member] = field(default_factory=list)
    #: Valori di tutti i parametri del registro (`pistak.params`), per chiave; riempiti da `compute_all`.
    values: dict[str, float] = field(default_factory=dict)
    index: int = -1
    fa: int = 0
    fb: int = 0
    #: Vero se almeno una delle due molecole ha più di un frammento π (allora conta quali sono fa e fb).
    fragmented: bool = False
    #: Area comune dei contorni proiettati / area del contorno minore (0–1): la sovrapposizione del criterio
    #: `overlap_measure = "area"`. Si misura sempre, così cambiare criterio non richiede di rimisurare.
    overlap_area: float = math.nan
    #: Vero se una delle due molecole è una rete estesa (polimero di coordinazione, catena covalente): i frammenti π
    #: sono anelli della rete, il contatto è quello tra i due frammenti e la coppia non si esporta come dimero.
    polymer: bool = False
    #: Vero se i due frammenti hanno un atomo in comune o un legame covalente tra loro (frammenti della stessa molecola o
    #: della stessa rete estesa: due anelli chelati che condividono il metallo, un biarile). Non è mai un π-stack: la
    #: coppia resta parallel o tilted, con tutti i suoi parametri.
    bonded: bool = False
    #: Vero se i due frammenti, senza essere legati tra loro, sono legati entrambi a uno stesso atomo (il metallo di un
    #: metallocene o di un complesso bis-arene, il metallo di un chelato e di un legante): come `bonded`, non è mai un
    #: π-stack (revisione 2, voce 58).
    bridged: bool = False
    #: Angolo (0–90°) tra gli assi lunghi proiettati dei due frammenti; NaN se uno dei due assi non è definito.
    twist: float = math.nan
    #: Vero se d·n_A e d·n_B hanno segni opposti: la retta d'intersezione dei due piani passa tra i centroidi e la
    #: «distanza interplanare» d⊥ = |h_A + h_B| / (2 cos θ/2) è più piccola di entrambe le distanze perpendicolari.
    planes_cross: bool = False
    #: Slittamento longitudinale e laterale nel sistema di ciascun frammento (revisione 4, M1.2): la posizione del
    #: centroide dell'altro frammento proiettata sul piano di A (normale n_A) e divisa lungo l'asse lungo di A e la
    #: direzione n_A × l_A; lo stesso per B. Sono le distanze di «pitch» e «roll» di Curtis et al. (2004) viste da una
    #: molecola. Non dipendono dalla rotazione tra i due assi: sono definiti anche quando L, la bisettrice, non lo è
    #: (`max_axis_twist`). NaN se l'asse lungo di quel frammento non è definito. A rotazione e angolo nulli coincidono
    #: con `long_slip` e `lat_slip`.
    long_slip_a: float = math.nan
    lat_slip_a: float = math.nan
    long_slip_b: float = math.nan
    lat_slip_b: float = math.nan
    #: Da dove viene l'asse lungo di ciascun frammento: "pca" (il frammento stesso), "atoms" (due atomi scelti),
    #: "molecule" o "substituents" (ripiego per i frammenti isotropi), "" se non è definito.
    axis_a: str = ""
    axis_b: str = ""
    #: Vero se i due frammenti stanno nella stessa molecola o nella stessa rete estesa e la distanza interplanare è
    #: sotto `min_interplanar`: la separazione è imposta dall'impalcatura che li unisce (revisione 4, M6.4). Come
    #: `bonded` e `bridged` non è mai un π-stack; a differenza di quelli dipende da una soglia dello stadio `pairs`.
    clamped: bool = False
    #: Vero se i due frammenti appartengono allo stesso insieme di atomi legati: la stessa molecola discreta o la
    #: stessa rete estesa, anche in due unità ripetute diverse (revisione 5, N1). Per una molecola discreta coincide con
    #: `intramolecular`; per una rete non dipende da come è tagliata l'unità ripetuta, che `intramolecular` invece
    #: segue. Decide la bandiera `clamped` e la categoria della coppia nelle statistiche.
    same_network: bool = False

    def get(self, key: str, default: float = float("nan")) -> float:
        """Valore di un parametro del registro per questa coppia."""
        return self.values.get(RENAMED_KEYS.get(key, key), default)

    @property
    def intramolecular(self) -> bool:
        return self.a == self.b and self.shift == (0, 0, 0)

    @property
    def per_molecule(self) -> int:
        """Quanti contatti di questa classe ha la molecola A (il frammento `fa` della molecola `a`): `count` è per
        cella e dipende da Z, questo no (revisione 4, T16). Un contatto con una copia traslata della stessa molecola
        conta due volte, perché la molecola ha quel vicino da tutte e due le parti (una colonna: 2). Intramolecolare: 1.
        Nel modo del sito (coppie raccolte per valori) è `count`."""
        if self.intramolecular or not self.members or len(self.members[0]) < 5:
            return 1 if self.intramolecular else self.count
        me = (self.a, self.fa)
        return sum((m[0], m[3]) == me for m in self.members) + sum((m[1], m[4]) == me for m in self.members)

    @property
    def key(self) -> Member:
        return (self.a, self.b, self.shift, self.fa, self.fb)


def as_fragments(cores: Sequence) -> Fragments:
    """Accetta sia l'elenco dei nuclei principali (`Core | None` per molecola) sia quello dei frammenti."""
    return [list(c) if isinstance(c, (list, tuple)) else ([c] if c is not None else []) for c in cores]


def area_overlap(hull_a: np.ndarray, hull_b: np.ndarray) -> float:
    """Area comune dei due contorni divisa per l'area del più piccolo (0–1)."""
    if len(hull_a) < 3 or len(hull_b) < 3:
        return 0.0
    smaller = min(polygon_area(hull_a), polygon_area(hull_b))
    return polygon_area(clip_convex(hull_a, hull_b)) / smaller if smaller > 0 else 0.0


def classify(dihedral: float, interplanar: float, overlap: float, slip: float, limits: Limits, bonded: bool = False,
             tethered: bool = False) -> str:
    """stack | parallel | tilted secondo le soglie. `overlap` è già misurata come chiede `limits.overlap_measure`.
    Due frammenti legati tra loro o a uno stesso atomo (`bonded`: legati o a ponte) non sono mai uno stack; nemmeno due
    frammenti della stessa molecola o rete (`tethered`) più vicini di `min_interplanar` (`is_clamped`)."""
    if dihedral > limits.max_dihedral:
        return "tilted"
    if (not bonded and not (tethered and interplanar < getattr(limits, "min_interplanar", 0.0)) and interplanar <= limits.max_interplanar and overlap > limits.min_overlap
            and (limits.max_slip <= 0 or slip <= limits.max_slip)):
        return "stack"
    return "parallel"


def stack_overlap(p: Pair, limits: Limits) -> float:
    """La sovrapposizione della coppia misurata come chiede il criterio di stack."""
    return p.overlap_area if limits.overlap_measure == "area" else p.overlap


def tethered(p: Pair) -> bool:
    """Vero se i due frammenti stanno nella stessa molecola o nella stessa rete estesa: un'impalcatura di legami li
    tiene alla distanza a cui sono. Per una rete conta la rete, non l'unità ripetuta, il cui taglio è arbitrario: due
    unità traslate della stessa catena sono tenute dalla catena come due anelli della stessa unità; due catene
    parallele no (revisione 5, N1; prima contava solo l'unità ripetuta)."""
    return p.same_network or p.intramolecular


def is_clamped(p: Pair, limits: Limits) -> bool:
    """La coppia sarebbe un π-stack per angolo, distanza e sovrapposizione, ma i frammenti sono nella stessa molecola o
    rete e più vicini di `min_interplanar`: coppia «vincolata» (`Pair.clamped`)."""
    if p.bonded or p.bridged or not tethered(p) or p.interplanar >= getattr(limits, "min_interplanar", 0.0):
        return False
    return classify(p.dihedral, p.interplanar, stack_overlap(p, limits), p.slip, limits) == "stack"


def reclassify(p: Pair, limits: Limits) -> str:
    """Il tipo di una coppia già misurata secondo altre soglie (stadio `pairs`: nessuna geometria da rifare)."""
    return classify(p.dihedral, p.interplanar, stack_overlap(p, limits), p.slip, limits, p.bonded or p.bridged,
                    tethered(p))


#: Due posizioni più vicine di così sono lo stesso atomo (frammenti che lo condividono), come `d² ≥ 0,16` di `structure`.
SAME_ATOM = 0.4


def fragments_link(s: Structure, ca: Core, cb: Core, a: int, b: int, v: np.ndarray, limits: Limits) -> str:
    """"bonded" se i frammenti condividono un atomo o hanno un legame tra loro, "bridged" se un altro atomo (non H) è
    legato a entrambi, "" altrimenti. Stessa regola dei legami di `structure`. Solo per frammenti della stessa molecola o
    di reti estese."""
    if fragments_bonded(s, ca, cb, a, b, v, limits):
        return "bonded"
    ma, mb = s.molecules[a], s.molecules[b]
    if ma.extended or mb.extended:
        # rete estesa: l'atomo che lega i due frammenti può stare in un'altra copia dell'unità ripetuta, il cui taglio è
        # arbitrario; si cerca tra i vicini legati nella rete infinita (`Molecule.links`), e solo se i due frammenti
        # sono nella stessa rete: tra reti diverse, o tra una rete e una molecola, non c'è nessun legame (revisione 5, N1)
        shift = tuple(int(x) for x in np.rint(np.linalg.solve(s.vectors.T, v)))
        if a != b or not same_network(s, a, b, shift) or not ma.links:
            return ""
        own_a, own_b = {int(i) for i in ca.atoms}, {int(i) for i in cb.atoms}
        zero = (0, 0, 0)
        near_a = {(kk, t) for i in own_a for kk, t in ma.links[i]
                  if ma.el[kk] != "H" and not (t == zero and kk in own_a) and not (t == shift and kk in own_b)}
        for j in own_b:
            for kk, t in ma.links[j]:
                x = (kk, (t[0] + shift[0], t[1] + shift[1], t[2] + shift[2]))
                if x in near_a:
                    return "bridged"
        return ""
    pa, pb = ma.pos[ca.atoms], mb.pos[cb.atoms] + v
    ra = np.array([covalent_radius(ma.el[i]) for i in ca.atoms])
    rb = np.array([covalent_radius(mb.el[i]) for i in cb.atoms])
    own_a, own_b = set(int(i) for i in ca.atoms), set(int(i) for i in cb.atoms)
    cand = [(ma.pos[i], ma.el[i]) for i in range(len(ma)) if ma.el[i] != "H" and i not in own_a and not (a == b and i in own_b)]
    if a != b or v.any():
        cand += [(mb.pos[i] + v, mb.el[i]) for i in range(len(mb)) if mb.el[i] != "H" and i not in own_b]
    for x, el in cand:
        rx = covalent_radius(el)
        if (np.linalg.norm(pa - x, axis=1) < ra + rx + limits.bond_tolerance).any() and \
                (np.linalg.norm(pb - x, axis=1) < rb + rx + limits.bond_tolerance).any():
            return "bridged"
    return ""


def fragments_bonded(s: Structure, ca: Core, cb: Core, a: int, b: int, v: np.ndarray, limits: Limits) -> bool:
    """Vero se i frammenti `ca` (molecola a) e `cb` (molecola b traslata di v) condividono un atomo o hanno un legame
    tra loro, con la stessa regola dei legami di `structure` (somma dei raggi covalenti + `bond_tolerance`). Serve solo
    per frammenti della stessa molecola o di reti estese: molecole diverse non sono legate per costruzione."""
    ma, mb = s.molecules[a], s.molecules[b]
    d = cdist(ma.pos[ca.atoms], mb.pos[cb.atoms] + v)
    if (d < SAME_ATOM).any():
        return True
    ra = np.array([covalent_radius(ma.el[i]) for i in ca.atoms])
    rb = np.array([covalent_radius(mb.el[i]) for i in cb.atoms])
    return bool((d < ra[:, None] + rb[None, :] + limits.bond_tolerance).any())


def _heavy_contact(s: Structure, a: int, b: int, v: np.ndarray) -> float:
    A, B = s.molecules[a], s.molecules[b]
    ha, hb = A.heavy, B.heavy
    return float(cdist(A.pos[ha], B.pos[hb] + v).min()) if ha.any() and hb.any() else float("inf")


def _fragment_distance(s: Structure, ca: Core, cb: Core, a: int, b: int, v: np.ndarray) -> float:
    return float(cdist(s.molecules[a].pos[ca.atoms], s.molecules[b].pos[cb.atoms] + v).min())


def fragment_polygons(core: Core, q: np.ndarray, hull_only: bool = False) -> list[np.ndarray]:
    """I poligoni proiettati (coordinate 2D `q` degli atomi della molecola) che descrivono l'ingombro del frammento: uno
    per anello, ciascuno convesso; la loro unione segue la forma del frammento (le baie di fenantrene e crisene, la
    cavità di un macrociclo restano fuori). Senza anelli (frammento «molecola») o con `hull_only` (modo del sito):
    l'involucro convesso di tutti gli atomi del frammento."""
    if hull_only or not core.rings:
        return [convex_hull(q[core.atoms])]
    small = [convex_hull(q[list(r)]) for r in core.rings if len(r) <= 8]
    if _self_overlapping(small):
        # anelli dello stesso frammento che si sovrappongono in proiezione (reti dense di anelli a 4, revisione 2): la
        # somma delle intersezioni a due a due non è più l'area dell'unione; si usa l'involucro convesso del frammento
        return [convex_hull(q[core.atoms])]
    out = list(small)
    for r in core.rings:
        if len(r) <= 8:
            continue
        # macrociclo (anello interno di porfirine e ftalocianine, `fragment_merge = "conjugated"`): non è convesso e il
        # suo involucro coprirebbe parte degli anelli pirrolici; il poligono vero, nell'ordine dell'anello, diviso in
        # triangoli dal centroide (è stellato rispetto al centroide). Si aggiungono solo i triangoli che non coprono
        # area già presente: due cammini equivalenti dello stesso macrociclo non contano due volte
        pts = q[list(r)]
        c = pts.mean(axis=0)
        for i in range(len(pts)):
            tri = convex_hull(np.array([c, pts[i], pts[(i + 1) % len(pts)]]))  # antiorario
            if polygon_area(tri) > 1e-9 and not any(polygon_area(clip_convex(tri, o)) > 1e-6 for o in out):
                out.append(tri)
    return out


def _self_overlapping(polys: list[np.ndarray], rel: float = 0.02) -> bool:
    """Vero se due poligoni dello stesso frammento hanno in comune più del 2 % dell'area del minore: non solo un lato
    (anelli condensati, che in proiezione da coordinate non perfettamente planari si sovrappongono per una scheggia)."""
    boxes = [(p[:, 0].min(), p[:, 0].max(), p[:, 1].min(), p[:, 1].max()) for p in polys]
    areas = [polygon_area(p) for p in polys]
    for i, j in itertools.combinations(range(len(polys)), 2):
        a, b = boxes[i], boxes[j]
        if a[1] < b[0] or b[1] < a[0] or a[3] < b[2] or b[3] < a[2]:
            continue
        if polygon_area(clip_convex(polys[i], polys[j])) > rel * min(areas[i], areas[j]):
            return True
    return False


def _union_area(polys: list[np.ndarray]) -> float:
    """Area dell'unione dei poligoni degli anelli di un frammento: gli anelli di un sistema planare si toccano solo sui
    lati, quindi è la somma delle aree."""
    return float(sum(polygon_area(pg) for pg in polys))


def polygons_overlap(pa: list[np.ndarray], pb: list[np.ndarray]) -> float:
    """Area dell'intersezione delle due unioni di poligoni convessi (anelli che si toccano solo sui lati: somma delle
    intersezioni a due a due, ciascuna calcolata esattamente con il ritaglio di Sutherland–Hodgman)."""
    return float(sum(polygon_area(clip_convex(a, b)) for a in pa for b in pb if len(a) >= 3 and len(b) >= 3))


def _inside_any(polys: list[np.ndarray], pts: np.ndarray) -> np.ndarray:
    mask = np.zeros(len(pts), dtype=bool)
    for pg in polys:
        mask |= points_inside(pg, pts)
    return mask


def axis_twist(ca: Core, cb: Core, n: np.ndarray) -> float:
    """Angolo (0–90°) tra gli assi lunghi dei due frammenti proiettati sul piano medio di normale `n`."""
    la = ca.long - n * np.dot(ca.long, n)
    lb = cb.long - n * np.dot(cb.long, n)
    c = abs(float(np.dot(la, lb))) / (float(np.linalg.norm(la) * np.linalg.norm(lb)) or 1.0)
    return float(np.degrees(np.arccos(min(1.0, c))))


def axis_defined(core: Core, limits: Limits) -> bool:
    """L'asse lungo è definito se è dato da due atomi o se il frammento non è quasi isotropo nel piano (un asse Cn con
    n ≥ 3 normale al piano rende isotropa la dispersione: benzene, trifenilene, coronene, porfirine)."""
    return core.axis_source in ("atoms", *FALLBACK_SOURCES) or core.plane.aspect <= limits.max_axis_isotropy


def own_slips(ca: Core, cb: Core) -> tuple[float, float]:
    """Slittamento longitudinale e laterale del centroide di `cb` nel sistema di `ca`: d = c_B − c_A proiettato sul
    piano di A, poi |s·l_A| e |s·(n_A × l_A)|."""
    d = cb.centroid - ca.centroid
    n = ca.normal
    sv = d - n * float(np.dot(d, n))
    long = unit(ca.long - n * float(np.dot(ca.long, n)))
    return abs(float(np.dot(sv, long))), abs(float(np.dot(sv, np.cross(n, long))))


def _measure(s: Structure, ca: Core, cb0: Core, a: int, b: int, shift: Shift, limits: Limits, contact: float,
             fa: int, fb: int, fragmented: bool, bonded: bool = False, bridged: bool = False) -> Pair:
    v = s.shift_vector(shift)
    pa, pb = s.molecules[a].pos, s.molecules[b].pos + v
    cb = cb0.moved(v)
    f = frame_for(ca, cb, limits)
    d = cb.centroid - ca.centroid
    z = float(np.dot(d, f.n))
    sv = d - f.n * z
    # atan2(|a×b|, |a·b|) invece di arccos(|a·b|): stessa grandezza, ma precisa anche vicino a 0° (arccos lì perde
    # metà delle cifre: 1e-6° su una macchina, 0 su un'altra)
    dihedral = float(np.degrees(np.arctan2(np.linalg.norm(np.cross(ca.normal, cb.normal)), abs(np.dot(ca.normal, cb.normal)))))
    site = limits.pair_mode == "site"
    qa = f.project(pa)[:, :2]
    qb = f.project(pb)[:, :2]
    polys_a = fragment_polygons(ca, qa, hull_only=site)
    polys_b = fragment_polygons(cb, qb, hull_only=site)
    overlap = float(max(_inside_any(polys_a, qb[cb.atoms]).mean(), _inside_any(polys_b, qa[ca.atoms]).mean()))
    smaller = min(_union_area(polys_a), _union_area(polys_b))
    area = float(polygons_overlap(polys_a, polys_b) / smaller) if smaller > 0 else 0.0
    interplanar = abs(z)
    slip = float(np.linalg.norm(sv))
    long_slip, lat_slip = abs(float(np.dot(sv, f.L))), abs(float(np.dot(sv, f.W)))
    twist = axis_twist(ca, cb, f.n)
    def_a, def_b = axis_defined(ca, limits), axis_defined(cb, limits)
    defined = def_a and def_b
    if not site and (not defined or twist > limits.max_axis_twist):
        long_slip = lat_slip = math.nan  # asse lungo non definito o assi incrociati: la divisione di s in L e W è arbitraria
    own_a = own_slips(ca, cb) if def_a and not site else (math.nan, math.nan)
    own_b = own_slips(cb, ca) if def_b and not site else (math.nan, math.nan)
    nb = cb.normal if np.dot(ca.normal, cb.normal) >= 0 else -cb.normal
    # la retta d'intersezione passa tra i centroidi; altezze sotto 10⁻³ Å (frammenti complanari) contano come zero
    planes_cross = float(np.dot(d, ca.normal)) * float(np.dot(d, nb)) < -1e-6
    p = Pair(
        a=a, b=b, shift=shift, kind="", interplanar=interplanar,
        long_slip=long_slip, lat_slip=lat_slip, slip=slip,
        dihedral=dihedral, centroid_distance=float(np.linalg.norm(d)), contact=contact, overlap=overlap,
        fa=fa, fb=fb, fragmented=fragmented, overlap_area=area, bonded=bonded, bridged=bridged,
        twist=twist if defined else math.nan, planes_cross=planes_cross,
        long_slip_a=own_a[0], lat_slip_a=own_a[1], long_slip_b=own_b[0], lat_slip_b=own_b[1],
        axis_a=ca.axis_source if def_a else "", axis_b=cb.axis_source if def_b else "",
        same_network=same_network(s, a, b, shift),
    )
    p.kind, p.clamped = reclassify(p, limits), is_clamped(p, limits)
    return p


def measure_pair(s: Structure, cores: Sequence, a: int, b: int, shift: Shift, limits: Limits,
                 max_contact: float | None = None, fa: int = 0, fb: int = 0) -> Pair | None:
    """Geometria della coppia (frammento fa di a nella cella 0, frammento fb di b traslata di `shift`). None se manca
    un frammento o se il contatto più corto tra atomi pesanti supera `max_contact`. `cores` può essere l'elenco dei
    nuclei principali o quello dei frammenti di ogni molecola."""
    frags = as_fragments(cores)
    if fa >= len(frags[a]) or fb >= len(frags[b]):
        return None
    ca, cb0 = frags[a][fa], frags[b][fb]
    v = s.shift_vector(shift)
    intra = a == b and shift == (0, 0, 0)
    if intra:
        contact = _fragment_distance(s, ca, cb0, a, b, v)
    else:
        contact = _heavy_contact(s, a, b, v)
    if max_contact is not None and contact > max_contact:
        return None
    link = fragments_link(s, ca, cb0, a, b, v, limits) if (intra or s.molecules[a].extended or s.molecules[b].extended) else ""
    return _measure(s, ca, cb0, a, b, shift, limits, contact, fa, fb, len(frags[a]) > 1 or len(frags[b]) > 1,
                    link == "bonded", link == "bridged")


def remeasure(s: Structure, frags: Fragments, pair: Pair, limits: Limits) -> Pair:
    """La stessa coppia (stesse molecole, frammenti e traslazione) misurata su un'altra geometria: serve al Monte Carlo
    delle incertezze, dove gli atomi si spostano di poco ma la topologia resta quella dell'analisi."""
    ca, cb0 = frags[pair.a][pair.fa], frags[pair.b][pair.fb]
    v = s.shift_vector(pair.shift)
    by_fragment = pair.intramolecular or pair.polymer
    contact = _fragment_distance(s, ca, cb0, pair.a, pair.b, v) if by_fragment else _heavy_contact(s, pair.a, pair.b, v)
    q = _measure(s, ca, cb0, pair.a, pair.b, pair.shift, limits, contact, pair.fa, pair.fb, pair.fragmented, pair.bonded,
                 pair.bridged)
    q.index, q.count, q.members, q.polymer = pair.index, pair.count, pair.members, pair.polymer
    return q


def shells_needed(s: Structure, cores: Sequence, limits: Limits) -> tuple[int, int, int]:
    """Celle da esplorare (±) lungo a, b, c perché due molecole qualunque a contatto siano trovate, ma almeno
    `neighbour_shells` e senza il limite `MAX_SHELLS`. Nel modo del sito sempre `neighbour_shells`.

    Due molecole si toccano solo se |Δbaricentri| ≤ D = E_a + E_b + contact_cutoff (E = estensione); la
    componente frazionaria i di un vettore g è al più |g|·|a*_i| e i baricentri stanno nella cella [0, 1), quindi
    basta esplorare |k_i| ≤ D·|a*_i| + 1."""
    R = limits.neighbour_shells
    if limits.pair_mode == "site":
        return (R, R, R)
    frags = as_fragments(cores)
    ext = [m.extent for m, f in zip(s.molecules, frags) if f]  # per una rete estesa: l'unità ripetuta
    if not ext:
        return (R, R, R)
    D = 2 * max(ext) + limits.contact_cutoff
    recip = np.linalg.norm(np.linalg.inv(s.vectors), axis=0)  # |a*_i|: colonne di M⁻¹
    need = np.floor(D * recip + 1).astype(int)
    return tuple(int(max(R, n)) for n in need)  # type: ignore[return-value]


def neighbour_shells(s: Structure, cores: Sequence, limits: Limits) -> tuple[int, int, int]:
    """Le celle esplorate davvero: `shells_needed` limitato a `MAX_SHELLS` (se il limite taglia, i controlli rapidi
    lo segnalano con l'avviso `shellsTruncated`)."""
    need = shells_needed(s, cores, limits)
    if limits.pair_mode == "site":
        return need
    return tuple(min(MAX_SHELLS, n) for n in need)  # type: ignore[return-value]


def _neg(sh: Shift) -> Shift:
    return (-sh[0], -sh[1], -sh[2])


def _canon(a: int, b: int, sh: Shift) -> tuple[int, int, Shift]:
    """Lo stesso contatto tra molecole scritto in un solo modo: (a, b, s) e (b, a, −s) sono la stessa coppia."""
    return min((a, b, sh), (b, a, _neg(sh)))


def _canon5(a: int, b: int, sh: Shift, fa: int, fb: int) -> Member:
    """Come `_canon`, con i frammenti: (a, b, s, fa, fb) e (b, a, −s, fb, fa) sono lo stesso contatto."""
    return min((a, b, sh, fa, fb), (b, a, _neg(sh), fb, fa))


def iter_neighbour_pairs(s: Structure, cores: Sequence, limits: Limits) -> Iterator[Pair]:
    """Tutte le coppie di frammenti vicini, una per traslazione, senza raccogliere quelle equivalenti. Le molecole
    devono toccarsi (atomi pesanti entro `contact_cutoff`); oltre alla coppia dei frammenti principali (0–0) si
    misurano quelle di frammenti entro `contact_cutoff + FRAGMENT_MARGIN`. In modalità site ogni contatto compare nei
    due versi, altrimenti in uno.

    Reti estese (polimeri): la molecola è l'unità ripetuta, srotolata, con gli anelli interi; ogni anello del cristallo
    è un suo frammento più una traslazione. Due unità della stessa catena sono legate, quindi il contatto tra molecole
    non dice nulla: si misurano le coppie di frammenti entro `contact_cutoff + FRAGMENT_MARGIN` (come i frammenti
    secondari delle molecole), e il contatto è quello tra i due frammenti."""
    frags = as_fragments(cores)
    na, nb, nc = neighbour_shells(s, frags, limits)
    shifts = np.array(list(itertools.product(range(-na, na + 1), range(-nb, nb + 1), range(-nc, nc + 1))), dtype=int)
    shift_vec = shifts @ s.vectors
    extent = [m.extent for m in s.molecules]
    usable = [i for i, m in enumerate(s.molecules) if frags[i]]
    one_way = limits.pair_mode != "site"
    reach = limits.contact_cutoff + FRAGMENT_MARGIN
    for a in usable:
        for b in usable:
            if one_way and b < a:
                continue
            # scarto veloce: baricentri troppo lontani perché le molecole si tocchino
            gap = np.linalg.norm(s.molecules[b].centroid + shift_vec - s.molecules[a].centroid, axis=1)
            near = gap <= extent[a] + extent[b] + limits.contact_cutoff
            fragmented = len(frags[a]) > 1 or len(frags[b]) > 1
            for k in np.flatnonzero(near):
                sh = tuple(int(x) for x in shifts[k])
                if a == b and sh == (0, 0, 0):
                    continue
                if one_way and _canon(a, b, sh) != (a, b, sh):
                    continue
                v = s.shift_vector(sh)
                if s.molecules[a].extended or s.molecules[b].extended:
                    for fa, fb in itertools.product(range(len(frags[a])), range(len(frags[b]))):
                        dist = _fragment_distance(s, frags[a][fa], frags[b][fb], a, b, v)
                        if dist <= reach:  # come i frammenti secondari delle molecole: contact_cutoff + FRAGMENT_MARGIN
                            link = fragments_link(s, frags[a][fa], frags[b][fb], a, b, v, limits)
                            p = _measure(s, frags[a][fa], frags[b][fb], a, b, sh, limits, dist, fa, fb, fragmented,
                                         link == "bonded", link == "bridged")
                            p.polymer = True
                            yield p
                    continue
                contact = _heavy_contact(s, a, b, v)
                if contact > limits.contact_cutoff:
                    continue
                for fa, fb in itertools.product(range(len(frags[a])), range(len(frags[b]))):
                    if (fa or fb) and _fragment_distance(s, frags[a][fa], frags[b][fb], a, b, v) > reach:
                        continue
                    yield _measure(s, frags[a][fa], frags[b][fb], a, b, sh, limits, contact, fa, fb, fragmented)
    # coppie di frammenti della stessa molecola (molecole con più sistemi π, p. es. ripiegate su se stesse)
    zero = np.zeros(3)
    for a in usable:
        for fa, fb in itertools.permutations(range(len(frags[a])), 2):
            if one_way and fb < fa:
                continue
            dist = _fragment_distance(s, frags[a][fa], frags[a][fb], a, a, zero)
            if dist <= reach:
                link = fragments_link(s, frags[a][fa], frags[a][fb], a, a, zero, limits)
                p = _measure(s, frags[a][fa], frags[a][fb], a, a, (0, 0, 0), limits, dist, fa, fb, True,
                             link == "bonded", link == "bridged")
                p.polymer = s.molecules[a].extended
                yield p


#: Un arco di stack (a, b, shift, coppia, fa, fb): il frammento fb della molecola b, traslata di `shift`, è impilato
#: sul frammento fa di a. I frammenti servono alle reti estese, dove i nodi del grafo sono i frammenti (`stacks`).
StackEdge = tuple[int, int, Shift, "Pair", int, int]


def collect_pairs(s: Structure, cores: Sequence, limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
    """Coppie di frammenti vicini, con quelle equivalenti raccolte in una sola, e, a parte, tutti gli archi `stack`
    tra molecole diverse (nei due versi: servono a ricostruire le colonne)."""
    frags = as_fragments(cores)
    contacts = list(iter_neighbour_pairs(s, frags, limits))
    return group_pairs(contacts, contact_keys(s, frags, contacts, limits), limits)


def contact_keys(s: Structure, frags: Fragments, contacts: list[Pair], limits: Limits) -> list[Member] | None:
    """La chiave di simmetria di ogni coppia misurata (modalità canonical); None nel modo del sito, dove le coppie
    equivalenti si riconoscono dai valori e quindi dipendono dalla classificazione."""
    if limits.pair_mode == "site":
        return None
    # sempre, anche con un solo frammento per molecola: un'operazione può portare la molecola in sé stessa o in una sua
    # copia senza portare il frammento in un frammento (due anelli equivalenti di cui uno solo è un frammento, COD
    # 7044129; atomi disordinati tolti in modo non simmetrico, COD 4326386), e allora non va usata (voce 86)
    fmaps = fragment_maps(s, frags)
    fimgs = fragment_images(s, frags) if any(p.polymer for p in contacts) else None
    orbits = [symmetry_orbit(s, p.a, p.b, p.shift, p.fa, p.fb, fmaps, fimgs) for p in contacts]
    keys = [min(o) for o in orbits]
    _complete_classes(s, frags, contacts, keys, orbits, limits)
    return keys


def _complete_classes(s: Structure, frags: Fragments, contacts: list[Pair], keys: list[Member],
                      orbits: list[set[Member]], limits: Limits) -> None:
    """Aggiunge a `contacts` (e a `keys`) i membri di ogni classe di simmetria che la ricerca non ha misurato.

    La coppia dei frammenti principali (0–0) di due molecole a contatto si misura sempre; le altre coppie di frammenti
    solo se i frammenti sono entro `contact_cutoff + FRAGMENT_MARGIN`. Se un'operazione porta il frammento 0 di una
    molecola nel frammento 1 della sua immagine (frammenti della stessa taglia), l'immagine di una coppia 0–0 lontana è
    una coppia 1–1 che la ricerca scartava: la classe aveva meno membri del dovuto (molteplicità sottostimata) e il suo
    rappresentante cambiava con `contact_cutoff` (revisione 4, S5). Ogni classe ha ora tutta la sua orbita e il
    rappresentante è sempre il membro più piccolo (`symmetry_key`), qualunque siano le soglie."""
    have = {p.key for p in contacts}
    by_class: dict[Member, tuple[Pair, set[Member]]] = {}
    for p, key, orbit in zip(list(contacts), list(keys), orbits):
        ref, members = by_class.setdefault(key, (p, set()))
        members |= orbit
    for key, (ref, members) in by_class.items():
        for a, b, sh, fa, fb in sorted(members - have):
            if fa >= len(frags[a]) or fb >= len(frags[b]):
                continue
            q = _measure(s, frags[a][fa], frags[b][fb], a, b, sh, limits, ref.contact, fa, fb,
                         len(frags[a]) > 1 or len(frags[b]) > 1, ref.bonded, ref.bridged)
            q.polymer = ref.polymer
            contacts.append(q)
            keys.append(key)
            have.add(q.key)


def group_pairs(contacts: list[Pair], keys: list[Member] | None, limits: Limits) -> tuple[list[Pair], list[StackEdge]]:
    """Classifica con `limits` copie fresche delle coppie misurate, raccoglie quelle equivalenti (per chiave di
    simmetria, o per valori arrotondati se `keys` è None) e le ordina. Le coppie misurate non vengono toccate: si
    può ripetere con altre soglie dello stadio `pairs`."""
    fresh = [replace(p, kind=reclassify(p, limits), clamped=is_clamped(p, limits), count=1, members=[], values={},
                     index=-1) for p in contacts]
    pairs, stack_edges = _group_site(fresh) if keys is None else _group_canonical(fresh, keys)
    order = {k: i for i, k in enumerate(KINDS)}
    # distanze arrotondate a 1e-6 Å e spareggio su grandezze e indici: due coppie con la stessa distanza (2.3655 contro
    # 2.3655 − 1 ulp) non si scambiano di posto da una macchina all'altra (REGISTRO, voce 66)
    pairs.sort(key=lambda p: (order[p.kind], round(p.interplanar, 6), round(p.contact, 6),
                              round(p.centroid_distance, 6), p.a, p.b, p.fa, p.fb, tuple(p.shift)))
    for i, p in enumerate(pairs):
        p.index = i
    return pairs, stack_edges


def _group_site(contacts: list[Pair]) -> tuple[list[Pair], list[StackEdge]]:
    found: dict[str, Pair] = {}
    stack_edges: list[StackEdge] = []
    for p in contacts:
        if p.kind == "stack" and not p.intramolecular:
            stack_edges.append((p.a, p.b, p.shift, p, p.fa, p.fb))
        key = "|".join([p.kind, f"{p.interplanar:.2f}", f"{p.long_slip:.2f}", f"{p.lat_slip:.2f}",
                        f"{p.dihedral:.1f}", f"{p.contact:.2f}"])
        prev = found.get(key)
        if prev is not None:
            prev.count += 1
            prev.members.append(p.key)
        else:
            p.members = [p.key]
            found[key] = p
    return list(found.values()), stack_edges


def fragment_maps(s: Structure, frags: Fragments) -> list[list[list[int | None] | None]]:
    """Per ogni operazione e molecola: in quale frammento della molecola immagine va ciascun frammento (None se non
    corrisponde a un frammento)."""
    sets = [[frozenset(int(i) for i in f.atoms) for f in fs] for fs in frags]
    out = []
    for row_img, row_map in zip(s.molecule_images(), s.atom_maps()):
        row = []
        for m, (img, amap) in enumerate(zip(row_img, row_map)):
            if img is None or amap is None:
                row.append(None)
                continue
            target = {fs: j for j, fs in enumerate(sets[img[0]])}
            row.append([target.get(frozenset(int(amap[i]) for i in fs)) for fs in sets[m]])
        out.append(row)
    return out


FragmentImage = tuple[int, int, Shift]  # (molecola, frammento, traslazione in celle)


def fragment_images(s: Structure, frags: Fragments, tol: float = 0.3) -> list[list[list[FragmentImage | None]]]:
    """Per ogni operazione, molecola e frammento: in quale frammento di quale molecola va, e con quale traslazione.

    Serve alle reti estese: una catena infinita non ha un'immagine come molecola intera (l'unità ripetuta è un taglio
    arbitrario della catena), ma ogni frammento π è locale e va in un frammento più una traslazione."""
    owner: dict[int, tuple[int, int]] = {}  # atomo della cella → (molecola, posizione nella molecola)
    for mi, m in enumerate(s.molecules):
        for p, ci in enumerate(m.cell_index):
            owner[int(ci)] = (mi, p)
    # un atomo può stare in più frammenti (un metallo in due anelli chelati): si provano tutti, non solo l'ultimo; con
    # uno solo l'immagine si cercava a volte nel frammento sbagliato e l'operazione si scartava, secondo il taglio
    # dell'unità ripetuta (revisione 5, N1)
    in_frag: list[dict[int, list[int]]] = []
    for fs in frags:
        d: dict[int, list[int]] = {}
        for f, c in enumerate(fs):
            for i in c.atoms:
                d.setdefault(int(i), []).append(f)
        in_frag.append(d)
    by_el = {e: np.array([i for i, x in enumerate(s.cell_el) if x == e]) for e in set(s.cell_el)}
    out = []
    for op in s.symops:
        row = []
        for m, fs in zip(s.molecules, frags):
            imgs: list[FragmentImage | None] = []
            for c in fs:
                idx = np.asarray(c.atoms)
                src = s.cell_frac[m.cell_index[idx]] + m.cell_shift[idx]
                img = src @ op.r.T + op.t
                cand = by_el[m.el[int(idx[0])]]
                d = img[0] - s.cell_frac[cand]
                d -= np.round(d)
                k = int(np.argmin(np.linalg.norm(d @ s.vectors, axis=1)))
                j = int(cand[k])
                hit = owner.get(j)
                found: FragmentImage | None = None
                if hit is not None:
                    m2 = s.molecules[hit[0]]
                    t = np.rint(img[0] - s.cell_frac[j] - m2.cell_shift[hit[1]]).astype(int)
                    for f2 in in_frag[hit[0]].get(hit[1], []):
                        t2 = np.asarray(frags[hit[0]][f2].atoms)
                        if len(t2) != len(idx):
                            continue
                        tgt = s.cell_frac[m2.cell_index[t2]] + m2.cell_shift[t2] + t
                        dist = np.linalg.norm((img[:, None, :] - tgt[None, :, :]) @ s.vectors, axis=2)
                        if (dist.min(axis=1) < tol).all():
                            found = (hit[0], f2, (int(t[0]), int(t[1]), int(t[2])))
                            break
                imgs.append(found)
            row.append(imgs)
        out.append(row)
    return out


def symmetry_key(s: Structure, a: int, b: int, sh: Shift, fa: int = 0, fb: int = 0,
                 fmaps: list[list[list[int | None] | None]] | None = None,
                 fimgs: list[list[list[FragmentImage | None]]] | None = None) -> Member:
    """Rappresentante dell'orbita del contatto sotto il gruppo spaziale: la più piccola delle sue immagini (scritte
    con `_canon5`). Due contatti sono equivalenti per simmetria se e solo se hanno la stessa chiave; le operazioni che
    non mappano molecole e frammenti non si usano, così non si uniscono mai contatti diversi. Con `fimgs` (reti
    estese) le immagini si prendono frammento per frammento (`fragment_images`)."""
    return min(symmetry_orbit(s, a, b, sh, fa, fb, fmaps, fimgs))


def symmetry_orbit(s: Structure, a: int, b: int, sh: Shift, fa: int = 0, fb: int = 0,
                   fmaps: list[list[list[int | None] | None]] | None = None,
                   fimgs: list[list[list[FragmentImage | None]]] | None = None) -> set[Member]:
    """Tutte le immagini del contatto sotto le operazioni del gruppo spaziale (scritte con `_canon5`), lui compreso."""
    out = {_canon5(a, b, sh, fa, fb)}
    if fimgs is not None and (s.molecules[a].extended or s.molecules[b].extended):
        for k, row in enumerate(fimgs):
            ia, ib = row[a][fa], row[b][fb]
            if ia is None or ib is None:
                continue
            (a2, fa2, ta), (b2, fb2, tb) = ia, ib
            rs = s.rotation(k) @ np.array(sh)
            new = (tb[0] + int(rs[0]) - ta[0], tb[1] + int(rs[1]) - ta[1], tb[2] + int(rs[2]) - ta[2])
            if a2 == b2 and new == (0, 0, 0) and fa2 == fb2:
                continue
            out.add(_canon5(a2, b2, new, fa2, fb2))
        return out
    for k, row in enumerate(s.molecule_images()):
        ia, ib = row[a], row[b]
        if ia is None or ib is None:
            continue
        fa2, fb2 = fa, fb
        if fmaps is not None:
            ma, mb = fmaps[k][a], fmaps[k][b]
            if ma is None or mb is None or ma[fa] is None or mb[fb] is None:
                continue
            fa2, fb2 = ma[fa], mb[fb]
        (a2, ta), (b2, tb) = ia, ib
        rs = s.rotation(k) @ np.array(sh)
        new = (tb[0] + int(rs[0]) - ta[0], tb[1] + int(rs[1]) - ta[1], tb[2] + int(rs[2]) - ta[2])
        if a2 == b2 and new == (0, 0, 0) and fa2 == fb2:
            continue
        out.add(_canon5(a2, b2, new, fa2, fb2))
    return out


def _group_canonical(contacts: list[Pair], keys: list[Member]) -> tuple[list[Pair], list[StackEdge]]:
    classes: dict[tuple, list[Pair]] = {}
    stack_edges: list[StackEdge] = []
    for p, key in zip(contacts, keys, strict=True):
        # tra molecole discrete un arco unisce due molecole; in una rete estesa i nodi sono i frammenti, e anche due
        # frammenti della stessa unità ripetuta impilati tra loro sono un arco (revisione 5, N1)
        if p.kind == "stack" and (not p.intramolecular or p.polymer):
            stack_edges += [(p.a, p.b, p.shift, p, p.fa, p.fb), (p.b, p.a, _neg(p.shift), p, p.fb, p.fa)]
        classes.setdefault(key, []).append(p)
    pairs = []
    for key, members in classes.items():
        members.sort(key=lambda p: p.key)
        rep = next((p for p in members if p.key == key), members[0])
        rep.count = len(members)
        rep.members = [rep.key] + [p.key for p in members if p is not rep]
        pairs.append(rep)
    return pairs, stack_edges


def molecule_key(s: Structure, pair: Pair, limits: Limits) -> tuple:
    """Chiave del contatto tra le due molecole della coppia (senza i frammenti): le coppie di frammenti con la stessa
    chiave appartengono allo stesso contatto tra molecole, a meno della simmetria in modalità canonical."""
    if pair.intramolecular:
        if limits.pair_mode == "site":
            return ("intra", pair.a)
        images = [row[pair.a][0] for row in s.molecule_images() if row[pair.a] is not None]
        return ("intra", min([pair.a, *images]))
    if limits.pair_mode == "site":
        return ("inter", *_canon(pair.a, pair.b, pair.shift))
    return ("inter", *symmetry_key(s, pair.a, pair.b, pair.shift)[:3])


def class_spread(s: Structure, cores: Sequence, limits: Limits, pair: Pair) -> float:
    """Differenza massima di distanza interplanare, slittamenti e contatto tra i membri della classe di `pair` (Å):
    ~0 se l'equivalenza per simmetria è esatta. Serve ai controlli."""
    vals = []
    for a, b, sh, *f in pair.members:
        fa, fb = f if f else (0, 0)
        q = measure_pair(s, cores, a, b, sh, limits, fa=fa, fb=fb)
        if q is not None:
            vals.append((q.interplanar, q.long_slip, q.lat_slip, q.contact))
    if len(vals) < 2:
        return 0.0
    arr = np.array(vals)
    nan = np.isnan(arr)
    if (nan.any(axis=0) & ~nan.all(axis=0)).any():  # NaN in un membro e non in un altro: classe non equivalente
        return float("nan")
    arr = arr[:, ~nan.all(axis=0)]  # slittamenti non definiti (asse isotropo) per tutti i membri: non contano
    return float(np.max(arr.max(axis=0) - arr.min(axis=0))) if math.isfinite(arr.sum()) else float("nan")


def neighbour_pairs(s: Structure, cores: Sequence, limits: Limits) -> list[Pair]:
    return collect_pairs(s, cores, limits)[0]
