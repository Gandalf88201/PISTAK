"""Descrizione geometrica della scena 3D di una coppia, indipendente dalla libreria di rendering.

`build_scene` produce solo array numpy e piccole strutture: la vista 3D (`views/view3d.py`) le disegna con
PyVista, ma si potrebbe scrivere un altro disegnatore (VTK puro, matplotlib, esportazione in un file) senza
toccare il calcolo. `build_measure_scene` mostra invece una misura della scheda Geometria (un contatto, un legame,
un angolo…) con le molecole in cui stanno i suoi atomi; `pick_measure` misura gli atomi scelti con il mouse.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial.distance import cdist

from pistak.core.analysis import PairContext
from pistak.core.elements import CYAN_CARBON, color, covalent_radius, vdw_radius
from pistak.core.geometry import convex_hull
from pistak.core.intramolecular import GeometryItem, angle, symmetry_codes
from pistak.core.structure import Structure

RGB = tuple[int, int, int]

COLOR_L: RGB = (214, 69, 65)  # slittamento longitudinale
COLOR_W: RGB = (46, 160, 67)  # slittamento laterale
COLOR_N: RGB = (47, 111, 214)  # distanza interplanare
COLOR_PLANE_A: RGB = (240, 165, 60)
COLOR_PLANE_B: RGB = (64, 190, 220)
COLOR_MEASURE: RGB = (225, 150, 0)  # misure: della scheda Geometria e scelte con il mouse


@dataclass
class Polygon:
    points: np.ndarray  # K×3, in ordine lungo il contorno
    color: RGB
    opacity: float = 0.35


@dataclass
class Path:
    points: np.ndarray  # K×3
    color: RGB
    radius: float = 0.05
    arrow: bool = False  # punta di freccia sull'ultimo punto
    dashed: bool = False  # tratteggiata (le misure)


@dataclass
class Label:
    pos: np.ndarray
    text: str
    color: RGB
    size: int = 13


@dataclass
class Marker:
    pos: np.ndarray
    color: RGB
    radius: float
    label: str = ""


@dataclass
class SceneFrame:
    """Assi per inquadrare una scena senza coppia (come `PairFrame`: L, W, n ortonormali, origine)."""

    origin: np.ndarray
    L: np.ndarray
    W: np.ndarray
    n: np.ndarray

    def project(self, points: np.ndarray) -> np.ndarray:
        return (np.asarray(points) - self.origin) @ np.column_stack([self.L, self.W, self.n])

    @classmethod
    def of(cls, points: np.ndarray, along: np.ndarray | None = None) -> SceneFrame:
        """Assi principali dei punti; L lungo `along` se dato (per esempio la linea di un contatto)."""
        c = points.mean(axis=0)
        _, _, vt = np.linalg.svd(points - c) if len(points) >= 3 else (None, None, np.eye(3))
        long = vt[0] if along is None or np.linalg.norm(along) < 1e-9 else along / np.linalg.norm(along)
        n = vt[2] - (vt[2] @ long) * long
        if np.linalg.norm(n) < 1e-6:
            n = np.cross(long, vt[1])
        n /= np.linalg.norm(n)
        return cls(c, long, np.cross(n, long), n)


@dataclass
class Scene:
    pos: np.ndarray  # N×3 atomi
    el: list[str]
    molecule: np.ndarray  # N: 0 = prima molecola, 1 = seconda
    in_core: np.ndarray  # N bool
    colors: np.ndarray  # N×3 uint8
    radii: np.ndarray  # N (raggio covalente, per le sfere)
    bonds: np.ndarray  # M×2, indici globali
    polygons: list[Polygon] = field(default_factory=list)
    paths: list[Path] = field(default_factory=list)
    markers: list[Marker] = field(default_factory=list)
    axes: list[Path] = field(default_factory=list)  # assi L, W, n della coppia, dal baricentro del primo nucleo
    labels: list[Label] = field(default_factory=list)  # nomi degli assi e valori misurati lungo ciascuno
    ghost_pos: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))  # altre molecole della colonna
    ghost_bonds: np.ndarray = field(default_factory=lambda: np.zeros((0, 2), dtype=int))
    box: np.ndarray = field(default_factory=lambda: np.zeros((0, 2, 3)))  # E×2×3: spigoli della cella
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))
    size: float = 10.0
    #: Etichetta di ogni atomo, per le misure scelte con il mouse ("C3", "C3 (B)", "O4#2").
    names: list[str] = field(default_factory=list)
    #: Una misura da mostrare sempre: atomi evidenziati, linee tratteggiate ed etichette con il valore.
    highlight: list[int] = field(default_factory=list)
    measure_paths: list[Path] = field(default_factory=list)
    measure_labels: list[Label] = field(default_factory=list)
    #: Assi per l'inquadratura quando la scena non è di una coppia.
    frame: SceneFrame | None = None
    #: Il contatto più corto tra atomi pesanti delle due molecole (la colonna «Contatto»): (atomo A, atomo B, Å).
    contact: tuple[int, int, float] | None = None


def _plane_polygon(points: np.ndarray, centroid: np.ndarray, normal: np.ndarray, long: np.ndarray, margin: float = 0.5) -> np.ndarray:
    """Contorno (convesso, allargato di `margin`) degli atomi del nucleo, sul loro piano medio."""
    w = np.cross(normal, long)
    q = np.column_stack([(points - centroid) @ long, (points - centroid) @ w])
    hull = convex_hull(q)
    if len(hull) < 3:
        return np.zeros((0, 3))
    c2 = hull.mean(axis=0)
    d = hull - c2
    hull = c2 + d * (1 + margin / np.maximum(np.linalg.norm(d, axis=1, keepdims=True), 1e-6))
    return centroid + hull[:, :1] * long + hull[:, 1:] * w


def column_ghosts(ctx: PairContext, column, reps: int = 2, show_h: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Atomi e legami delle altre molecole della colonna di stack, `reps` ripetizioni per lato (rango 1)."""
    base = column.positions[ctx.pair.a]
    ks = range(-reps, reps + 1) if column.rank == 1 else [0]
    t = np.array(column.translations[0]) if column.rank == 1 else np.zeros(3, int)
    pts: list[np.ndarray] = []
    bonds: list[np.ndarray] = []
    n = 0
    for m in column.molecules:
        mol = ctx.structure.molecules[m]
        keep = np.ones(len(mol), bool) if show_h else mol.heavy
        idx = np.cumsum(keep) - 1
        for k in ks:
            shift = tuple(int(x) for x in np.array(column.positions[m]) - np.array(base) + k * t)
            if (m == ctx.pair.a and shift == (0, 0, 0)) or (m == ctx.pair.b and shift == tuple(ctx.pair.shift)):
                continue
            pts.append(mol.pos[keep] + ctx.structure.shift_vector(shift))
            b = [(idx[i] + n, idx[j] + n) for i, j in mol.bonds if keep[i] and keep[j]]
            if b:
                bonds.append(np.array(b))
            n += int(keep.sum())
    if not pts:
        return np.zeros((0, 3)), np.zeros((0, 2), dtype=int)
    return np.vstack(pts), (np.vstack(bonds) if bonds else np.zeros((0, 2), dtype=int))


def build_scene(ctx: PairContext, show_h: bool = True, show_cell: bool = False, show_column: bool = False) -> Scene:
    A, B = ctx.A, ctx.B
    keep_a = np.ones(len(A), bool) if show_h else A.heavy
    keep_b = np.ones(len(B), bool) if show_h else B.heavy
    pos = np.vstack([ctx.pos_a[keep_a], ctx.pos_b[keep_b]])
    el = [e for e, k in zip(A.el, keep_a) if k] + [e for e, k in zip(B.el, keep_b) if k]
    mol = np.r_[np.zeros(keep_a.sum(), int), np.ones(keep_b.sum(), int)]
    core = np.zeros(len(pos), bool)
    ia = np.cumsum(keep_a) - 1  # indice vecchio → nuovo
    ib = np.cumsum(keep_b) - 1 + keep_a.sum()
    core[ia[ctx.core_a.atoms[keep_a[ctx.core_a.atoms]]]] = True
    core[ib[ctx.core_b.atoms[keep_b[ctx.core_b.atoms]]]] = True
    colors = np.array([CYAN_CARBON if (m == 1 and e == "C") else color(e) for e, m in zip(el, mol)], dtype=np.uint8)
    radii = np.array([covalent_radius(e) for e in el])

    bonds = []
    for mol_obj, keep, idx in ((A, keep_a, ia), (B, keep_b, ib)):
        for i, j in mol_obj.bonds:
            if keep[i] and keep[j]:
                bonds.append((idx[i], idx[j]))
    bonds_arr = np.array(bonds, dtype=int).reshape(-1, 2)

    ca, cb = ctx.core_a, ctx.core_b
    scene = Scene(pos, el, mol, core, colors, radii, bonds_arr)
    scene.names = [lab for lab, k in zip(A.labels, keep_a) if k] + [f"{lab} (B)" for lab, k in zip(B.labels, keep_b) if k]
    for c, pts, col in ((ca, ctx.pos_a, COLOR_PLANE_A), (cb, ctx.pos_b, COLOR_PLANE_B)):
        poly = _plane_polygon(pts[c.atoms], c.centroid, c.normal, c.long)
        if len(poly):
            scene.polygons.append(Polygon(poly, col))
    # scomposizione di d: slittamento longitudinale, laterale e distanza interplanare, a gradini
    f = ctx.frame
    d = ctx.d
    dl, dw, dn = float(d @ f.L), float(d @ f.W), float(d @ f.n)
    p0 = ca.centroid
    p1 = p0 + f.L * dl
    p2 = p1 + f.W * dw
    p3 = p2 + f.n * dn
    scene.paths += [Path(np.array([p0, p1]), COLOR_L), Path(np.array([p1, p2]), COLOR_W), Path(np.array([p2, p3]), COLOR_N)]
    scene.paths.append(Path(np.array([p0, cb.centroid]), (90, 90, 90), 0.025))
    # le misure proiettate sugli assi: L, W, n passano per Cg A e ogni gradino porta il suo valore (quelli della tabella)
    reach = max(abs(dl), abs(dw), abs(dn), 2.0) + 1.5
    # (il nome di n è scostato lungo -W: guardando dall'alto la punta di n cade su Cg A)
    for axis, col, name, aside in ((f.L, COLOR_L, "L", 0.0), (f.W, COLOR_W, "W", 0.0), (f.n, COLOR_N, "n", -0.7)):
        scene.axes.append(Path(np.array([p0 - axis * reach, p0 + axis * reach]), col, 0.018, arrow=True))
        scene.labels.append(Label(p0 + axis * (reach + 0.45) + f.W * aside, name, col, 16))
    # etichette scostate dal proprio gradino, dalla parte opposta agli altri (n: oltre Cg B, che dall'alto la coprirebbe)
    sl, sw = (1.0 if dl >= 0 else -1.0), (1.0 if dw >= 0 else -1.0)
    for pos, value, col, name in (((p0 + p1) / 2 - f.W * sw * 0.5, dl, COLOR_L, "L"), ((p1 + p2) / 2 + f.L * sl * 0.3, dw, COLOR_W, "W"),
                                  ((p2 + p3) / 2 + f.W * sw * 0.7, dn, COLOR_N, "n")):
        scene.labels.append(Label(pos, f"{name} {abs(value):.2f} Å", col))
    scene.markers += [Marker(ca.centroid, COLOR_PLANE_A, 0.16, "Cg A"), Marker(cb.centroid, COLOR_PLANE_B, 0.16, "Cg B")]
    # dove sta il contatto più corto tra le due molecole (gli atomi pesanti, come `Pair.contact`)
    ha, hb = np.flatnonzero(A.heavy), np.flatnonzero(B.heavy)
    if not ctx.pair.intramolecular and len(ha) and len(hb):
        dist = cdist(ctx.pos_a[ha], ctx.pos_b[hb])
        i, j = np.unravel_index(int(np.argmin(dist)), dist.shape)
        scene.contact = (int(ia[ha[i]]), int(ib[hb[j]]), float(dist[i, j]))

    if show_column:
        col = ctx.analysis.column_for(ctx.pair)
        if col is not None:
            scene.ghost_pos, scene.ghost_bonds = column_ghosts(ctx, col)

    if show_cell:
        v = ctx.structure.vectors
        corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=float) @ v
        edges = [(a, b) for a in range(8) for b in range(a + 1, 8) if bin(a ^ b).count("1") == 1]
        scene.box = np.array([[corners[a], corners[b]] for a, b in edges])

    allpts = np.vstack([pos, ca.centroid[None], cb.centroid[None]])
    scene.center = allpts.mean(axis=0)
    scene.size = float(np.linalg.norm(allpts.max(axis=0) - allpts.min(axis=0))) or 10.0
    return scene


def build_measure_scene(structure: Structure, item: GeometryItem, show_h: bool = True, show_cell: bool = False) -> Scene:
    """La scena di una misura della scheda Geometria: le molecole (intere, nella loro posizione) che contengono i suoi
    atomi, gli atomi evidenziati, la misura tratteggiata con il valore. Il carbonio delle molecole dopo la prima è
    ciano, come la seconda molecola di una coppia."""
    places = list(item.places)
    instances = list(dict.fromkeys((m, tuple(t)) for m, t, _ in places if m >= 0))
    labels, _legend = symmetry_codes([item])
    pos_parts, el, mol, names, bonds = [], [], [], [], []
    where: dict[tuple, int] = {}  # (istanza, indice nella molecola) → indice globale
    n = 0
    for k, (m, t) in enumerate(instances):
        molecule = structure.molecules[m]
        keep = np.ones(len(molecule), bool) if show_h else molecule.heavy.copy()
        for mm, tt, i in places:  # gli atomi misurati restano anche senza idrogeni
            if (mm, tuple(tt)) == (m, t):
                keep[i] = True
        idx = np.cumsum(keep) - 1 + n
        pos_parts.append(molecule.pos[keep] + structure.shift_vector(t))
        el += [e for e, kk in zip(molecule.el, keep) if kk]
        mol += [min(k, 1)] * int(keep.sum())
        names += [lab if k == 0 else f"{lab} ({k + 1})" for lab, kk in zip(molecule.labels, keep) if kk]
        bonds += [(idx[i], idx[j]) for i, j in molecule.bonds if keep[i] and keep[j]]
        for i in np.flatnonzero(keep):
            where[(k, int(i))] = int(idx[i])
        n += int(keep.sum())
    pos = np.vstack(pos_parts) if pos_parts else np.zeros((0, 3))
    colors = np.array([CYAN_CARBON if (m == 1 and e == "C") else color(e) for e, m in zip(el, mol)], dtype=np.uint8).reshape(-1, 3)
    scene = Scene(pos, el, np.array(mol, int), np.zeros(len(pos), bool), colors,
                  np.array([covalent_radius(e) for e in el]), np.array(bonds, dtype=int).reshape(-1, 2))
    atoms = [where[(instances.index((m, tuple(t))), i)] for m, t, i in places if m >= 0]
    scene.names = names
    for a_, name in zip(atoms, labels[0]):  # gli atomi misurati con il codice di simmetria, come nella tabella
        scene.names[a_] = name
    scene.highlight = atoms
    pts = pos[atoms]
    if len(pts) >= 2:
        unit = {"angle": "°", "torsion": "°"}.get(item.kind, " Å")
        digits = 1 if unit == "°" else 3
        scene.measure_paths.append(Path(pts, COLOR_MEASURE, 0.035, dashed=item.kind == "contact"))
        at = pts.mean(axis=0) if item.kind in ("bond", "contact") else pts[1] if item.kind == "angle" else (pts[1] + pts[2]) / 2
        scene.measure_labels.append(Label(at, f"{item.value:.{digits}f}{unit}", COLOR_MEASURE, 15))
    if show_cell:
        v = structure.vectors
        corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=float) @ v
        edges = [(a, b) for a in range(8) for b in range(a + 1, 8) if bin(a ^ b).count("1") == 1]
        scene.box = np.array([[corners[a], corners[b]] for a, b in edges])
    scene.center = pos.mean(axis=0) if len(pos) else np.zeros(3)
    scene.size = float(np.linalg.norm(pos.max(axis=0) - pos.min(axis=0))) if len(pos) else 10.0
    scene.frame = SceneFrame.of(pos, pts[-1] - pts[0] if len(pts) >= 2 else None) if len(pos) else None
    return scene


def contact_measure(scene: Scene) -> tuple[list[int], list[Path], list[Label], str]:
    """Il contatto più corto della scena di una coppia come misura: atomi, linea tratteggiata, valore e testo."""
    if scene.contact is None:
        return [], [], [], ""
    i, j, d = scene.contact
    pts = scene.pos[[i, j]]
    text = f"{scene.names[i]}···{scene.names[j]} {d:.2f} Å" if scene.names else f"{d:.2f} Å"
    return [i, j], [Path(pts, COLOR_MEASURE, 0.03, dashed=True)], [Label(pts.mean(axis=0), f"{d:.2f} Å", COLOR_MEASURE, 15)], text


def scene_frame(scene: Scene) -> SceneFrame:
    """Gli assi di una scena (quelli calcolati, o gli assi principali dei suoi atomi)."""
    return scene.frame or SceneFrame.of(scene.pos if len(scene.pos) else np.eye(3))


def nearest_atom(scene: Scene, point, slack: float = 0.35) -> int | None:
    """L'atomo della scena sotto un punto scelto con il mouse (il punto sta sulla superficie della sfera disegnata):
    il più vicino, se il punto cade entro il suo raggio di van der Waals più `slack`."""
    if not len(scene.pos):
        return None
    d = np.linalg.norm(scene.pos - np.asarray(point, float), axis=1)
    i = int(np.argmin(d))
    return i if d[i] <= vdw_radius(scene.el[i]) + slack else None


def pick_measure(scene: Scene, picks: list[int]) -> tuple[list[Path], list[Label], str]:
    """Le misure tra gli atomi scelti con il mouse: con due la distanza, con tre anche l'angolo nel secondo.
    Restituisce linee, etichette e il testo per la barra della vista."""
    if not picks:
        return [], [], ""
    pts = scene.pos[picks]
    names = [scene.names[i] if i < len(scene.names) else scene.el[i] for i in picks]
    paths, labels, text = [], [], [names[0]]
    for k in range(1, len(picks)):
        d = float(np.linalg.norm(pts[k] - pts[k - 1]))
        paths.append(Path(pts[k - 1:k + 1], COLOR_MEASURE, 0.035, dashed=True))
        labels.append(Label((pts[k] + pts[k - 1]) / 2, f"{d:.3f} Å", COLOR_MEASURE, 15))
        text = [f"{names[k - 1]}···{names[k]} {d:.3f} Å"] if k == 1 else text + [f"{names[k - 1]}···{names[k]} {d:.3f} Å"]
    if len(picks) == 3:
        ang = angle(pts)
        labels.append(Label(pts[1], f"{ang:.1f}°", COLOR_MEASURE, 15))
        text.append(f"{names[0]}–{names[1]}–{names[2]} {ang:.1f}°")
    return paths, labels, " · ".join(text)
