"""Descrizione geometrica della scena 3D di una coppia, indipendente dalla libreria di rendering.

`build_scene` produce solo array numpy e piccole strutture: la vista 3D (`views/view3d.py`) le disegna con
PyVista, ma si potrebbe scrivere un altro disegnatore (VTK puro, matplotlib, esportazione in un file) senza
toccare il calcolo.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from pistak.core.analysis import PairContext
from pistak.core.elements import CYAN_CARBON, color, covalent_radius
from pistak.core.geometry import convex_hull

RGB = tuple[int, int, int]

COLOR_L: RGB = (214, 69, 65)  # slittamento longitudinale
COLOR_W: RGB = (46, 160, 67)  # slittamento laterale
COLOR_N: RGB = (47, 111, 214)  # distanza interplanare
COLOR_PLANE_A: RGB = (240, 165, 60)
COLOR_PLANE_B: RGB = (64, 190, 220)


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


@dataclass
class Marker:
    pos: np.ndarray
    color: RGB
    radius: float
    label: str = ""


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
    ghost_pos: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))  # altre molecole della colonna
    ghost_bonds: np.ndarray = field(default_factory=lambda: np.zeros((0, 2), dtype=int))
    box: np.ndarray = field(default_factory=lambda: np.zeros((0, 2, 3)))  # E×2×3: spigoli della cella
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))
    size: float = 10.0


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
    scene.markers += [Marker(ca.centroid, COLOR_PLANE_A, 0.16, "Cg A"), Marker(cb.centroid, COLOR_PLANE_B, 0.16, "Cg B")]

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
