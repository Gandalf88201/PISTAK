"""Vista 3D della coppia (PyVista/VTK dentro Qt).

`populate()` traduce una `Scene` (array numpy) in mesh PyVista su un qualunque `pv.Plotter`: è separata dal
widget Qt perché si possa provare anche fuori schermo (vedi tests/test_scene.py) e riusare per esportare immagini.
"""
from __future__ import annotations

import os
from dataclasses import asdict, dataclass, fields

os.environ.setdefault("QT_API", "pyside6")  # pyvistaqt/QtPy: usare PySide6

import numpy as np
import pyvista as pv
from PySide6.QtWidgets import QCheckBox, QComboBox, QFileDialog, QPushButton

from pistak.core.config import load_section
from pistak.core.elements import vdw_radius
from pistak.gui import theme
from pistak.gui.scene import Scene, build_scene
from pistak.gui.session import Session
from pistak.gui.views.base import PairView
from pistak.i18n import tf, tr

STYLES = ["Sfere e bastoncini", "Bastoncini", "Spazio occupato (vdW)"]


@dataclass
class Options:
    style: int = 0
    show_h: bool = True
    show_planes: bool = True
    show_vectors: bool = True
    show_cell: bool = False
    show_column: bool = True
    highlight_core: bool = True
    labels: bool = True


def _sticks(scene: Scene, radius: float) -> pv.PolyData | None:
    """Legami come cilindri, ciascuna metà del colore del suo atomo."""
    if not len(scene.bonds):
        return None
    pi, pj = scene.pos[scene.bonds[:, 0]], scene.pos[scene.bonds[:, 1]]
    mid = (pi + pj) / 2
    n = len(scene.bonds)
    pts = np.empty((4 * n, 3))
    pts[0::4], pts[1::4], pts[2::4], pts[3::4] = pi, mid, mid, pj
    rgb = np.empty((4 * n, 3), dtype=np.uint8)
    ci, cj = scene.colors[scene.bonds[:, 0]], scene.colors[scene.bonds[:, 1]]
    rgb[0::4] = rgb[1::4] = ci
    rgb[2::4] = rgb[3::4] = cj
    lines = np.empty((2 * n, 3), dtype=np.int64)
    base = np.arange(n) * 4
    lines[0::2] = np.column_stack([np.full(n, 2), base, base + 1])
    lines[1::2] = np.column_stack([np.full(n, 2), base + 2, base + 3])
    mesh = pv.PolyData(pts, lines=lines.ravel())
    mesh["rgb"] = rgb
    return mesh.tube(radius=radius, n_sides=14)


def populate(plotter: pv.Plotter, scene: Scene, opts: Options, sphere_scale: float = 0.28, stick_radius: float = 0.11) -> None:
    """Disegna la scena (senza svuotare il plotter né toccare la camera)."""
    if opts.style == 2:
        radii = np.array([vdw_radius(e) for e in scene.el])
    elif opts.style == 1:
        radii = np.full(len(scene.el), stick_radius * 1.05)
    else:
        radii = scene.radii * sphere_scale
    if opts.highlight_core and opts.style != 2:
        radii = np.where(scene.in_core, radii * 1.3, radii)
    atoms = pv.PolyData(scene.pos)
    atoms["rgb"] = scene.colors
    atoms["radius"] = radii
    glyphs = atoms.glyph(scale="radius", factor=1.0, orient=False, geom=pv.Sphere(radius=1.0, theta_resolution=24, phi_resolution=24))
    plotter.add_mesh(glyphs, scalars="rgb", rgb=True, smooth_shading=True, specular=0.4)
    if opts.style != 2:
        sticks = _sticks(scene, stick_radius)
        if sticks is not None:
            plotter.add_mesh(sticks, scalars="rgb", rgb=True, smooth_shading=True)
    if opts.show_planes:
        for poly in scene.polygons:
            n = len(poly.points)
            face = pv.PolyData(poly.points, faces=np.r_[n, np.arange(n)])
            plotter.add_mesh(face, color=tuple(c / 255 for c in poly.color), opacity=poly.opacity, show_edges=True,
                             edge_color=tuple(c / 255 for c in poly.color), line_width=1.5)
    if opts.show_vectors:
        for p in scene.paths:
            seg = pv.Line(p.points[0], p.points[1]).tube(radius=p.radius, n_sides=10)
            plotter.add_mesh(seg, color=tuple(c / 255 for c in p.color), smooth_shading=True)
        for m in scene.markers:
            plotter.add_mesh(pv.Sphere(radius=m.radius, center=m.pos), color=tuple(c / 255 for c in m.color))
        if opts.labels and scene.markers:
            plotter.add_point_labels(np.array([m.pos for m in scene.markers]), [m.label for m in scene.markers],
                                     font_size=14, shape=None, show_points=False, always_visible=True,
                                     text_color=tuple(c / 255 for c in theme.foreground()))
    if len(scene.ghost_pos):
        ghost = pv.PolyData(scene.ghost_pos)
        plotter.add_mesh(ghost.glyph(scale=False, factor=0.2, orient=False, geom=pv.Sphere(radius=1.0, theta_resolution=12, phi_resolution=12)),
                         color=(0.62, 0.62, 0.66), opacity=0.35)
        if len(scene.ghost_bonds):
            n = len(scene.ghost_bonds)
            lines = np.column_stack([np.full(n, 2), scene.ghost_bonds]).ravel()
            plotter.add_mesh(pv.PolyData(scene.ghost_pos, lines=lines).tube(radius=0.05, n_sides=8), color=(0.62, 0.62, 0.66), opacity=0.35)
    if len(scene.box):
        for a, b in scene.box:
            plotter.add_mesh(pv.Line(a, b), color="gray", line_width=1, opacity=0.6)


def view_directions(frame) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Direzione di sguardo e verso 'alto' per le tre viste standard della coppia."""
    return {
        "top": (-frame.n, frame.W),  # dall'alto, perpendicolare ai piani: si vede la sovrapposizione
        "side": (-frame.W, frame.n),  # di lato: L orizzontale, distanza interplanare verticale
        "end": (-frame.L, frame.n),  # lungo l'asse lungo
        # obliqua: dall'alto e di tre quarti, per leggere insieme sovrapposizione e distanza
        "oblique": (-(0.45 * frame.L - 1.0 * frame.W + 0.9 * frame.n), frame.n),
    }


class _NullPlotter:
    """Sostituto del plotter quando non c'è OpenGL (o con PISTAK_NO_3D=1): le altre viste funzionano lo stesso."""

    camera_position = None

    def __init__(self, parent, reason: str):
        from PySide6.QtWidgets import QLabel

        self.interactor = QLabel(tf("Vista 3D non disponibile: {reason}", reason=reason), parent)

    def __getattr__(self, name):
        return lambda *a, **k: None


class View3D(PairView):
    TITLE = "3D"

    def __init__(self, session: Session, group: str | None = "A", parent=None):
        super().__init__(session, group, parent)
        from pyvistaqt import QtInteractor

        cfg = load_section("gui")
        self.sphere_scale = float(cfg.get("sphere_scale", 0.28))
        self.stick_radius = float(cfg.get("stick_radius", 0.11))
        self.opts = Options(show_h=bool(cfg.get("show_hydrogens", True)))
        self.scene: Scene | None = None
        self.ctx = None
        #: Da che cosa è stata costruita la scena corrente (vedi `on_pair`).
        self._scene_structure = None
        self._scene_key: tuple | None = None

        self.style_combo = QComboBox()
        self.style_combo.addItems([tr(s) for s in STYLES])
        self.style_combo.currentIndexChanged.connect(lambda i: self._set("style", i))
        self.toolbar.addWidget(self.style_combo)
        self.boxes: dict[str, QCheckBox] = {}
        for attr, text, tip in (("show_h", "H", tr("Mostra gli idrogeni")), ("show_planes", tr("Piani"), tr("Piani medi dei nuclei π")),
                                ("show_vectors", tr("Vettori"), tr("Scomposizione dello slittamento: L, W e distanza interplanare")),
                                ("highlight_core", tr("Nucleo"), tr("Ingrandisce gli atomi del nucleo π")),
                                ("show_column", tr("Colonna"), tr("Mostra le altre molecole della colonna di stack (grigio)")),
                                ("show_cell", tr("Cella"), tr("Spigoli della cella elementare"))):
            cb = QCheckBox(text)
            cb.setToolTip(tip)
            cb.setChecked(getattr(self.opts, attr))
            cb.toggled.connect(lambda on, a=attr: self._set(a, on))
            self.boxes[attr] = cb
            self.toolbar.addWidget(cb)
        self.lock_cam = QCheckBox(tr("Blocca camera"))
        self.lock_cam.setToolTip(tr("Non riporta la camera sulla coppia quando cambia la selezione"))
        self.toolbar.addWidget(self.lock_cam)
        for text, key in ((tr("Obliqua"), "oblique"), (tr("Dall'alto"), "top"), (tr("Di lato"), "side"), (tr("Lungo L"), "end")):
            b = QPushButton(text)
            b.setToolTip({"top": tr("Guarda perpendicolarmente ai piani"), "side": tr("Guarda di lato: si vede la distanza interplanare"),
                          "end": tr("Guarda lungo l'asse lungo della molecola"), "oblique": tr("Vista di tre quarti")}[key])
            b.clicked.connect(lambda _=False, k=key: self.look(k))
            self.toolbar.addWidget(b)
        shot = QPushButton(tr("PNG…"))
        shot.clicked.connect(self._save_png)
        self.toolbar.addWidget(shot)

        self._pending_camera = None
        try:
            if os.environ.get("PISTAK_NO_3D"):
                raise RuntimeError(tr("disattivata con PISTAK_NO_3D"))
            self.plotter = QtInteractor(self)
        except Exception as e:  # noqa: BLE001  (nessun contesto OpenGL)
            self.plotter = _NullPlotter(self, str(e))
        self.plotter.set_background(theme.background())
        self.body.addWidget(self.plotter.interactor)

    # ---- opzioni -------------------------------------------------------------------
    def _set(self, attr: str, value) -> None:
        setattr(self.opts, attr, value)
        if attr in ("show_h", "show_cell", "show_column"):
            self._rebuild(keep_camera=True)
        else:
            self._draw(keep_camera=True)

    # ---- disegno -------------------------------------------------------------------
    def on_analysis(self) -> None:
        if not self.analysis or not self.analysis.pairs:
            self.scene = self.ctx = None
            self._draw(False)

    def _shown_as(self) -> tuple | None:
        """Ciò da cui dipende la scena: la coppia (contatto e tipo) e la sua colonna, sulla stessa struttura."""
        a, p = self.analysis, self.pair
        if not (a and p):
            return None
        col = a.column_for(p)
        return (p.key, p.kind, tuple(col.molecules) if col else None)

    def on_pair(self, index: int) -> None:
        # con le soglie in tempo reale la stessa coppia torna a ogni scatto del cursore: se la scena non cambia
        # (stessa struttura, stesso contatto, stesso tipo e colonna) non si ricostruisce
        if self.scene is not None and self.analysis is not None and self._scene_structure is self.analysis.structure \
                and self._scene_key is not None and self._scene_key == self._shown_as():
            self.ctx = self.analysis.pair_context(self.pair)
            return
        self._rebuild(keep_camera=self.lock_cam.isChecked())

    def _rebuild(self, keep_camera: bool) -> None:
        pair = self.pair
        self._scene_structure = self.analysis.structure if self.analysis else None
        self._scene_key = self._shown_as()
        if not (self.analysis and pair):
            self.scene = self.ctx = None
        else:
            self.ctx = self.analysis.pair_context(pair)
            self.scene = build_scene(self.ctx, self.opts.show_h, self.opts.show_cell, self.opts.show_column)
        self._draw(keep_camera)

    def _draw(self, keep_camera: bool) -> None:
        pl = self.plotter
        cam = pl.camera_position if keep_camera and self.scene is not None else None
        pl.clear()
        pl.set_background(theme.background())
        if self.scene is not None:
            populate(pl, self.scene, self.opts, self.sphere_scale, self.stick_radius)
            if self._pending_camera is not None:  # camera ripristinata da una sessione
                pl.camera_position, self._pending_camera = self._pending_camera, None
            elif cam is not None:
                pl.camera_position = cam
            else:
                self.look("oblique", render=False)
        pl.render()

    def look(self, key: str, render: bool = True) -> None:
        if self.scene is None or self.ctx is None:
            return
        direction, up = view_directions(self.ctx.frame)[key]
        direction = direction / np.linalg.norm(direction)
        center = self.scene.center
        dist = max(self.scene.size, 4.0) * 1.9
        self.plotter.camera_position = [tuple(center - direction * dist), tuple(center), tuple(up)]
        self.plotter.reset_camera_clipping_range()
        if render:
            self.plotter.render()

    def _save_png(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, tr("Salva immagine"), "vista3d.png", "PNG (*.png)")
        if path:
            self.screenshot(path)

    def screenshot(self, path: str) -> None:
        self.plotter.screenshot(path)

    def extra_state(self) -> dict:
        cam = self.plotter.camera_position
        return {"opts": asdict(self.opts), "lock_camera": self.lock_cam.isChecked(),
                "camera": [list(map(float, p)) for p in cam] if cam is not None else None}

    def apply_state(self, st: dict) -> None:
        names = {f.name for f in fields(Options)}
        for k, v in st.get("opts", {}).items():
            if k in names:
                setattr(self.opts, k, v)
        self.style_combo.blockSignals(True)
        self.style_combo.setCurrentIndex(self.opts.style)
        self.style_combo.blockSignals(False)
        for attr, cb in self.boxes.items():
            cb.blockSignals(True)
            cb.setChecked(bool(getattr(self.opts, attr)))
            cb.blockSignals(False)
        self.lock_cam.setChecked(bool(st.get("lock_camera", False)))
        cam = st.get("camera")
        self._pending_camera = [tuple(p) for p in cam] if cam else None

    def closeEvent(self, event) -> None:  # noqa: N802
        try:
            self.plotter.close()
        finally:
            super().closeEvent(event)
