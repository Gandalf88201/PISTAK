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
    show_measures: bool = True  # assi L, W, n e valori misurati lungo ciascuno
    show_cell: bool = False
    show_column: bool = True
    highlight_core: bool = True
    labels: bool = True
    ortho: bool = True  # proiezione ortogonale: dall'alto la vista coincide con la Proiezione


def _sticks(scene: Scene, colors: np.ndarray, radius: float, sides: int = 14) -> pv.PolyData | None:
    """Legami come cilindri, ciascuna metà del colore del suo atomo."""
    if not len(scene.bonds):
        return None
    pi, pj = scene.pos[scene.bonds[:, 0]], scene.pos[scene.bonds[:, 1]]
    mid = (pi + pj) / 2
    n = len(scene.bonds)
    pts = np.empty((4 * n, 3))
    pts[0::4], pts[1::4], pts[2::4], pts[3::4] = pi, mid, mid, pj
    rgb = np.empty((4 * n, 3), dtype=np.uint8)
    ci, cj = colors[scene.bonds[:, 0]], colors[scene.bonds[:, 1]]
    rgb[0::4] = rgb[1::4] = ci
    rgb[2::4] = rgb[3::4] = cj
    lines = np.empty((2 * n, 3), dtype=np.int64)
    base = np.arange(n) * 4
    lines[0::2] = np.column_stack([np.full(n, 2), base, base + 1])
    lines[1::2] = np.column_stack([np.full(n, 2), base + 2, base + 3])
    mesh = pv.PolyData(pts, lines=lines.ravel())
    mesh["rgb"] = rgb
    return mesh.tube(radius=radius, n_sides=sides)


def populate(plotter: pv.Plotter, scene: Scene, opts: Options, sphere_scale: float = 0.28, stick_radius: float = 0.11,
             px: float = 1.0) -> None:
    """Disegna la scena (senza svuotare il plotter né toccare la camera). `px` moltiplica ciò che si misura in pixel
    (testi, spessori delle linee) e affina sfere e cilindri: per le immagini esportate più grandi della finestra."""
    res = 24 if px <= 1 else 48
    if opts.style == 2:
        radii = np.array([vdw_radius(e) for e in scene.el])
    elif opts.style == 1:
        radii = np.full(len(scene.el), stick_radius * 1.05)
    else:
        radii = scene.radii * sphere_scale
    if opts.highlight_core and opts.style != 2:
        radii = np.where(scene.in_core, radii * 1.3, radii)
    colors = scene.colors
    if not theme.is_dark():  # idrogeni bianchi su fondo bianco: un grigio chiaro
        colors = np.where((colors == 255).all(axis=1, keepdims=True), np.uint8(200), colors).astype(np.uint8)
    atoms = pv.PolyData(scene.pos)
    atoms["rgb"] = colors
    atoms["radius"] = radii
    glyphs = atoms.glyph(scale="radius", factor=1.0, orient=False, geom=pv.Sphere(radius=1.0, theta_resolution=res, phi_resolution=res))
    plotter.add_mesh(glyphs, scalars="rgb", rgb=True, smooth_shading=True, specular=0.4)
    if opts.style != 2:
        sticks = _sticks(scene, colors, stick_radius, res // 2 + 2)
        if sticks is not None:
            plotter.add_mesh(sticks, scalars="rgb", rgb=True, smooth_shading=True)
    if opts.show_planes:
        for poly in scene.polygons:
            n = len(poly.points)
            face = pv.PolyData(poly.points, faces=np.r_[n, np.arange(n)])
            plotter.add_mesh(face, color=tuple(c / 255 for c in poly.color), opacity=poly.opacity, show_edges=True,
                             edge_color=tuple(c / 255 for c in poly.color), line_width=1.5 * px)
    if opts.show_vectors:
        for p in scene.paths:
            _add_path(plotter, p, sides=res // 2)
        for m in scene.markers:
            plotter.add_mesh(pv.Sphere(radius=m.radius, center=m.pos, theta_resolution=res, phi_resolution=res), color=tuple(c / 255 for c in m.color))
        if opts.labels and scene.markers:
            plotter.add_point_labels(np.array([m.pos for m in scene.markers]), [m.label for m in scene.markers],
                                     font_size=round(14 * px), shape=None, show_points=False, always_visible=True,
                                     text_color=tuple(c / 255 for c in theme.foreground()))
    if opts.show_measures:
        for p in scene.axes:
            _add_path(plotter, p, opacity=0.8, sides=res // 2)
        for color_, group in _by_color(scene.labels).items():
            for size, labels in group.items():
                plotter.add_point_labels(np.array([lb.pos for lb in labels]), [lb.text for lb in labels], font_size=round(size * px), bold=True,
                                         shape="rounded_rect", shape_color=tuple(c / 255 for c in theme.background()), shape_opacity=0.7,
                                         margin=round(2 * px), show_points=False, always_visible=True, text_color=tuple(c / 255 for c in color_))
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
            plotter.add_mesh(pv.Line(a, b), color="gray", line_width=px, opacity=0.6)


def _add_path(plotter: pv.Plotter, p, opacity: float = 1.0, sides: int = 12) -> None:
    color_ = tuple(c / 255 for c in p.color)
    a, b = p.points[0], p.points[-1]
    if p.arrow:  # la punta occupa l'ultimo tratto del segmento
        axis = (b - a) / np.linalg.norm(b - a)
        head = min(0.45, 0.3 * float(np.linalg.norm(b - a)))
        plotter.add_mesh(pv.Cone(center=b - axis * head / 2, direction=axis, height=head, radius=max(4 * p.radius, 0.09), resolution=2 * sides),
                         color=color_, opacity=opacity)
        b = b - axis * head
    plotter.add_mesh(pv.Line(a, b).tube(radius=p.radius, n_sides=sides), color=color_, smooth_shading=True, opacity=opacity)


def _by_color(labels) -> dict:
    """Etichette raggruppate per colore e dimensione: una chiamata a add_point_labels per gruppo."""
    out: dict = {}
    for lb in labels:
        out.setdefault(tuple(lb.color), {}).setdefault(lb.size, []).append(lb)
    return out


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
        #: Le immagini PNG si disegnano tante volte più grandi del pannello (testi e linee in proporzione).
        self.export_scale = max(1, int(cfg.get("export_scale", 3)))
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
                                ("show_measures", tr("Misure"), tr("Assi L, W, n della coppia con i valori misurati lungo ciascuno")),
                                ("highlight_core", tr("Nucleo"), tr("Ingrandisce gli atomi del nucleo π")),
                                ("show_column", tr("Colonna"), tr("Mostra le altre molecole della colonna di stack (grigio)")),
                                ("show_cell", tr("Cella"), tr("Spigoli della cella elementare")),
                                ("ortho", tr("Ortogonale"), tr("Proiezione ortogonale, senza prospettiva: vista dall'alto coincide con la Proiezione"))):
            cb = QCheckBox(text)
            cb.setToolTip(tip)
            cb.setChecked(getattr(self.opts, attr))
            cb.toggled.connect(lambda on, a=attr: self._set(a, on))
            self.boxes[attr] = cb
            self.toolbar.addWidget(cb)
        self.lock_cam = QCheckBox(tr("Mantieni la vista"))
        self.lock_cam.setToolTip(tr("Cambiando coppia la camera conserva angolazione e zoom rispetto alla coppia (assi L, W, n) "
                                    "invece di tornare alla vista dall'alto: per confrontare coppie diverse dallo stesso punto di vista"))
        self.toolbar.addWidget(self.lock_cam)
        for text, key in ((tr("Dall'alto"), "top"), (tr("Obliqua"), "oblique"), (tr("Di lato"), "side"), (tr("Lungo L"), "end")):
            b = QPushButton(text)
            b.setToolTip({"top": tr("Guarda perpendicolarmente ai piani, orientata come la Proiezione (L a destra, W in alto)"), "side": tr("Guarda di lato: si vede la distanza interplanare"),
                          "end": tr("Guarda lungo l'asse lungo della molecola"), "oblique": tr("Vista di tre quarti")}[key])
            b.clicked.connect(lambda _=False, k=key: self.look(k))
            self.toolbar.addWidget(b)
        shot = QPushButton(tr("PNG…"))
        shot.clicked.connect(self._save_png)
        self.toolbar.addWidget(shot)

        self._pending_camera = None
        self._pending_scale: float | None = None
        #: L'inquadratura data da `look` (vista, camera): finché l'utente non muove la camera, segue le dimensioni del pannello.
        self._auto_view: tuple[str, tuple] | None = None
        try:
            if os.environ.get("PISTAK_NO_3D"):
                raise RuntimeError(tr("disattivata con PISTAK_NO_3D"))
            self.plotter = QtInteractor(self)
            # su macOS il multisampling dentro Qt non agisce e i bordi escono a gradini: supersampling
            self.plotter.enable_anti_aliasing("ssaa")
        except Exception as e:  # noqa: BLE001  (nessun contesto OpenGL)
            self.plotter = _NullPlotter(self, str(e))
        self.plotter.set_background(theme.background())
        self.body.addWidget(self.plotter.interactor)

    # ---- opzioni -------------------------------------------------------------------
    def _set(self, attr: str, value) -> None:
        setattr(self.opts, attr, value)
        if attr == "ortho":
            self._apply_projection()
            self.plotter.render()
        elif attr in ("show_h", "show_cell", "show_column"):
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
        # la camera si conserva rispetto alla coppia (assi L, W, n, origine in Cg A), non nello spazio del cristallo:
        # così vale anche passando a un'altra coppia, che può stare altrove nella cella. Un'inquadratura ancora automatica
        # (un pulsante, mai toccata) si rifà con lo stesso pulsante sulla coppia nuova.
        auto = self._auto_view[0] if keep_camera and self._still_auto() else None
        relative = self._relative_camera() if keep_camera and auto is None else None
        pair = self.pair
        self._scene_structure = self.analysis.structure if self.analysis else None
        self._scene_key = self._shown_as()
        if not (self.analysis and pair):
            self.scene = self.ctx = None
        else:
            self.ctx = self.analysis.pair_context(pair)
            self.scene = build_scene(self.ctx, self.opts.show_h, self.opts.show_cell, self.opts.show_column)
        self._draw(False, relative=relative, view=auto or "top")

    def _draw(self, keep_camera: bool, relative: dict | None = None, view: str = "top") -> None:
        pl = self.plotter
        cam = pl.camera_position if keep_camera and self.scene is not None else None
        scale = self._parallel_scale() if keep_camera else None
        pl.clear()
        pl.set_background(theme.background())
        if self.scene is not None:
            populate(pl, self.scene, self.opts, self.sphere_scale, self.stick_radius)
            self._apply_projection()
            if self._pending_camera is not None:  # camera ripristinata da una sessione
                pl.camera_position, self._pending_camera = self._pending_camera, None
                self._set_parallel_scale(self._pending_scale)
                self._pending_scale = None
            elif relative is not None:
                self._apply_relative_camera(relative)
            elif cam is not None:
                pl.camera_position = cam
                self._set_parallel_scale(scale)
            else:
                self.look(view, render=False)  # di norma «top»: si parte come la Proiezione
        pl.render()

    # ---- camera --------------------------------------------------------------------
    def _camera(self):
        return None if isinstance(self.plotter, _NullPlotter) else self.plotter.camera

    def _apply_projection(self) -> None:
        cam = self._camera()
        if cam is not None and bool(cam.GetParallelProjection()) != self.opts.ortho:
            cam.SetParallelProjection(self.opts.ortho)

    def _parallel_scale(self) -> float | None:
        cam = self._camera()
        return float(cam.GetParallelScale()) if cam is not None else None

    def _set_parallel_scale(self, scale: float | None) -> None:
        cam = self._camera()
        if cam is not None and scale:
            cam.SetParallelScale(scale)

    def _framing(self, key: str) -> tuple[np.ndarray, float]:
        """Centro e semi-altezza (Å) dell'inquadratura. Dall'alto: lo stesso riquadro della Proiezione, cioè gli
        atomi nel piano (L, W) della coppia, centrati come lì."""
        if key != "top":
            return self.scene.center, max(self.scene.size, 4.0) * 0.51
        f = self.ctx.frame
        q = f.project(self.scene.pos)
        lo, hi = q[:, :2].min(axis=0), q[:, :2].max(axis=0)
        mid = (lo + hi) / 2
        center = f.origin + f.L * mid[0] + f.W * mid[1] + f.n * float(q[:, 2].mean())
        w = self.plotter.interactor.width() or 1
        h = self.plotter.interactor.height() or 1
        half = max((hi[1] - lo[1]) / 2, (hi[0] - lo[0]) / 2 * h / w)
        return center, max(half, 2.0) * 1.08

    def look(self, key: str, render: bool = True) -> None:
        if self.scene is None or self.ctx is None:
            return
        direction, up = view_directions(self.ctx.frame)[key]
        direction = direction / np.linalg.norm(direction)
        center, half = self._framing(key)
        # in prospettiva, la distanza a cui `half` riempie metà dell'angolo di vista (30°)
        dist = max(half / np.tan(np.radians(15)), max(self.scene.size, 4.0))
        self._apply_projection()
        self.plotter.camera_position = [tuple(center - direction * dist), tuple(center), tuple(up)]
        self._set_parallel_scale(half)
        self.plotter.reset_camera_clipping_range()
        self._auto_view = (key, self._camera_state())
        if render:
            self.plotter.render()

    def _relative_camera(self) -> dict | None:
        """La camera nelle coordinate (L, W, n) della coppia mostrata, con origine in Cg A."""
        if self.ctx is None or self._camera() is None or self.plotter.camera_position is None:
            return None
        f = self.ctx.frame
        axes = np.column_stack([f.L, f.W, f.n])
        pos, focal, up = (np.array(list(map(float, p))) for p in self.plotter.camera_position)
        return {"pos": (pos - f.origin) @ axes, "focal": (focal - f.origin) @ axes, "up": up @ axes, "scale": self._parallel_scale()}

    def _apply_relative_camera(self, rel: dict) -> None:
        f = self.ctx.frame
        axes = np.column_stack([f.L, f.W, f.n])
        self.plotter.camera_position = [tuple(f.origin + axes @ rel["pos"]), tuple(f.origin + axes @ rel["focal"]), tuple(axes @ rel["up"])]
        self._set_parallel_scale(rel["scale"])
        self.plotter.reset_camera_clipping_range()

    def _camera_state(self) -> tuple:
        cam = self.plotter.camera_position
        return (np.array([list(map(float, p)) for p in cam]) if cam is not None else None, self._parallel_scale())

    def _still_auto(self) -> bool:
        if self._auto_view is None or self.scene is None:
            return False
        (pos, scale), (pos0, scale0) = self._camera_state(), self._auto_view[1]
        return pos is not None and pos0 is not None and np.allclose(pos, pos0) and scale == scale0

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if self._still_auto():  # il pannello si è assestato dopo il primo disegno: si rifà la stessa inquadratura
            self.look(self._auto_view[0])

    def _save_png(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, tr("Salva immagine"), "vista3d.png", "PNG (*.png)")
        if path:
            self.screenshot(path)

    def screenshot(self, path: str, scale: int | None = None) -> None:
        """Salva la vista in PNG ad alta risoluzione: la scena si ridisegna fuori schermo `scale` volte più grande del
        pannello, con la stessa camera, testi e linee in proporzione e supersampling."""
        scale = scale or self.export_scale
        cam = self._camera()
        if self.scene is None or cam is None:
            self.plotter.screenshot(path)
            return
        w, h = max(self.plotter.interactor.width(), 100), max(self.plotter.interactor.height(), 100)
        pl = pv.Plotter(off_screen=True, window_size=(w * scale, h * scale))
        try:
            pl.set_background(theme.background())
            populate(pl, self.scene, self.opts, self.sphere_scale, self.stick_radius, px=scale)
            pl.camera_position = self.plotter.camera_position
            pl.camera.SetParallelProjection(cam.GetParallelProjection())
            pl.camera.SetParallelScale(cam.GetParallelScale())
            pl.camera.SetViewAngle(cam.GetViewAngle())
            pl.reset_camera_clipping_range()
            pl.enable_anti_aliasing("ssaa")
            pl.screenshot(path)
        finally:
            pl.close()

    def extra_state(self) -> dict:
        cam = self.plotter.camera_position
        return {"opts": asdict(self.opts), "lock_camera": self.lock_cam.isChecked(),
                "camera": [list(map(float, p)) for p in cam] if cam is not None else None, "scale": self._parallel_scale()}

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
        self._pending_scale = st.get("scale")

    def closeEvent(self, event) -> None:  # noqa: N802
        try:
            self.plotter.close()
        finally:
            super().closeEvent(event)
