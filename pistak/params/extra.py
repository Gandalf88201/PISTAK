"""Parametri aggiuntivi rispetto al sito (fase 2). Si possono modificare liberamente o copiare in `plugins/`."""
from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist

from pistak.core.elements import vdw_radius
from pistak.core.geometry import acute_angle_deg, clip_convex, convex_hull, polygon_area
from pistak.params.registry import parameter

G = "angoli"


@parameter("slip_angle", "Direzione dello slittamento", "°", "Angolo di s nel piano rispetto all'asse lungo: 0° = slittamento tutto longitudinale, 90° = tutto laterale.", 1, G)
def slip_angle(ctx):
    return float(np.degrees(np.arctan2(ctx.pair.lat_slip, ctx.pair.long_slip)))


@parameter("beta_a", "β (Cg→normale A)", "°", "Angolo tra il vettore tra i centroidi e la normale del nucleo A (stile Janiak).", 1, G)
def beta_a(ctx):
    return acute_angle_deg(ctx.d, ctx.core_a.normal)


@parameter("beta_b", "γ (Cg→normale B)", "°", "Angolo tra il vettore tra i centroidi e la normale del nucleo B.", 1, G)
def beta_b(ctx):
    return acute_angle_deg(ctx.d, ctx.core_b.normal)


@parameter("height_b_over_a", "Quota di B sul piano A", "Å", "Distanza del centroide di B dal piano medio di A.", 2, "distanze")
def height_b_over_a(ctx):
    return abs(float(np.dot(ctx.d, ctx.core_a.normal)))


@parameter("height_a_over_b", "Quota di A sul piano B", "Å", "Distanza del centroide di A dal piano medio di B.", 2, "distanze")
def height_a_over_b(ctx):
    return abs(float(np.dot(ctx.d, ctx.core_b.normal)))


@parameter("mean_atom_height", "Quota media degli atomi di B", "Å", "Distanza media degli atomi del nucleo B dal piano medio di A: meno sensibile al baricentro.", 2, "distanze")
def mean_atom_height(ctx):
    z = ctx.proj_b[ctx.core_b.atoms, 2]
    return float(abs(z.mean()))


def _core_heavy(ctx, which: str):
    mol, core, pos = (ctx.A, ctx.core_a, ctx.pos_a) if which == "a" else (ctx.B, ctx.core_b, ctx.pos_b)
    return core.atoms, pos[core.atoms], [mol.el[i] for i in core.atoms]


@parameter("core_contact", "Contatto più corto tra i nuclei", "Å", "Distanza minima tra atomi dei due nuclei π.", 2, "contatti")
def core_contact(ctx):
    _, pa, _ = _core_heavy(ctx, "a")
    _, pb, _ = _core_heavy(ctx, "b")
    return float(cdist(pa, pb).min())


@parameter("cc_contact", "Contatto C···C più corto", "Å", "Distanza minima tra due atomi di carbonio delle due molecole (NaN se mancano).", 2, "contatti")
def cc_contact(ctx):
    ia = [i for i, e in enumerate(ctx.A.el) if e == "C"]
    ib = [i for i, e in enumerate(ctx.B.el) if e == "C"]
    return float(cdist(ctx.pos_a[ia], ctx.pos_b[ib]).min()) if ia and ib else float("nan")


@parameter("n_vdw_contacts", "Contatti sotto vdW", "", "Numero di coppie di atomi pesanti, uno per molecola, più vicine della somma dei raggi di van der Waals.", 0, "contatti")
def n_vdw_contacts(ctx):
    ha, hb = np.flatnonzero(ctx.A.heavy), np.flatnonzero(ctx.B.heavy)
    ra = np.array([vdw_radius(ctx.A.el[i]) for i in ha])
    rb = np.array([vdw_radius(ctx.B.el[i]) for i in hb])
    d = cdist(ctx.pos_a[ha], ctx.pos_b[hb])
    return float((d < ra[:, None] + rb[None, :]).sum())


def _hulls(ctx):
    qa = ctx.proj_a[ctx.core_a.atoms, :2]
    qb = ctx.proj_b[ctx.core_b.atoms, :2]
    return convex_hull(qa), convex_hull(qb)


@parameter("overlap_area", "Area di sovrapposizione", "Å²", "Area dell'intersezione dei contorni dei due nuclei, proiettati sul piano medio.", 2, "sovrapposizione")
def overlap_area(ctx):
    ha, hb = _hulls(ctx)
    if len(ha) < 3 or len(hb) < 3:
        return 0.0
    return polygon_area(clip_convex(ha, hb))


@parameter("overlap_ratio", "Sovrapposizione / nucleo minore", "", "Area di sovrapposizione divisa per l'area del contorno più piccolo (0–1).", 2, "sovrapposizione")
def overlap_ratio(ctx):
    ha, hb = _hulls(ctx)
    smaller = min(polygon_area(ha), polygon_area(hb))
    return overlap_area(ctx) / smaller if smaller > 0 else 0.0


@parameter("core_size_a", "Atomi nel nucleo A", "", "", 0, "nuclei")
def core_size_a(ctx):
    return len(ctx.core_a.atoms)


@parameter("core_size_b", "Atomi nel nucleo B", "", "", 0, "nuclei")
def core_size_b(ctx):
    return len(ctx.core_b.atoms)


@parameter("core_rms_a", "Planarità nucleo A", "Å", "Scarto quadratico medio degli atomi del nucleo dal piano medio.", 3, "nuclei")
def core_rms_a(ctx):
    return ctx.core_a.rms


@parameter("core_rms_b", "Planarità nucleo B", "Å", "", 3, "nuclei")
def core_rms_b(ctx):
    return ctx.core_b.rms
