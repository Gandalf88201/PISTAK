"""Parametri aggiuntivi rispetto al sito (fase 2). Si possono modificare liberamente o copiare in `plugins/`.

Dipendono tutti solo dalla geometria della coppia (`geometric=True`): cambiando le soglie del criterio di stack si
riusano invece di ricalcolarli. Chi ne copia uno in un plugin e lo fa dipendere dal tipo della coppia o dalle soglie
deve togliere `geometric`."""
from __future__ import annotations

import math
from functools import partial

import numpy as np
from scipy.spatial.distance import cdist

from pistak.core.elements import vdw_radius
from pistak.core.geometry import acute_angle_deg, polygon_area
from pistak.core.pairs import fragment_polygons, polygons_overlap
from pistak.params.registry import parameter as _parameter

parameter = partial(_parameter, geometric=True)

G = "angoli"


@parameter("slip_angle", "Direzione dello slittamento", "°", "Angolo di s nel piano rispetto all'asse lungo: 0° = slittamento tutto longitudinale, 90° = tutto laterale.", 1, G)
def slip_angle(ctx):
    return float(np.degrees(np.arctan2(ctx.pair.lat_slip, ctx.pair.long_slip)))


@parameter("long_slip_signed", "Slittamento longitudinale con segno", "Å", "Componente di s (da A verso B) lungo L. Il verso di L è quello dell'asse definito dagli atomi (Soglie ▸ Asse lungo da due atomi), altrimenti una convenzione riproducibile; scambiando A e B cambia segno.", 2, G)
def long_slip_signed(ctx):
    if math.isnan(ctx.pair.long_slip):  # asse lungo non definito o assi incrociati (vedi Pair.long_slip)
        return math.nan
    return float(ctx.slip_vector @ ctx.frame.L)


@parameter("lat_slip_signed", "Slittamento laterale con segno", "Å", "Componente di s (da A verso B) lungo W = n × L.", 2, G)
def lat_slip_signed(ctx):
    if math.isnan(ctx.pair.lat_slip):
        return math.nan
    return float(ctx.slip_vector @ ctx.frame.W)


@parameter("slip_direction", "Direzione dello slittamento con segno", "°", "Angolo di s nel piano, da L verso W (−180…180°): distingue i quattro quadranti che «Direzione dello slittamento» confonde.", 1, G)
def slip_direction(ctx):
    if math.isnan(ctx.pair.long_slip):
        return math.nan
    a = float(np.degrees(np.arctan2(ctx.slip_vector @ ctx.frame.W, ctx.slip_vector @ ctx.frame.L)))
    return a + 360.0 if a <= -180.0 + 1e-6 else a  # in (−180, 180]: lo zero con segno non decide tra −180 e 180


@parameter("beta", "β (Cg→normale A)", "°", "Angolo tra il vettore tra i centroidi e la normale del nucleo A (β di PLATON e Janiak).", 1, G)
def beta(ctx):
    return acute_angle_deg(ctx.d, ctx.core_a.normal)


@parameter("gamma", "γ (Cg→normale B)", "°", "Angolo tra il vettore tra i centroidi e la normale del nucleo B (γ di PLATON e Janiak).", 1, G)
def gamma(ctx):
    return acute_angle_deg(ctx.d, ctx.core_b.normal)


@parameter("height_b_over_a", "Quota di B sul piano A", "Å", "Distanza del centroide di B dal piano medio di A.", 2, "distanze")
def height_b_over_a(ctx):
    return abs(float(np.dot(ctx.d, ctx.core_a.normal)))


@parameter("height_a_over_b", "Quota di A sul piano B", "Å", "Distanza del centroide di A dal piano medio di B.", 2, "distanze")
def height_a_over_b(ctx):
    return abs(float(np.dot(ctx.d, ctx.core_b.normal)))


@parameter("mean_atom_height", "Quota media degli atomi di B", "Å", "Distanza media degli atomi del nucleo B dal piano medio di A (media delle distanze, non distanza della media): con piani inclinati non coincide con la quota del baricentro.", 2, "distanze")
def mean_atom_height(ctx):
    # la media delle quote con segno sarebbe la quota del baricentro di B (height_b_over_a): servono i valori assoluti
    z = (ctx.pos_b[ctx.core_b.atoms] - ctx.core_a.centroid) @ ctx.core_a.normal
    return float(np.abs(z).mean())


def _core_heavy(ctx, which: str):
    mol, core, pos = (ctx.A, ctx.core_a, ctx.pos_a) if which == "a" else (ctx.B, ctx.core_b, ctx.pos_b)
    return core.atoms, pos[core.atoms], [mol.el[i] for i in core.atoms]


@parameter("core_contact", "Contatto più corto tra i nuclei", "Å", "Distanza minima tra atomi dei due nuclei π.", 2, "contatti")
def core_contact(ctx):
    _, pa, _ = _core_heavy(ctx, "a")
    _, pb, _ = _core_heavy(ctx, "b")
    return float(cdist(pa, pb).min())


def _contact_atoms(ctx, which: str) -> list[int]:
    """Gli atomi tra cui cercare i contatti: tutta la molecola, o solo il frammento π se le due «molecole» sono
    legate tra loro (coppia nella stessa molecola, o unità di una rete estesa), dove gli altri atomi sono legami."""
    mol, core = (ctx.A, ctx.core_a) if which == "a" else (ctx.B, ctx.core_b)
    if ctx.pair.intramolecular or ctx.pair.polymer:
        return [int(i) for i in core.atoms]
    return list(range(len(mol.el)))


@parameter("cc_contact", "Contatto C···C più corto", "Å", "Distanza minima tra due atomi di carbonio delle due molecole (NaN se mancano; nella stessa molecola o in un polimero, tra i due frammenti π).", 2, "contatti")
def cc_contact(ctx):
    ia = [i for i in _contact_atoms(ctx, "a") if ctx.A.el[i] == "C"]
    ib = [i for i in _contact_atoms(ctx, "b") if ctx.B.el[i] == "C"]
    return float(cdist(ctx.pos_a[ia], ctx.pos_b[ib]).min()) if ia and ib else float("nan")


@parameter("n_vdw_contacts", "Contatti sotto vdW", "", "Numero di coppie di atomi pesanti, uno per molecola, più vicine della somma dei raggi di van der Waals (nella stessa molecola o in un polimero, tra i due frammenti π).", 0, "contatti")
def n_vdw_contacts(ctx):
    ha = np.array([i for i in _contact_atoms(ctx, "a") if ctx.A.el[i] != "H"], dtype=int)
    hb = np.array([i for i in _contact_atoms(ctx, "b") if ctx.B.el[i] != "H"], dtype=int)
    ra = np.array([vdw_radius(ctx.A.el[i]) for i in ha])
    rb = np.array([vdw_radius(ctx.B.el[i]) for i in hb])
    d = cdist(ctx.pos_a[ha], ctx.pos_b[hb])
    return float((d < ra[:, None] + rb[None, :]).sum())


def _polygons(ctx):
    hull_only = ctx.limits.pair_mode == "site"
    qa = ctx.proj_a[:, :2]
    qb = ctx.proj_b[:, :2]
    return fragment_polygons(ctx.core_a, qa, hull_only), fragment_polygons(ctx.core_b, qb, hull_only)


@parameter("overlap_area", "Area di sovrapposizione", "Å²", "Area comune dei due frammenti proiettati sul piano medio: unione dei poligoni dei singoli anelli (nel modo del sito, involucro convesso).", 2, "sovrapposizione")
def overlap_area(ctx):
    pa, pb = _polygons(ctx)
    return polygons_overlap(pa, pb)


@parameter("overlap_ratio", "Sovrapposizione / nucleo minore", "", "Area di sovrapposizione divisa per l'area del frammento più piccolo (0–1): la misura del criterio di stack predefinito.", 2, "sovrapposizione")
def overlap_ratio(ctx):
    pa, pb = _polygons(ctx)
    smaller = min(sum(polygon_area(p) for p in pa), sum(polygon_area(p) for p in pb))
    return polygons_overlap(pa, pb) / smaller if smaller > 0 else 0.0


@parameter("twist", "Rotazione tra gli assi lunghi", "°", "Angolo (0–90°) tra gli assi lunghi dei due frammenti proiettati sul piano medio; NaN se uno dei due assi non è definito (nucleo quasi isotropo).", 1, G)
def twist(ctx):
    return ctx.pair.twist


@parameter("planes_cross", "Piani che si intersecano tra i centroidi", "", "1 se la retta d'intersezione dei due piani passa tra i due centroidi (d·n_A e d·n_B di segno opposto): la distanza interplanare è allora più piccola delle due distanze perpendicolari e va letta con cautela.", 0, "distanze")
def planes_cross(ctx):
    return float(ctx.pair.planes_cross)


@parameter("local_interplanar", "Distanza interplanare locale", "Å", "Media delle distanze degli atomi di ciascun frammento che cadono, in proiezione, sull'altro frammento dal piano medio dell'altro: la separazione nella zona di sovrapposizione (NaN se non c'è sovrapposizione).", 2, "distanze")
def local_interplanar(ctx):
    pa, pb = _polygons(ctx)
    qa = ctx.proj_a[ctx.core_a.atoms, :2]
    qb = ctx.proj_b[ctx.core_b.atoms, :2]
    from pistak.core.geometry import points_inside

    in_a = np.zeros(len(qa), dtype=bool)
    for pg in pb:
        in_a |= points_inside(pg, qa)
    in_b = np.zeros(len(qb), dtype=bool)
    for pg in pa:
        in_b |= points_inside(pg, qb)
    h = []
    if in_b.any():
        h += list(np.abs((ctx.pos_b[ctx.core_b.atoms][in_b] - ctx.core_a.centroid) @ ctx.core_a.normal))
    if in_a.any():
        h += list(np.abs((ctx.pos_a[ctx.core_a.atoms][in_a] - ctx.core_b.centroid) @ ctx.core_b.normal))
    return float(np.mean(h)) if h else math.nan


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


@parameter("fragment_a", "Frammento π di A", "", "Numero del frammento π della molecola A (1 = il sistema di anelli più grande).", 0, "nuclei")
def fragment_a(ctx):
    return ctx.pair.fa + 1


@parameter("fragment_b", "Frammento π di B", "", "", 0, "nuclei")
def fragment_b(ctx):
    return ctx.pair.fb + 1


@parameter("axis_aspect_a", "Isotropia del nucleo A", "", "Dispersione lungo l'asse corto / lungo l'asse lungo, nel piano (0–1). Vicino a 1 (coronene, porfirina) l'asse lungo automatico è mal definito e la divisione dello slittamento in L e W è arbitraria: conviene definire l'asse con due atomi.", 2, "nuclei")
def axis_aspect_a(ctx):
    return ctx.core_a.plane.aspect


@parameter("axis_aspect_b", "Isotropia del nucleo B", "", "", 2, "nuclei")
def axis_aspect_b(ctx):
    return ctx.core_b.plane.aspect
