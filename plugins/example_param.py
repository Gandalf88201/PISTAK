"""Esempio di plugin: due parametri nuovi. Copia questo file, cambialo e premi Analisi ▸ Ricarica i plugin.

Un parametro è una funzione che riceve il `PairContext` della coppia e restituisce un numero. Il decoratore
@parameter lo registra: da quel momento compare nella tabella (Colonne ▾), negli assi dello scatter e
dell'istogramma, nel pannello Coppia e negli export CSV/JSON. Se la funzione solleva un'eccezione la coppia
vale NaN e l'errore si legge nel pannello Dettaglio ▸ Avvisi: un plugin sbagliato non ferma l'analisi.

Cosa trovi in `ctx` (vedi pistak/core/analysis.py, classe PairContext):
    ctx.pair                la coppia (kind, interplanar, long_slip, lat_slip, dihedral, contact, ...)
    ctx.A, ctx.B            le due molecole: .el (elementi), .pos (N×3, Å), .bonds, .heavy (maschera non-H)
    ctx.core_a, ctx.core_b  nuclei π: .atoms (indici), .centroid, .normal, .long, .rms
    ctx.pos_a, ctx.pos_b    posizioni degli atomi, con B già nella cella giusta
    ctx.frame               sistema di riferimento della coppia: .L, .W, .n (versori), .origin
    ctx.d                   vettore tra i baricentri dei nuclei; ctx.proj_a/proj_b: atomi nel sistema (L, W, n)
    ctx.limits              le soglie in uso; ctx.analysis: tutto il risultato (celle, altre coppie, colonne)
"""
import numpy as np

from pistak.params import parameter


@parameter("slip_over_distance", "Slittamento / distanza centroidi", unit="", group="plugin", digits=3,
           description="Quanto del vettore tra i centroidi sta nel piano: 0 = sovrapposti, 1 = affiancati.")
def slip_over_distance(ctx):
    return ctx.pair.slip / ctx.pair.centroid_distance


@parameter("heavy_atoms_in_contact", "Atomi pesanti entro 3.6 Å", unit="", group="plugin", digits=0,
           description="Quanti atomi pesanti di A hanno almeno un atomo pesante di B più vicino di 3.6 Å.")
def heavy_atoms_in_contact(ctx):
    from scipy.spatial.distance import cdist

    ha, hb = ctx.A.heavy, ctx.B.heavy
    d = cdist(ctx.pos_a[ha], ctx.pos_b[hb])
    return float(np.count_nonzero((d < 3.6).any(axis=1)))
