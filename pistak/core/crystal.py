"""Grandezze dell'intero cristallo: contenuto della cella, formula per unità formula, densità, Z′ e i valori
dell'affinamento dichiarati nel CIF.

Il contenuto della cella si conta sul modello del CIF con tutte le componenti disordinate (anche quelle che la regola
`disorder` ha scartato per l'analisi delle coppie), ciascun atomo pesato con la sua occupanza e moltiplicato per il
numero di posizioni distinte che genera la simmetria (le posizioni speciali ne generano meno).
"""
from __future__ import annotations

import math
import re

import numpy as np

from pistak.core.cif import Cif, CifAtom, cif_number
from pistak.core.elements import atomic_mass
from pistak.core.structure import cell_vectors

#: Da u/Å³ a g/cm³ (1 u = 1,66054·10⁻²⁴ g, 1 Å³ = 10⁻²⁴ cm³).
U_PER_A3_TO_G_PER_CM3 = 1.66053907

#: Valori dell'affinamento letti dal CIF: chiave fissa → etichette CIF in ordine di preferenza.
REFINEMENT_TAGS: dict[str, tuple[str, ...]] = {
    "r1": ("_refine_ls_r_factor_gt", "_refine_ls_r_factor_obs"),
    "wr2": ("_refine_ls_wr_factor_ref", "_refine_ls_wr_factor_all"),
    "goof": ("_refine_ls_goodness_of_fit_ref", "_refine_ls_goodness_of_fit_all"),
    "rho_max": ("_refine_diff_density_max",),
    "rho_min": ("_refine_diff_density_min",),
    "temperature": ("_diffrn_ambient_temperature", "_cell_measurement_temperature"),
    "density_cif": ("_exptl_crystal_density_diffrn",),
    "z": ("_cell_formula_units_z",),
}


def cif_value(cif: Cif, key: str) -> float:
    """Un valore numerico di `REFINEMENT_TAGS` (NaN se il CIF non lo dà)."""
    for tag in REFINEMENT_TAGS[key]:
        v = cif_number(cif.items.get(tag))
        if math.isfinite(v):
            return v
    return math.nan


def cell_volume(cif: Cif) -> float:
    return float(abs(np.linalg.det(cell_vectors(cif.cell))))


def site_multiplicity(cif: Cif, atom: CifAtom, vectors: np.ndarray | None = None, tol: float = 0.1) -> int:
    """Quante posizioni distinte della cella genera la simmetria a partire da `atom`."""
    vectors = cell_vectors(cif.cell) if vectors is None else vectors
    seen: list[np.ndarray] = []
    for op in cif.symops:
        f = op.apply(atom.frac)
        f = f - np.floor(f + 1e-9)
        if not any(np.linalg.norm(((f - g) - np.round(f - g)) @ vectors) < tol for g in seen):
            seen.append(f)
    return len(seen)


def cell_contents(cif: Cif, include_omitted: bool = True, tol: float = 0.1) -> dict[str, float]:
    """Elemento → numero di atomi nella cella (pesati con l'occupanza). Con `include_omitted` conta anche le componenti
    disordinate scartate dalla regola di disordine dell'analisi. Una posizione già occupata da un altro atomo dello
    stesso elemento (entro `tol` Å: un CIF che elenca anche atomi equivalenti per simmetria) non si conta due volte;
    le componenti disordinate stanno sempre più lontane di così."""
    vectors = cell_vectors(cif.cell)
    out: dict[str, float] = {}
    placed: dict[str, np.ndarray] = {}
    for atom in cif.atoms + (cif.omitted if include_omitted else []):
        prev = placed.get(atom.el, np.zeros((0, 3)))
        new = []
        for op in cif.symops:
            f = op.apply(atom.frac)
            f = f - np.floor(f + 1e-9)
            others = np.vstack([prev, *new]) if new else prev
            d = f - others
            d -= np.round(d)
            if not len(d) or (np.linalg.norm(d @ vectors, axis=1) >= tol).all():
                new.append(f[None])
        if new:
            placed[atom.el] = np.vstack([prev, *new])
        out[atom.el] = out.get(atom.el, 0.0) + atom.occ * len(new)
    return out


def parse_formula(text: str | None) -> dict[str, float]:
    """'C14 H16 N6 O2' → {'C': 14, 'H': 16, 'N': 6, 'O': 2}; numeri decimali ammessi."""
    out: dict[str, float] = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d+(?:\.\d+)?|\.\d+)?", text or ""):
        if el == "D":
            el = "H"
        out[el] = out.get(el, 0.0) + (float(n) if n else 1.0)
    return out


def reported_formula(cif: Cif) -> dict[str, float]:
    """La formula per unità formula dichiarata nel CIF (`_chemical_formula_sum`), vuota se manca."""
    return parse_formula(cif.items.get("_chemical_formula_sum"))


def formula_text(counts: dict[str, float]) -> str:
    """{'C': 14, 'H': 16} → 'C14 H16', in ordine di Hill (C, H, poi alfabetico); 2 decimali se non interi."""
    order = sorted(counts, key=lambda e: (e != "C", e != "H", e)) if "C" in counts else sorted(counts)

    def num(n: float) -> str:
        return "" if abs(n - 1) < 1e-6 else (f"{n:.0f}" if abs(n - round(n)) < 0.005 else f"{n:.2f}")

    return " ".join(f"{e}{num(counts[e])}" for e in order if counts[e] > 1e-6)


def density(contents: dict[str, float], volume: float) -> float:
    """Densità in g/cm³ dal contenuto della cella (NaN se un elemento non ha massa nota)."""
    mass = sum(n * atomic_mass(e) for e, n in contents.items())
    return mass * U_PER_A3_TO_G_PER_CM3 / volume if volume > 0 else math.nan


def z_prime(cif: Cif) -> float:
    """Z′ = Z / numero di posizioni generali (le operazioni di simmetria della cella); NaN senza Z."""
    z = cif_value(cif, "z")
    return z / len(cif.symops) if math.isfinite(z) and cif.symops else math.nan
