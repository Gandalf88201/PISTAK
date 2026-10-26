"""Costruzione di cristalli di prova (NON sperimentali) in formato CIF.

Si dà una molecola in coordinate cartesiane, una cella e le operazioni di simmetria: il CIF contiene
l'intera molecola come unità asimmetrica e la simmetria genera il resto. Serve a fabbricare casi dalla
geometria nota (π-stack, herringbone, reti infinite…) senza dipendere da file esterni.
"""
from __future__ import annotations

import math

import numpy as np

from pistak.core.cif import Cell
from pistak.core.structure import cell_vectors


def ring_system(n_rings: int, cc: float = 1.4, ch: float = 1.08, with_h: bool = True) -> list[tuple[str, np.ndarray]]:
    """Acene lineare (1 = benzene, 2 = naftalene, 3 = antracene…) nel piano xy, asse lungo x, centrato nell'origine."""
    h = cc * math.cos(math.pi / 6)
    centers = [(k - (n_rings - 1) / 2) * 2 * h for k in range(n_rings)]
    sites: list[list] = []  # [posizione, angolo di uscita dell'idrogeno, quante volte compare]
    for cx in centers:
        for k in range(6):
            ang = math.pi / 6 + k * math.pi / 3
            p = np.array([cx + cc * math.cos(ang), cc * math.sin(ang), 0.0])
            for s in sites:
                if np.linalg.norm(s[0] - p) < 1e-6:
                    s[2] += 1
                    break
            else:
                sites.append([p, ang, 1])
    out: list[tuple[str, np.ndarray]] = []
    for p, ang, n in sites:
        out.append(("C", p))
        if with_h and n == 1:  # gli atomi di fusione (in due anelli) non hanno idrogeno
            out.append(("H", np.array([p[0] + ch * math.cos(ang), p[1] + ch * math.sin(ang), 0.0])))
    return out


def rotation(axis: str, deg: float) -> np.ndarray:
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    return {"x": np.array([[1, 0, 0], [0, c, -s], [0, s, c]]), "y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]),
            "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[axis]


def make_cif(atoms: list[tuple[str, np.ndarray]], cell: Cell, symops: list[str], name: str = "test",
             space_group: str = "?", origin_frac=(0.0, 0.0, 0.0)) -> str:
    """CIF con la molecola data (coordinate cartesiane, già posizionata) come unità asimmetrica."""
    m = cell_vectors(cell)
    inv = np.linalg.inv(m)
    counts: dict[str, int] = {}
    rows = []
    for el, p in atoms:
        counts[el] = counts.get(el, 0) + 1
        f = p @ inv + np.asarray(origin_frac)
        rows.append(f"{el}{counts[el]} {el} {f[0]:.6f} {f[1]:.6f} {f[2]:.6f}")
    sym = "\n".join(f"'{s}'" for s in symops)
    return "\n".join([
        f"data_{name}", f"_symmetry_space_group_name_H-M '{space_group}'",
        f"_cell_length_a {cell.a:.5f}", f"_cell_length_b {cell.b:.5f}", f"_cell_length_c {cell.c:.5f}",
        f"_cell_angle_alpha {cell.alpha:.4f}", f"_cell_angle_beta {cell.beta:.4f}", f"_cell_angle_gamma {cell.gamma:.4f}",
        "loop_", "_symmetry_equiv_pos_as_xyz", sym,
        "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        *rows, "",
    ])


P21C = ["x, y, z", "-x, y+1/2, -z+1/2", "-x, -y, -z", "x, -y+1/2, z+1/2"]


def _rot(ax: float, ay: float, az: float) -> np.ndarray:
    return rotation("z", az) @ rotation("y", ay) @ rotation("x", ax)


def anthracene_p21a_like() -> str:
    """Antracene in P21/c sull'origine (centro di inversione): due orientazioni, una coppia π-stack lungo
    una direzione e coppie quasi perpendicolari (lisca di pesce). Parametri trovati con una ricerca casuale."""
    R = _rot(11.525658880859515, 136.48144433162713, 106.39792492763716)
    mol = [(el, R @ p) for el, p in ring_system(3)]
    cell = Cell(8.602535319031425, 5.5620235029395, 12.596581710936837, 90.0, 111.81872721663667, 90.0)
    return make_cif(mol, cell, P21C, "herringbone", "P 21/c")


def pyrene_like_dimers() -> str:
    """Naftalene ruotato in P-1: due molecole sovrapposte con slittamento, più una seconda coppia in cella."""
    mol = [(el, rotation("y", 12) @ rotation("z", 40) @ p) for el, p in ring_system(2)]
    cell = Cell(6.2, 7.9, 9.4, 83.0, 97.0, 101.0)
    return make_cif(mol, cell, ["x, y, z", "-x, -y, -z"], "dimers", "P -1", origin_frac=(0.1, 0.2, 0.15))


def tilted_sandwich() -> str:
    """Seconda disposizione di antracene in P21/c, con coppie a ~45° e un contatto di stack molto stretto."""
    R = _rot(55.74060962574102, 13.854726846974142, 107.95730558338813)
    mol = [(el, R @ p) for el, p in ring_system(3)]
    cell = Cell(8.062755524350635, 5.894769712856839, 11.223808406850972, 90.0, 115.26167807418354, 90.0)
    return make_cif(mol, cell, P21C, "tilted", "P 21/c")


def infinite_chain() -> str:
    """Un atomo di carbonio per cella, a 1.6 Å dai suoi traslati lungo a: catena infinita, nessuna molecola."""
    cell = Cell(1.6, 8.0, 8.0, 90.0, 90.0, 90.0)
    return make_cif([("C", np.zeros(3))], cell, ["x, y, z"], "chain", "P 1", origin_frac=(0.0, 0.5, 0.5))


def chain_molecule() -> str:
    """Esano planare a zig-zag senza anelli: il nucleo è l'intera molecola (kind = molecule)."""
    pts = [(el, np.array([i * 1.25, 0.7 * (i % 2), 0.0])) for i, el in enumerate("CCCCCC")]
    cell = Cell(20.0, 6.0, 4.0, 90.0, 90.0, 90.0)
    return make_cif(pts, cell, ["x, y, z", "-x, -y, -z"], "zigzag", "P -1", origin_frac=(0.1, 0.3, 0.1))


ALL = {
    "herringbone": anthracene_p21a_like,
    "dimers": pyrene_like_dimers,
    "tilted": tilted_sandwich,
    "chain_molecule": chain_molecule,
    "infinite_chain": infinite_chain,
}
