"""Cristallo dimostrativo, costruito a tavolino e NON sperimentale: un naftalene ideale
(C–C 1,40 Å, C–H 1,08 Å) in una cella P-1 scelta in modo che le molecole impilate lungo a
abbiano distanza interplanare 3,40 Å, slittamento longitudinale 1,50 Å e laterale 0,80 Å.
Serve a provare lo strumento senza un file proprio e come caso di prova nei test.

Porting di `src/lib/crystal/demo.ts`.
"""
from __future__ import annotations

import math

import numpy as np

DEMO_STACK = {"interplanar": 3.4, "long_slip": 1.5, "lat_slip": 0.8}


def naphthalene() -> list[tuple[str, np.ndarray]]:
    """Atomi del naftalene nel piano xy, asse lungo x, centro di inversione nell'origine."""
    cc, ch = 1.4, 1.08
    h = cc * math.cos(math.pi / 6)  # distanza dal centro dell'anello al lato condiviso
    out: list[tuple[str, np.ndarray]] = []
    for side in (1, -1):
        cx = side * h
        for k in range(6):
            ang = math.pi / 6 + k * math.pi / 3
            p = np.array([cx + cc * math.cos(ang), cc * math.sin(ang), 0.0])
            shared = abs(p[0]) < 1e-9
            if shared and side == -1:
                continue  # gli atomi del lato comune una volta sola
            out.append(("C", p))
            if not shared:
                out.append(("H", np.array([p[0] + ch * math.cos(ang), p[1] + ch * math.sin(ang), 0.0])))
    return out


def _deg(u: np.ndarray, v: np.ndarray) -> float:
    return math.degrees(math.acos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))))


def demo_cif() -> str:
    # vettori di cella nel sistema della molecola (L = x, W = y, normale = z)
    a = np.array([DEMO_STACK["long_slip"], DEMO_STACK["lat_slip"], DEMO_STACK["interplanar"]])
    b = np.array([0.0, 7.5, 0.0])
    c = np.array([-11.0, 0.0, 0.0])  # verso scelto perché la terna a, b, c sia destrorsa
    vol = np.dot(a, np.cross(b, c))
    inv = [np.cross(b, c) / vol, np.cross(c, a) / vol, np.cross(a, b) / vol]
    f = lambda x: f"{x:.5f}"  # noqa: E731
    # unità asimmetrica: metà molecola, l'altra metà la genera il centro di inversione
    asym = [(el, p) for el, p in naphthalene() if p[0] > 1e-9 or (abs(p[0]) < 1e-9 and p[1] > 0)]
    n_c = n_h = 0
    rows = []
    for el, p in asym:
        if el == "C":
            n_c += 1
            label = f"C{n_c}"
        else:
            n_h += 1
            label = f"H{n_h}"
        rows.append(f"{label} {el} " + " ".join(f(np.dot(r, p)) for r in inv))
    return "\n".join([
        "# Struttura dimostrativa costruita a mano: non è una misura sperimentale.",
        "data_demo_naphthalene_stack",
        "_symmetry_space_group_name_H-M 'P -1'",
        f"_cell_length_a {f(np.linalg.norm(a))}",
        f"_cell_length_b {f(np.linalg.norm(b))}",
        f"_cell_length_c {f(np.linalg.norm(c))}",
        f"_cell_angle_alpha {f(_deg(b, c))}",
        f"_cell_angle_beta {f(_deg(a, c))}",
        f"_cell_angle_gamma {f(_deg(a, b))}",
        "loop_",
        "_symmetry_equiv_pos_as_xyz",
        "'x, y, z'",
        "'-x, -y, -z'",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
        *rows,
        "",
    ])
