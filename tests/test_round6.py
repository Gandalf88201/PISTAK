"""v1.4.0 (referee report of 13 November 2026, issue 1): a ring enters a π fragment only if its non-metal atoms can be
sp² or sp (`Limits.ring_atoms = "sp2"`, `cores.sp2_ring`); `"any"` keeps the planarity-only rule of v1.3.1."""
from __future__ import annotations

import numpy as np

from cifgen import make_cif, ring_system
from pistak.core.analysis import analyse_text
from pistak.core.cif import Cell
from pistak.core.config import Limits
from pistak.core.cores import MIN_SP2_ANGLE_SUM

CELL = Cell(20.0, 20.0, 20.0, 90.0, 90.0, 90.0)


def cif(atoms, name):
    return make_cif(atoms, CELL, ["x, y, z"], name, "P 1", origin_frac=(0.5, 0.5, 0.5))


def planar_cyclobutane() -> list:
    """A square C₄ ring (planar, as on an inversion centre) with two H atoms on each carbon, above and below."""
    atoms = []
    for x, y in ((0.775, 0.775), (-0.775, 0.775), (-0.775, -0.775), (0.775, -0.775)):
        c = np.array([x, y, 0.0])
        out = c / np.linalg.norm(c) * 0.63
        atoms += [("C", c), ("H", c + out + [0, 0, 0.89]), ("H", c + out - [0, 0, 0.89])]
    return atoms


def benzocyclobutene() -> list:
    """Benzene fused to a planar four-membered ring whose two outer carbons are CH₂ (sp³)."""
    benz = ring_system(1)
    p1, p2 = np.array([1.2124, -0.7, 0.0]), np.array([1.2124, 0.7, 0.0])  # the fused bond (angles 330° and 30°)
    atoms = [(e, p) for e, p in benz if not (e == "H" and min(np.linalg.norm(p - p1), np.linalg.norm(p - p2)) < 1.2)]
    for y in (0.7, -0.7):
        c = np.array([2.7624, y, 0.0])
        out = np.array([0.54, np.sign(y) * 0.4, 0.0])
        atoms += [("C", c), ("H", c + out + [0, 0, 0.89]), ("H", c + out - [0, 0, 0.89])]
    return atoms


def test_planar_saturated_ring_is_not_a_pi_fragment():
    text = cif(planar_cyclobutane(), "cb")
    assert analyse_text(text, "cb", Limits()).fragments[0] == []
    old = analyse_text(text, "cb", Limits(ring_atoms="any")).fragments[0]
    assert len(old) == 1 and len(old[0].atoms) == 4


def test_aromatic_ring_fused_to_a_saturated_ring_keeps_only_the_aromatic_ring():
    text = cif(benzocyclobutene(), "bcb")
    (f,) = analyse_text(text, "bcb", Limits()).fragments[0]
    assert len(f.atoms) == 6 and len(f.rings) == 1
    (g,) = analyse_text(text, "bcb", Limits(ring_atoms="any")).fragments[0]
    assert len(g.atoms) == 8 and len(g.rings) == 2


def test_aromatic_rings_are_kept():
    for n in (1, 2, 3):
        (f,) = analyse_text(cif(ring_system(n), f"acene{n}"), "a", Limits()).fragments[0]
        assert len(f.atoms) == 4 * n + 2


def test_pyramidal_ring_atom_excludes_the_ring():
    # benzene with one ring carbon bonded to a substituent 50° out of the ring plane: angle sum about 318°
    atoms = [(e, p) for e, p in ring_system(1) if not (e == "H" and p[0] > 2.0)]
    c = np.array([1.4, 0.0, 0.0])
    atoms.append(("Cl", c + 1.74 * np.array([np.cos(np.radians(50)), 0.0, np.sin(np.radians(50))])))
    assert analyse_text(cif(atoms, "pyr"), "pyr", Limits()).fragments[0] == []
    assert MIN_SP2_ANGLE_SUM == 345.0


def test_site_mode_keeps_planarity_only():
    assert Limits.site().ring_atoms == "any" and Limits().ring_atoms == "sp2"
