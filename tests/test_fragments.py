"""Phase 3: several π fragments per molecule, fragment–fragment pairs, per-molecule summary, user-defined long axis,
signed slips."""
import json
import math
from pathlib import Path

import numpy as np
import pytest

import cifgen
from pistak.core.analysis import analyse_file, analyse_text
from pistak.core.cif import Cell
from pistak.core.config import Limits
from pistak.core.cores import pi_core, pi_fragments
from pistak.core.pairs import class_spread
from pistak.core.structure import cell_vectors
from pistak.labels import relation
from pistak.params import compute_all

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"


def naphthyl_benzene(twist: float = 50.0) -> list[tuple[str, np.ndarray]]:
    """2-phenylnaphthalene: naphthalene in the xy plane plus a benzene bonded to C2 and turned by `twist` degrees about
    the bond. Two separate (not fused) ring systems: two π fragments, of 10 and 6 atoms."""
    naph = cifgen.ring_system(2)
    carbons = [p for el, p in naph if el == "C"]
    ipso_n = max(carbons, key=lambda p: (round(p[0], 6), p[1]))  # C2 (largest x, positive y)
    u = np.array([math.cos(math.pi / 6), math.sin(math.pi / 6), 0.0])
    # remove the C2 hydrogen: the benzene goes in its place
    atoms = [(el, p) for el, p in naph if not (el == "H" and np.linalg.norm(p - ipso_n - 1.08 * u) < 1e-6)]
    v = np.array([-u[1], u[0], 0.0])
    w = np.cross(u, v)
    t = math.radians(twist)
    v2 = math.cos(t) * v + math.sin(t) * w
    center = ipso_n + (1.48 + 1.4) * u
    for k in range(6):
        ang = math.pi + k * math.pi / 3  # k = 0: the ipso carbon, towards the naphthalene
        d = math.cos(ang) * u + math.sin(ang) * v2
        atoms.append(("C", center + 1.4 * d))
        if k:
            atoms.append(("H", center + 2.48 * d))
    return atoms


def stacked_crystal(twist: float = 50.0) -> str:
    """P-1 with the molecule stacked by translation along a (≈3.75 Å, nearly perpendicular to the naphthalene)."""
    mine = np.array([[1.2, 0.6, 3.5], [0.0, 16.0, 0.0], [-30.0, 0.0, 0.0]])
    ln = np.linalg.norm(mine, axis=1)
    ang = lambda i, j: math.degrees(math.acos(mine[i] @ mine[j] / (ln[i] * ln[j])))  # noqa: E731
    cell = Cell(ln[0], ln[1], ln[2], ang(1, 2), ang(0, 2), ang(0, 1))
    to_cell = np.linalg.inv(mine) @ cell_vectors(cell)  # same fractional coordinates in the cell_vectors frame
    if np.linalg.det(to_cell) < 0:
        raise AssertionError("left-handed axes")
    mol = [(el, p @ to_cell) for el, p in naphthyl_benzene(twist)]
    return cifgen.make_cif(mol, cell, ["x, y, z", "-x, -y, -z"], "phenylnaphthalene", "P -1", origin_frac=(0.25, 0.25, 0.25))


@pytest.fixture(scope="module")
def crystal():
    return stacked_crystal()


# ---- fragments -------------------------------------------------------------------------------------------
def test_fragments_of_a_molecule(crystal):
    a = analyse_text(crystal, "fn", Limits())
    m = a.structure.molecules[0]
    assert [len(f.atoms) for f in a.fragments[0]] == [10, 6]
    assert pi_core(m, a.limits).atoms.tolist() == a.fragments[0][0].atoms.tolist()  # the main one is the usual core
    assert len(pi_fragments(m, a.limits.with_values(fragments="largest"))) == 1
    assert len(pi_fragments(m, a.limits.with_values(min_fragment_atoms=7))) == 1
    g = a.molecule_groups()[0]
    assert g.fragments == 2 and g.with_core == 2


def test_fragment_fragment_pairs(crystal):
    full = compute_all(analyse_text(crystal, "fn", Limits()))
    largest = analyse_text(crystal, "fn", Limits(fragments="largest"))
    assert all((p.fa, p.fb) == (0, 0) for p in largest.pairs)
    combos = {(p.fa, p.fb) for p in full.pairs if not p.intramolecular}
    assert {(0, 0), (1, 1)} <= combos
    # the naphthalene–naphthalene stack along a is there in both modes, with the same values
    n_full = next(p for p in full.pairs if (p.fa, p.fb) == (0, 0) and p.kind == "stack")
    n_big = next(p for p in largest.pairs if p.kind == "stack")
    assert n_full.interplanar == pytest.approx(n_big.interplanar) == pytest.approx(3.5, abs=1e-4)
    # phenyl–phenyl along a: parallel planes (same molecule translated)
    ph = next(p for p in full.pairs if (p.fa, p.fb) == (1, 1) and p.dihedral < 1)
    assert ph.get("fragment_a") == 2 and ph.get("fragment_b") == 2 and "π2–π2" in relation(ph)
    for p in full.pairs:
        assert class_spread(full.structure, full.fragments, full.limits, p) < 1e-6


def test_intramolecular_pairs(crystal):
    a = analyse_text(crystal, "fn", Limits())
    intra = [p for p in a.pairs if p.intramolecular]
    assert len(intra) == 1  # π1–π2 in the same molecule; the two molecules of the cell are equivalent
    p = intra[0]
    assert (p.fa, p.fb) == (0, 1) and p.count == 2 and p.kind == "tilted"
    assert p.dihedral == pytest.approx(50.0, abs=1e-3) and relation(p).startswith("intramolecolare · π1–π2")


def test_summary_per_molecule(crystal):
    a = analyse_text(crystal, "fn", Limits())
    contacts = a.molecule_contacts()
    assert sum(len(c.pairs) for c in contacts) == len(a.pairs)
    along_a = next(c for c in contacts if not c.intramolecular and c.a == c.b and tuple(abs(x) for x in c.shift) == (1, 0, 0))
    assert len(along_a.pairs) >= 2 and along_a.count == 2  # one per molecule of the cell (Z = 2)
    assert along_a.best.kind == "stack" and along_a.best.fa == along_a.best.fb == 0
    intra = [c for c in contacts if c.intramolecular]
    assert len(intra) == 1 and intra[0].count == 2


def test_json_with_fragments_and_contacts(crystal, tmp_path):
    from pistak.io import export

    a = compute_all(analyse_text(crystal, "fn", Limits()))
    data = json.loads(export.write_json(a, tmp_path / "a.json").read_text())
    assert data["version"] == 5
    assert [(f["molecule"], f["fragment"], f["atoms"]) for f in data["fragments"]] == [(1, 1, 10), (1, 2, 6), (2, 1, 10), (2, 2, 6)]
    assert {p["fragment_a"] for p in data["pairs"]} == {1, 2}
    assert sum(len(c["pairs"]) for c in data["molecule_contacts"]) == len(data["pairs"])


def test_single_fragment_examples_unchanged():
    for cif in sorted(p for p in EX.glob("*.cif") if not p.name.startswith(".")):
        full, largest = analyse_file(cif, Limits()), analyse_file(cif, Limits(fragments="largest"))
        assert [(p.kind, round(p.interplanar, 9), p.count) for p in full.pairs] == \
               [(p.kind, round(p.interplanar, 9), p.count) for p in largest.pairs]


# ---- user-defined long axis and signed slips --------------------------------------------------------------
def _short_axis_labels(a) -> tuple[str, str]:
    """Two core atoms aligned along the short axis (same L coordinate, opposite W)."""
    m, core = a.structure.molecules[0], a.cores[0]
    W = np.cross(core.normal, core.long)
    best = None
    for i in core.atoms:
        for j in core.atoms:
            d = m.pos[j] - m.pos[i]
            score = abs(d @ W) - 10 * abs(d @ core.long)
            if best is None or score > best[0]:
                best = (score, i, j)
    return m.labels[best[1]], m.labels[best[2]]


def test_long_axis_from_atoms():
    cif = EX / "herringbone_anthracene.cif"
    auto = compute_all(analyse_file(cif, Limits()))
    i, j = _short_axis_labels(auto)
    turned = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{i},{j}")))
    assert all(c.axis_source == "atoms" for c in turned.cores)
    p0, p1 = auto.pairs[0], turned.pairs[0]
    assert p0.kind == p1.kind == "stack" and p1.interplanar == pytest.approx(p0.interplanar)
    # the «long» axis is now the short one: the two slips swap
    assert p1.long_slip == pytest.approx(p0.lat_slip, abs=1e-4) and p1.lat_slip == pytest.approx(p0.long_slip, abs=1e-4)
    # labels that do not exist: the automatic axis stays
    assert all(c.axis_source == "pca" for c in analyse_file(cif, Limits(long_axis_atoms="X98,X99")).cores)


def test_swapping_axis_atoms_flips_sign():
    # molecule off inversion centres (each label once per molecule): the axis direction is defined
    cif = EX / "dimers_naphthalene.cif"
    auto = analyse_file(cif, Limits())
    i, j = _short_axis_labels(auto)
    fwd = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{i},{j}")))
    rev = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{j},{i}")))
    for p, q in zip(fwd.pairs, rev.pairs):
        assert q.get("long_slip_signed") == pytest.approx(-p.get("long_slip_signed"), abs=1e-9)
        assert q.long_slip == pytest.approx(p.long_slip, abs=1e-9)


def test_signed_slips():
    a = compute_all(analyse_file(EX / "dimers_naphthalene.cif", Limits()))
    for p in a.pairs:
        sl, sw = p.get("long_slip_signed"), p.get("lat_slip_signed")
        assert abs(sl) == pytest.approx(p.long_slip, abs=1e-9) and abs(sw) == pytest.approx(p.lat_slip, abs=1e-9)
        assert p.get("slip_direction") == pytest.approx(math.degrees(math.atan2(sw, sl)), abs=1e-9)
        assert 0 <= p.get("axis_aspect_a") <= 1


def test_core_isotropy():
    naph = analyse_file(EX / "demo_naphthalene.cif", Limits()).cores[0].plane.aspect
    anth = analyse_file(EX / "herringbone_anthracene.cif", Limits()).cores[0].plane.aspect
    assert 0 < anth < naph < 1  # anthracene is more elongated than naphthalene


def test_axis_from_thresholds():
    assert Limits.from_dict({"long_axis_atoms": " C1,C9 "}).long_axis_atoms == "C1,C9"
    assert Limits.from_dict({"long_axis_atoms": ["C1"]}).long_axis_atoms == ""
