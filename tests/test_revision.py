"""Changes made for the referee report on the paper (REGISTRO items 50–56): π fragments of organometallics and
phthalocyanines, the conjugated-merge option, overlap on the union of ring polygons, undefined long axes, twist,
planes crossing between the centroids, local interplanar distance, and the linear fallback.

The CIFs in tests/data/cod are from the Crystallography Open Database (public domain)."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from pistak import api
from pistak.core.config import Limits
from pistak.core.cores import Core, pi_fragments
from pistak.core.geometry import fit_plane, polygon_area
from pistak.core.pairs import axis_defined, fragment_polygons, polygons_overlap
from pistak.params import compute_all

COD = Path(__file__).resolve().parent / "data" / "cod"


def _frags(cid: str, **kw):
    a = api.analyse_file(COD / f"{cid}.cif", Limits(**kw))
    m = a.structure.molecules
    return a, [[sorted(m[i].el[j] for j in f.atoms) for f in fs] for i, fs in enumerate(a.fragments)]


@pytest.mark.parametrize("cid, ring", [("7033930", 5), ("4103814", 6), ("2101898", 6)])
def test_metal_bound_rings_are_found(cid, ring):
    """Ferrocene-type metallocene (Cp), Cr(η⁶-arene)(CO)₃ and Ag(C₆H₆)ClO₄: the M–C–C triangles hid the aromatic
    ring from the shortest-cycle search; now the ring is the fragment and the metal is not in it."""
    _, frags = _frags(cid)
    sizes = {len(f) for fs in frags for f in fs}
    assert ring in sizes
    assert all(set(f) == {"C"} for fs in frags for f in fs if len(f) == ring)


def test_phthalocyanines_free_base_and_metal():
    """β-H₂Pc: four isoindoles with the default (fused rings), one 40-atom fragment with `fragment_merge =
    "conjugated"`, stacked along b (4.73 Å); CuPc (P1̄, Z = 1): one 41-atom fragment with Cu either way (the chelate rings fuse
    the macrocycle), stacked by translation (3.77 Å)."""
    a, f = _frags("2008099")
    assert sorted(len(x) for x in f[0]) == [9, 9, 9, 9]
    a, f = _frags("2008099", fragment_merge="conjugated")
    assert [len(x) for x in f[0]] == [40]
    st = [p for p in a.pairs if p.kind == "stack"]
    assert len(st) == 1 and st[0].centroid_distance == pytest.approx(4.73, abs=0.01) and [c.rank for c in a.stacks] == [1, 1]
    assert st[0].overlap_area > 0.3  # the macrocycle interior counts (16-membered inner ring)
    for merge in ("fused", "conjugated"):
        a, f = _frags("2103883", fragment_merge=merge)
        assert [len(x) for x in f[0]] == [41] and "Cu" in f[0][0]
        st = [p for p in a.pairs if p.kind == "stack"]
        assert len(st) == 1 and st[0].centroid_distance == pytest.approx(3.77, abs=0.01)


def _flat(atoms, bonds):
    """A Molecule-like object with positions in the plane z = 0 (or as given)."""
    from pistak.core.structure import Molecule

    el = [e for e, *_ in atoms]
    pos = np.array([p for _, *p in atoms], dtype=float)
    return Molecule(el=el, labels=[f"{e}{k + 1}" for k, e in enumerate(el)], pos=pos,
                    bonds=np.array(bonds, dtype=int).reshape(-1, 2), centroid=pos.mean(axis=0), extended=False)


def _ring(n, r, cx=0.0, cy=0.0, z=0.0, phase=0.0):
    return [(cx + r * math.cos(phase + 2 * math.pi * k / n), cy + r * math.sin(phase + 2 * math.pi * k / n), z) for k in range(n)]


def _biphenyl(twist_deg: float):
    """Two benzene rings joined by a 1.48 Å C–C bond along x, the second ring rotated about x by `twist_deg`."""
    a = _ring(6, 1.39)
    t = math.radians(twist_deg)
    b = [(x + 2 * 1.39 + 1.48, y * math.cos(t), y * math.sin(t)) for x, y, _ in a]
    atoms = [("C", *p) for p in a + b]
    bonds = [(k, (k + 1) % 6) for k in range(6)] + [(6 + k, 6 + (k + 1) % 6) for k in range(6)] + [(0, 9)]
    return _flat(atoms, bonds)


def test_conjugated_merge_coplanar_and_twisted_biaryl():
    lim = Limits()
    flat, twisted = _biphenyl(0.0), _biphenyl(40.0)
    assert [len(f.atoms) for f in pi_fragments(flat, lim.with_values(fragments="all"))] == [6, 6]
    conj = lim.with_values(fragments="all", fragment_merge="conjugated")
    assert [len(f.atoms) for f in pi_fragments(flat, conj)] == [12]
    assert [len(f.atoms) for f in pi_fragments(twisted, conj)] == [6, 6]  # a rotated biaryl stays in two fragments


def test_linear_molecule_has_no_fallback_plane():
    """A linear heavy-atom skeleton (no rings) has no defined normal: no 'molecule' fragment."""
    m = _flat([("C", 1.2 * k, 0.0, 0.0) for k in range(6)], [(k, k + 1) for k in range(5)])
    on = Limits(ringless_fragments="planar")
    assert pi_fragments(m, on) == []
    bent = _flat([("C", 1.2 * k, 0.6 * (k % 2), 0.0) for k in range(6)], [(k, k + 1) for k in range(5)])
    assert [f.kind for f in pi_fragments(bent, on)] == ["molecule"]
    assert pi_fragments(bent, Limits()) == []  # off by default since round 4


def test_union_of_rings_versus_hull():
    """Naphthalene: the union of the two hexagons is smaller than the convex hull by the two notches at the ring
    fusion; the overlap of a fragment with itself is its own area."""
    h = 1.4 * math.sqrt(3) / 2
    a = _ring(6, 1.4, -h, 0.0, phase=math.pi / 6)
    b = _ring(6, 1.4, h, 0.0, phase=math.pi / 6)
    pts = np.array(a + b)
    # merge the two shared atoms
    keep, idx = [], []
    for p in pts:
        j = next((k for k, q in enumerate(keep) if np.linalg.norm(p - q) < 1e-3), None)
        if j is None:
            keep.append(p)
            j = len(keep) - 1
        idx.append(j)
    keep = np.array(keep)
    rings = (tuple(idx[:6]), tuple(idx[6:]))
    core = Core(fit_plane(keep), np.arange(len(keep)), "rings", rings=rings)
    q = keep[:, :2]
    union = fragment_polygons(core, q)
    hull = fragment_polygons(core, q, hull_only=True)
    hexagon = 3 * math.sqrt(3) / 2 * 1.4 ** 2
    assert sum(polygon_area(p) for p in union) == pytest.approx(2 * hexagon, rel=1e-6)
    assert polygon_area(hull[0]) > 2 * hexagon * 1.05
    assert polygons_overlap(union, union) == pytest.approx(2 * hexagon, rel=1e-6)


def test_isotropic_core_has_no_long_axis_and_slips_are_nan():
    """Benzene, coronene: no long axis (aspect ≈ 1), so s_L and s_W are NaN; the total slip is still defined."""
    a = api.analyse_file(COD.parent.parent.parent / "examples" / "dimers_naphthalene.cif", Limits())
    assert all(axis_defined(f, a.limits) for fs in a.fragments for f in fs)  # naphthalene: defined
    benz = Core(fit_plane(np.array(_ring(6, 1.39))), np.arange(6), "rings", rings=(tuple(range(6)),))
    assert not axis_defined(benz, Limits())
    a = compute_all(api.analyse_file(COD / "2101898.cif", Limits(overlap_measure="atoms", min_overlap=0.0)))
    for p in a.pairs:
        assert math.isnan(p.long_slip) and math.isnan(p.lat_slip) and math.isfinite(p.slip)
        assert math.isnan(p.twist)


def test_twist_planes_cross_and_local_interplanar():
    a = compute_all(api.analyse_file(COD / "2008099.cif", Limits()))  # H₂Pc, isoindole fragments: anisotropic
    st = next(p for p in a.pairs if p.kind == "stack")
    # the stacked isoindoles are different fragments (fa ≠ fb), crossed at ~87°: above max_axis_twist (75°) the
    # long-axis frame is meaningless and s_L, s_W are NaN; the total slip stays
    assert st.fa != st.fb and 80 < st.twist <= 90 and st.get("twist") == pytest.approx(st.twist)
    assert math.isnan(st.long_slip) and math.isnan(st.lat_slip) and math.isfinite(st.slip)
    same = [p for p in a.pairs if p.fa == p.fb and p.dihedral < 1 and not p.intramolecular]
    assert same and all(p.twist == pytest.approx(0.0, abs=1e-6) for p in same)  # translation-related: parallel axes
    a = compute_all(api.analyse_file(COD / "2103883.cif", Limits()))  # CuPc: D4h, no long axis
    st = next(p for p in a.pairs if p.kind == "stack")
    assert math.isnan(st.twist) and math.isnan(st.long_slip)
    assert st.get("planes_cross") == 0.0
    assert st.get("local_interplanar") == pytest.approx(st.interplanar, abs=0.05)
    for p in a.pairs:
        assert p.planes_cross == (p.get("planes_cross") == 1.0)


# ---- second round of review (REGISTRO items 58–62) --------------------------------------------------------------
def test_metallocene_rings_are_bridged_not_stacked():
    """Ferrocene: the two Cp rings are not bonded to each other but both to Fe: 'bridged', never a π-stack, although
    they are parallel, 3.3 Å apart and overlap almost completely."""
    a = api.analyse_file(COD / "7033930.cif", Limits())
    intra = [p for p in a.pairs if p.intramolecular]
    assert len(intra) == 1
    p = intra[0]
    assert p.bridged and not p.bonded and p.kind == "parallel"
    assert p.interplanar == pytest.approx(3.31, abs=0.01) and p.overlap_area > 0.8


def test_structures_without_organic_bonds_have_no_fragments():
    """No C–C or C–H bond and no planar non-metal ring with B, C or N (here a planar Si₃O₃ ring, as in a
    cyclosilicate): no π fragment by default; pi_scope = 'all' restores the search (and Limits.site()). A B₃N₃ ring
    (borazine type) is in scope since round 3 (N1): see `test_triazines_and_borazines_are_in_scope`."""
    from pistak.core.analysis import has_organic_bonds, in_pi_scope

    text = """data_k
_cell_length_a 7
_cell_length_b 7
_cell_length_c 7
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
"""
    # a planar Si3O3 ring (Si–O 1.64 Å, cell edge 7 Å: radius 1.64 Å = 0.2343 in fractions): planar ring, no B, C or N
    import math as _m
    for k in range(6):
        el = "Si" if k % 2 == 0 else "O"
        text += f"{el}{k + 1} {el} {0.5 + 0.2343 * _m.cos(_m.pi * k / 3):.5f} {0.5 + 0.2343 * _m.sin(_m.pi * k / 3):.5f} 0.5\n"
    a = api.analyse_text(text, "sio", Limits())
    assert not has_organic_bonds(a.structure) and not in_pi_scope(a.structure, a.limits)
    assert all(not f for f in a.fragments)
    b = api.analyse_text(text, "sio", Limits(pi_scope="all"))
    assert any(f for f in b.fragments)


def test_split_positions_are_not_merged():
    """Two partially occupied positions of the same atom 0.2 Å apart are kept (split site); symmetry images of one
    site are still merged."""
    base = """data_s
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P -1'
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
-x,-y,-z
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_disorder_group
Br1A Br 0.30 0.30 0.30 0.5 1
Br1B Br 0.32 0.30 0.30 0.5 2
Cu1 Cu 0 0 0 1 .
"""
    a = api.analyse_text(base, "split", Limits(disorder="all"))
    assert sorted(a.structure.cell_el).count("Br") == 4  # two sites × inversion
    assert sorted(a.structure.cell_el).count("Cu") == 1  # on the inversion centre: images merged
    # the same atom written twice at half occupancy (its inversion image listed as another atom, as in COD 2100449):
    # the two coincide and are merged, whatever their occupancy
    dup = base.replace("Cu1 Cu 0 0 0 1 .\n", "Cu1 Cu 0 0 0 1 .\nS1 S 0.10 0.20 0.30 0.5 .\nS1i S -0.10 -0.20 -0.30 0.5 .\n")
    a = api.analyse_text(dup, "dup", Limits(disorder="all"))
    assert sorted(a.structure.cell_el).count("S") == 2


def test_monte_carlo_keeps_special_position_constraints():
    """An atom on a mirror (x, y, 0) stays on it in every perturbed copy: the displacement is projected onto the
    site symmetry."""
    from pistak.core.uncertainty import perturbed, site_projectors

    text = """data_m
_cell_length_a 7
_cell_length_b 7
_cell_length_c 7
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
x,y,-z
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.1000(3) 0.2000(3) 0.0000(3)
C2 C 0.3000(3) 0.2000(3) 0.1500(3)
"""
    a = api.analyse_text(text, "m", Limits())
    proj = site_projectors(a.cif)
    assert set(proj) == {0} and np.allclose(proj[0], np.diag([1, 1, 0]))
    rng = np.random.default_rng(3)
    for _ in range(20):
        a2 = perturbed(a, rng)
        z = [m.pos[i][2] for m in a2.structure.molecules for i in range(len(m)) if m.labels[i] == "C1"]
        assert z and all(abs(v) < 1e-9 for v in z)


# ---- REGISTRO items 63–64 ---------------------------------------------------------------------------------------
def test_self_overlapping_rings_fall_back_to_the_hull():
    """Two rings of one fragment that overlap each other in projection (cages, ring stacks inside a network): the sum
    of pairwise clips would count the common area twice and give ratios above 1, so the outline is the hull."""
    sq = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0), (0.0, 2.0, 0.0)]
    pts = np.array(sq + [(x + 1.0, y + 1.0, 0.05) for x, y, _ in sq])
    core = Core(fit_plane(pts), np.arange(8), "rings", rings=((0, 1, 2, 3), (4, 5, 6, 7)))
    polys = fragment_polygons(core, pts[:, :2])
    assert len(polys) == 1 and polygon_area(polys[0]) == pytest.approx(8.0, rel=1e-6)  # hull of the two squares (3×3 − 2 corners)
    assert polygons_overlap(polys, polys) == pytest.approx(polygon_area(polys[0]), rel=1e-6)
    apart = np.array(sq + [(x + 2.0, y, 0.0) for x, y, _ in sq])  # touching along an edge: not an overlap
    core = Core(fit_plane(apart), np.arange(8), "rings", rings=((0, 1, 2, 3), (4, 5, 6, 7)))
    assert len(fragment_polygons(core, apart[:, :2])) == 2


def test_conjugated_macrocycle_rings_are_not_duplicated():
    """H₂Pc with fragment_merge = "conjugated": one 40-atom fragment with its four isoindoles and the 16-membered inner
    ring, each ring once; the overlap of the column stack is 0.351 (exact union computed independently with shapely,
    `paper/scripts/geometry_oracle.py`), not the 0.66 of the hull fallback triggered by the duplicated rings."""
    a = compute_all(api.analyse_file(COD / "2008099.cif", Limits(fragment_merge="conjugated")))
    for f in (f for fs in a.fragments for f in fs):
        assert sorted(len(r) for r in f.rings) == [5, 5, 5, 5, 6, 6, 6, 6, 16]
        assert len(set(map(frozenset, f.rings))) == len(f.rings)
    st = next(p for p in a.pairs if p.kind == "stack")
    assert st.overlap_area == pytest.approx(0.3511, abs=5e-4)


def test_bowl_shaped_polycycle_has_no_fragment():
    """Corannulene (COD 4020270): all rings fuse into one bowl-shaped system that fails the planarity test, and a
    non-planar fused system is not split into planar parts (only rings with a metal are dropped and retried)."""
    a = api.analyse_file(COD / "4020270.cif", Limits())
    assert a.structure.molecules and not any(a.fragments)
    assert not a.pairs
