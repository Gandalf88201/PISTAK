"""Third round of the pre-submission review of the paper (REGISTRO items 67–71).

* N1: `pi_scope = "organic"` also keeps structures with no C–C or C–H bond when they have a planar ring of non-metal
  atoms containing B, C or N (melamine, cyanuric acid, cyanuric chloride, borazines), and still skips carbides,
  silicate-type rings and networks whose rings contain metals.
* N2: in the Monte Carlo, an atom on a special position gets one deviate per free parameter, with the s.u. of that
  parameter; the site constraints give the other coordinates. The plain projection onto the site symmetry averaged
  two deviates and gave σ/√2 along x = y.

The CIFs in tests/data/cod are from the Crystallography Open Database (public domain)."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

import cifgen
from pistak import api
from pistak.core.analysis import has_organic_bonds, in_pi_scope
from pistak.core.config import Limits
from pistak.core.uncertainty import constraint_matrix, perturbed, site_constraints, site_projectors

COD = Path(__file__).resolve().parent / "data" / "cod"


# ---- N1: pi_scope ------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("kind", ["melamine", "cyanuric_acid", "cyanuric_chloride"])
def test_triazines_are_in_scope_and_stack(kind):
    """Melamine, cyanuric acid and cyanuric chloride have no C–C or C–H bond. Stacked along c at 3.40 Å with a
    1.17 Å offset, each triazine ring is a π fragment with the defaults and the pair with its translate is a π-stack."""
    a = api.analyse_text(cifgen.triazine_stack(kind), kind, Limits())
    assert not has_organic_bonds(a.structure) and in_pi_scope(a.structure, a.limits)
    frags = [f for fs in a.fragments for f in fs]
    assert len(frags) == 1 and frags[0].kind == "rings"
    m = a.structure.molecules[0]
    assert sorted(m.el[i] for i in frags[0].atoms) == ["C", "C", "C", "N", "N", "N"]
    stacks = [p for p in a.pairs if p.kind == "stack"]
    assert len(stacks) == 1
    p = stacks[0]
    assert p.interplanar == pytest.approx(3.40, abs=1e-4)
    assert p.dihedral == pytest.approx(0.0, abs=1e-3)
    assert p.slip == pytest.approx(math.hypot(1.0, 0.6), abs=1e-4)
    assert p.overlap > 0.10


def test_triazines_and_borazines_are_in_scope():
    """A planar B₃N₃ ring (borazine type, no H) has no carbon at all: in scope since round 3 (it was skipped)."""
    text = cifgen.make_cif([("B" if k % 2 == 0 else "N", 1.435 * np.array([math.cos(k * math.pi / 3), math.sin(k * math.pi / 3), 0.0]))
                            for k in range(6)], cifgen.Cell(5.0, 5.0, 5.0, 90.0, 90.0, 90.0), ["x, y, z"], "bn", "P 1",
                           origin_frac=(0.5, 0.5, 0.5))
    a = api.analyse_text(text, "bn", Limits())
    assert not has_organic_bonds(a.structure) and in_pi_scope(a.structure, a.limits)
    assert any(f for f in a.fragments)


def test_melaminium_salt_is_in_scope():
    """COD 7206203, (H₂mela)₂[CuCl₅]Cl: no C–C or C–H bond, so the old criterion skipped it. With the ring criterion
    the four melaminium rings are fragments and the six pairs of `pi_scope = "all"` are found (none is a π-stack:
    the parallel pairs at 3.23 and 3.32 Å do not overlap)."""
    a = api.analyse_file(COD / "7206203.cif", Limits())
    b = api.analyse_file(COD / "7206203.cif", Limits(pi_scope="all"))
    assert not has_organic_bonds(a.structure) and in_pi_scope(a.structure, a.limits)
    m = a.structure.molecules
    rings = [sorted(m[i].el[j] for j in f.atoms) for i, fs in enumerate(a.fragments) for f in fs]
    assert len(rings) == 4 and all(r == ["C", "C", "C", "N", "N", "N"] for r in rings)
    assert len(a.pairs) == len(b.pairs) == 6
    assert not any(p.kind == "stack" for p in a.pairs)


def test_carbide_is_still_skipped():
    """COD 1538535 (Ta₂VC₂): the Ta/V/C rings of the network contain metals, so no π fragment by default (they were
    π-stacks at 0.00 and 2.14 Å before REGISTRO item 49); `pi_scope = "all"` still finds them."""
    a = api.analyse_file(COD / "1538535.cif", Limits())
    assert not in_pi_scope(a.structure, a.limits)
    assert all(not f for f in a.fragments) and not a.pairs
    b = api.analyse_file(COD / "1538535.cif", Limits(pi_scope="all"))
    assert any(f for f in b.fragments)


def test_carbonate_and_cyanide_are_still_skipped():
    """A bidentate carbonate (planar Cu–O–C–O ring) and a Cu–C≡N–Cu–N≡C square: their rings have two adjacent
    non-metals but contain the metal, so they are not π rings."""
    carb = [("Cu", np.array([0.0, 0.0, 0.0])), ("O", np.array([1.95, 0.0, 0.0])), ("O", np.array([0.0, 1.95, 0.0])),
            ("C", np.array([1.6, 1.6, 0.0])), ("O", np.array([2.45, 2.45, 0.0]))]
    cyan = [("Cu", np.array([0.0, 0.0, 0.0])), ("C", np.array([1.9, 0.0, 0.0])), ("N", np.array([3.05, 0.0, 0.0])),
            ("Cu", np.array([4.95, 0.0, 0.0])), ("N", np.array([0.0, 1.9, 0.0])), ("C", np.array([0.0, 3.05, 0.0]))]
    for name, atoms in (("carbonate", carb), ("cyanide", cyan)):
        text = cifgen.make_cif(atoms, cifgen.Cell(12.0, 12.0, 12.0, 90.0, 90.0, 90.0), ["x, y, z"], name, "P 1",
                               origin_frac=(0.3, 0.3, 0.5))
        a = api.analyse_text(text, name, Limits())
        assert not has_organic_bonds(a.structure) and not in_pi_scope(a.structure, a.limits), name
        assert all(not f for f in a.fragments)


# ---- N2: special positions in the Monte Carlo --------------------------------------------------------------------
P6MM = ["x, y, z", "-y, x-y, z", "-x+y, -x, z", "-x, -y, z", "y, -x+y, z", "x-y, x, z",
        "-y, -x, z", "-x+y, y, z", "x, x-y, z", "y, x, z", "x-y, -y, z", "-x, -x+y, z"]


def _cif(symops: list[str], cell: tuple[float, ...], atoms: list[str]) -> str:
    lines = ["data_sp", *(f"_cell_length_{k} {v}" for k, v in zip("abc", cell[:3])),
             *(f"_cell_angle_{k} {v}" for k, v in zip(("alpha", "beta", "gamma"), cell[3:])),
             "loop_", "_symmetry_equiv_pos_as_xyz", *(f"'{s}'" for s in symops),
             "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
             "_atom_site_fract_z", *atoms, ""]
    return "\n".join(lines)


def _deltas(text: str, label: str, n: int = 4000, seed: int = 11) -> tuple[np.ndarray, np.ndarray]:
    """Fractional displacements of the CIF position of `label` in `n` perturbed copies (cell s.u. not used), and the
    s.u. of the CIF."""
    a = api.analyse_text(text, "sp", Limits())
    i = next(k for k, at in enumerate(a.cif.atoms) if at.label == label)
    x0, su = a.cif.atoms[i].frac, a.cif.atoms[i].su
    inv = np.linalg.inv(a.structure.vectors)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        s2 = perturbed(a, rng, include_cell=False).structure
        frac = np.vstack([m.pos for m in s2.molecules]) @ inv
        d = frac - x0
        d -= np.round(d)
        k = int(np.argmin(np.linalg.norm(d, axis=1)))
        out.append(d[k])
    return np.array(out), np.asarray(su)


def test_constraint_matrices_of_linear_constraints():
    """(x, x, z), (x, 2x, z), (x, −x, z) and a mirror z = 0: one column per free parameter, dependent and fixed
    coordinates follow; every matrix maps the allowed subspace onto itself (C·v = v)."""
    su = np.array([3e-4, 3e-4, 4e-4])
    swap = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    cases = {
        "x,x,z": (0.5 * (np.eye(3) + swap), [[1, 0, 0], [1, 0, 0], [0, 0, 1]]),
        "x,2x,z": (0.5 * (np.eye(3) + np.array([[-1, 1, 0], [0, 1, 0], [0, 0, 1]])), [[1, 0, 0], [2, 0, 0], [0, 0, 1]]),
        "x,-x,z": (0.5 * (np.eye(3) + np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])), [[1, 0, 0], [-1, 0, 0], [0, 0, 1]]),
        "x,y,0": (np.diag([1.0, 1.0, 0.0]), np.diag([1.0, 1.0, 0.0])),
    }
    for name, (proj, expected) in cases.items():
        c = constraint_matrix(proj, su)
        assert np.allclose(c, expected), name
        assert np.allclose(c @ c, c) and np.allclose(c @ proj, proj), name
    # the s.u. only on y: y is the free parameter and x follows
    assert np.allclose(constraint_matrix(cases["x,x,z"][0], np.array([0.0, 3e-4, 4e-4])), [[0, 1, 0], [0, 1, 0], [0, 0, 1]])


def test_monte_carlo_special_position_keeps_the_su_of_the_free_parameter():
    """C1 at (x, x, z) on the mirror y, x, z of a tetragonal cell: the sample s.d. of x and of y is the CIF s.u. of x
    (not s.u./√2 as with the averaging projector) and x − y stays 0; z keeps its own s.u."""
    text = _cif(["x, y, z", "y, x, z"], (7.0, 7.0, 9.0, 90, 90, 90),
                ["C1 C 0.1000(3) 0.1000(3) 0.2000(4)", "C2 C 0.3000(3) 0.1500(3) 0.2500(4)"])
    d, su = _deltas(text, "C1")
    assert np.max(np.abs(d[:, 0] - d[:, 1])) < 1e-12
    sd = d.std(axis=0, ddof=1)
    rel = 3 / math.sqrt(2 * (len(d) - 1))  # three standard errors of a sample s.d.
    assert sd[0] == pytest.approx(su[0], rel=rel) and sd[1] == pytest.approx(su[1], rel=rel)
    assert sd[2] == pytest.approx(su[2], rel=rel)
    assert sd[0] > 0.9 * su[0]  # the old projection gave su/√2 = 0.71 su


def test_monte_carlo_hexagonal_x_2x_z():
    """P6mm, C1 at (x, 2x, z) on the mirror −x+y, y, z: s.d. of x is s.u.(x), y moves by exactly 2δx (s.d. 2 s.u.(x),
    the s.u. a refinement program prints for y), z free; a second atom at (x, −x, z) keeps y = −x."""
    text = _cif(P6MM, (8.0, 8.0, 5.0, 90, 90, 120),
                ["C1 C 0.1000(3) 0.2000(6) 0.3000(4)", "C2 C 0.3500(2) -0.3500(2) 0.1000(5)"])
    a = api.analyse_text(text, "sp", Limits())
    assert set(site_projectors(a.cif)) == {0, 1}
    assert np.allclose(site_constraints(a.cif)[0], [[1, 0, 0], [2, 0, 0], [0, 0, 1]])
    assert np.allclose(site_constraints(a.cif)[1], [[1, 0, 0], [-1, 0, 0], [0, 0, 1]])
    d, su = _deltas(text, "C1", n=3000)
    assert np.max(np.abs(d[:, 1] - 2 * d[:, 0])) < 1e-12
    sd = d.std(axis=0, ddof=1)
    rel = 3 / math.sqrt(2 * (len(d) - 1))
    assert sd[0] == pytest.approx(su[0], rel=rel)
    assert sd[1] == pytest.approx(su[1], rel=rel)  # 2 × 0.0003 = 0.0006, as in the CIF
    assert sd[2] == pytest.approx(su[2], rel=rel)
    d2, _ = _deltas(text, "C2", n=200)
    assert np.max(np.abs(d2[:, 0] + d2[:, 1])) < 1e-12


def test_monte_carlo_atom_fixed_by_symmetry_does_not_move():
    """An atom on the sixfold axis of P6mm (0, 0, z) moves only along z, even though the CIF gives s.u. for x and y
    (the mirror z = 0 is `test_monte_carlo_keeps_special_position_constraints`, REGISTRO item 61)."""
    text = _cif(P6MM, (8.0, 8.0, 5.0, 90, 90, 120), ["C1 C 0.0000(3) 0.0000(3) 0.3000(4)"])
    d, _ = _deltas(text, "C1", n=300)
    assert np.max(np.abs(d[:, :2])) < 1e-12 and d[:, 2].std() > 0


# ---- N5: the planar-skeleton fallback is limited to sp²/sp skeletons ---------------------------------------------
def _dichloroethane() -> list[tuple[str, np.ndarray]]:
    """anti-1,2-Dichloroethane (C–C 1.513, C–Cl 1.79, C–H 1.09 Å, tetrahedral angles): Cl–C–C–Cl is planar."""
    c1, c2 = np.zeros(3), np.array([1.513, 0.0, 0.0])
    t = math.radians(109.5)
    u1 = np.array([math.cos(t), math.sin(t), 0.0])
    u2 = -u1
    atoms = [("C", c1), ("C", c2), ("Cl", c1 + 1.79 * u1), ("Cl", c2 + 1.79 * u2)]
    for c, a, b in ((c1, np.array([1.0, 0, 0]), u1), (c2, np.array([-1.0, 0, 0]), u2)):
        s = -(a + b) / np.linalg.norm(a + b)
        for z in (1, -1):
            atoms.append(("H", c + 1.09 * (s * math.cos(math.radians(54.75)) + np.array([0, 0, z]) * math.sin(math.radians(54.75)))))
    return atoms


def _trigonal(centre: str, arm: str, d: float, h: bool) -> list[tuple[str, np.ndarray]]:
    """Nitrate (N, O, 1.25 Å) or guanidinium (C, N, 1.33 Å, two N–H per arm in the plane)."""
    atoms = [(centre, np.zeros(3))]
    for k in range(3):
        a = 2 * math.pi * k / 3
        u = np.array([math.cos(a), math.sin(a), 0.0])
        n = d * u
        atoms.append((arm, n))
        if h:
            for s in (1, -1):
                b = a + s * math.pi / 3
                atoms.append(("H", n + 1.01 * np.array([math.cos(b), math.sin(b), 0.0])))
    return atoms


def _fallback_fragments(atoms, **kw):
    text = cifgen.make_cif(atoms, cifgen.Cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0), ["x, y, z"], "fb", "P 1",
                           origin_frac=(0.4, 0.4, 0.5))
    a = api.analyse_text(text, "fb", Limits(**kw))
    return [f for fs in a.fragments for f in fs]


def test_fallback_keeps_sp2_skeletons():
    """Nitrate and guanidinium (no ring, planar, every atom with at most three neighbours) are fallback fragments.
    Neither has a C–C or C–H bond, so `pi_scope = "all"` is needed to reach the fallback at all."""
    for atoms in (_trigonal("N", "O", 1.25, False), _trigonal("C", "N", 1.33, True)):
        frags = _fallback_fragments(atoms, pi_scope="all")
        assert len(frags) == 1 and frags[0].kind == "molecule" and len(frags[0].atoms) == 4


def test_fallback_rejects_planar_sp3_chains():
    """anti-1,2-Dichloroethane: the heavy-atom skeleton is planar and not linear, but both carbons have four
    neighbours (C, Cl, 2 H), so it is not a π fragment (it was one before round 3, N5). Without its hydrogens the
    carbons look two-coordinated and the fallback accepts it: the documented limitation."""
    assert _fallback_fragments(_dichloroethane()) == []
    no_h = [(e, p) for e, p in _dichloroethane() if e != "H"]
    frags = _fallback_fragments(no_h, pi_scope="all")
    assert len(frags) == 1 and frags[0].kind == "molecule"


# ---- cell parameters in the Monte Carlo: only the free ones are drawn -------------------------------------------
P4 = ["x, y, z", "-y, x, z", "-x, -y, z", "y, -x, z"]


def test_cell_constraints_from_symmetry():
    """Constraints from G = RᵀGR, not from a crystal-system name: tetragonal b = a and fixed angles, hexagonal b = a
    and γ fixed, monoclinic α and γ fixed, triclinic all free; a cell inconsistent with its operations gives None."""
    from pistak.core.uncertainty import cell_constraints

    def cons(symops, cell):
        a = api.analyse_text(_cif(symops, cell, ["C1 C 0.1000(3) 0.2000(3) 0.3000(4)"]), "c", Limits())
        return cell_constraints(a.cif)

    tet = cons(P4, ("7.000(2)", "7.000(2)", "9.000(3)", 90, 90, 90))
    assert np.allclose(tet, [[1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]] + [[0] * 6] * 3)
    hexa = cons(P6MM, ("8.000(2)", "8.000(2)", "5.000(1)", 90, 90, 120))
    assert np.allclose(hexa[:3, :3], [[1, 0, 0], [1, 0, 0], [0, 0, 1]]) and np.allclose(hexa[3:], 0)
    mono = cons(["x, y, z", "-x, y+1/2, -z+1/2"], ("7.000(2)", "8.000(2)", "9.000(3)", 90, "100.00(2)", 90))
    assert np.allclose(mono, np.diag([1, 1, 1, 0, 1, 0]))
    tri = cons(["x, y, z", "-x, -y, -z"], ("7.000(2)", "8.000(2)", "9.000(3)", "80.00(2)", "100.00(2)", "95.00(2)"))
    assert np.allclose(tri, np.eye(6))
    assert cons(P4, ("7.000(2)", "7.100(2)", "9.000(3)", 90, 90, 90)) is None  # a ≠ b with a fourfold axis


def test_tetragonal_cell_keeps_a_equal_b_and_equivalent_pairs():
    """P4 with two benzenes per asymmetric unit stacked along c: a and b stay equal in every sample with σ(a) = s.u.(a),
    and the four copies of a contact related by the fourfold axis keep identical descriptors in every sample."""
    from pistak.core.uncertainty import perturbed

    mol = [(el, p + np.array([4.0, 2.0, 0.0])) for el, p in cifgen.ring_system(1)]
    text = cifgen.make_cif(mol, cifgen.Cell(14.0, 14.0, 3.6, 90.0, 90.0, 90.0), P4, "p4", "P 4")
    text = text.replace("_cell_length_a 14.00000", "_cell_length_a 14.0000(20)").replace(
        "_cell_length_b 14.00000", "_cell_length_b 14.0000(20)").replace("_cell_length_c 3.60000", "_cell_length_c 3.6000(10)")
    a = api.analyse_text(text, "p4", Limits())
    assert a.cif.cell_su[0] == pytest.approx(0.002) and a.cif.cell_su[1] == pytest.approx(0.002)
    assert len(a.structure.molecules) == 4
    shifts = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], dtype=float)

    def contacts(s) -> list[float]:
        """Shortest centroid distance from molecule k to the one generated from it by the fourfold axis, k = 0…3."""
        c = [m.centroid for m in s.molecules]
        order = sorted(range(4), key=lambda k: int(s.cell_origin[s.molecules[k].cell_index[0], 1]))
        return [min(np.linalg.norm(c[order[(k + 1) % 4]] + t @ s.vectors - c[order[k]]) for t in shifts) for k in range(4)]

    rng = np.random.default_rng(5)
    lengths = []
    for _ in range(400):
        s2 = perturbed(a, rng).structure
        v = s2.vectors
        la, lb = np.linalg.norm(v[0]), np.linalg.norm(v[1])
        assert la == pytest.approx(lb, abs=1e-12)
        assert abs(np.dot(v[0], v[1])) < 1e-9 and abs(np.dot(v[0], v[2])) < 1e-9
        d = contacts(s2)
        assert max(d) - min(d) < 1e-9  # equivalent contacts stay equivalent
        lengths.append(la)
    assert np.std(lengths, ddof=1) == pytest.approx(0.002, rel=3 / math.sqrt(2 * 399))
