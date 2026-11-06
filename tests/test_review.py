"""Defects found by the code review (validation/REGISTRO.md), each with an oracle that does not use pistak.core:
known topology of polycycles, normals recomputed with numpy, theoretical limit of the cells to explore."""
import itertools
import math
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import cifgen
from pistak import api
from pistak.api import analyse_text, finish, reanalyse
from pistak.core.cif import Cell
from pistak.core.config import Limits
from pistak.core.cores import small_rings
from pistak.core.pairs import MAX_SHELLS

EX = Path(__file__).resolve().parents[1] / "examples"

ROOT = Path(__file__).resolve().parents[1]


# ---- tokenize: whitespace other than space and tab ---------------------------------------------------------------
@pytest.mark.parametrize("space", ["\xa0", "\x0c", "\x0b", " ", "\x1c"])
def test_tokenize_does_not_hang_on_unicode_spaces(space):
    # in a separate process: the defect was an infinite loop that fills the memory
    code = f"from pistak.core.cif import tokenize; print(tokenize('_cell_length_a{space}5.0 # nota'))"
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, timeout=30)
    assert out.stdout.strip() == "['_cell_length_a', '5.0']"


# ---- small_rings: rings whose bonds are all shared ----------------------------------------------------------------
def honeycomb(centres):
    """Atoms, bonds and expected hexagons of a benzenoid given by the centres of its rings (axial coordinates)."""
    cc, h = 1.4, 1.4 * math.sqrt(3) / 2
    atoms, rings = [], []
    for q, r in centres:
        c = np.array([q * 2 * h + r * h, r * 1.5 * cc])
        ring = []
        for k in range(6):
            p = c + cc * np.array([math.cos(math.pi / 6 + k * math.pi / 3), math.sin(math.pi / 6 + k * math.pi / 3)])
            i = next((i for i, a in enumerate(atoms) if np.linalg.norm(a - p) < 1e-6), None)
            if i is None:
                atoms.append(p)
                i = len(atoms) - 1
            ring.append(i)
        rings.append(frozenset(ring))
    bonds = [(i, j) for i, j in itertools.combinations(range(len(atoms)), 2) if np.linalg.norm(atoms[i] - atoms[j]) < 1.45]
    return len(atoms), bonds, rings


CORONENE = [(0, 0), (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
HEXABENZOCORONENE = CORONENE + [(2, -1), (1, 1), (-1, 2), (-2, 1), (-1, -1), (1, -2)]


@pytest.mark.parametrize("centres", [CORONENE, HEXABENZOCORONENE], ids=["coronene", "hexabenzocoronene"])
def test_small_rings_finds_inner_rings_in_any_order(centres):
    n, bonds, expected = honeycomb(centres)
    assert len(bonds) - n + 1 == len(expected)  # cyclomatic number = number of hexagons
    rng = random.Random(0)
    for _ in range(300):
        perm = list(range(n))
        rng.shuffle(perm)
        inv = {old: new for new, old in enumerate(perm)}
        b = np.array(sorted(tuple(sorted((inv[i], inv[j]))) for i, j in bonds))
        found = {frozenset(perm[i] for i in r) for r in small_rings(n, b, lambda i: True, 8)}
        assert found == set(expected)


# ---- mean_atom_height: mean distance from the plane of A ----------------------------------------------------------
def test_mean_height_from_plane_of_a():
    an = analyse_text((ROOT / "examples/tilted_anthracene.cif").read_text(encoding="utf-8"), "tilted", Limits())
    finish(an)
    tilted = 0
    for p in an.pairs:
        ctx = an.pair_context(p)
        pa, pb = ctx.pos_a[ctx.core_a.atoms], ctx.pos_b[ctx.core_b.atoms]
        ca = pa.mean(axis=0)
        na = np.linalg.svd(pa - ca)[2][-1]
        assert p.get("mean_atom_height") == pytest.approx(np.abs((pb - ca) @ na).mean(), abs=1e-9)
        if p.dihedral > 30:  # no longer the interplanar distance under another name
            tilted += 1
            assert abs(p.get("mean_atom_height") - p.interplanar) > 0.1
    assert tilted


# ---- MAX_SHELLS: a truncated search is reported -------------------------------------------------------------------
def long_acene_column(a: float = 3.5) -> str:
    """Undecacene (extent ~14 Å) stacked along a short a axis: more than MAX_SHELLS cells are needed."""
    mol = [(el, cifgen.rotation("y", 90) @ p) for el, p in cifgen.ring_system(11)]
    return cifgen.make_cif(mol, Cell(a, 8.0, 34.0, 90, 90, 90), ["x, y, z"], "undecacene", "P 1", origin_frac=(0.5, 0.5, 0.5))


def _truncated(an):
    return [x for x in an.alerts if x.code == "shellsTruncated"]


def test_truncated_shells_give_an_alert():
    an = finish(analyse_text(long_acene_column(), "undecacene", Limits()))
    extent = an.structure.molecules[0].extent
    need = lambda cutoff: math.floor((2 * extent + cutoff) / 3.5 + 1)  # theoretical limit along a (orthogonal cell)
    assert need(4.0) > MAX_SHELLS
    alerts = _truncated(an)
    assert len(alerts) == 1 and alerts[0].level == "B" and alerts[0].values["needed"].startswith(f"±{need(4.0)},")
    assert any(p.kind == "stack" for p in an.pairs)
    # with another maximum contact other cells are needed: the alert reused by reanalyse is recomputed
    smaller = reanalyse(an, an.limits.with_values(contact_cutoff=0.5))
    assert need(0.5) < need(4.0) and _truncated(smaller)[0].values["needed"].startswith(f"±{need(0.5)},")


def test_no_alert_for_normal_cells():
    an = finish(analyse_text(long_acene_column(a=6.0), "undecacene", Limits()))
    assert not _truncated(an)


# ---- contact_cutoff: every stack within max_interplanar is measured -----------------------------------------------
def naphthalene_column(interplanar: float, slip: float) -> str:
    """Infinite naphthalene column with exactly parallel planes: height and slip set by the monoclinic cell."""
    mol = cifgen.ring_system(2)
    c = math.hypot(interplanar, slip)
    beta = math.degrees(math.atan2(interplanar, slip))
    return cifgen.make_cif(mol, Cell(20.0, 20.0, c, 90, beta, 90), ["x, y, z"], "colonna", "P 1", origin_frac=(0.5, 0.5, 0.0))


@pytest.mark.parametrize("slip", [0.0, 1.4])
@pytest.mark.parametrize("d", [3.5, 4.05, 4.1, 4.19])
def test_stack_between_4_0_and_4_2_angstrom(d, slip):
    an = analyse_text(naphthalene_column(d, slip), "colonna", Limits())
    stacks = [p for p in an.pairs if p.kind == "stack"]
    assert len(stacks) == 1 and stacks[0].interplanar == pytest.approx(d, abs=1e-4)


def test_beyond_max_interplanar_pair_is_parallel():
    an = analyse_text(naphthalene_column(4.3, 0.0), "colonna", Limits())
    assert [p.kind for p in an.pairs] == ["parallel"]


def test_website_mode_keeps_4_angstrom():
    assert Limits.site().contact_cutoff == 4.0 and Limits.site(contact_cutoff=5.0).contact_cutoff == 5.0
    an = analyse_text(naphthalene_column(4.1, 0.0), "colonna", Limits.site())
    assert not an.pairs


# ---- beta_a / beta_b → beta / gamma -------------------------------------------------------------------------------
def test_old_names_of_beta_and_gamma():
    from pistak.params import REGISTRY

    an = finish(analyse_text((ROOT / "examples/tilted_anthracene.cif").read_text(encoding="utf-8"), "tilted", Limits()))
    p = an.pairs[0]
    assert "beta" in p.values and "gamma" in p.values and "beta_a" not in p.values
    assert p.get("beta_a") == p.get("beta") and p.get("beta_b") == p.get("gamma")
    assert REGISTRY.get("beta_b").key == "gamma" and REGISTRY.get("gamma").label.startswith("γ")


# ---- X–H distances by hybridisation ----------------------------------------------------------------------------
def test_hydrogens_by_hybridisation():
    from pistak.io.qc import normalize_hydrogens

    # C0 sp³ (CH3–), C4 sp² (=CH–), N7 with four neighbours (N⁺H3), P11 with one H: all X–H at 0.95 Å along different axes
    el = ["C", "H", "H", "H", "C", "H", "C", "N", "H", "H", "H", "P", "H"]
    bonds = np.array([(0, 1), (0, 2), (0, 3), (0, 4), (4, 5), (4, 6), (6, 7), (7, 8), (7, 9), (7, 10), (11, 12)])
    pos = np.zeros((len(el), 3))
    rng = np.random.default_rng(3)
    for x, h in bonds:
        if el[h] == "H":
            v = rng.normal(size=3)
            pos[h] = pos[x] + 0.95 * v / np.linalg.norm(v)
    out, moved = normalize_hydrogens(el, pos, bonds)
    length = lambda x, h: float(np.linalg.norm(out[h] - out[x]))
    # values from Allen & Bruno, Acta Cryst. B66, 380 (2010), means at T ≤ 60 K
    assert [length(0, h) for h in (1, 2, 3)] == pytest.approx([1.095] * 3, abs=1e-9)
    assert length(4, 5) == pytest.approx(1.085, abs=1e-9)
    assert [length(7, h) for h in (8, 9, 10)] == pytest.approx([1.040] * 3, abs=1e-9)
    assert length(11, 12) == pytest.approx(0.95, abs=1e-9) and moved == 7  # P–H without a source: left as it is


# ---- H-M symbols with an explicit setting (:1/:2, :H/:R) ----------------------------------------------------------
def test_setting_suffixes():
    from pistak.core.spacegroups import operations

    rhombo, source = operations(hm="R -3 :R")
    assert source == "table" and len(rhombo) == 6 and "z,x,y" in rhombo  # rhombohedral axes: no R centring
    assert len(operations(hm="R -3")[0]) == len(operations(hm="R -3 :H")[0]) == 18
    # origin 2 (default): inversion centre at the origin; origin 1: not
    assert "-x,-y,-z" in operations(hm="F d d d")[0] and "-x,-y,-z" in operations(hm="F d d d :2")[0]
    assert "-x,-y,-z" not in operations(hm="F d d d :1")[0]
    assert "-x,-y,-z" not in operations(hm="I 41/a :1")[0] and "-x,-y,-z" in operations(hm="I 41/a")[0]


# ---- element_of: ambiguous labels and dummy sites -----------------------------------------------------------------
@pytest.mark.parametrize("symbol, label, expected", [
    (None, "Cl1", "Cl"), (None, "CL1", "Cl"), (None, "BR2", "Br"), (None, "CA1", "C"), (None, "CD1", "C"),
    (None, "Co1", "Co"), (None, "D1", "H"), (None, "Cg1", ""), (None, "Q1", ""), ("Q", "Q1", ""),
    ("Fe2+", "Fe1", "Fe"), ("C", "Co1", "C"), ("Cl", "CL1", "Cl"),
])
def test_element_from_label(symbol, label, expected):
    from pistak.core.cif import element_of

    assert element_of(symbol, label) == expected


def test_dummy_sites_discarded():
    from pistak.core.cif import parse_cif

    text = "\n".join([
        "data_t", "_cell_length_a 10", "_cell_length_b 10", "_cell_length_c 10", "_cell_angle_alpha 90",
        "_cell_angle_beta 90", "_cell_angle_gamma 90", "loop_", "_symmetry_equiv_pos_as_xyz", "'x, y, z'",
        "loop_", "_atom_site_label", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        "_atom_site_calc_flag", "C1 0.1 0.1 0.1 d", "CL1 0.2 0.1 0.1 .", "X1 0.3 0.1 0.1 dum", "Cg1 0.4 0.1 0.1 .",
        "Q1 0.5 0.1 0.1 .", ""])
    assert [(a.label, a.el) for a in parse_cif(text).atoms] == [("C1", "C"), ("CL1", "Cl")]


# ---- listed operations against the group symbol -----------------------------------------------------------------
def test_operations_inconsistent_with_symbol():
    ok = cifgen.make_cif(cifgen.ring_system(2), Cell(8, 6, 9, 90, 100, 90), cifgen.P21C, "ok", "P 21/c", origin_frac=(0.3, 0.3, 0.3))
    wrong = cifgen.make_cif(cifgen.ring_system(2), Cell(8, 6, 9, 90, 100, 90), ["x, y, z", "-x, -y, -z"], "no", "P 21/c",
                            origin_frac=(0.3, 0.3, 0.3))
    assert not [x for x in finish(analyse_text(ok, "ok", Limits())).alerts if x.code == "symopsMismatch"]
    (alert,) = [x for x in finish(analyse_text(wrong, "no", Limits())).alerts if x.code == "symopsMismatch"]
    assert alert.level == "C" and alert.values == {"declared": "P 21/c", "n_cif": 2, "n_symbol": 4}


# ---- cell contents: atoms repeated by the CIF and special positions -----------------------------------------------
def test_cell_contents_with_repeated_atoms():
    from pistak.core.cif import parse_cif
    from pistak.core.crystal import cell_contents

    text = "\n".join([
        "data_t", "_cell_length_a 10", "_cell_length_b 10", "_cell_length_c 10", "_cell_angle_alpha 90",
        "_cell_angle_beta 90", "_cell_angle_gamma 90", "loop_", "_symmetry_equiv_pos_as_xyz", "'x, y, z'", "'-x, -y, -z'",
        "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
        "_atom_site_fract_z", "_atom_site_occupancy",
        "C1 C 0.1 0.2 0.3 0.5", "C1i C -0.1 -0.2 -0.3 0.5",  # image written again, half occupancy: 2 C in all
        "N1 N 0.3 0.3 0.1 1", "N1i N -0.3 -0.3 -0.1 1",       # image written again, full occupancy: 2 N, not 4
        "O1 O 0 0 0 0.5",                                      # on the inversion centre, half occupied: 0.5 O
        ""])
    assert cell_contents(parse_cif(text, disorder="all")) == pytest.approx({"C": 2.0, "N": 2.0, "O": 0.5})


def _cif_p1(atoms: list[tuple[str, str, tuple[float, float, float]]], a: float = 15.0) -> str:
    rows = "\n".join(f"{lab} {el} {x / a:.6f} {y / a:.6f} {z / a:.6f}" for lab, el, (x, y, z) in atoms)
    return (f"data_x\n_cell_length_a {a}\n_cell_length_b {a}\n_cell_length_c {a}\n_cell_angle_alpha 90\n"
            f"_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\nloop_\n"
            f"_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n{rows}\n")


@pytest.mark.parametrize("name, atoms, expected", [
    # terminal C triple-bonded to N or O: no missing hydrogen
    ("isonitrile", [("C1", "C", (7, 7, 7)), ("N1", "N", (8.42, 7, 7)), ("C2", "C", (9.58, 7, 7)),
                    ("H1", "H", (6.64, 8.03, 7)), ("H2", "H", (6.64, 6.49, 7.89)), ("H3", "H", (6.64, 6.49, 6.11))], set()),
    ("cyanide", [("C1", "C", (7, 7, 7)), ("N1", "N", (8.15, 7, 7))], set()),
    ("carbonyl", [("Cr1", "Cr", (7, 7, 7)), ("C1", "C", (8.91, 7, 7)), ("O1", "O", (10.05, 7, 7))], set()),
    # a C chain without hydrogens: the alert stays
    ("no H", [("C1", "C", (7, 7, 7)), ("C2", "C", (8.5, 7, 7)), ("C3", "C", (9.25, 8.3, 7))], {"noHydrogens"}),
])
def test_missing_hydrogens_without_false_positives(name, atoms, expected):
    from pistak.core.analysis import analyse_text
    from pistak.core.validation import validate

    a = analyse_text(_cif_p1(atoms), name)
    assert {x.code for x in validate(a)} & {"noHydrogens", "missingH"} == expected


def _benzene_p1(occ: float) -> str:
    """A benzene and its inversion image, written in P1: the P-1 symmetry is missed."""
    import math

    rows = []
    for k in range(6):
        t = 2 * math.pi * k / 6
        x, y, z = 3 + 1.39 * math.cos(t), 3 + 1.39 * math.sin(t), 3.0
        rows.append((f"C{k + 1}", (x / 12, y / 12, z / 12)))
        rows.append((f"C{k + 7}", (1 - x / 12, 1 - y / 12, 1 - z / 12)))
    atoms = "\n".join(f"{lab} C {f[0]:.6f} {f[1]:.6f} {f[2]:.6f} {occ}" for lab, f in rows)
    return ("data_x\n_cell_length_a 12\n_cell_length_b 12\n_cell_length_c 12\n_cell_angle_alpha 90\n_cell_angle_beta 90\n"
            "_cell_angle_gamma 90\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\nloop_\n_atom_site_label\n"
            "_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
            f"{atoms}\n")


def test_addsym_leaves_out_disordered_atoms():
    pytest.importorskip("spglib")
    from pistak.core.analysis import analyse_text
    from pistak.core.validation import validate

    full = {x.code: x for x in validate(analyse_text(_benzene_p1(1.0), "ordered"))}
    assert full["missedSymmetry"].level == "A"  # spglib finds Immm: the benzene written this way also has mirrors
    half = {x.code for x in validate(analyse_text(_benzene_p1(0.5), "disordered"))}
    assert "addsymDisorder" in half and "missedSymmetry" not in half


def test_missing_values_in_cif():
    # unknown density and Z («?»): no error, the formula is not checked
    from pistak.core.analysis import analyse_text
    from pistak.core.validation import validate

    text = _benzene_p1(1.0).replace("_cell_angle_gamma 90\n", "_cell_angle_gamma 90\n_exptl_crystal_density_diffrn ?\n"
                                    "_cell_formula_units_Z ?\n_chemical_formula_sum 'C6 H6'\n")
    codes = {x.code for x in validate(analyse_text(text, "missing"))}
    assert "formulaUnchecked" in codes and "density" not in codes and "checkFailed" not in codes


# ---- section 9: silent behaviour --------------------------------------------------------------------------------
def test_failing_spglib_not_mistaken_for_missing_spglib(monkeypatch):
    """If spglib is there but fails, the alert is checkFailed with the message, not «the spglib package is needed»."""
    import spglib

    from pistak.core.validation import validate

    def broken(*_a, **_k):
        return None  # like spglib < 2.5 when it fails

    monkeypatch.setattr(spglib, "get_symmetry_dataset", broken)
    a = api.analyse_file(EX / "demo_naphthalene.cif")
    codes = {x.code: x for x in validate(a)}
    assert "addsymUnavailable" not in codes
    assert codes["checkFailed"].values["check"] == "addsym" and "spglib" in codes["checkFailed"].values["error"]


def test_monte_carlo_counts_missing_samples():
    """A parameter missing in half the samples: σ comes from the others and `incomplete` says so."""
    import math

    from pistak.core.uncertainty import geometry_values, monte_carlo

    a = api.analyse_file(EX / "demo_naphthalene.cif")
    calls = {"n": 0}

    def flaky(analysis, pair):
        v = geometry_values(analysis, pair)
        calls["n"] += 1
        if calls["n"] > len(a.pairs) and calls["n"] % 2:  # the nominal values always, then one sample in two is missing
            v["interplanar"] = math.nan
        return v

    u = monte_carlo(a, samples=20, seed=1, evaluate=flaky)
    i = a.pairs[0].index
    assert 0 < u.incomplete[i]["interplanar"] < 20 == u.valid[i]["slip"]
    assert "slip" not in u.incomplete.get(i, {})


def test_untabulated_vdw_radius_is_reported():
    """Pm has no van der Waals radius in any of the tables (Mantina 2009, Alvarez 2013): the G alert says the mean value
    is used. Fe now has one (Alvarez 2013, table 1) and does not trigger the alert."""
    from pistak.core.elements import DEFAULT_VDW, VDW_RADII
    from pistak.core.validation import validate

    assert "Pm" not in VDW_RADII and VDW_RADII["Ra"] == 2.83  # Ra from Mantina 2009, table 12
    assert VDW_RADII["Fe"] == 2.44 and VDW_RADII["Cu"] == 2.38 and VDW_RADII["Pt"] == 2.29  # Alvarez 2013, table 1
    for el, alert in (("Pm", True), ("Fe", False)):
        cif = cifgen.make_cif([(el, np.array([1.0, 1.0, 1.0])), ("O", np.array([1.0, 1.0, 3.2]))],
                              Cell(10, 10, 10, 90, 90, 90), ["x, y, z"], "metal")
        alerts = {x.code: x for x in validate(analyse_text(cif, "metal"))}
        assert ("vdwDefault" in alerts) == alert
        if alert:
            assert alerts["vdwDefault"].values == {"elements": el, "r": DEFAULT_VDW}
    assert "vdwDefault" not in {x.code for x in validate(api.analyse_file(EX / "demo_naphthalene.cif"))}


# ---- item 14: disorder with negative groups ---------------------------------------------------------------------
def _disorder_cif(symops: list[str], rows: list[str]) -> str:
    return "\n".join(["data_dis", "_cell_length_a 10", "_cell_length_b 10", "_cell_length_c 10", "_cell_angle_alpha 90",
                      "_cell_angle_beta 90", "_cell_angle_gamma 90", "loop_", "_symmetry_equiv_pos_as_xyz",
                      *[f"'{s}'" for s in symops], "loop_", "_atom_site_label", "_atom_site_type_symbol",
                      "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy",
                      "_atom_site_disorder_assembly", "_atom_site_disorder_group", *rows, ""])


def test_groups_1_and_minus_1_alternatives_are_separated():
    """COD 4113834: 1 (0.61) and −1 (0.39) are two alternatives of the same site, not the same group on a symmetry
    element: the `major` rule keeps only part 1 (it used to keep both, atoms 0.4 Å apart)."""
    from pistak.core.cif import parse_cif

    cif = parse_cif(_disorder_cif(["x, y, z"], [
        "C1 C 0.50 0.50 0.50 1 . .",
        "C2A C 0.60 0.50 0.50 0.61 A 1", "C2B C 0.64 0.50 0.50 0.39 A -1",
    ]), disorder="major")
    assert [a.label for a in cif.atoms] == ["C1", "C2A"]


def test_negative_group_on_inversion_centre_one_orientation():
    """A −1 group on an inversion centre (occupancy 0.5): its image is the other orientation. Only one stays in the
    cell; an ordered fragment on the same centre instead keeps both halves."""
    from pistak.core.cif import parse_cif
    from pistak.core.structure import build_structure

    p1bar = ["x, y, z", "-x, -y, -z"]
    disordered = build_structure(parse_cif(_disorder_cif(p1bar, [
        "C1 C 0.03 0.00 0.00 0.5 . -1", "C2 C 0.03 0.15 0.00 0.5 . -1",  # image of C1 at 0.6 Å: the other orientation
    ]), disorder="major"))
    assert sorted(disordered.cell_labels) == ["C1", "C2"]
    ordered = build_structure(parse_cif(_disorder_cif(p1bar, [
        "C1 C 0.07 0.00 0.00 1 . .",  # C1···C1' 1.4 Å across the centre: a whole molecule
    ]), disorder="major"))
    assert sorted(ordered.cell_labels) == ["C1", "C1"]


def test_alternatives_in_different_assemblies_one_is_kept():
    """COD 7208500: the two orientations of an acetonitrile in different assemblies (A 0.65, B 0.35) overlap: they are
    alternatives, the major one is kept. Groups of different assemblies far from each other instead both stay."""
    from pistak.core.cif import parse_cif

    cif = parse_cif(_disorder_cif(["x, y, z"], [
        "C1 C 0.50 0.50 0.50 1 . .",
        "C2 C 0.60 0.50 0.50 0.65 A 1", "C3 C 0.63 0.50 0.50 0.35 B 2",  # 0.3 Å: alternatives
        "C4 C 0.10 0.10 0.10 0.70 C 1", "C5 C 0.30 0.30 0.30 0.60 D 1",  # far apart: independent disorders
    ]), disorder="major")
    assert [a.label for a in cif.atoms] == ["C1", "C2", "C4", "C5"]


# ---- extended networks: π pairs in polymers --------------------------------------------------------------------
def _polyphenylene(symops: list[str], centre: tuple[float, float, float], rings_per_cell: int = 1) -> str:
    """Poly(p-phenylene) chains along a (para C–C bond of 1.48 Å between rings), rings in the xz plane stacked along b
    every 3.6 Å: an infinite chain, with known π-stacks between neighbouring chains."""
    a, b, c = 4.28 * rings_per_cell, 3.6, 10.0
    ring = [(1.4, 0, 0), (-1.4, 0, 0), (0.7, 0, 1.2124), (-0.7, 0, 1.2124), (0.7, 0, -1.2124), (-0.7, 0, -1.2124)]
    rows = [f"C{k + 1} C {centre[0] + x / a:.6f} {centre[1] + y / b:.6f} {centre[2] + z / c:.6f}" for k, (x, y, z) in enumerate(ring)]
    return "\n".join(["data_ppp", f"_cell_length_a {a}", f"_cell_length_b {b}", f"_cell_length_c {c}",
                      "_cell_angle_alpha 90", "_cell_angle_beta 90", "_cell_angle_gamma 90", "loop_",
                      "_symmetry_equiv_pos_as_xyz", *[f"'{s}'" for s in symops], "loop_", "_atom_site_label",
                      "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", *rows, ""])


def test_polymer_pi_pairs_between_chains():
    """An infinite chain has no separate molecules, but its rings form π-stacks with the neighbouring chains: the pair
    at 3.6 Å along b is found, is marked as a polymer pair and is not exported as a dimer."""
    from pistak.core.diagnosis import why_no_pairs
    from pistak.io import qc

    a = analyse_text(_polyphenylene(["x, y, z"], (0.5, 0.5, 0.5)), "ppp")
    assert a.structure.molecules[0].extended
    stacks = [p for p in a.pairs if p.kind == "stack"]
    assert len(stacks) == 1 and stacks[0].polymer and stacks[0].shift in ((0, 1, 0), (0, -1, 0))
    assert abs(stacks[0].interplanar - 3.6) < 1e-4 and abs(stacks[0].dihedral) < 1e-4
    assert abs(stacks[0].contact - 3.6) < 1e-4  # between the two rings, not between bonded units of the chain
    assert why_no_pairs(a) is None and qc.exportable(a.pairs, None) == []


def test_polymer_symmetry_equivalent_pairs_merged(monkeypatch):
    """2₁ axis along the chain (two rings per cell): it maps the first ring onto the second, and the repeat unit onto
    another cut of the same chain, so the chain has no image as a whole molecule. With fragment images the two
    equivalent π-stacks become a single pair, ×2; without, two pairs remain."""
    import pistak.core.pairs as pairs_module

    text = _polyphenylene(["x, y, z", "x+1/2, -y, -z"], (0.25, 0.0, 0.5), rings_per_cell=2)
    stacks = [p for p in analyse_text(text, "ppp").pairs if p.kind == "stack"]
    assert len(stacks) == 1 and stacks[0].count == 2 and abs(stacks[0].interplanar - 3.6) < 1e-4
    monkeypatch.setattr(pairs_module, "fragment_images", lambda *_a, **_k: None)  # the oracle: without them, they are not merged
    assert len([p for p in analyse_text(text, "ppp").pairs if p.kind == "stack"]) == 2


def test_polymer_rings_whole_in_repeat_unit():
    """COD 7024607 (3D Mn(II) network): the breadth-first search reached one of the two N5 pyridines from two sides in
    different cells and left it split, so without a π fragment. Now the repeat unit keeps rings and ligands whole: there
    are four fragments and the ring pairs are the five of PLATON (CALC ALL, `validation/platon_polimeri/`), including
    N5···N5′ (−x, −y, −z) at 5.247 Å, β 58.0°."""
    from pistak.core.ringpairs import ring_pairs

    a = api.analyse_file(ROOT / "tests" / "data" / "cod" / "7024607.cif", Limits())
    m = a.structure.molecules[0]
    assert m.extended and len(a.fragments[0]) == 4
    for f in a.fragments[0]:  # whole rings: the six atoms within 3 Å of each other
        p = m.pos[list(f.atoms)]
        assert np.linalg.norm(p[:, None] - p[None], axis=2).max() < 3.0
    platon = [4.453, 4.717, 5.247, 5.630, 5.954]
    mine = sorted(rp.cg_distance for p in a.pairs for rp in ring_pairs(a.pair_context(p)))
    assert len(mine) == len(platon) and all(abs(x - y) < 1e-3 for x, y in zip(mine, platon))
    self_pair = [rp for p in a.pairs for rp in ring_pairs(a.pair_context(p)) if abs(rp.cg_distance - 5.247) < 1e-3]
    assert abs(self_pair[0].alpha) < 0.1 and abs(self_pair[0].beta - 58.0) < 0.1


def test_old_default_reported_and_user_toml_only_differences(tmp_path):
    """Item 31: contact_cutoff = 4.0 was the default until 3 November 2026. A user.toml saved earlier froze it: the
    `oldDefault` quick check says so (not in website mode, where 4.0 is the convention). Now user.toml gets only the
    thresholds that differ from the defaults, so a changed default also applies to the user."""
    import tomllib
    from dataclasses import replace

    from pistak.core.config import load_limits, save_user_config
    from pistak.core.validation import validate

    a = api.analyse_file(EX / "demo_naphthalene.cif", Limits())
    assert "oldDefault" not in {x.code for x in validate(a)}
    old = reanalyse(a, replace(a.limits, contact_cutoff=4.0))
    alert = {x.code: x for x in validate(old)}["oldDefault"]
    assert alert.level == "C" and alert.values == {"key": "contact_cutoff", "value": 4.0, "default": 4.5, "date": "2026-11-03"}
    assert "oldDefault" not in {x.code for x in validate(reanalyse(a, Limits.site()))}

    save_user_config(replace(Limits(), max_dihedral=15.0), config_dir=tmp_path)
    saved = tomllib.loads((tmp_path / "user.toml").read_text(encoding="utf-8"))["limits"]
    assert saved == {"max_dihedral": 15.0}
    assert load_limits(tmp_path) == replace(Limits(), max_dihedral=15.0)
