"""Round-4 review of the paper (PISTAK_review_round4.md): fallback long axis (M1.1), slips in each fragment's own
frame (M1.2), clamped pairs (M6.4), complete symmetry classes (S5), exporter warnings (M3.2–M3.3), reported value,
seed and covariance matrix in the Monte Carlo (M4.5–M4.6), terms of the listing (E5, S7)."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from cifgen import make_cif, ring_system, rotation
from pistak import api
from pistak.core.analysis import analyse_file, analyse_text
from pistak.core.cif import Cell
from pistak.core.config import Limits
from pistak.core.uncertainty import covariance_factor, monte_carlo, read_covariance

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
COD = ROOT / "tests" / "data" / "cod"


# ---- synthetic molecules ---------------------------------------------------------------------------------------
def para_disubstituted(x: str = "Cl", bond: float = 1.74) -> list[tuple[str, np.ndarray]]:
    """Benzene in the xy plane with the H atoms on the y axis (ring vertices at 90° and 270°) replaced by `x`: a
    1,4-disubstituted benzene whose long axis is y. The ring alone is isotropic."""
    out = []
    for el, p in ring_system(1):
        if el == "H" and abs(p[0]) < 1e-6:
            out.append((x, np.array([0.0, math.copysign(1.4 + bond, p[1]), 0.0])))
        else:
            out.append((el, p))
    return out


def translation_stack(mol, interplanar: float, slip: tuple[float, float], name: str = "stack") -> str:
    """One molecule per cell (P1, a = b = 14 Å along x and y) stacked along c = (slip_x, slip_y, interplanar)."""
    c = np.array([slip[0], slip[1], interplanar])
    lc = float(np.linalg.norm(c))
    cell = Cell(14.0, 14.0, lc, math.degrees(math.acos(slip[1] / lc)), math.degrees(math.acos(slip[0] / lc)), 90.0)
    return make_cif(mol, cell, ["x, y, z"], name, "P 1", origin_frac=(0.5, 0.5, 0.5))


def stack_pair(a):
    return next(p for p in a.pairs if p.kind == "stack")


# ---- M1.1: fallback long axis ----------------------------------------------------------------------------------
@pytest.mark.parametrize("mode", ["molecule", "substituents"])
def test_fallback_axis_of_a_para_disubstituted_benzene(mode):
    # slips imposed along the molecular axis (y, 1.30 Å) and across it (x, 0.70 Å)
    text = translation_stack(para_disubstituted(), 3.50, (0.70, 1.30))
    a = analyse_text(text, "pdcb", Limits(axis_fallback=mode))
    core = a.fragments[0][0]
    assert core.plane.aspect > 0.99 and core.axis_source == mode  # the ring is isotropic: the axis is the fallback
    assert abs(abs(core.long[1]) - 1.0) < 1e-9
    p = stack_pair(a)
    assert p.interplanar == pytest.approx(3.50, abs=1e-4)
    assert (p.long_slip, p.lat_slip) == pytest.approx((1.30, 0.70), abs=1e-4)
    assert (p.long_slip_a, p.lat_slip_a) == pytest.approx((1.30, 0.70), abs=1e-4)
    assert p.axis_a == p.axis_b == mode


def test_no_fallback_leaves_the_isotropic_ring_undefined():
    text = translation_stack(para_disubstituted(), 3.50, (0.70, 1.30))
    p = stack_pair(analyse_text(text, "pdcb", Limits(axis_fallback="none")))
    assert math.isnan(p.long_slip) and math.isnan(p.long_slip_a) and p.axis_a == "" and p.slip == pytest.approx(math.hypot(0.7, 1.3), abs=1e-4)


def test_benzene_has_no_long_axis_even_with_the_fallback():
    for mode in ("molecule", "substituents"):
        p = stack_pair(analyse_text(translation_stack(ring_system(1), 3.5, (0.7, 1.3)), "bz", Limits(axis_fallback=mode)))
        assert math.isnan(p.long_slip) and p.axis_a == ""


def test_elongated_cores_keep_their_own_axis():
    # regression on the existing benchmark: naphthalene, slips 1.50 / 0.80 Å by construction
    a = analyse_file(EX / "demo_naphthalene.cif", Limits())
    p = stack_pair(a)
    assert a.fragments[0][0].axis_source == "pca" and p.axis_a == "pca"
    assert (p.interplanar, p.long_slip, p.lat_slip) == pytest.approx((3.4, 1.5, 0.8), abs=1e-3)
    off = stack_pair(analyse_file(EX / "demo_naphthalene.cif", Limits(axis_fallback="none")))
    assert (off.long_slip, off.lat_slip) == (p.long_slip, p.lat_slip)


def test_site_mode_has_no_fallback_axis():
    a = analyse_text(translation_stack(para_disubstituted(), 3.5, (0.7, 1.3)), "pdcb", Limits.site())
    assert a.fragments[0][0].axis_source == "pca"


# ---- M1.2: slips in each fragment's own frame ------------------------------------------------------------------
def test_own_frame_slips_coincide_at_zero_twist():
    p = stack_pair(analyse_file(EX / "demo_naphthalene.cif", Limits()))
    assert p.twist == pytest.approx(0.0, abs=1e-6)
    for own in ((p.long_slip_a, p.lat_slip_a), (p.long_slip_b, p.lat_slip_b)):
        assert own == pytest.approx((p.long_slip, p.lat_slip), abs=1e-9)


def twisted_pair(tau: float, d=(1.10, 0.45, 3.40)) -> str:
    """Two independent naphthalenes in P1 (cell 24 Å): A in the xy plane with its long axis along x, B parallel to it,
    rotated by `tau` about z and displaced by `d`."""
    a = ring_system(2)
    b = [(el, rotation("z", tau) @ p + np.array(d)) for el, p in ring_system(2)]
    return make_cif(a + b, Cell(24.0, 24.0, 24.0, 90.0, 90.0, 90.0), ["x, y, z"], "twisted", "P 1", origin_frac=(0.4, 0.4, 0.4))


def test_own_frame_slips_of_a_twisted_pair():
    tau, (dx, dy, dz) = 50.0, (1.10, 0.45, 3.40)
    c, s = math.cos(math.radians(tau)), math.sin(math.radians(tau))
    from_a = (abs(dx), abs(dy))  # B seen from A: A's long axis is x
    from_b = (abs(dx * c + dy * s), abs(-dx * s + dy * c))  # A seen from B: B's long axis is at tau from x
    a = analyse_text(twisted_pair(tau), "tw", Limits())
    p = next(q for q in a.pairs if not q.intramolecular)
    assert p.kind == "stack" and p.interplanar == pytest.approx(dz, abs=1e-4)
    got = {(round(p.long_slip_a, 4), round(p.lat_slip_a, 4)), (round(p.long_slip_b, 4), round(p.lat_slip_b, 4))}
    assert got == {tuple(round(v, 4) for v in from_a), tuple(round(v, 4) for v in from_b)}
    # the pair frame is not reported above max_axis_twist (30°), the own-frame values are
    assert math.isnan(p.long_slip) and math.isnan(p.lat_slip)
    # with the old limit the bisector frame is defined: L at tau/2 from x
    wide = next(q for q in analyse_text(twisted_pair(tau), "tw", Limits(max_axis_twist=75)).pairs if not q.intramolecular)
    c2, s2 = math.cos(math.radians(tau / 2)), math.sin(math.radians(tau / 2))
    assert (wide.long_slip, wide.lat_slip) == pytest.approx((abs(dx * c2 + dy * s2), abs(-dx * s2 + dy * c2)), abs=1e-4)
    assert wide.twist == pytest.approx(tau, abs=2e-3)


def test_default_twist_limit():
    assert Limits().max_axis_twist == 30.0
    small = next(q for q in analyse_text(twisted_pair(20.0), "tw", Limits()).pairs if not q.intramolecular)
    assert math.isfinite(small.long_slip) and small.twist == pytest.approx(20.0, abs=2e-3)


def test_own_frame_slips_in_registry_and_monte_carlo():
    a = api.analyse_file(EX / "demo_naphthalene.cif", Limits())
    p = stack_pair(a)
    assert {"long_slip_a", "lat_slip_a", "long_slip_b", "lat_slip_b"} <= set(a.param_keys)
    assert p.get("long_slip_a") == p.long_slip_a
    u = monte_carlo(a, samples=5)
    assert "long_slip_a" in u.keys


# ---- M6.4: clamped pairs ---------------------------------------------------------------------------------------
def cyclophane(gap: float) -> str:
    """A [2.2]paracyclophane-like molecule: two benzene rings face to face `gap` Å apart, joined by two C–C bridges
    on the y axis (ring C – C – C – ring C: not bonded to each other and not bonded to a common atom)."""
    ring = [(el, p) for el, p in ring_system(1) if not (el == "H" and abs(p[0]) < 1e-6)]
    top = [(el, p + np.array([0.0, 0.0, gap])) for el, p in ring]
    bridges = []
    for sign in (+1, -1):
        bridges += [("C", np.array([0.0, sign * 2.5, gap / 2 - 0.75])), ("C", np.array([0.0, sign * 2.5, gap / 2 + 0.75]))]
    return make_cif(ring + top + bridges, Cell(20.0, 20.0, 20.0, 90.0, 90.0, 90.0), ["x, y, z"], "cyclophane", "P 1",
                    origin_frac=(0.5, 0.5, 0.4))


def test_clamped_pair_is_not_a_stack():
    a = analyse_text(cyclophane(3.00), "cp", Limits())
    assert len(a.structure.molecules) == 1 and len(a.fragments[0]) == 2
    (p,) = [q for q in a.pairs if q.intramolecular]
    assert not p.bonded and not p.bridged and p.clamped
    assert p.kind == "parallel" and p.interplanar == pytest.approx(3.0, abs=1e-4) and p.overlap_area > 0.9


def test_clamped_depends_on_min_interplanar_only():
    a = analyse_text(cyclophane(3.00), "cp", Limits())
    free = api.reanalyse(api.finish(a), Limits(min_interplanar=0.0))  # stage `pairs`: nothing is remeasured
    (p,) = [q for q in free.pairs if q.intramolecular]
    assert p.kind == "stack" and not p.clamped
    (q,) = [q for q in analyse_text(cyclophane(3.40), "cp", Limits()).pairs if q.intramolecular]
    assert q.kind == "stack" and not q.clamped  # above the limit a folded molecule can stack on itself


def test_short_intermolecular_stack_is_not_clamped():
    # between separate molecules a short distance is a real contact (TCNQ stacks at 3.17 Å)
    p = stack_pair(analyse_text(translation_stack(ring_system(2), 3.10, (1.2, 0.5)), "short", Limits()))
    assert p.kind == "stack" and not p.clamped


def test_cod_7239588_chelate_pair_is_clamped():
    a = analyse_file(COD / "7239588.cif", Limits())
    clamped = [p for p in a.pairs if p.clamped]
    assert len(clamped) == 1
    p = clamped[0]
    assert p.intramolecular and p.kind == "parallel" and p.interplanar == pytest.approx(3.103, abs=2e-3)
    assert sum(q.kind == "stack" for q in a.pairs) == 2  # the two intermolecular π-stacked pairs stay


def test_clamped_in_exports(tmp_path):
    from pistak.io import export
    from pistak.labels import relation

    a = api.analyse_text(cyclophane(3.00), "cp", Limits())
    data = json.loads(export.write_json(a, tmp_path / "a.json").read_text())
    intra = [d for d in data["pairs"] if d["clamped"]]
    assert len(intra) == 1 and data["version"] == 5 and "long_axis_a" in intra[0]
    assert "vincolata" in relation(next(p for p in a.pairs if p.clamped))


# ---- S5: complete symmetry classes -----------------------------------------------------------------------------
@pytest.mark.parametrize("cod", ["7123350", "7239588"])
def test_classes_do_not_depend_on_the_contact_cutoff(cod):
    base = analyse_file(COD / f"{cod}.cif", Limits())
    wide = analyse_file(COD / f"{cod}.cif", Limits(contact_cutoff=6.0))
    kb, kw = {p.key: p for p in base.pairs}, {p.key: p for p in wide.pairs}
    assert not kb.keys() - kw.keys()  # no pair is "lost" by a larger cut-off: the representative is canonical
    for k, p in kb.items():
        assert kw[k].count == p.count and kw[k].interplanar == pytest.approx(p.interplanar, abs=1e-9)
    # new − lost = change in the total (trivially, now with lost = 0)
    assert len(kw.keys() - kb.keys()) == len(wide.pairs) - len(base.pairs)


def test_representative_is_the_smallest_member_of_its_class():
    for name in ("7123350", "7239588", "2014747"):
        path = COD / f"{name}.cif"
        if not path.exists():
            continue
        for p in analyse_file(path, Limits()).pairs:
            assert p.key == min(p.members) and len(set(p.members)) == p.count


# ---- M3.2 / M3.3: exporter -------------------------------------------------------------------------------------
def test_exported_dimer_has_neutron_c_h_distances():
    from pistak.io import qc

    a = api.analyse_file(EX / "herringbone_anthracene.cif", Limits())
    p = stack_pair(a)
    d = qc.dimer(a, p)
    for mono, mol in ((d.a, a.structure.molecules[p.a]), (d.b, a.structure.molecules[p.b])):
        n = 0
        for i, j in mol.bonds:
            i, j = int(i), int(j)
            if {mono.el[i], mono.el[j]} == {"C", "H"}:
                assert np.linalg.norm(mono.pos[i] - mono.pos[j]) == pytest.approx(qc.XH_NEUTRON["Csp2"], abs=1e-9)
                n += 1
        assert n == 10  # every C–H of anthracene


def test_charge_and_multiplicity_warnings():
    from pistak.io import qc

    a = api.analyse_file(EX / "herringbone_anthracene.cif", Limits())
    p = stack_pair(a)
    ok = qc.dimer_input(a, p, qc.QCOptions())
    assert "is impossible" not in ok and "different chemical species" not in ok
    bad = qc.dimer_input(a, p, qc.QCOptions(charges=(1, 0)))  # a radical cation declared as a singlet
    assert "molecule A has 93 electrons with charge 1: multiplicity 1 is impossible" in bad
    fine = qc.dimer_input(a, p, qc.QCOptions(charges=(1, 0), multiplicities=(2, 1)))
    assert "is impossible" not in fine
    assert qc.electrons(["C", "H", "H", "H", "D"], 0) == 10


def test_metal_and_salt_warnings():
    from pistak.io import qc

    a = api.analyse_file(COD / "7239588.cif", Limits())
    p = next(q for q in a.pairs if q.kind == "stack" and not q.intramolecular and not q.polymer)
    text = qc.dimer_input(a, p, qc.QCOptions())
    assert "contains Zn" in text and "are not derived from the CIF" in text


# ---- M4.5 / M4.6: Monte Carlo ----------------------------------------------------------------------------------
def with_su(text: str, last_digit: int = 3) -> str:
    out, atoms = [], False
    for line in text.splitlines():
        parts = line.split()
        if atoms and len(parts) == 5:
            parts[2:] = [f"{float(v):.4f}({last_digit})" for v in parts[2:]]
            line = " ".join(parts)
        if line.strip() == "_atom_site_fract_z":
            atoms = True
        out.append(line)
    return "\n".join(out) + "\n"


def test_reported_value_is_nominal_and_seed_is_reproducible(tmp_path):
    from pistak.io import export

    text = with_su((EX / "herringbone_anthracene.cif").read_text())
    a = api.analyse_text(text, "h", Limits())
    p = stack_pair(a)
    before = (p.interplanar, p.dihedral, p.slip)
    u1 = api.uncertainty(a, samples=40, seed=11)
    assert (p.interplanar, p.dihedral, p.slip) == before and p.dihedral == 0.0  # θ = 0 by symmetry stays 0
    csv1 = export.write_csv(a, tmp_path / "1.csv").read_text()
    u2 = api.uncertainty(a, samples=40, seed=11)
    assert repr(u1.sigma) == repr(u2.sigma) and export.write_csv(a, tmp_path / "2.csv").read_text() == csv1
    data = json.loads(export.write_json(a, tmp_path / "a.json").read_text())
    assert data["summary"]["uncertainty"]["reported_value"] == "nominal" and data["summary"]["uncertainty"]["seed"] == 11
    # for a tilted pair the sample mean of a non-negative quantity near zero is biased upwards; the value reported is not
    tilted = next(q for q in a.pairs if q.kind == "tilted")
    assert u1.mean[tilted.index]["dihedral"] != tilted.dihedral


def test_covariance_matrix_diagonal_and_rigid(tmp_path):
    text = with_su(twisted_pair(10.0), 4)
    a = analyse_text(text, "tw", Limits())
    p = next(q for q in a.pairs if not q.intramolecular)
    su = np.array([at.su for at in a.cif.atoms]).ravel()
    independent = monte_carlo(a, samples=400, seed=2, include_cell=False)
    diagonal = monte_carlo(a, samples=400, seed=3, include_cell=False, covariance=np.diag(su**2))
    assert diagonal.covariance and not independent.covariance
    for k in ("interplanar", "slip", "dihedral"):
        assert diagonal.get(p, k) == pytest.approx(independent.get(p, k), rel=0.2)
    # perfectly correlated displacements (every atom moves by the same vector): nothing changes
    j = np.tile(np.eye(3), (len(a.cif.atoms), 1))  # 3N×3: the three rigid translations
    rigid = monte_carlo(a, samples=50, seed=1, include_cell=False, covariance=1e-6 * j @ j.T)
    assert independent.get(p, "interplanar") > 1e-4
    for k in ("interplanar", "slip", "dihedral", "centroid_distance"):
        assert rigid.get(p, k) < 1e-6
    f = covariance_factor(np.diag(su**2))
    assert np.allclose(f @ f.T, np.diag(su**2))
    # files: .npy and text; wrong shape rejected
    np.save(tmp_path / "c.npy", np.diag(su**2))
    np.savetxt(tmp_path / "c.txt", np.diag(su**2))
    assert np.allclose(read_covariance(tmp_path / "c.npy", a.cif), read_covariance(tmp_path / "c.txt", a.cif))
    np.save(tmp_path / "bad.npy", np.eye(5))
    with pytest.raises(ValueError):
        read_covariance(tmp_path / "bad.npy", a.cif)


def test_command_line_covariance_and_terms(tmp_path, capsys, monkeypatch):
    from pistak import i18n
    from pistak.cli import main

    cif = tmp_path / "tw.cif"
    cif.write_text(with_su(twisted_pair(10.0), 4))
    a = analyse_text(cif.read_text(), "tw", Limits())
    su = np.array([at.su for at in a.cif.atoms]).ravel()
    np.save(tmp_path / "c.npy", np.diag(su**2))
    lang = i18n.get_language()
    try:
        i18n.set_language("en")
        assert main(["analyze", str(cif), "--mc", "5", "--mc-covariance", str(tmp_path / "c.npy")]) == 0
        out = capsys.readouterr().out
        assert "π-stacked" in out and "angle" in out and "dihedral" not in out
        assert main(["analyze", str(EX / "herringbone_anthracene.cif"), "--rows", "50"]) == 0
        rows = [ln for ln in capsys.readouterr().out.splitlines() if " parallel " in ln or " tilted " in ln]
        assert rows and any(" – " in ln for ln in rows)  # no projected overlap: no "interplanar distance" printed
    finally:
        i18n.set_language(lang)


def test_copy_of_cod_files_is_listed():
    readme = (COD / "README.md").read_text()
    for name in ("7239588", "7123350"):
        assert (COD / f"{name}.cif").exists() and name in readme


# ---- T16: contacts per molecule --------------------------------------------------------------------------------
def test_contacts_per_molecule_do_not_depend_on_z():
    a = api.analyse_file(EX / "demo_naphthalene.cif", Limits())
    p = stack_pair(a)
    assert p.per_molecule == 2 and p.get("per_molecule") == 2  # a column: one neighbour above, one below
    h = api.analyse_file(EX / "herringbone_anthracene.cif", Limits())
    q = stack_pair(h)
    assert q.count == 2 and q.per_molecule == 2  # Z = 2: two contacts per cell, still two per molecule
    (c,) = [x for x in analyse_text(cyclophane(3.4), "cp", Limits()).pairs if x.intramolecular]
    assert c.per_molecule == 1


# ---- two elements on one site (mixed occupancy) ----------------------------------------------------------------
def test_no_angle_or_torsion_through_two_atoms_on_one_site():
    from pistak.core.intramolecular import molecule_geometry

    a = analyse_file(COD / "2108505.cif", Limits())  # Ho and Bi share a site (occupancies 0.963 and 0.037)
    items = list(molecule_geometry(a.structure, a.cif, kinds=("bond", "angle", "torsion")))
    assert any(g.kind == "bond" and "Ho1" in g.atoms for g in items)
    s = a.structure
    n = 0
    for g in items:
        if g.kind in ("angle", "torsion") and {"Ho1", "Bi1"} <= set(g.atoms):
            n += 1  # Ho and Bi of different sites may appear together; never the two on the same site
            pos = [s.molecules[m].pos[i] + s.shift_vector(sh) for m, sh, i in g.places]
            assert min(np.linalg.norm(p - q) for k, p in enumerate(pos) for q in pos[k + 1:]) > 0.4, g.atoms
    assert n > 0


# ---- symmetry classes hold only congruent contacts -------------------------------------------------------------
def test_class_members_are_congruent_when_an_operation_does_not_map_the_fragment():
    """COD 7044129 (C2/c, molecule on a twofold axis): the π fragment is one of two equivalent four-membered rings and
    the other is not a fragment, so half of the operations map the molecule but not the fragment. Those operations
    must not be used to build the classes: before, contacts with interatomic distances up to 5.8 Å apart were merged."""
    a = api.analyse_file(COD / "7044129.cif", Limits())
    s = a.structure

    def fingerprint(member):
        i, j, sh, fa, fb = member
        xa = s.molecules[i].pos[a.fragments[i][fa].atoms]
        xb = s.molecules[j].pos[a.fragments[j][fb].atoms] + s.shift_vector(sh)
        return np.sort(np.linalg.norm(xa[:, None, :] - xb[None, :, :], axis=2).ravel())

    assert all(len(f) == 1 for f in a.fragments)  # one fragment per molecule: the case that skipped the fragment maps
    n = 0
    for p in a.pairs:
        base = fingerprint(p.members[0])
        for m in p.members[1:]:
            n += 1
            assert np.abs(fingerprint(m) - base).max() < 1e-4, (p.index, m)
    assert n > 0
