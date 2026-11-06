"""Phase 4a: whole-structure quantities (structure registry) and quick checks (pre-checkCIF)."""
import csv
import math
from pathlib import Path

import pytest

from pistak import api, i18n
from pistak.core.config import Limits
from pistak.core.crystal import cell_contents, formula_text, parse_formula, site_multiplicity
from pistak.core.validation import counts, validate
from pistak.labels import ALERT_TEXT, alert_text
from pistak.params import StructureRegistry, compute_structure

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
COD = ROOT / "tests" / "data" / "cod"  # COD structures (public domain)
REAL = ROOT / "DATI_REALI"  # CIFs that cannot be redistributed: if missing, the tests that use them are skipped
DEMO = (EX / "demo_naphthalene.cif").read_text(encoding="utf-8")


def demo(extra: str = "", drop_h: bool = False, atoms: str = ""):
    """The demo naphthalene (P-1, molecule on the inversion centre) with extra CIF items."""
    text = DEMO.replace("loop_\n_symmetry_equiv_pos_as_xyz", extra + "\nloop_\n_symmetry_equiv_pos_as_xyz", 1)
    if drop_h:
        text = "\n".join(line for line in text.splitlines() if not line.startswith("H"))
    return api.analyse_text(text.rstrip() + "\n" + atoms, "demo.cif")


def codes(a) -> set[str]:
    return {al.code for al in a.alerts}


def real(name: str) -> Path:
    p = REAL / name
    if not p.is_file():
        pytest.skip(f"{p} is missing (real CIFs are not distributed)")
    return p


# ---- cell contents and formula --------------------------------------------------------------------------
def test_formula_read_from_cif():
    assert parse_formula("C14 H16 N6 O2") == {"C": 14, "H": 16, "N": 6, "O": 2}
    assert parse_formula("'C2.5 H4 Cl'") == {"C": 2.5, "H": 4, "Cl": 1}
    assert parse_formula(None) == {}
    assert formula_text({"H": 8, "C": 10, "N": 1}) == "C10 H8 N"


def test_cell_contents_with_site_multiplicity():
    a = demo()
    assert cell_contents(a.cif) == pytest.approx({"C": 10, "H": 8})  # 5 C and 4 H on the inversion centre: ×2
    assert all(site_multiplicity(a.cif, at) == 2 for at in a.cif.atoms)


@pytest.mark.parametrize("name", ["herringbone_anthracene.cif", "tilted_anthracene.cif", "dimers_naphthalene.cif",
                                  "zigzag_chain.cif"])
def test_contents_match_rebuilt_cell(name):
    """Without disorder the contents counted from the CIF are those of the rebuilt cell, even when the CIF also lists
    symmetry-equivalent atoms (like the synthetic examples)."""
    from collections import Counter

    a = api.analyse_file(EX / name)
    assert cell_contents(a.cif) == pytest.approx(dict(Counter(a.structure.cell_el)))


def test_structure_quantities():
    a = demo("_chemical_formula_sum 'C10 H8'\n_cell_formula_units_Z 1\n_exptl_crystal_density_diffrn 0.759\n"
             "_refine_ls_R_factor_gt 0.0412\n_refine_ls_wR_factor_ref 0.1050\n_refine_ls_goodness_of_fit_ref 1.03")
    v = a.structure_values
    assert v["z"] == 1 and v["z_prime"] == pytest.approx(0.5)
    assert v["density_calc"] == pytest.approx(128.17 * 1.66054 / v["volume"], rel=1e-3)
    assert (v["r1"], v["wr2"], v["goof"]) == (0.0412, 0.1050, 1.03)
    assert math.isnan(v["temperature"])
    assert codes(a) <= {"shortContact"} and counts(a.alerts)["B"] == 0


# ---- single checks --------------------------------------------------------------------------------------
def test_formula_differs_from_model():
    assert "formulaH" in codes(demo("_chemical_formula_sum 'C10 H10'\n_cell_formula_units_Z 1"))
    a = demo("_chemical_formula_sum 'C12 H8'\n_cell_formula_units_Z 1")
    al = next(x for x in a.alerts if x.code == "formulaHeavy")
    assert al.level == "B" and "C 10.00/12" in al.values["diff"] and al.values["model"] == "C10 H8"
    assert "formulaUnchecked" in codes(demo())


def test_density_differs():
    al = next(x for x in demo("_exptl_crystal_density_diffrn 0.90").alerts if x.code == "density")
    assert al.level == "B" and al.values["calc"] == pytest.approx(0.759, abs=1e-3)


def test_missing_hydrogens():
    a = demo(drop_h=True)
    assert "noHydrogens" in codes(a) and counts(a.alerts)["B"] >= 1


def test_overlapping_atoms():
    a = demo(atoms="N9 N 0.00000 0.10900 -0.22044\n")  # 0.12 Å from C1
    al = next(x for x in a.alerts if x.code == "overlap")
    assert al.level == "A" and "C1···N9" in al.values["atoms"]


def test_refinement():
    a = demo("_refine_ls_R_factor_gt 0.12\n_refine_ls_wR_factor_ref 0.22\n_refine_ls_goodness_of_fit_ref 0.5\n"
             "_refine_diff_density_max 1.5\n_refine_diff_density_min -0.4")
    levels = {al.code: al.level for al in a.alerts}
    assert levels["highR1"] == "B" and levels["highWR2"] == "C" and levels["goof"] == "B"
    assert levels["residualDensity"] == "C"
    assert "noRefinement" in codes(demo())


def test_missed_symmetry():
    pytest.importorskip("spglib")
    a = api.analyse_file(EX / "infinite_chain.cif")  # declared P1, actually tetragonal
    al = next(x for x in a.alerts if x.code == "missedSymmetry")
    assert al.level == "A" and al.values["ops_found"] > al.values["ops_declared"] == 1
    assert not {"missedSymmetry", "pseudoSymmetry"} & codes(api.analyse_file(EX / "herringbone_anthracene.cif"))


def test_without_spglib(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def no_spglib(name, *args, **kw):
        if name == "spglib":
            raise ImportError(name)
        return real_import(name, *args, **kw)

    monkeypatch.setattr(builtins, "__import__", no_spglib)
    a = api.analyse_file(EX / "infinite_chain.cif")
    assert "addsymUnavailable" in codes(a) and "missedSymmetry" not in codes(a)


def test_failing_check_does_not_stop_the_others(monkeypatch):
    from pistak.core import validation

    def broken(_a):
        raise RuntimeError("broken")

    monkeypatch.setitem(validation.CHECKS, "broken", broken)
    alerts = validate(api.analyse_file(EX / "demo_naphthalene.cif"))
    failed = [x for x in alerts if x.code == "checkFailed"]
    assert failed and failed[0].values["check"] == "broken" and "noRefinement" in {x.code for x in alerts}
    assert [x.rank for x in alerts] == sorted(x.rank for x in alerts)


def test_alert_texts_in_two_languages():
    a = demo("_chemical_formula_sum 'C12 H8'\n_cell_formula_units_Z 1\n_refine_ls_R_factor_gt 0.12", drop_h=True)
    assert all(al.code in ALERT_TEXT for al in a.alerts)
    it = [alert_text(al) for al in a.alerts]
    try:
        i18n.set_language("en")
        en = [alert_text(al) for al in a.alerts]
    finally:
        i18n.set_language("it")
    assert "R1 alto: 0.12." in it and "High R1: 0.12." in en
    assert all("{" not in t for t in it + en)


# ---- registry, reuse, export ----------------------------------------------------------------------------
def test_own_structure_registry_and_errors():
    reg = StructureRegistry()

    @reg.parameter("n_mol", "Molecules", digits=0)
    def n_mol(a):
        return len(a.structure.molecules)

    @reg.parameter("broken", "Broken")
    def broken(_a):
        raise ValueError("no")

    a = compute_structure(api.analyse_file(EX / "demo_naphthalene.cif"), reg)
    assert a.structure_values["n_mol"] == 1 and math.isnan(a.structure_values["broken"])
    assert "structure:broken" in a.param_errors


def test_stack_thresholds_reuse_the_checks():
    a = api.analyse_file(EX / "herringbone_anthracene.cif", Limits())
    b = api.reanalyse(a, Limits(max_dihedral=30))
    assert b is not a and b.alerts is a.alerts and b.structure_values.keys() == a.structure_values.keys()
    c = api.reanalyse(a, Limits(bond_tolerance=0.3))  # structure redone: checks redone
    assert c.alerts is not a.alerts


def test_json_and_batch_summary(tmp_path):
    from pistak.io import export
    from pistak.io.batch import STRUCTURE_KEYS, run_batch, write_summary_csv

    a = api.analyse_file(EX / "infinite_chain.cif")
    data = export.summary_dict(a)
    assert set(data["structure"]) == set(a.structure_values) and data["structure"]["r1"] is None
    assert {x["code"] for x in data["alerts"]} == codes(a)
    rows = list(csv.DictReader(export.write_alerts_csv(a, tmp_path / "al.csv").open(encoding="utf-8")))
    assert {r["code"] for r in rows} == codes(a)
    entries = run_batch([EX / "demo_naphthalene.cif", EX / "infinite_chain.cif"], Limits(), workers=1)
    out = list(csv.DictReader(write_summary_csv(entries, tmp_path / "s.csv").open(encoding="utf-8")))
    assert all(k in out[0] for k in STRUCTURE_KEYS) and "alerts" in out[0]
    assert "G:" not in out[0]["alerts"]  # information-only alerts do not go into the summary


def test_cli_check(capsys, tmp_path):
    from pistak.cli import main

    assert main(["check", str(EX / "demo_naphthalene.cif"), "--csv", str(tmp_path / "a.csv")]) == 0
    out = capsys.readouterr().out
    assert "non sostituiscono checkCIF" in out and "[G]" in out and (tmp_path / "a.csv").exists()
    try:
        import spglib  # noqa: F401
    except ImportError:
        return
    assert main(["check", str(EX / "infinite_chain.cif")]) == 2  # A alerts: exit code for scripts


# ---- real structures ----------------------------------------------------------------------------------------
def test_real_cif_z_prime_and_density():
    a = api.analyse_file(COD / "4500022.cif")
    v = a.structure_values
    assert v["z_prime"] == 3 and v["density_calc"] == pytest.approx(v["density_cif"], abs=0.002)
    assert counts(a.alerts)["A"] == counts(a.alerts)["B"] == 0
    assert formula_text({e: n / 12 for e, n in cell_contents(a.cif).items()}) == "C14 H16 N6 O2"


@pytest.mark.parametrize("name, z_prime, alerts", [
    ("2311088.cif", 0.5, {"goof"}),  # naphthalene: declared GoF 0.25
    ("2012157.cif", 1.0, set()),  # pentacene
    ("2100417.cif", 0.25, set()),  # rubrene, Cmca (= Cmce) with the molecule on an mm site
    ("1550133.cif", 1.0, set()),
])
def test_cod_structures(name, z_prime, alerts):
    a = api.analyse_file(COD / name)
    v = a.structure_values
    assert v["z_prime"] == pytest.approx(z_prime)
    assert v["density_calc"] == pytest.approx(v["density_cif"], abs=0.003)
    assert {al.code for al in a.alerts if al.level in "AB"} == alerts
    assert "missedSymmetry" not in codes(a)


# ---- real structures that cannot be redistributed (DATI_REALI/ is not in the repository) ------------------
def test_real_cif_with_disorder():
    a = api.analyse_file(real("UNITCELL.cif"))
    assert {"highR1", "highWR2", "partialOccupancy", "disorderSkipped"} <= codes(a)
    assert not {"formulaHeavy", "density", "overlap"} & codes(a)  # contents with all the disordered components


# ---- short contacts: fixed thresholds or van der Waals radii sum -------------------------------------------------
def test_contacts_with_van_der_waals_radii():
    fixed = api.analyse_file(EX / "herringbone_anthracene.cif")
    vdw = api.reanalyse(fixed, fixed.limits.with_values(contact_check="vdw"))
    assert "shortContact" in codes(fixed) and "vdwContact" not in codes(fixed)
    pairs = {al.values["pair"] for al in vdw.alerts if al.code == "vdwContact"}
    assert "shortContact" not in codes(vdw) and {"C···C", "C···H", "H···H"} <= pairs  # also C···H, which the fixed thresholds do not see
    strict = api.reanalyse(fixed, fixed.limits.with_values(contact_check="vdw", contact_check_margin=1.0))
    assert sum(al.values["n"] for al in strict.alerts if al.code == "vdwContact") < \
        sum(al.values["n"] for al in vdw.alerts if al.code == "vdwContact")
    assert api.reanalyse(vdw, vdw.limits.with_values(max_dihedral=25)).alerts is vdw.alerts  # same checks: reused


def test_vdw_contacts_without_hydrogen_bonds():
    # 4500022: the H···N and H···O of hydrogen bonds are well below the radii sum, but they are not alerts
    a = api.analyse_file(COD / "4500022.cif", Limits(contact_check="vdw"))
    assert "vdwContact" not in codes(a)
    a = api.analyse_file(COD / "4500022.cif", Limits(contact_check="vdw", contact_check_margin=0.0))
    assert all("H" not in al.values["pair"] or "O" not in al.values["pair"] for al in a.alerts if al.code == "vdwContact")


def test_check_thresholds_from_toml():
    lim = Limits.from_dict({"contact_check": "VDW", "contact_check_margin": "0.3"})
    assert (lim.contact_check, lim.contact_check_margin) == ("vdw", 0.3)
    assert Limits.from_dict({"contact_check": "boh"}).contact_check == "fixed"
