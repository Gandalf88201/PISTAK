"""Code structure (Phase 1): separate layers, exports with fixed keys, the registry as an object, user configuration
outside the sources, errors with stable codes."""
import ast
import csv
import json
import math
from pathlib import Path

import pytest

from pistak import i18n
from pistak.core.analysis import analyse_file as core_analyse
from pistak.core.config import Limits

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
PKG = ROOT / "pistak"


def _imports(folder: Path) -> dict[str, set[str]]:
    """File → imported `pistak.*` modules (also inside functions)."""
    out: dict[str, set[str]] = {}
    for f in sorted(folder.rglob("*.py")):
        if f.name.startswith("._"):
            continue
        mods: set[str] = set()
        for n in ast.walk(ast.parse(f.read_text(encoding="utf-8"))):
            if isinstance(n, ast.ImportFrom) and n.module and n.module.startswith("pistak"):
                mods.add(n.module)
            elif isinstance(n, ast.Import):
                mods |= {a.name for a in n.names if a.name.startswith("pistak")}
        out[str(f.relative_to(ROOT))] = mods
    return out


def test_core_imports_only_itself():
    bad = {f: sorted(m for m in mods if not m.startswith("pistak.core")) for f, mods in _imports(PKG / "core").items()}
    assert {f: m for f, m in bad.items() if m} == {}


@pytest.mark.parametrize("layer", ["io", "params"])
def test_io_and_params_depend_on_no_interface_or_translation(layer):
    forbidden = ("pistak.cli", "pistak.gui", "pistak.i18n", "pistak.labels")
    bad = {f: sorted(m for m in mods if m.startswith(forbidden)) for f, mods in _imports(PKG / layer).items()}
    # params translates the parameter labels (ParamSpec.name): the only exception allowed
    if layer == "params":
        bad = {f: [m for m in ms if m != "pistak.i18n"] for f, ms in bad.items()}
    assert {f: m for f, m in bad.items() if m} == {}


# ---- exports with fixed keys ---------------------------------------------------------------------------
def test_exports_identical_in_italian_and_english(tmp_path):
    from pistak import analyse_file
    from pistak.io import export
    from pistak.io.batch import run_batch, write_summary_csv
    from pistak.io.compare import compare, write_csv
    from pistak.io.stats import PairTable, describe, write_describe_csv

    a = analyse_file(EX / "herringbone_anthracene.cif")
    entries = run_batch([EX / "demo_naphthalene.cif", EX / "herringbone_anthracene.cif"], Limits(), workers=1)

    def write_all(d: Path) -> dict[str, str]:
        d.mkdir()
        t = PairTable.from_entries(entries)
        files = {
            "pairs": export.write_csv(a, d / "p.csv"), "json": export.write_json(a, d / "a.json"),
            "xyz": export.write_pair_xyz(a, a.pairs[0], d / "p.xyz"), "summary": write_summary_csv(entries, d / "s.csv"),
            "compare": write_csv(compare(entries, 0), d / "c.csv"), "table": t.write_csv(d / "t.csv"),
            "describe": write_describe_csv(describe(t, ["interplanar"]), d / "d.csv"),
        }
        return {k: p.read_text(encoding="utf-8") for k, p in files.items()}

    try:
        i18n.set_language("en")
        en = write_all(tmp_path / "en")
    finally:
        i18n.set_language("it")
    it = write_all(tmp_path / "it")
    assert en == it
    heads = {k: next(csv.reader([v.splitlines()[0]])) for k, v in it.items() if k not in ("json", "xyz")}
    assert heads["pairs"][:5] == ["n", "kind", "molecule_a", "molecule_b", "shift"] and "interplanar" in heads["pairs"]
    assert heads["summary"][:3] == ["file", "error", "space_group"] and heads["compare"][-1] == "motif_distance"
    assert heads["table"][:4] == ["file", "n", "kind", "count"]
    assert heads["describe"] == ["parameter", "n", "pairs", "mean", "std", "min", "q25", "median", "q75", "max"]
    data = json.loads(it["json"])
    assert data["format"] == "pistak-analysis" and data["version"] == 4 and data["pairs"][0]["kind"] == "stack"
    assert "interplanar=" in it["xyz"].splitlines()[1]


def test_summary_columns_have_stable_keys(tmp_path):
    from pistak.io.batch import run_batch, write_summary_csv

    entries = run_batch([EX / "herringbone_anthracene.cif"], Limits(), workers=1)
    row = list(csv.DictReader(write_summary_csv(entries, tmp_path / "s.csv").open(encoding="utf-8")))[0]
    assert row["columns"] == "column; column"


# ---- errors with codes ---------------------------------------------------------------------------------
def test_errors_with_code_and_message(tmp_path, monkeypatch):
    from pistak.core import analysis
    from pistak.core.errors import CifError, InputError, PistakError
    from pistak.errors import describe_error, error_code
    from pistak.io.batch import analyse_one

    bad = tmp_path / "b.cif"
    bad.write_text("data_x\n_cell_length_a 5\n")
    with pytest.raises(CifError) as e:
        core_analyse(bad)
    assert isinstance(e.value, PistakError) and error_code(e.value) == "noCell"
    assert describe_error(e.value) == "mancano i parametri di cella"
    entry = analyse_one(str(bad), Limits())
    assert entry.error_code == "noCell" and "cella" in entry.error

    monkeypatch.setattr(analysis, "MAX_BYTES", 10)
    with pytest.raises(InputError) as e:
        analysis.read_cif_file(EX / "demo_naphthalene.cif")
    assert e.value.code == "fileTooLarge" and describe_error(e.value).startswith("file troppo grande (> 0 MB)")


# ---- parameter registry as an object -------------------------------------------------------------------
def test_independent_registry_and_per_analysis_errors():
    from pistak.params import REGISTRY, Registry, compute_all

    reg = Registry()

    @reg.parameter("double_interplanar", "Double interplanar", unit="Å")
    def _double(ctx):
        return 2 * ctx.pair.interplanar

    @reg.parameter("broken", "Failing parameter")
    def _broken(ctx):
        raise ValueError("on purpose")

    assert "double_interplanar" not in REGISTRY and len(reg) == 2
    a = compute_all(core_analyse(EX / "demo_naphthalene.cif", Limits()), reg)
    p = a.pairs[0]
    assert a.param_keys == ["double_interplanar", "broken"]
    assert p.get("double_interplanar") == pytest.approx(2 * p.interplanar) and math.isnan(p.get("broken"))
    assert list(a.param_errors) == ["broken"] and "on purpose" in a.param_errors["broken"]
    # the default registry knows nothing of this one's errors: the next analysis starts clean
    b = compute_all(core_analyse(EX / "demo_naphthalene.cif", Limits()))
    assert b.param_errors == {} and "broken" not in b.param_keys and "interplanar" in b.param_keys


# ---- user configuration --------------------------------------------------------------------------------
def test_configuration_in_user_folder(tmp_path, monkeypatch):
    from pistak.core import config

    monkeypatch.setenv("PISTAK_CONFIG_DIR", str(tmp_path / "cfg"))
    assert config.user_config_dir() == tmp_path / "cfg"
    assert config.load_limits() == Limits()  # only the default values
    path = config.save_user_config(Limits(max_dihedral=15.0), {"gui": {"language": "en"}})
    assert path == tmp_path / "cfg" / "user.toml" and not str(path).startswith(str(ROOT / "pistak"))
    assert config.load_limits().max_dihedral == 15.0 and config.load_section("gui")["language"] == "en"
    assert config.load_section("gui")["background"] == "auto"  # from the package default.toml


def test_old_settings_copied_once(tmp_path, monkeypatch):
    from pistak.core import config

    legacy = tmp_path / "project" / "config"
    legacy.mkdir(parents=True)
    (legacy / "user.toml").write_text('[limits]\nmax_interplanar = 3.9\n', encoding="utf-8")
    (legacy / "last.pistak").write_text("{}", encoding="utf-8")
    target = tmp_path / "system"
    monkeypatch.delenv("PISTAK_CONFIG_DIR")
    monkeypatch.setattr(config, "LEGACY_CONFIG_DIR", legacy)
    monkeypatch.setattr(config, "_system_config_dir", lambda: target)
    assert config.user_config_dir() == target
    assert config.load_limits().max_interplanar == 3.9 and (target / "last.pistak").exists()
    (legacy / "user.toml").write_text('[limits]\nmax_interplanar = 3.0\n', encoding="utf-8")
    assert config.load_limits().max_interplanar == 3.9  # the new copy is no longer overwritten


def test_config_command(capsys):
    from pistak.cli import main

    assert main(["config"]) == 0
    out = capsys.readouterr().out
    assert "user.toml" in out and "default.toml" in out and "plugins" in out
