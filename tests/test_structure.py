"""Struttura del codice (Fase 1): strati separati, export con chiavi fisse, registro come oggetto, configurazione
dell'utente fuori dai sorgenti, errori con codici stabili."""
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
    """File → moduli `pistak.*` importati (anche dentro le funzioni)."""
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


def test_il_nucleo_importa_solo_se_stesso():
    bad = {f: sorted(m for m in mods if not m.startswith("pistak.core")) for f, mods in _imports(PKG / "core").items()}
    assert {f: m for f, m in bad.items() if m} == {}


@pytest.mark.parametrize("layer", ["io", "params"])
def test_io_e_params_non_dipendono_da_interfacce_ne_traduzioni(layer):
    forbidden = ("pistak.cli", "pistak.gui", "pistak.i18n", "pistak.labels")
    bad = {f: sorted(m for m in mods if m.startswith(forbidden)) for f, mods in _imports(PKG / layer).items()}
    # params traduce le etichette dei parametri (ParamSpec.name): l'unica eccezione ammessa
    if layer == "params":
        bad = {f: [m for m in ms if m != "pistak.i18n"] for f, ms in bad.items()}
    assert {f: m for f, m in bad.items() if m} == {}


# ---- export con chiavi fisse ---------------------------------------------------------------------------
def test_export_identici_in_italiano_e_in_inglese(tmp_path):
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
    assert data["format"] == "pistak-analysis" and data["version"] == 2 and data["pairs"][0]["kind"] == "stack"
    assert "interplanar=" in it["xyz"].splitlines()[1]


def test_colonne_nel_riepilogo_con_chiavi_stabili(tmp_path):
    from pistak.io.batch import run_batch, write_summary_csv

    entries = run_batch([EX / "herringbone_anthracene.cif"], Limits(), workers=1)
    row = list(csv.DictReader(write_summary_csv(entries, tmp_path / "s.csv").open(encoding="utf-8")))[0]
    assert row["columns"] == "column; column"


# ---- errori con codici ---------------------------------------------------------------------------------
def test_errori_con_codice_e_messaggio(tmp_path, monkeypatch):
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


# ---- registro dei parametri come oggetto ---------------------------------------------------------------
def test_registro_indipendente_e_errori_per_analisi():
    from pistak.params import REGISTRY, Registry, compute_all

    reg = Registry()

    @reg.parameter("doppia_interplanare", "Doppia interplanare", unit="Å")
    def _doppia(ctx):
        return 2 * ctx.pair.interplanar

    @reg.parameter("rotto", "Parametro che fallisce")
    def _rotto(ctx):
        raise ValueError("apposta")

    assert "doppia_interplanare" not in REGISTRY and len(reg) == 2
    a = compute_all(core_analyse(EX / "demo_naphthalene.cif", Limits()), reg)
    p = a.pairs[0]
    assert a.param_keys == ["doppia_interplanare", "rotto"]
    assert p.get("doppia_interplanare") == pytest.approx(2 * p.interplanar) and math.isnan(p.get("rotto"))
    assert list(a.param_errors) == ["rotto"] and "apposta" in a.param_errors["rotto"]
    # il registro predefinito non sa niente degli errori di questo: la prossima analisi parte pulita
    b = compute_all(core_analyse(EX / "demo_naphthalene.cif", Limits()))
    assert b.param_errors == {} and "rotto" not in b.param_keys and "interplanar" in b.param_keys


# ---- configurazione dell'utente ------------------------------------------------------------------------
def test_configurazione_nella_cartella_dell_utente(tmp_path, monkeypatch):
    from pistak.core import config

    monkeypatch.setenv("PISTAK_CONFIG_DIR", str(tmp_path / "cfg"))
    assert config.user_config_dir() == tmp_path / "cfg"
    assert config.load_limits() == Limits()  # solo i valori predefiniti
    path = config.save_user_config(Limits(max_dihedral=15.0), {"gui": {"language": "en"}})
    assert path == tmp_path / "cfg" / "user.toml" and not str(path).startswith(str(ROOT / "pistak"))
    assert config.load_limits().max_dihedral == 15.0 and config.load_section("gui")["language"] == "en"
    assert config.load_section("gui")["background"] == "auto"  # dal default.toml del pacchetto


def test_le_vecchie_impostazioni_si_copiano_una_volta(tmp_path, monkeypatch):
    from pistak.core import config

    legacy = tmp_path / "progetto" / "config"
    legacy.mkdir(parents=True)
    (legacy / "user.toml").write_text('[limits]\nmax_interplanar = 3.9\n', encoding="utf-8")
    (legacy / "last.pistak").write_text("{}", encoding="utf-8")
    target = tmp_path / "sistema"
    monkeypatch.delenv("PISTAK_CONFIG_DIR")
    monkeypatch.setattr(config, "LEGACY_CONFIG_DIR", legacy)
    monkeypatch.setattr(config, "_system_config_dir", lambda: target)
    assert config.user_config_dir() == target
    assert config.load_limits().max_interplanar == 3.9 and (target / "last.pistak").exists()
    (legacy / "user.toml").write_text('[limits]\nmax_interplanar = 3.0\n', encoding="utf-8")
    assert config.load_limits().max_interplanar == 3.9  # la copia nuova non viene più sovrascritta


def test_comando_config(capsys):
    from pistak.cli import main

    assert main(["config"]) == 0
    out = capsys.readouterr().out
    assert "user.toml" in out and "default.toml" in out and "plugins" in out
