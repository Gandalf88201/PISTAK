import os
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

os.environ["PISTAK_LANG"] = "it"  # i test controllano i testi italiani; quelli inglesi cambiano lingua da soli
# impostazioni personali, ultima sessione e plugin dell'utente: in una cartella temporanea, mai quella vera
os.environ["PISTAK_CONFIG_DIR"] = tempfile.mkdtemp(prefix="pistak-test-config-")

# Marcatori assegnati per modulo (vedi [tool.pytest.ini_options] in pyproject.toml):
#   pytest -m "not gui"     solo il calcolo, senza finestre (molto più veloce)
#   pytest -m snapshot      solo il confronto con i risultati salvati
GUI_MODULES = {"test_gui", "test_layout", "test_session"}
MODULE_MARKS = {"test_parity": "parity", "test_snapshots": "snapshot"}


def pytest_addoption(parser):
    parser.addoption("--update-snapshots", action="store_true", default=False,
                     help="riscrive tests/snapshots/*.json con i risultati attuali invece di confrontarli")


def pytest_collection_modifyitems(config, items):
    for item in items:
        module = item.module.__name__.rsplit(".", 1)[-1]
        # in test_i18n solo i test che aprono finestre (fixture `app`) sono della GUI
        if module in GUI_MODULES or "app" in getattr(item, "fixturenames", ()):
            item.add_marker(pytest.mark.gui)
        if module in MODULE_MARKS:
            item.add_marker(getattr(pytest.mark, MODULE_MARKS[module]))
