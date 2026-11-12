import gc
import os
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

os.environ["PISTAK_LANG"] = "it"  # the tests check the Italian texts; the English ones switch language themselves
# the user's personal settings, last session and plugins: in a temporary folder, never the real one
os.environ["PISTAK_CONFIG_DIR"] = tempfile.mkdtemp(prefix="pistak-test-config-")

# Markers assigned per module (see [tool.pytest.ini_options] in pyproject.toml):
#   pytest -m "not gui"     the computation only, no windows (much faster)
#   pytest -m snapshot      only the comparison with the saved results
GUI_MODULES = {"test_gui", "test_layout", "test_session"}
MODULE_MARKS = {"test_parity": "parity", "test_snapshots": "snapshot"}


def pytest_addoption(parser):
    parser.addoption("--update-snapshots", action="store_true", default=False,
                     help="rewrite tests/snapshots/*.json with the current results instead of comparing them")


def pytest_collection_modifyitems(config, items):
    for item in items:
        module = item.module.__name__.rsplit(".", 1)[-1]
        # in test_i18n only the tests that open windows (fixture `app`) belong to the GUI
        if module in GUI_MODULES or "app" in getattr(item, "fixturenames", ()):
            item.add_marker(pytest.mark.gui)
        if module in MODULE_MARKS:
            item.add_marker(getattr(pytest.mark, MODULE_MARKS[module]))


@pytest.fixture(autouse=True)
def _collect_after_gui_test(request):
    """Dopo ogni test con finestre, raccoglie i cicli irraggiungibili nel thread principale: le finestre lasciate dal
    test non devono essere distrutte dal garbage collector in un thread secondario del test seguente (registro, voce 91)."""
    yield
    if request.node.get_closest_marker("gui") is not None:
        gc.collect()
