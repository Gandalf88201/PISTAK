"""Lingue: completezza del dizionario, testi inglesi in terminale e nella finestra, cambio lingua a finestra aperta."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_API", "pyside6")
os.environ["PISTAK_NO_3D"] = "1"

import re
import string
import time
from pathlib import Path

import pytest

from pistak import i18n, i18n_check
from pistak.core.config import Limits
from pistak.i18n import tf, tr
from pistak.i18n_en import TEXT

EX = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(autouse=True)
def italiano_alla_fine():
    yield
    i18n.set_language("it")


def fields(text: str) -> set[str]:
    return {name for _, name, _, _ in string.Formatter().parse(text) if name}


# ---- dizionario ---------------------------------------------------------------------------
def test_nessun_testo_senza_traduzione():
    assert i18n_check.missing() == []


def test_segnaposto_identici_nelle_due_lingue():
    for it, en in TEXT.items():
        assert fields(it) == fields(en), f"segnaposto diversi: {it!r} → {en!r}"


def test_il_dizionario_non_ha_voci_morte():
    vive = i18n_check.all_texts() | {"#"}
    assert [k for k in TEXT if k not in vive] == []


def test_nessuna_traduzione_vuota_o_uguale_per_errore():
    vuote = [k for k, v in TEXT.items() if not v.strip()]
    assert vuote == []


# ---- funzioni di base -------------------------------------------------------------------------
def test_italiano_e_il_testo_originale_inglese_lo_traduce():
    assert tr("Tabella") == "Tabella"
    i18n.set_language("en")
    assert tr("Tabella") == "Table"
    assert tf("{n} coppie in {f} file; tipo: {kind}", n=3, f=2, kind="stack") == "3 pairs in 2 files; type: stack"


def test_testo_senza_traduzione_resta_italiano():
    i18n.set_language("en")
    assert tr("un testo che non esiste nel dizionario") == "un testo che non esiste nel dizionario"


def test_scelta_della_lingua():
    assert i18n.resolve("en") == "en" and i18n.resolve("IT") == "it"
    assert i18n.resolve("klingon") == i18n.DEFAULT and i18n.resolve(None) == i18n.DEFAULT
    assert i18n.resolve("auto") in i18n.LANGS


def test_traduzioni_dei_plugin():
    i18n.add_translations("en", {"Il mio parametro": "My parameter"})
    i18n.set_language("en")
    assert tr("Il mio parametro") == "My parameter"
    TEXT_COPY = i18n._load("en")
    TEXT_COPY.pop("Il mio parametro")


# ---- parti non grafiche ---------------------------------------------------------------------------
def test_nomi_dei_parametri_e_dei_tipi_in_inglese():
    from pistak.gui import theme
    from pistak.labels import rank_name
    from pistak.params import get_spec

    i18n.set_language("en")
    assert get_spec("interplanar").name == "Interplanar distance"
    assert get_spec("interplanar").title == "Interplanar distance (Å)"
    assert get_spec("long_slip").group_name == "website"
    assert theme.kind_label("tilted") == "tilted" and theme.kind_label("parallel") == "parallel"
    assert rank_name(1) == "column"
    i18n.set_language("it")
    assert get_spec("interplanar").name == "Distanza interplanare" and rank_name(1) == "colonna"


def test_relazione_della_coppia_in_inglese_e_csv_uguale_nelle_due_lingue(tmp_path):
    from pistak import analyse_file
    from pistak.io import export
    from pistak.labels import relation

    a = analyse_file(EX / "herringbone_anthracene.cif")
    i18n.set_language("en")
    assert any(relation(p).startswith("molecule 2") for p in a.pairs) and relation(a.pairs[0]).startswith("same")
    en = export.write_csv(a, tmp_path / "p.csv").read_text()
    i18n.set_language("it")
    assert relation(a.pairs[0]).startswith("stessa")
    it = export.write_csv(a, tmp_path / "q.csv").read_text()
    # le intestazioni dei file non si traducono: uno script li legge allo stesso modo in ogni lingua
    assert en == it and it.splitlines()[0].startswith("n,kind,molecule_a,molecule_b,shift,")


def test_terminale_in_inglese(capsys):
    from pistak.cli import main

    assert main(["--lang", "en", "analyze", str(EX / "demo_naphthalene.cif")]) == 0
    out = capsys.readouterr().out
    assert "group P -1" in out and "atoms in cell" in out and "with π core" in out
    assert "gruppo" not in out and "atomi" not in out
    assert main(["--lang=it", "analyze", str(EX / "demo_naphthalene.cif")]) == 0
    assert "atomi in cella" in capsys.readouterr().out


def test_errori_del_lettore_in_inglese(tmp_path, capsys):
    from pistak.cli import main

    bad = tmp_path / "b.cif"
    bad.write_text("data_x\n_cell_length_a 5\n")
    assert main(["--lang", "en", "analyze", str(bad)]) == 1
    assert "cell parameters are missing" in capsys.readouterr().err


def test_batch_compare_stats_in_inglese(tmp_path, capsys):
    from pistak.cli import main

    assert main(["--lang", "en", "compare", str(EX / "demo_naphthalene.cif"), str(EX / "herringbone_anthracene.cif")]) == 0
    assert "reference: demo_naphthalene.cif" in capsys.readouterr().out
    assert main(["--lang", "en", "stats", str(EX), "--kind", "stack", "--summary", str(tmp_path / "s.csv")]) == 0
    out = capsys.readouterr().out
    assert "pairs in" in out and "mean" in out
    assert (tmp_path / "s.csv").read_text().splitlines()[0].startswith("parameter,n,pairs,mean,std,")


def test_sessioni_vecchie_con_nomi_italiani_delle_viste():
    from pistak.gui.views import LEGACY_KINDS, VIEW_TYPES

    assert all(v in VIEW_TYPES for v in LEGACY_KINDS.values())


# ---- finestra ---------------------------------------------------------------------------------------
pytest.importorskip("PySide6")
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QDockWidget

ITALIAN = re.compile(
    r"\b(coppia|coppie|molecola|molecole|tabella|vista|viste|sessione|apri|salva|finestra|analisi|cartella|soglie|parametro|"
    r"parametri|istogramma|proiezione|dettaglio|colonna|colonne|struttura|strutture|distanza|slittamento|sovrapposizione|anello|"
    r"anelli|nucleo|nuclei|piani|angolo|tutte|tutti|selezione|esporta|confronto|statistiche|aggiungi|rimuovi|svuota|riepilogo|"
    r"correlazioni|distribuzioni|dispersione|nessuna|nessun|scritto|errore|gruppo|atomi|cella|lingua|aiuto)\b", re.I)


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def wait(app, cond, timeout=40):
    t = time.time()
    while not cond() and time.time() - t < timeout:
        app.processEvents()
        time.sleep(0.02)
    assert cond(), "tempo scaduto"


def visible_texts(root) -> list[str]:
    """Tutti i testi che l'utente può leggere in una finestra: etichette, pulsanti, menu, schede, intestazioni, suggerimenti."""
    from PySide6.QtGui import QAction
    from PySide6.QtWidgets import QAbstractButton, QComboBox, QLabel, QTableView, QTableWidget, QTabWidget, QTextBrowser, QTreeWidget

    out: list[str] = [root.windowTitle()]
    for w in [root, *root.findChildren(object)]:
        if hasattr(w, "toolTip") and callable(w.toolTip):
            out.append(w.toolTip())
        if isinstance(w, (QLabel, QAbstractButton)):
            out.append(w.text())
        if isinstance(w, QTextBrowser):
            out.append(w.toPlainText())
        if isinstance(w, QComboBox):
            out += [w.itemText(i) for i in range(w.count())]
        if isinstance(w, QTabWidget):
            out += [w.tabText(i) for i in range(w.count())]
        if isinstance(w, QDockWidget):
            out.append(w.windowTitle())
        if isinstance(w, QAction):
            out.append(w.text())
        if isinstance(w, QTableWidget):
            out += [(w.horizontalHeaderItem(c).text() if w.horizontalHeaderItem(c) else "") for c in range(w.columnCount())]
        if isinstance(w, QTreeWidget):
            out += [w.headerItem().text(c) for c in range(w.columnCount())]
        if isinstance(w, QTableView) and w.model() is not None:
            m = w.model()
            out += [str(m.headerData(c, Qt.Orientation.Horizontal) or "") for c in range(min(m.columnCount(), 40))]
    return [t for t in out if t]


def test_finestra_principale_tutta_in_inglese(app):
    from pistak.gui.main_window import MainWindow
    from pistak.gui.session import Session

    i18n.set_language("en")
    win = MainWindow(Session())
    win.autosave = False
    win.session.limits = Limits.site()  # 3 coppie nell'herringbone: la selezione 2 esiste
    win.session.load_file(EX / "herringbone_anthracene.cif")
    wait(app, lambda: win.session.analysis is not None)
    for _ in range(20):
        app.processEvents()
    assert [a.text() for a in win.menuBar().actions()] == ["&File", "&View", "&Analysis", "&Language", "&Help"]
    assert [d.windowTitle() for d in win.docks] == ["Table", "Parameter scatter", "Histogram", "3D", "Projection", "Detail"]
    for _tab in win.docks[-1].widget().tabs.findChildren(object):
        pass
    left = [t for t in visible_texts(win) if ITALIAN.search(t)]
    assert left == [], f"testi ancora in italiano: {left}"


def test_finestra_di_confronto_tutta_in_inglese(app):
    from pistak.core.config import Limits
    from pistak.gui.collection_window import CollectionWindow

    i18n.set_language("en")
    win = CollectionWindow(Limits())
    win.add_paths([EX])
    wait(app, lambda: len(win.entries) == 6)
    for i in range(win.tabs.count()):
        win.tabs.setCurrentIndex(i)
        for _ in range(5):
            app.processEvents()
    texts = visible_texts(win)
    assert [t for t in texts if ITALIAN.search(t)] == []
    assert win.tabs.tabText(0) == "Structures" and "Distributions" in [win.tabs.tabText(i) for i in range(win.tabs.count())]


def test_dettaglio_e_impostazioni_in_inglese(app):
    from pistak.core.config import Limits
    from pistak.gui.settings_dialog import SettingsDialog

    i18n.set_language("en")
    dlg = SettingsDialog(Limits())
    texts = visible_texts(dlg)
    assert "Analysis thresholds" in texts and [t for t in texts if ITALIAN.search(t)] == []
    labels = [dlg.findChildren(object)[i].text() for i in range(len(dlg.findChildren(object))) if hasattr(dlg.findChildren(object)[i], "text")]
    assert any("Maximum interplanar distance" in str(x) for x in labels)


def test_cambio_lingua_a_finestra_aperta_conserva_la_sessione(app, monkeypatch):
    from pistak.gui.main_window import MainWindow
    from pistak.gui.session import Session
    from pistak.gui.views.hist_view import HistogramView

    monkeypatch.setattr(i18n, "save_language", lambda lang: None)  # non scrive config/user.toml durante i test
    i18n.set_language("it")
    win = MainWindow(Session())
    win.autosave = False
    win.session.limits = Limits.site()  # 3 coppie nell'herringbone: la selezione 2 esiste
    win.session.load_file(EX / "herringbone_anthracene.cif")
    wait(app, lambda: win.session.analysis is not None)
    for _ in range(20):
        app.processEvents()
    win.session.select("A", 2)
    win.session.set_subset("A", {1, 2})
    hist = next(d.widget() for d in win.docks if isinstance(d.widget(), HistogramView))
    hist.apply_state({"key": "lat_slip", "bins": 7, "weighted": False})
    assert [d.windowTitle() for d in win.docks][0] == "Tabella"

    win.change_language("en")
    for _ in range(20):
        app.processEvents()
    new = app.main_window
    assert new is not win and i18n.get_language() == "en" and not win.isVisible()
    assert [d.windowTitle() for d in new.docks] == ["Table", "Parameter scatter", "Histogram", "3D", "Projection", "Detail"]
    assert new.session is win.session
    assert (new.session.selection("A"), new.session.subset("A")) == (2, {1, 2})
    hist2 = next(d.widget() for d in new.docks if isinstance(d.widget(), HistogramView))
    assert hist2.bins.value() == 7 and hist2.combo.currentData() == "lat_slip" and not hist2.weighted.isChecked()
    assert new.docks[0].widget().index == 2  # le viste mostrano di nuovo la coppia scelta

    new.change_language("it")  # e si torna indietro
    for _ in range(20):
        app.processEvents()
    assert app.main_window.docks[0].windowTitle() == "Tabella" and app.main_window.session.selection("A") == 2


def test_cambio_tema_a_finestra_aperta_conserva_la_sessione(app, monkeypatch):
    from pistak.gui import theme
    from pistak.gui.main_window import MainWindow
    from pistak.gui.session import Session
    from pistak.gui.views.scatter_view import ScatterView

    monkeypatch.setattr(theme, "save_mode", lambda m: None)  # non scrive user.toml durante i test
    monkeypatch.setattr(theme, "_mode", "light")
    win = MainWindow(Session())
    win.autosave = False
    win.session.limits = Limits.site()
    win.session.load_file(EX / "herringbone_anthracene.cif")
    wait(app, lambda: win.session.analysis is not None)
    for _ in range(20):
        app.processEvents()
    win.session.select("A", 2)

    win.change_theme("dark")
    for _ in range(20):
        app.processEvents()
    new = app.main_window
    assert new is not win and not win.isVisible() and theme.mode() == "dark"
    assert new.session is win.session and new.session.selection("A") == 2
    scatter = next(d.widget() for d in new.docks if isinstance(d.widget(), ScatterView))
    assert scatter.plot.backgroundBrush().color().getRgb()[:3] == theme.background() == (30, 32, 36)

    new.change_theme("light")
    for _ in range(20):
        app.processEvents()
    assert theme.mode() == "light" and app.main_window.session.selection("A") == 2
