"""Finestre ridimensionabili e dentro lo schermo, anche piccolo."""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_API", "pyside6")
os.environ["PISTAK_NO_3D"] = "1"

import pytest

pytest.importorskip("PySide6")
from PySide6.QtCore import QRect
from PySide6.QtWidgets import QApplication, QPushButton, QWidget

from pistak.core.config import Limits
from pistak.gui.collection_window import CollectionWindow
from pistak.gui.flow_layout import FlowLayout
from pistak.gui.main_window import MainWindow
from pistak.gui.screen import clamp_into_screen, fit_to_screen
from pistak.gui.session import Session
from pistak.gui.settings_dialog import SettingsDialog


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def test_fit_to_screen_riduce_e_centra(app):
    w = QWidget()
    avail = QRect(0, 25, 1024, 600)
    fit_to_screen(w, 1500, 920, avail=avail)
    assert w.width() <= 1024 and w.height() <= 600
    assert w.minimumWidth() <= w.width() and w.minimumHeight() <= w.height()
    assert avail.contains(w.geometry().topLeft()) and w.x() + w.width() <= 1024 and w.y() + w.height() <= 625
    big = QWidget()  # schermo grande: si usa la dimensione preferita
    fit_to_screen(big, 1500, 920, avail=QRect(0, 0, 2560, 1440))
    assert (big.width(), big.height()) == (1500, 920)


def test_clamp_into_screen_riporta_dentro_una_finestra_troppo_grande(app):
    w = QWidget()
    w.setMinimumSize(100, 100)
    w.resize(3000, 2000)
    w.move(2500, 1800)
    clamp_into_screen(w, avail=QRect(0, 0, 1280, 720))
    assert w.width() <= 1280 and w.height() <= 720
    assert 0 <= w.x() and w.x() + w.width() <= 1280 and 0 <= w.y() and w.y() + w.height() <= 720


def test_flow_layout_va_a_capo(app):
    host = QWidget()
    lay = FlowLayout(host)
    for i in range(6):
        b = QPushButton(f"Pulsante molto lungo numero {i}")
        lay.addWidget(b)
    wide = lay.heightForWidth(4000)
    narrow = lay.heightForWidth(250)
    assert narrow > 3 * wide  # stretto: una riga per pulsante
    assert lay.minimumSize().width() < 300  # la larghezza minima è quella del pulsante più largo, non della somma


def test_la_finestra_principale_puo_stare_in_uno_schermo_piccolo(app):
    win = MainWindow(Session())
    win.show()
    for _ in range(5):
        app.processEvents()
    hint = win.minimumSizeHint()
    assert hint.width() <= 700 and hint.height() <= 500, f"finestra troppo rigida: {hint}"
    win.resize(720, 480)
    for _ in range(5):
        app.processEvents()
    assert (win.width(), win.height()) == (720, 480)  # il layout non impone più di così
    win.resize(1800, 1100)  # e si può ingrandire
    app.processEvents()
    assert win.width() == 1800


def test_finestra_principale_non_supera_lo_schermo(app):
    avail = app.primaryScreen().availableGeometry()
    win = MainWindow(Session())
    assert win.width() <= avail.width() and win.height() <= avail.height()


def test_finestra_di_confronto_e_impostazioni_si_adattano(app):
    avail = app.primaryScreen().availableGeometry()
    for w in (CollectionWindow(Limits()), SettingsDialog(Limits())):
        w.show()
        app.processEvents()
        assert w.width() <= avail.width() and w.height() <= avail.height()
        assert w.minimumSizeHint().width() < 700 and w.minimumSizeHint().height() < 500
