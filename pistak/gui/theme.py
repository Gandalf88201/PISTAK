"""Colori condivisi dalle viste e tema chiaro/scuro.

Il tema si sceglie da Vista ▸ Tema (o con `background` nella sezione [gui] di user.toml): «auto» segue il sistema,
«light» e «dark» lo impongono a tutta l'applicazione, finestre e grafici. `apply()` lo applica a Qt e a pyqtgraph.
"""
from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QGuiApplication, QPalette, QPen

from pistak.core.config import load_section, save_user_config
from pistak.i18n import tr

KIND_COLORS = {"stack": (46, 160, 67), "parallel": (224, 160, 30), "tilted": (170, 100, 210)}
KIND_LABELS = {"stack": "π-stack", "parallel": "parallela", "tilted": "inclinata"}  # testi italiani, da passare a tr()
MODES = {"auto": "Automatico (come il sistema)", "light": "Chiaro", "dark": "Scuro"}  # testi italiani, da passare a tr()

_mode: str | None = None  # letto da config alla prima richiesta


def kind_label(kind: str) -> str:
    """Il nome del tipo di coppia (stack | parallel | tilted) nella lingua corrente."""
    return tr(KIND_LABELS[kind])


def mode() -> str:
    """Il tema scelto: auto | light | dark."""
    global _mode
    if _mode is None:
        m = str(load_section("gui").get("background", "auto")).strip().lower()
        _mode = m if m in MODES else "auto"
    return _mode


def set_mode(m: str) -> None:
    """Cambia il tema per questa esecuzione e lo applica (non lo salva: vedi `save_mode`)."""
    global _mode
    if m not in MODES:
        raise ValueError(m)
    _mode = m
    apply()


def save_mode(m: str) -> None:
    """Ricorda il tema in user.toml per le prossime esecuzioni."""
    save_user_config(sections={"gui": {"background": m}})


def is_dark() -> bool:
    m = mode()
    if m != "auto":
        return m == "dark"
    app = QGuiApplication.instance()
    return bool(app and app.palette().color(QPalette.ColorRole.Window).lightness() < 128)


def background() -> tuple[int, int, int]:
    return (30, 32, 36) if is_dark() else (255, 255, 255)


def foreground() -> tuple[int, int, int]:
    return (220, 222, 226) if is_dark() else (30, 32, 36)


def qcolor(rgb: tuple[int, int, int], alpha: int = 255) -> QColor:
    return QColor(rgb[0], rgb[1], rgb[2], alpha)


def point_pen(rgb: tuple[int, int, int], alpha: int = 255) -> QPen:
    """Bordo dei punti degli scatter: un pixel, una tonalità più scura del riempimento. Penna cosmetica: una penna
    normale pyqtgraph la ingrandisce insieme al simbolo e dietro ogni punto resta un quadrato grigio."""
    pen = QPen(qcolor(rgb, alpha).darker(150))
    pen.setWidthF(1.0)
    pen.setCosmetic(True)
    return pen


def apply() -> None:
    """Porta il tema scelto su finestre (Qt) e grafici (pyqtgraph). I grafici già creati non cambiano: la finestra
    principale si ricostruisce dopo un cambio di tema."""
    import pyqtgraph as pg

    app = QGuiApplication.instance()
    m = mode()
    if app is not None:
        hints = app.styleHints()
        if hasattr(hints, "setColorScheme"):  # Qt ≥ 6.8: anche i controlli nativi (macOS) cambiano aspetto
            if m == "auto":
                hints.unsetColorScheme()
            else:
                hints.setColorScheme(Qt.ColorScheme.Dark if m == "dark" else Qt.ColorScheme.Light)
        elif m != "auto":  # Qt più vecchio: stile Fusion con una tavolozza propria
            app.setStyle("Fusion")
            app.setPalette(_palette(m == "dark"))
    pg.setConfigOptions(background=background(), foreground=foreground())


def _palette(dark: bool) -> QPalette:
    if not dark:
        return QPalette(QColor(239, 239, 239))
    p = QPalette()
    R = QPalette.ColorRole
    for role, rgb in ((R.Window, (45, 47, 52)), (R.WindowText, (220, 222, 226)), (R.Base, (30, 32, 36)),
                      (R.AlternateBase, (40, 42, 47)), (R.Text, (220, 222, 226)), (R.Button, (55, 57, 62)),
                      (R.ButtonText, (220, 222, 226)), (R.ToolTipBase, (55, 57, 62)), (R.ToolTipText, (220, 222, 226)),
                      (R.Highlight, (47, 111, 214)), (R.HighlightedText, (255, 255, 255)), (R.PlaceholderText, (140, 142, 146))):
        p.setColor(role, QColor(*rgb))
    return p


def plot_widget():
    """Un grafico pyqtgraph con lo sfondo del tema. Il menu contestuale tiene le voci della vista (View All, assi,
    Export) ma non «Plot Options»: trasformazioni, FFT, medie e trasparenza valgono per le curve e sui punti e sulle
    barre di PISTAK danno solo errori."""
    import pyqtgraph as pg

    w = pg.PlotWidget(background=background())
    w.getPlotItem().setMenuEnabled(False, enableViewBoxMenu=True)
    return w
