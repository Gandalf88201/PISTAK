"""Colori condivisi dalle viste: seguono il tema chiaro/scuro del sistema (o la scelta in config)."""
from __future__ import annotations

from PySide6.QtGui import QColor, QGuiApplication, QPalette

from pistak.core.config import load_section
from pistak.i18n import tr

KIND_COLORS = {"stack": (46, 160, 67), "parallel": (224, 160, 30), "tilted": (170, 100, 210)}
KIND_LABELS = {"stack": "π-stack", "parallel": "parallela", "tilted": "inclinata"}  # testi italiani, da passare a tr()


def kind_label(kind: str) -> str:
    """Il nome del tipo di coppia (stack | parallel | tilted) nella lingua corrente."""
    return tr(KIND_LABELS[kind])


def is_dark() -> bool:
    mode = str(load_section("gui").get("background", "auto"))
    if mode in ("light", "dark"):
        return mode == "dark"
    app = QGuiApplication.instance()
    return bool(app and app.palette().color(QPalette.ColorRole.Window).lightness() < 128)


def background() -> tuple[int, int, int]:
    return (30, 32, 36) if is_dark() else (255, 255, 255)


def foreground() -> tuple[int, int, int]:
    return (220, 222, 226) if is_dark() else (30, 32, 36)


def qcolor(rgb: tuple[int, int, int], alpha: int = 255) -> QColor:
    return QColor(rgb[0], rgb[1], rgb[2], alpha)
