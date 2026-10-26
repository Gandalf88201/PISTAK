"""Dimensioni delle finestre rispetto allo schermo: mai più grandi dell'area disponibile, sempre ridimensionabili."""
from __future__ import annotations

from PySide6.QtCore import QRect
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QWidget


def available_rect(widget: QWidget | None = None) -> QRect:
    """L'area utilizzabile dello schermo su cui sta il widget (esclusi barra dei menu e Dock)."""
    screen = (widget.screen() if widget is not None else None) or QGuiApplication.primaryScreen()
    return screen.availableGeometry() if screen is not None else QRect(0, 0, 1280, 720)


def fit_to_screen(widget: QWidget, width: int, height: int, min_width: int = 640, min_height: int = 420,
                  avail: QRect | None = None, fraction: float = 0.94) -> tuple[int, int]:
    """Dimensiona il widget a `width`×`height` ma al massimo `fraction` dello schermo, lo centra e gli dà una
    dimensione minima piccola: l'utente può sempre rimpicciolirlo o ingrandirlo."""
    a = avail or available_rect(widget)
    w, h = min(width, int(a.width() * fraction)), min(height, int(a.height() * fraction))
    widget.setMinimumSize(min(min_width, w), min(min_height, h))
    widget.resize(w, h)
    widget.move(a.x() + (a.width() - w) // 2, a.y() + (a.height() - h) // 2)
    return w, h


def clamp_into_screen(widget: QWidget, avail: QRect | None = None) -> None:
    """Se la finestra (per esempio ripristinata da una sessione salvata su uno schermo più grande) esce dallo
    schermo, la riduce e la sposta dentro."""
    a = avail or available_rect(widget)
    frame = widget.frameGeometry()
    w, h = min(frame.width(), a.width()), min(frame.height(), a.height())
    x = min(max(frame.x(), a.x()), a.x() + a.width() - w)
    y = min(max(frame.y(), a.y()), a.y() + a.height() - h)
    if (w, h) != (frame.width(), frame.height()):
        widget.resize(max(widget.minimumWidth(), w - (frame.width() - widget.width())),
                      max(widget.minimumHeight(), h - (frame.height() - widget.height())))
    widget.move(x, y)
