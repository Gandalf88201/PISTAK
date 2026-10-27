"""Controllo delle traduzioni: elenca i testi italiani del programma e quelli che mancano in inglese.

    python -m pistak.i18n_check            stampa i testi senza traduzione (esce con 1 se ce ne sono)
    python -m pistak.i18n_check --all      stampa tutti i testi traducibili

Serve a chi aggiunge testi: scrivili in italiano dentro `tr("...")` / `tf("...", ...)`, poi aggiungi qui la coppia in
`pistak/i18n_en.py`. I parametri dei plugin (etichette, descrizioni, gruppi) si traducono con `add_translations`.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def literal_texts() -> dict[str, list[str]]:
    """Testo → elenco dei punti del codice dove compare come primo argomento di tr()/tf()."""
    found: dict[str, list[str]] = {}
    for f in sorted(ROOT.rglob("*.py")):
        if f.name.startswith("._") or f.name == "i18n_en.py":
            continue
        for n in ast.walk(ast.parse(f.read_text(encoding="utf-8"))):
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in ("tr", "tf") and n.args:
                a = n.args[0]
                if isinstance(a, ast.Constant) and isinstance(a.value, str):
                    found.setdefault(a.value, []).append(f"{f.relative_to(ROOT.parent)}:{n.lineno}")
    return found


def constant_texts() -> set[str]:
    """I testi tenuti in costanti e tabelle e tradotti solo al momento dell'uso."""
    from dataclasses import fields

    from pistak import params
    from pistak.api import ensure_plugins
    from pistak.core.config import Limits
    from pistak.core.ringpairs import HEADERS
    from pistak.core.stacks import RANK_NAMES
    from pistak.errors import ERROR_TEXT
    from pistak.gui import collection_window, theme
    from pistak.gui.views import VIEW_TYPES
    from pistak.gui.views.table_view import FIXED
    from pistak.gui.views.view3d import STYLES

    ensure_plugins()
    out: set[str] = set(ERROR_TEXT.values()) | set(HEADERS.values()) | set(RANK_NAMES.values()) | set(theme.KIND_LABELS.values())
    out |= set(theme.MODES.values())
    out |= set(FIXED) | set(STYLES) | {label for label, _ in collection_window.KIND_FILTER}
    out |= {f.metadata["label"] for f in fields(Limits)} | {cls.TITLE for cls in VIEW_TYPES.values()}
    out |= {text for f in fields(Limits) for text in f.metadata.get("choices", {}).values()}
    from pistak.gui import periodic_dialog, qc_dialog
    from pistak.gui.views.detail_view import GEOMETRY_KINDS
    from pistak.labels import ALERT_LEVELS, ALERT_TEXT, NO_PAIRS_HINT, NO_PAIRS_TEXT

    out |= set(qc_dialog.METHOD_TEXT.values()) | set(qc_dialog.DISP_TEXT.values())
    out |= set(periodic_dialog.TASK_TEXT.values()) | set(periodic_dialog.DISP_TEXT.values())
    out |= set(ALERT_TEXT.values()) | set(ALERT_LEVELS.values()) | set(NO_PAIRS_TEXT.values()) | {NO_PAIRS_HINT} | {title for title, _u, _d in GEOMETRY_KINDS.values()}
    for s in [*params.specs(), *params.structure_specs()]:
        out |= {t for t in (s.label, s.description, s.group) if t}
    return out


def all_texts() -> set[str]:
    return set(literal_texts()) | constant_texts()


def missing() -> list[str]:
    from pistak.i18n_en import TEXT

    return sorted(t for t in all_texts() if t not in TEXT)


def main() -> int:
    if "--all" in sys.argv:
        print("\n".join(sorted(all_texts())))
        return 0
    miss = missing()
    for t in miss:
        print(repr(t))
    print(f"{len(miss)} testi senza traduzione inglese su {len(all_texts())}", file=sys.stderr)
    return 1 if miss else 0


if __name__ == "__main__":
    sys.exit(main())
