"""Caricamento dei plugin dell'utente.

Ogni file `.py` nella cartella `plugins/` (e in quelle elencate in $PISTAK_PLUGINS, separate da ':')
viene eseguito all'avvio. Basta che usi `@parameter(...)` per aggiungere parametri. I file che iniziano
con '_' o '.' sono ignorati (su volumi exFAT macOS crea file '._nome' che non sono codice).
"""
from __future__ import annotations

import importlib.util
import os
import sys
import traceback
from pathlib import Path

from pistak.core.config import ROOT
from pistak.params import registry

PLUGIN_DIRS = [ROOT / "plugins"]
#: Problemi incontrati caricando i plugin: nome file → messaggio.
load_errors: dict[str, str] = {}


def plugin_dirs() -> list[Path]:
    extra = [Path(p).expanduser() for p in os.environ.get("PISTAK_PLUGINS", "").split(os.pathsep) if p]
    return PLUGIN_DIRS + extra


def plugin_files() -> list[Path]:
    files: list[Path] = []
    for d in plugin_dirs():
        if d.is_dir():
            files += sorted(p for p in d.glob("*.py") if not p.name.startswith(("_", ".")))
    return files


def load_plugins() -> list[str]:
    """(Ri)carica tutti i plugin. Ritorna i nomi dei file caricati con successo."""
    ok: list[str] = []
    load_errors.clear()
    for path in plugin_files():
        name = path.name
        registry.remove_source(name)  # ricaricando, i parametri vecchi del file spariscono
        prev = registry.set_source(name)
        try:
            spec = importlib.util.spec_from_file_location(f"pistak_plugin_{path.stem}", path)
            assert spec and spec.loader
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            ok.append(name)
        except Exception:
            registry.remove_source(name)
            load_errors[name] = traceback.format_exc(limit=4)
        finally:
            registry.set_source(prev)
    return ok
