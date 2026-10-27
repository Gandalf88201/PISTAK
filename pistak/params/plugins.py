"""Caricamento dei plugin dell'utente.

Ogni file `.py` nelle cartelle dei plugin viene eseguito all'avvio; basta che usi `@parameter(...)` per aggiungere
parametri al registro predefinito. Le cartelle, in ordine:

* `plugins/` del progetto, se si lavora dai sorgenti;
* `plugins/` nella cartella di configurazione dell'utente (`pistak.core.config.user_config_dir()`);
* quelle elencate in $PISTAK_PLUGINS (separate da ':', ';' su Windows).

I file che iniziano con '_' o '.' sono ignorati (su volumi exFAT macOS crea file '._nome' che non sono codice).
"""
from __future__ import annotations

import importlib.util
import os
import sys
import traceback
from pathlib import Path

from pistak.core.config import ROOT, user_config_dir
from pistak.params.registry import REGISTRY, Registry
from pistak.params.structure import STRUCTURE_REGISTRY

#: Problemi incontrati caricando i plugin nel registro predefinito: nome file → messaggio.
load_errors: dict[str, str] = REGISTRY.load_errors


def plugin_dirs() -> list[Path]:
    dirs = [d for d in (ROOT / "plugins",) if d.is_dir()]
    dirs.append(user_config_dir() / "plugins")
    dirs += [Path(p).expanduser() for p in os.environ.get("PISTAK_PLUGINS", "").split(os.pathsep) if p]
    return dirs


def plugin_files() -> list[Path]:
    files: list[Path] = []
    for d in plugin_dirs():
        if d.is_dir():
            files += sorted(p for p in d.glob("*.py") if not p.name.startswith(("_", ".")))
    return files


def load_plugins(registry: Registry = REGISTRY) -> list[str]:
    """(Ri)carica tutti i plugin. Ritorna i nomi dei file caricati con successo.

    I plugin registrano con `@parameter`, cioè nel registro predefinito: per questo, durante il caricamento, la
    sorgente si imposta su quello. Un registro diverso riceve solo gli errori di caricamento."""
    ok: list[str] = []
    registry.load_errors.clear()
    for path in plugin_files():
        name = path.name
        REGISTRY.remove_source(name)  # ricaricando, i parametri vecchi del file spariscono
        STRUCTURE_REGISTRY.remove_source(name)
        prev = REGISTRY.set_source(name)
        prev_structure = STRUCTURE_REGISTRY.set_source(name)
        try:
            spec = importlib.util.spec_from_file_location(f"pistak_plugin_{path.stem}", path)
            assert spec and spec.loader
            module = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = module
            spec.loader.exec_module(module)
            ok.append(name)
        except Exception:
            REGISTRY.remove_source(name)
            STRUCTURE_REGISTRY.remove_source(name)
            registry.load_errors[name] = traceback.format_exc(limit=4)
        finally:
            REGISTRY.set_source(prev)
            STRUCTURE_REGISTRY.set_source(prev_structure)
    registry.plugins_loaded = True
    return ok
