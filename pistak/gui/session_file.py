"""Formato del file di sessione (.pistak): JSON leggibile, con il CIF incorporato.

    {"version": 1, "app": "PISTAK", "saved": "<data>",
     "source":  {"name": ..., "path": ..., "text": <CIF, se non troppo grande>},
     "limits":  {...soglie...},
     "selection": {"A": 0, ...}, "subset": {"A": [..] | null, ...},
     "views":   [{"kind": "3D", "object": "3D-1", "title": ..., "state": {...}}, ...],
     "window":  {"geometry": <base64>, "dock_state": <base64>},
     "collection": {"paths": [...], "reference": 0}}

Si può aprire e correggere a mano. Le chiavi sconosciute vengono ignorate, così un file scritto da una versione
più nuova si apre comunque (con meno dettagli).
"""
from __future__ import annotations

import base64
import json
from datetime import datetime
from pathlib import Path

from pistak.i18n import tf, tr

VERSION = 1
SUFFIX = ".pistak"
#: Oltre questa dimensione il CIF non viene incorporato (resta solo il percorso).
MAX_EMBED = 2 * 1024 * 1024


def b64(data: bytes) -> str:
    return base64.b64encode(data).decode("ascii")


def unb64(text: str | None) -> bytes:
    return base64.b64decode(text) if text else b""


class SessionFileError(Exception):
    pass


def write(path: str | Path, state: dict) -> Path:
    path = Path(path)
    state = {"version": VERSION, "app": "PISTAK", "saved": datetime.now().isoformat(timespec="seconds"), **state}
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(state, indent=1, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)  # scrittura atomica: un errore a metà non rovina la sessione precedente
    return path


def read(path: str | Path) -> dict:
    try:
        st = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as e:
        raise SessionFileError(tf("file di sessione illeggibile: {e}", e=e)) from e
    if not isinstance(st, dict) or st.get("app") != "PISTAK":
        raise SessionFileError(tr("non è un file di sessione di PISTAK"))
    if int(st.get("version", 0)) > VERSION:
        raise SessionFileError(tf("sessione scritta da una versione più nuova (formato {v}, questa legge fino a {max})", v=st["version"], max=VERSION))
    return st


def embed_source(name: str, path: str | None, text: str) -> dict:
    src: dict = {"name": name, "path": path}
    if len(text.encode("utf-8")) <= MAX_EMBED:
        src["text"] = text
    return src
