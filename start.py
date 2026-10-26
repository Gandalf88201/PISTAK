#!/usr/bin/env python3
"""Avvio di PISTAK: prepara l'ambiente virtuale e lancia l'applicazione.

    ./start.py                      apre la finestra grafica
    ./start.py esempio.cif          apre la GUI già con quel file
    ./start.py analyze a.cif --csv out.csv     analisi da riga di comando (nessuna GUI)
    ./start.py --reinstall          reinstalla i pacchetti nell'ambiente virtuale
    ./start.py --setup-only         prepara l'ambiente e basta
    ./start.py --dev                installa anche gli strumenti di sviluppo (pytest)
    ./start.py --python /path/py    usa un interprete preciso per creare .venv
    ./start.py --lang en            interfaccia in inglese / interface in English (anche PISTAK_LANG=en)

Usa solo la libreria standard: può essere lanciato con qualunque Python 3.8+.
I pacchetti vengono scaricati in `.venv/` dentro questa cartella, mai nel Python di sistema.
Per aggiungere pacchetti propri scrivili in `requirements-user.txt` (uno per riga): al
prossimo avvio vengono installati da soli.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV = ROOT / ".venv"
STAMP = VENV / ".pistak-stamp"
MIN_PY, MAX_PY = (3, 11), (3, 13)  # intervallo per cui esistono le ruote di PySide6 e VTK
REQUIREMENTS = ["requirements.txt", "requirements-user.txt"]
DEV_REQUIREMENTS = ["requirements-dev.txt"]
IS_WINDOWS = os.name == "nt"


MESSAGES = {  # id: (italiano, inglese)
    "bad_interpreter": ("interprete non utilizzabile: {p}", "unusable interpreter: {p}"),
    "out_of_range": ("attenzione: Python {v} è fuori dall'intervallo provato ({lo}–{hi}); procedo comunque.",
                     "warning: Python {v} is outside the tested range ({lo}–{hi}); continuing anyway."),
    "no_python": ("non trovo un Python {lo}–{hi} (questo è {v}). Installane uno oppure usa --python /percorso/python3.12",
                  "cannot find a Python {lo}–{hi} (this one is {v}). Install one or use --python /path/to/python3.12"),
    "recreate": ("ricreo l'ambiente virtuale…", "recreating the virtual environment…"),
    "create": ("creo l'ambiente virtuale in {d}/ con Python {v}", "creating the virtual environment in {d}/ with Python {v}"),
    "venv_failed": ("impossibile creare l'ambiente virtuale", "could not create the virtual environment"),
    "installing": ("installo i pacchetti (la prima volta ci vuole qualche minuto)…", "installing packages (the first time takes a few minutes)…"),
    "install_failed": ("installazione fallita per {f}: correggi il file e rilancia ./start.py",
                       "installation failed for {f}: fix the file and run ./start.py again"),
    "purged": ("rimossi {n} file di metadati '._*' creati dal volume", "removed {n} '._*' metadata files created by the volume"),
    "ready": ("pacchetti pronti.", "packages ready."),
    "env_ready": ("ambiente pronto: {p}", "environment ready: {p}"),
}


def _config_files() -> list[Path]:
    """Dove cercare la lingua, in ordine: user.toml dell'utente (stesse regole di pistak/core/config.py), quello
    vecchio in config/ del progetto, poi i valori predefiniti del pacchetto."""
    env = os.environ.get("PISTAK_CONFIG_DIR")
    if env:
        user = Path(env).expanduser()
    elif sys.platform == "darwin":
        user = Path.home() / "Library" / "Application Support" / "PISTAK"
    elif IS_WINDOWS:
        user = Path(os.environ.get("APPDATA") or Path.home() / "AppData" / "Roaming") / "PISTAK"
    else:
        user = Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config") / "pistak"
    return [user / "user.toml", ROOT / "config" / "user.toml", ROOT / "pistak" / "default.toml"]


def _config_language() -> str:
    """La lingua scelta nei file di configurazione (lettura semplice: start.py non dipende da tomllib)."""
    for f in _config_files():
        if f.exists():
            m = re.search(r'^\s*language\s*=\s*"(\w+)"', f.read_text(encoding="utf-8"), re.M)
            if m:
                return m.group(1)
    return "it"


def resolve_language(explicit: str | None = None) -> str:
    value = (explicit or os.environ.get("PISTAK_LANG") or _config_language()).lower()
    if value == "auto":
        import locale

        value = "it" if (locale.getlocale()[0] or os.environ.get("LANG", "")).lower().startswith("it") else "en"
    return value if value in ("it", "en") else "it"


LANG = resolve_language()


def t(key: str, **kw) -> str:
    it, en = MESSAGES[key]
    return (en if LANG == "en" else it).format(**kw)


def say(msg: str) -> None:
    print(f"[pistak] {msg}", flush=True)


def probe(cmd: list[str]) -> tuple[int, int] | None:
    """Versione (major, minor) dell'interprete `cmd`, o None se non parte."""
    try:
        out = subprocess.run(
            [*cmd, "-c", "import sys; print(sys.version_info[0], sys.version_info[1])"],
            capture_output=True, text=True, timeout=30,
        )
        if out.returncode != 0:
            return None
        major, minor = out.stdout.split()
        return int(major), int(minor)
    except (OSError, ValueError, subprocess.SubprocessError):
        return None


def supported(v: tuple[int, int] | None) -> bool:
    return v is not None and MIN_PY <= v <= MAX_PY


def find_python(explicit: str | None) -> list[str]:
    """Il primo interprete adatto: --python, $PISTAK_PYTHON, poi i python3.x più comuni."""
    wanted = explicit or os.environ.get("PISTAK_PYTHON")
    if wanted:
        v = probe([wanted])
        if v is None:
            sys.exit("[pistak] " + t("bad_interpreter", p=wanted))
        if not supported(v):
            say(t("out_of_range", v=f"{v[0]}.{v[1]}", lo=f"{MIN_PY[0]}.{MIN_PY[1]}", hi=f"{MAX_PY[0]}.{MAX_PY[1]}"))
        return [wanted]
    candidates: list[list[str]] = []
    for minor in range(MAX_PY[1], MIN_PY[1] - 1, -1):
        name = f"python3.{minor}"
        found = shutil.which(name)
        if found:
            candidates.append([found])
        for d in ("/opt/homebrew/bin", "/usr/local/bin", "/opt/local/bin", "/usr/bin"):
            p = Path(d) / name
            if p.exists():
                candidates.append([str(p)])
        if IS_WINDOWS:
            candidates.append(["py", f"-3.{minor}"])
    candidates.append([sys.executable])
    for c in candidates:
        if supported(probe(c)):
            return c
    v = probe([sys.executable])
    sys.exit("[pistak] " + t("no_python", lo=f"{MIN_PY[0]}.{MIN_PY[1]}", hi=f"{MAX_PY[0]}.{MAX_PY[1]}", v=f"{v[0]}.{v[1]}"))


def venv_python() -> Path:
    return VENV / ("Scripts/python.exe" if IS_WINDOWS else "bin/python")


def ensure_venv(py: list[str], recreate: bool) -> Path:
    vpy = venv_python()
    if vpy.exists() and not recreate and supported(probe([str(vpy)])):
        return vpy
    if VENV.exists():
        say(t("recreate"))
        shutil.rmtree(VENV, ignore_errors=True)
    v = probe(py)
    say(t("create", d=VENV.name, v=f"{v[0]}.{v[1]}"))
    for extra in ([], ["--copies"]):  # su volumi exFAT/FAT i collegamenti simbolici possono mancare
        r = subprocess.run([*py, "-m", "venv", *extra, str(VENV)])
        if r.returncode == 0 and vpy.exists():
            return vpy
        shutil.rmtree(VENV, ignore_errors=True)
    sys.exit("[pistak] " + t("venv_failed"))


def requirement_files(dev: bool) -> list[Path]:
    names = REQUIREMENTS + (DEV_REQUIREMENTS if dev else [])
    return [ROOT / n for n in names if (ROOT / n).exists()]


def fingerprint(vpy: Path, files: list[Path]) -> str:
    h = hashlib.sha256()
    h.update(str(probe([str(vpy)])).encode())
    for f in files:
        h.update(f.name.encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def purge_appledouble(root: Path) -> int:
    """Su volumi exFAT/FAT macOS crea file '._nome' accanto a ogni file con attributi estesi. Dentro
    site-packages rompono librerie che leggono intere cartelle (matplotlib: 'cannot decode configuration
    file ... .__classic_test_patch.mplstyle'). Si cancellano dopo ogni installazione."""
    n = 0
    for dirpath, _dirs, names in os.walk(root):
        for name in names:
            if name.startswith("._"):
                try:
                    os.unlink(os.path.join(dirpath, name))
                    n += 1
                except OSError:
                    pass
    return n


def install(vpy: Path, files: list[Path], force: bool) -> None:
    fp = fingerprint(vpy, files)
    if not force and STAMP.exists() and STAMP.read_text().strip() == fp:
        return
    say(t("installing"))
    pip = [str(vpy), "-m", "pip", "--disable-pip-version-check"]
    subprocess.run([*pip, "install", "--quiet", "--upgrade", "pip"], check=False)
    for f in files:
        say(f"  {f.name}")
        r = subprocess.run([*pip, "install", "-r", str(f)])
        if r.returncode != 0:
            sys.exit("[pistak] " + t("install_failed", f=f.name))
    purged = purge_appledouble(VENV)
    if purged:
        say(t("purged", n=purged))
    STAMP.write_text(fp)
    say(t("ready"))


def main() -> None:
    ap = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    ap.add_argument("--reinstall", action="store_true")
    ap.add_argument("--setup-only", action="store_true")
    ap.add_argument("--dev", action="store_true")
    ap.add_argument("--python")
    ap.add_argument("--lang", choices=["it", "en"])
    ap.add_argument("--start-help", action="store_true")
    own, rest = ap.parse_known_args()
    global LANG
    LANG = resolve_language(own.lang)
    if own.start_help:
        print(__doc__)
        return

    # ambiente già pronto: niente ricerca dell'interprete a ogni avvio
    vpy = venv_python()
    if own.reinstall or own.python or not vpy.exists():
        vpy = ensure_venv(find_python(own.python), recreate=bool(own.python))
    install(vpy, requirement_files(own.dev), force=own.reinstall)
    if own.setup_only:
        say(t("env_ready", p=vpy))
        return

    os.chdir(ROOT)
    cmd = [str(vpy), "-m", "pistak", *(["--lang", own.lang] if own.lang else []), *rest]
    if IS_WINDOWS:
        sys.exit(subprocess.call(cmd))
    os.execv(cmd[0], cmd)


if __name__ == "__main__":
    main()
