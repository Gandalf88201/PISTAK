"""Lingue dell'interfaccia: italiano (la lingua dei testi nel codice) e inglese.

Come funziona: nel codice ogni testo mostrato all'utente è scritto in italiano dentro `tr("...")`. Con la lingua
italiana `tr` restituisce il testo così com'è; con l'inglese lo cerca nel dizionario `pistak/i18n_en.py` e, se manca,
lo lascia in italiano (mai un errore). I testi con numeri usano i segnaposto di `str.format`:

    tf("{n} coppie in {f} file", n=3, f=2)       →  "3 pairs in 2 files"

Scegliere la lingua, in ordine di priorità:
    1. variabile d'ambiente PISTAK_LANG (it | en)       2. `./start.py --lang en`
    3. user.toml nella cartella dell'utente (`pistak config`), sezione [gui], language = "it" | "en" | "auto"
       (auto = la lingua del sistema)
    4. italiano

Si cambia anche dal menu (Lingua / Language): la finestra si ricostruisce da sola conservando la sessione.

Per aggiungere una lingua nuova: inserirla in LANGS, creare `pistak/i18n_<codice>.py` con un dizionario `TEXT`
e aggiungerla a `_load`. Per tradurre i parametri di un plugin: `add_translations("en", {"Il mio parametro": "My parameter"})`.
"""
from __future__ import annotations

import locale
import os

#: Lingue disponibili: codice → nome (nella lingua stessa, come si usa nei selettori).
LANGS = {"it": "Italiano", "en": "English"}
DEFAULT = "it"

_lang: str | None = None  # None finché non viene scelta (la prima tr() la sceglie)
_catalogs: dict[str, dict[str, str]] = {}


def _load(lang: str) -> dict[str, str]:
    if lang not in _catalogs:
        if lang == "en":
            from pistak.i18n_en import TEXT
        else:
            TEXT = {}
        _catalogs[lang] = dict(TEXT)
    return _catalogs[lang]


def add_translations(lang: str, texts: dict[str, str]) -> None:
    """Aggiunge (o sostituisce) traduzioni: utile ai plugin. `texts` va dal testo italiano a quello nella lingua `lang`."""
    _load(lang).update(texts)


def system_language() -> str:
    """La lingua del sistema se è tra quelle disponibili, altrimenti l'inglese."""
    try:
        code = (locale.getlocale()[0] or os.environ.get("LANG", "") or "").lower()
    except (ValueError, TypeError):
        code = ""
    for lang in LANGS:
        if code.startswith(lang):
            return lang
    return "en"


def resolve(value: str | None) -> str:
    """Da un valore di configurazione ("it", "en", "auto", None) a un codice lingua valido."""
    value = (value or "").strip().lower()
    if value == "auto":
        return system_language()
    return value if value in LANGS else DEFAULT


def configured_language() -> str:
    """La lingua da usare all'avvio: ambiente, poi configurazione, poi italiano."""
    env = os.environ.get("PISTAK_LANG")
    if env:
        return resolve(env)
    try:
        from pistak.core.config import load_section

        return resolve(str(load_section("gui").get("language", DEFAULT)))
    except Exception:  # noqa: BLE001  (configurazione illeggibile: si parte in italiano)
        return DEFAULT


def set_language(lang: str) -> str:
    global _lang
    _lang = resolve(lang)
    _load(_lang)
    return _lang


def get_language() -> str:
    global _lang
    if _lang is None:
        set_language(configured_language())
    return _lang  # type: ignore[return-value]


def tr(text: str) -> str:
    """Il testo nella lingua corrente. Se non c'è una traduzione resta quello originale."""
    lang = get_language()
    return text if lang == DEFAULT else _load(lang).get(text, text)


def tf(text: str, **values) -> str:
    """`tr` seguito da `str.format`: per i testi con numeri o nomi."""
    return tr(text).format(**values)


def save_language(lang: str) -> None:
    """Salva la scelta in user.toml (cartella dell'utente), così vale anche al prossimo avvio."""
    from pistak.core.config import save_user_config

    save_user_config(sections={"gui": {"language": resolve(lang)}})
