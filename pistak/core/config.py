"""Soglie dell'analisi e preferenze, caricabili da TOML.

Due file, il secondo sovrascrive il primo:

* `pistak/default.toml`: i valori di partenza (quelli del sito), distribuiti con il pacchetto;
* `user.toml` nella cartella di configurazione dell'utente (`user_config_dir()`), scritto anche dalla finestra
  «Impostazioni». La cartella si sceglie con la variabile PISTAK_CONFIG_DIR; altrimenti è quella del sistema
  (macOS: ~/Library/Application Support/PISTAK, Windows: %APPDATA%\\PISTAK, Linux: ~/.config/pistak).

Lì stanno anche l'ultima sessione (`last.pistak`) e una cartella `plugins/`. Le versioni precedenti tenevano
`user.toml` e `last.pistak` in `config/` dentro il progetto: alla prima lettura vengono copiati nella nuova cartella.

Per usare soglie diverse da codice::

    analyse_file("a.cif", Limits(max_dihedral=15, max_interplanar=4.0))
"""
from __future__ import annotations

import os
import shutil
import sys
import tomllib
from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path

#: Cartella del progetto quando si lavora dai sorgenti (esempi, plugins/ del repository).
ROOT = Path(__file__).resolve().parents[2]
DEFAULTS_FILE = Path(__file__).resolve().parents[1] / "default.toml"
USER_FILE = "user.toml"
#: Dove le versioni precedenti tenevano le impostazioni personali.
LEGACY_CONFIG_DIR = ROOT / "config"
_LEGACY_FILES = ("user.toml", "last.pistak")


def _system_config_dir() -> Path:
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "PISTAK"
    if os.name == "nt":
        return Path(os.environ.get("APPDATA") or Path.home() / "AppData" / "Roaming") / "PISTAK"
    return Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config") / "pistak"


def _migrate_legacy(target: Path) -> None:
    """Copia user.toml e last.pistak dalla vecchia cartella config/ del progetto, se nella nuova mancano."""
    for name in _LEGACY_FILES:
        old, new = LEGACY_CONFIG_DIR / name, target / name
        if old.is_file() and not new.exists():
            try:
                target.mkdir(parents=True, exist_ok=True)
                shutil.copy2(old, new)
            except OSError:
                pass  # cartella non scrivibile: si continua con i valori predefiniti


def user_config_dir() -> Path:
    """Cartella delle impostazioni personali (può non esistere ancora: la crea chi scrive)."""
    env = os.environ.get("PISTAK_CONFIG_DIR")
    if env:
        return Path(env).expanduser()
    d = _system_config_dir()
    _migrate_legacy(d)
    return d


def user_config_file() -> Path:
    return user_config_dir() / USER_FILE


@dataclass(frozen=True)
class Limits:
    """Soglie della classificazione e della ricostruzione delle molecole. Unità: Å e gradi."""

    max_dihedral: float = field(default=20.0, metadata={"label": "Angolo massimo tra i piani (°)", "min": 0, "max": 90, "step": 1})
    max_interplanar: float = field(default=4.2, metadata={"label": "Distanza interplanare massima (Å)", "min": 1, "max": 10, "step": 0.1})
    #: 4,5 Å: con piani paralleli il contatto più corto non scende sotto la distanza interplanare, e un atomo di B dentro
    #: il contorno di A dista al più ~1,4 Å (raggio di un anello) in proiezione dall'atomo di A più vicino; √(4,2² +
    #: 1,4²) ≈ 4,43 Å, quindi ogni stack entro `max_interplanar` = 4,2 Å viene misurato. Il sito usava 4,0 Å.
    contact_cutoff: float = field(default=4.5, metadata={"label": "Contatto massimo tra atomi pesanti (Å)", "min": 2, "max": 10, "step": 0.1})
    ring_rms: float = field(default=0.10, metadata={"label": "Planarità massima di un anello (Å rms)", "min": 0.01, "max": 1, "step": 0.01})
    molecule_rms: float = field(default=0.15, metadata={"label": "Planarità massima del nucleo (Å rms)", "min": 0.01, "max": 1, "step": 0.01})
    bond_tolerance: float = field(default=0.4, metadata={"label": "Tolleranza dei legami (Å)", "min": 0, "max": 1.5, "step": 0.05})
    same_site: float = field(default=0.3, metadata={"label": "Distanza minima tra atomi distinti (Å)", "min": 0.01, "max": 1, "step": 0.05})
    max_ring_size: int = field(default=8, metadata={"label": "Atomi massimi in un anello", "min": 3, "max": 20, "step": 1})
    neighbour_shells: int = field(default=2, metadata={"label": "Celle vicine esplorate (± minimo)", "min": 1, "max": 4, "step": 1})
    ring_pair_max_distance: float = field(default=6.0, metadata={"label": "Cg···Cg massimo per le coppie di anelli (Å)", "min": 3, "max": 12, "step": 0.5})
    max_cell_atoms: int = field(default=6000, metadata={"label": "Atomi massimi nella cella", "min": 100, "max": 50000, "step": 100})
    # ---- criterio di stack (oltre ad angolo e distanza interplanare) --------------------------------------
    #: Predefinito: area comune dei poligoni degli anelli / area del frammento minore, sopra 0,10 (revisione del paper,
    #: voce 50). Il sito usava «un atomo dentro l'involucro convesso dell'altro nucleo» (`Limits.site()`), che accetta
    #: coppie con meno dell'1 % di area in comune.
    overlap_measure: str = field(default="area", metadata={"label": "Sovrapposizione misurata come", "choices": {
        "atoms": "frazione di atomi dentro il contorno dell'altro nucleo",
        "area": "area comune / area del nucleo minore"}})
    min_overlap: float = field(default=0.10, metadata={"label": "Sovrapposizione minima per uno stack (sopra questo valore)", "min": 0, "max": 1, "step": 0.05})
    max_slip: float = field(default=0.0, metadata={"label": "Slittamento massimo per uno stack (Å, 0 = nessun limite)", "min": 0, "max": 20, "step": 0.1})
    # ---- coppie e disordine ------------------------------------------------------------------------------
    pair_mode: str = field(default="canonical", metadata={"label": "Coppie", "choices": {
        "canonical": "una per contatto, equivalenti per simmetria del gruppo",
        "site": "come il sito (a→b e b→a, equivalenza per valori arrotondati)"}})
    disorder: str = field(default="major", metadata={"label": "Atomi disordinati", "choices": {
        "major": "componente con occupanza maggiore",
        "first": "gruppo 1 (come il sito)",
        "all": "tutti gli atomi"}, "numeric": True})
    # ---- frammenti π e assi ------------------------------------------------------------------------------
    fragments: str = field(default="all", metadata={"label": "Frammenti π per molecola", "choices": {
        "all": "tutti i sistemi di anelli planari",
        "largest": "solo il più grande (come il sito)"}})
    min_fragment_atoms: int = field(default=5, metadata={"label": "Atomi minimi in un frammento π secondario", "min": 3, "max": 60, "step": 1})
    fragment_merge: str = field(default="fused", metadata={"label": "Frammenti π", "choices": {
        "fused": "sistemi di anelli condensati (che condividono un legame)",
        "conjugated": "anche sistemi complanari collegati da legami o atomi sp/sp² (biarili, porfirine)"}})
    #: Strutture inorganiche (reti, carburi, carbonati, cianuri metallici): per impostazione nessun frammento π, perché i
    #: loro anelli planari non sono sistemi π e la loro unità ripetuta dipende dalla cella (revisione 2). «organic»: un
    #: legame C–C o C–H, oppure un anello planare di soli non metalli con B, C o N (melammina, acido cianurico,
    #: borazine: revisione 3, N1; `analysis.in_pi_scope`).
    pi_scope: str = field(default="organic", metadata={"label": "Strutture in cui cercare frammenti π", "choices": {
        "organic": "con legami C–C o C–H o un anello planare di non metalli con B, C o N",
        "all": "tutte, anche le inorganiche"}})
    long_axis_atoms: str = field(default="", metadata={"label": "Asse lungo da due atomi (etichette, es. C1,C9; vuoto = automatico)", "free": True})
    #: Sopra questa isotropia ((λ2/λ3)^½ nel piano, la radice del rapporto degli autovalori) l'asse lungo automatico non è definito (benzene, coronene, porfirine:
    #: un asse Cn con n ≥ 3 rende isotropa la dispersione nel piano) e gli slittamenti L/W valgono NaN.
    max_axis_isotropy: float = field(default=0.9, metadata={"label": "Isotropia massima per l'asse lungo automatico", "min": 0.5, "max": 1, "step": 0.01})
    #: Sopra questo angolo tra gli assi lunghi proiettati la bisettrice salta di 90° per variazioni infinitesime: L/W NaN.
    max_axis_twist: float = field(default=75.0, metadata={"label": "Angolo massimo tra gli assi lunghi per L/W (°)", "min": 0, "max": 90, "step": 1})
    # ---- controlli rapidi --------------------------------------------------------------------------------
    contact_check: str = field(default="fixed", metadata={"label": "Controllo dei contatti corti", "choices": {
        "fixed": "soglie fisse (C···C sotto 3,0 Å, H···H sotto 1,8 Å)",
        "vdw": "somma dei raggi di van der Waals meno un margine, per tutti gli elementi"}})
    contact_check_margin: float = field(default=0.5, metadata={"label": "Margine sotto la somma dei raggi di van der Waals (Å)", "min": 0, "max": 1.5, "step": 0.05})

    def with_values(self, **kw) -> Limits:
        return replace(self, **kw)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def site(cls, **kw) -> Limits:
        """Le soglie e le convenzioni dell'analizzatore del sito (per i confronti di parità)."""
        return cls(**{"pair_mode": "site", "disorder": "first", "fragments": "largest", "contact_cutoff": 4.0,
                      "overlap_measure": "atoms", "min_overlap": 0.0, "pi_scope": "all", **kw})

    @classmethod
    def from_dict(cls, data: dict, base: Limits | None = None) -> Limits:
        """Soglie da un dizionario (TOML, sessione): le chiavi sconosciute e i valori non validi si ignorano."""
        out = asdict(base or cls())
        for f in fields(cls):
            if f.name not in data:
                continue
            try:
                out[f.name] = _coerce(f, data[f.name])
            except (TypeError, ValueError):
                pass
        return cls(**out)


def _coerce(f, v):
    """Un valore di configurazione convertito nel tipo del campo; ValueError se non è accettabile."""
    if f.type == "str" and f.metadata.get("free"):
        if isinstance(v, (str, int, float)) and not isinstance(v, bool):
            return str(v).strip()
        raise TypeError(f.name)
    if f.type == "str":
        v = str(v).strip().lower()
        if v in f.metadata.get("choices", {}) or (f.metadata.get("numeric") and v.isdigit()):
            return v
        raise ValueError(f"{f.name}: {v!r}")
    if isinstance(v, bool) or not isinstance(v, (int, float, str)):
        raise TypeError(f.name)
    return int(v) if f.type == "int" else float(v)


def _read(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("rb") as f:
        return tomllib.load(f)


#: Valori che erano predefiniti in versioni precedenti: (vecchio valore, data del cambio). Un user.toml salvato prima
#: li ha congelati (la finestra «Impostazioni» scriveva tutte le soglie); il controllo rapido `oldDefault` lo segnala.
SUPERSEDED_DEFAULTS: dict[str, tuple[object, str]] = {
    "contact_cutoff": (4.0, "2026-11-03"),  # 4,0 → 4,5 Å (registro, voce 6); il modo del sito resta a 4,0
}


def superseded_values(limits: Limits) -> list[tuple[str, object, object, str]]:
    """Soglie uguali a un vecchio valore predefinito: (chiave, valore, predefinito attuale, data del cambio). Nel modo
    del sito `contact_cutoff` = 4,0 è la sua convenzione, non un valore rimasto indietro."""
    current = default_limits().to_dict()
    out = []
    for key, (old, when) in SUPERSEDED_DEFAULTS.items():
        value = getattr(limits, key)
        if value == old != current[key] and not (key == "contact_cutoff" and limits.pair_mode == "site"):
            out.append((key, value, current[key], when))
    return out


def default_limits() -> Limits:
    """Le soglie predefinite: quelle del codice con sopra default.toml, senza il file dell'utente."""
    return Limits.from_dict(dict(_read(DEFAULTS_FILE).get("limits", {})))


def _merged(section: str, config_dir: Path | None) -> dict:
    merged: dict = dict(_read(DEFAULTS_FILE).get(section, {}))
    merged.update(_read((config_dir or user_config_dir()) / USER_FILE).get(section, {}))
    return merged


def load_limits(config_dir: Path | None = None) -> Limits:
    """Soglie del sito, poi quelle di default.toml, poi quelle dell'utente."""
    return Limits.from_dict(_merged("limits", config_dir))


def load_section(name: str, config_dir: Path | None = None) -> dict:
    """Una sezione qualunque di default.toml/user.toml (per le preferenze della GUI)."""
    return _merged(name, config_dir)


def _toml_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    return '"' + str(v).replace("\\", "\\\\").replace('"', '\\"') + '"'


def save_user_config(limits: Limits | None = None, sections: dict[str, dict] | None = None,
                     config_dir: Path | None = None) -> Path:
    """Scrive user.toml nella cartella dell'utente, conservando le altre sezioni già presenti. Delle soglie si
    scrivono solo quelle diverse dai valori predefiniti."""
    config_dir = config_dir or user_config_dir()
    path = config_dir / USER_FILE
    data = _read(path)
    if limits is not None:  # solo ciò che differisce dai predefiniti: un predefinito cambiato in seguito vale anche qui
        base = default_limits().to_dict()
        data["limits"] = {k: v for k, v in limits.to_dict().items() if v != base.get(k)}
    for k, v in (sections or {}).items():
        data.setdefault(k, {}).update(v)
    config_dir.mkdir(parents=True, exist_ok=True)
    lines = ["# Impostazioni personali / Personal settings: sovrascrivono / override pistak/default.toml. Modificabile a mano / editable.\n"]
    for section, values in data.items():
        lines.append(f"\n[{section}]\n")
        lines.extend(f"{k} = {_toml_value(v)}\n" for k, v in values.items())
    path.write_text("".join(lines), encoding="utf-8")
    return path
