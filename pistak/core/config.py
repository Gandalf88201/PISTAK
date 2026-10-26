"""Soglie dell'analisi, caricabili da TOML.

`config/default.toml` contiene i valori del sito; `config/user.toml` (se esiste) li sovrascrive.
Anche la finestra "Impostazioni" scrive su `user.toml`. Per usare soglie diverse da codice::

    analyse_file("a.cif", Limits(max_dihedral=15, max_interplanar=4.0))
"""
from __future__ import annotations

import tomllib
from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = ROOT / "config"


@dataclass(frozen=True)
class Limits:
    """Soglie della classificazione e della ricostruzione delle molecole. Unità: Å e gradi."""

    max_dihedral: float = field(default=20.0, metadata={"label": "Angolo massimo tra i piani (°)", "min": 0, "max": 90, "step": 1})
    max_interplanar: float = field(default=4.2, metadata={"label": "Distanza interplanare massima (Å)", "min": 1, "max": 10, "step": 0.1})
    contact_cutoff: float = field(default=4.0, metadata={"label": "Contatto massimo tra atomi pesanti (Å)", "min": 2, "max": 10, "step": 0.1})
    ring_rms: float = field(default=0.10, metadata={"label": "Planarità massima di un anello (Å rms)", "min": 0.01, "max": 1, "step": 0.01})
    molecule_rms: float = field(default=0.15, metadata={"label": "Planarità massima del nucleo (Å rms)", "min": 0.01, "max": 1, "step": 0.01})
    bond_tolerance: float = field(default=0.4, metadata={"label": "Tolleranza dei legami (Å)", "min": 0, "max": 1.5, "step": 0.05})
    same_site: float = field(default=0.3, metadata={"label": "Distanza minima tra atomi distinti (Å)", "min": 0.01, "max": 1, "step": 0.05})
    max_ring_size: int = field(default=8, metadata={"label": "Atomi massimi in un anello", "min": 3, "max": 20, "step": 1})
    neighbour_shells: int = field(default=2, metadata={"label": "Celle vicine esplorate (±)", "min": 1, "max": 4, "step": 1})
    ring_pair_max_distance: float = field(default=6.0, metadata={"label": "Cg···Cg massimo per le coppie di anelli (Å)", "min": 3, "max": 12, "step": 0.5})
    max_cell_atoms: int = field(default=6000, metadata={"label": "Atomi massimi nella cella", "min": 100, "max": 50000, "step": 100})

    def with_values(self, **kw) -> Limits:
        return replace(self, **kw)

    def to_dict(self) -> dict:
        return asdict(self)


def _read(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("rb") as f:
        return tomllib.load(f)


def load_limits(config_dir: Path = CONFIG_DIR) -> Limits:
    """Soglie del sito, poi quelle di default.toml, poi quelle dell'utente."""
    merged: dict = {}
    for name in ("default.toml", "user.toml"):
        merged.update(_read(config_dir / name).get("limits", {}))
    known = {f.name: f.type for f in fields(Limits)}
    clean = {k: (int(v) if known[k] == "int" else float(v)) for k, v in merged.items() if k in known}
    return Limits(**clean)


def load_section(name: str, config_dir: Path = CONFIG_DIR) -> dict:
    """Una sezione qualunque di default.toml/user.toml (per le preferenze della GUI)."""
    merged: dict = {}
    for fn in ("default.toml", "user.toml"):
        merged.update(_read(config_dir / fn).get(name, {}))
    return merged


def _toml_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    return '"' + str(v).replace("\\", "\\\\").replace('"', '\\"') + '"'


def save_user_config(limits: Limits | None = None, sections: dict[str, dict] | None = None,
                     config_dir: Path = CONFIG_DIR) -> Path:
    """Scrive config/user.toml, conservando le altre sezioni già presenti."""
    path = config_dir / "user.toml"
    data = _read(path)
    if limits is not None:
        data["limits"] = limits.to_dict()
    for k, v in (sections or {}).items():
        data.setdefault(k, {}).update(v)
    config_dir.mkdir(parents=True, exist_ok=True)
    lines = ["# Impostazioni personali / Personal settings: sovrascrivono / override config/default.toml. Modificabile a mano / editable.\n"]
    for section, values in data.items():
        lines.append(f"\n[{section}]\n")
        lines.extend(f"{k} = {_toml_value(v)}\n" for k, v in values.items())
    path.write_text("".join(lines), encoding="utf-8")
    return path
