"""Confronto tra strutture: per ciascuna il motivo di stack migliore, le differenze rispetto a una struttura di
riferimento e la distanza geometrica tra i motivi di stack più simili."""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from pistak.i18n import tr
from pistak.io.batch import PairRow, StructureEntry

#: Parametri della tabella di confronto (chiavi del registro).
COMPARE_KEYS = ["interplanar", "long_slip", "lat_slip", "dihedral", "centroid_distance", "contact", "overlap_area", "beta_a"]
#: Parametri che descrivono la geometria del motivo di stack (per la distanza tra due motivi).
MOTIF_KEYS = ("interplanar", "long_slip", "lat_slip")


@dataclass
class CompareRow:
    name: str
    error: str | None
    formulas: str
    space_group: str
    n_stack: int
    n_columns: int
    best: dict[str, float] = field(default_factory=dict)
    #: differenza (questa − riferimento) per ogni parametro del motivo migliore
    delta: dict[str, float] = field(default_factory=dict)
    #: distanza (Å) tra i due motivi di stack più simili tra questa struttura e il riferimento
    motif_distance: float = float("nan")


def _motif(p: PairRow) -> np.ndarray:
    return np.array([p.values.get(k, np.nan) for k in MOTIF_KEYS])


def motif_distance(a: StructureEntry, b: StructureEntry) -> float:
    """Distanza euclidea (Å) in (interplanare, slittamento L, slittamento W) tra le due coppie di stack più simili,
    una per struttura. NaN se una delle due non ha stack."""
    ma = [_motif(p) for p in a.pairs if p.kind == "stack"]
    mb = [_motif(p) for p in b.pairs if p.kind == "stack"]
    if not ma or not mb:
        return float("nan")
    return float(min(np.linalg.norm(x - y) for x in ma for y in mb))


def compare(entries: list[StructureEntry], reference: int = 0, keys: list[str] | None = None) -> list[CompareRow]:
    keys = keys or COMPARE_KEYS
    ref = entries[reference] if entries else None
    ref_best = ref.best if ref else {}
    rows = []
    for e in entries:
        best = e.best
        delta = {k: best[k] - ref_best[k] for k in keys if k in best and k in ref_best and math.isfinite(best[k] - ref_best[k])}
        rows.append(CompareRow(e.name, e.error, e.formulas, e.space_group, e.counts.get("stack", 0), len(e.columns),
                               {k: best[k] for k in keys if k in best}, delta,
                               0.0 if e is ref else (motif_distance(e, ref) if ref else float("nan"))))
    return rows


def write_csv(rows: list[CompareRow], path: str | Path, keys: list[str] | None = None) -> Path:
    keys = keys or COMPARE_KEYS
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["file", tr("errore"), tr("gruppo"), tr("molecole"), "n_stack", tr("colonne"), *keys, *[f"delta_{k}" for k in keys],
                    tr("distanza_motivo")])
        fmt = lambda v: "" if v is None or not math.isfinite(v) else f"{v:.4g}"  # noqa: E731
        for r in rows:
            w.writerow([r.name, r.error or "", r.space_group, r.formulas, r.n_stack, r.n_columns,
                        *[fmt(r.best.get(k)) for k in keys], *[fmt(r.delta.get(k)) for k in keys], fmt(r.motif_distance)])
    return path
