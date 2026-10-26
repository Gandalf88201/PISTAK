"""Statistiche su molte coppie, anche di file diversi: tabella unica, descrizione, correlazioni, istogrammi.

Solo numpy: i risultati sono array e dizionari, disegnabili da qualunque libreria (la GUI usa pyqtgraph).
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pistak.io.batch import StructureEntry

KINDS = ("stack", "parallel", "tilted")


@dataclass
class PairTable:
    """Tutte le coppie di un insieme di strutture in una matrice file × parametro."""

    files: list[str]  # nome del file di ogni riga
    paths: list[str]
    index: np.ndarray  # indice della coppia nel suo file
    kinds: np.ndarray  # str
    counts: np.ndarray
    keys: list[str]
    data: np.ndarray  # N × K

    def __len__(self) -> int:
        return len(self.files)

    @classmethod
    def from_entries(cls, entries: list[StructureEntry]) -> PairTable:
        keys: list[str] = []
        for e in entries:
            for p in e.pairs:
                for k in p.values:
                    if k not in keys:
                        keys.append(k)
        rows = [(e, p) for e in entries for p in e.pairs]
        data = np.full((len(rows), len(keys)), np.nan)
        for i, (_e, p) in enumerate(rows):
            for j, k in enumerate(keys):
                data[i, j] = p.values.get(k, np.nan)
        return cls([e.name for e, _ in rows], [e.path for e, _ in rows], np.array([p.index for _, p in rows], dtype=int),
                   np.array([p.kind for _, p in rows], dtype=str), np.array([p.count for _, p in rows], dtype=int), keys, data)

    def mask(self, kind: str | None = None, files: set[str] | None = None) -> np.ndarray:
        m = np.ones(len(self), bool)
        if kind:
            m &= self.kinds == kind
        if files is not None:
            m &= np.isin(np.array(self.files, dtype=str), list(files))
        return m

    def column(self, key: str) -> np.ndarray:
        return self.data[:, self.keys.index(key)]

    def weights(self, weighted: bool) -> np.ndarray:
        return self.counts.astype(float) if weighted else np.ones(len(self))

    def write_csv(self, path: str | Path) -> Path:
        path = Path(path)
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["file", "n", "kind", "count", *self.keys])
            for i in range(len(self)):
                w.writerow([self.files[i], self.index[i] + 1, self.kinds[i], self.counts[i],
                            *[("" if not math.isfinite(v) else f"{v:.5g}") for v in self.data[i]]])
        return path


def _wquantile(x: np.ndarray, w: np.ndarray, q: float) -> float:
    order = np.argsort(x)
    x, w = x[order], w[order]
    cum = (np.cumsum(w) - 0.5 * w) / w.sum()
    return float(np.interp(q, cum, x))


def describe(table: PairTable, keys: list[str] | None = None, kind: str | None = None,
             files: set[str] | None = None, weighted: bool = True) -> list[dict]:
    """Media, deviazione standard, minimo, quartili e massimo di ogni parametro. Con `weighted` ogni coppia
    pesa quanto le sue coppie equivalenti per simmetria (come nell'istogramma)."""
    m = table.mask(kind, files)
    out = []
    for k in keys or table.keys:
        x, w = table.column(k)[m], table.weights(weighted)[m]
        ok = np.isfinite(x)
        x, w = x[ok], w[ok]
        if not len(x):
            out.append({"key": k, "n": 0})
            continue
        mean = float(np.average(x, weights=w))
        std = float(math.sqrt(np.average((x - mean) ** 2, weights=w)))
        out.append({"key": k, "n": int(w.sum()), "pairs": len(x), "mean": mean, "std": std, "min": float(x.min()),
                    "q25": _wquantile(x, w, 0.25), "median": _wquantile(x, w, 0.5), "q75": _wquantile(x, w, 0.75),
                    "max": float(x.max())})
    return out


DESCRIBE_FIELDS = ("mean", "std", "min", "q25", "median", "q75", "max")


def write_describe_csv(rows: list[dict], path: str | Path) -> Path:
    """Scrive il risultato di `describe` (una riga per parametro, intestazioni fisse)."""
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "n", "pairs", *DESCRIBE_FIELDS])
        for r in rows:
            w.writerow([r["key"], r["n"], r.get("pairs", ""), *[f"{r[k]:.5g}" if k in r else "" for k in DESCRIBE_FIELDS]])
    return path


def correlation(table: PairTable, keys: list[str], kind: str | None = None, files: set[str] | None = None) -> np.ndarray:
    """Matrice di correlazione di Pearson tra i parametri scelti; NaN dove non è definita (parametro costante o
    meno di 3 coppie)."""
    m = table.mask(kind, files)
    X = np.column_stack([table.column(k)[m] for k in keys])
    X = X[np.isfinite(X).all(axis=1)]
    n = len(keys)
    out = np.full((n, n), np.nan)
    if len(X) < 3:
        return out
    sd = X.std(axis=0)
    for i in range(n):
        for j in range(n):
            if sd[i] > 0 and sd[j] > 0:
                out[i, j] = float(np.corrcoef(X[:, i], X[:, j])[0, 1])
    return out


def histogram(table: PairTable, key: str, bins: int = 15, weighted: bool = True,
              files: set[str] | None = None) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Estremi degli intervalli e conteggi per tipo di coppia (stessi estremi per tutti i tipi)."""
    m = table.mask(None, files)
    x = table.column(key)
    ok = m & np.isfinite(x)
    if not ok.any():
        return np.zeros(0), {}
    edges = np.histogram_bin_edges(x[ok], bins=bins)
    w = table.weights(weighted)
    return edges, {k: np.histogram(x[ok & (table.kinds == k)], bins=edges, weights=w[ok & (table.kinds == k)])[0] for k in KINDS}


def summary_rows(entries: list[StructureEntry]) -> list[dict]:
    """Una riga per file: quante coppie di ogni tipo e le medie dei parametri degli stack."""
    rows = []
    for e in entries:
        stacks = [p for p in e.pairs if p.kind == "stack"]
        row = {"file": e.name, "error": e.error, "pairs": e.n_pairs, **{f"n_{k}": e.counts.get(k, 0) for k in KINDS},
               "columns": len(e.columns)}
        for k in ("interplanar", "long_slip", "lat_slip", "dihedral"):
            vals = [p.values.get(k, np.nan) for p in stacks]
            row[f"mean_stack_{k}"] = float(np.nanmean(vals)) if vals and np.isfinite(vals).any() else float("nan")
        rows.append(row)
    return rows
