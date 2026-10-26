"""Analisi di molti file CIF in parallelo.

Ogni file diventa un `StructureEntry`: riepilogo del cristallo, una `PairRow` per ogni coppia (con tutti i
parametri del registro) e le colonne di stack. Sono i dati su cui lavorano `stats.py` e `compare.py` e la
finestra «Confronto e statistiche».
"""
from __future__ import annotations

import csv
import math
import os
from collections.abc import Callable
from concurrent.futures import BrokenExecutor, ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from pistak.core.config import Limits

#: Parametri riportati per la coppia di π-stack migliore (la prima: la più vicina tra gli stack).
SUMMARY_KEYS = ["interplanar", "long_slip", "lat_slip", "dihedral", "centroid_distance", "contact", "overlap_area", "beta_a"]


@dataclass
class PairRow:
    file: str
    index: int  # posizione nell'elenco delle coppie di quel file (come nella GUI)
    kind: str
    count: int
    values: dict[str, float]


@dataclass
class StructureEntry:
    path: str
    name: str
    #: Messaggio leggibile, nella lingua del processo che ha analizzato il file (il codice stabile è `error_code`).
    error: str | None = None
    space_group: str = ""
    cell: tuple[float, ...] = ()
    formulas: str = ""
    n_pairs: int = 0
    counts: dict[str, int] = field(default_factory=dict)
    pairs: list[PairRow] = field(default_factory=list)
    columns: list[dict] = field(default_factory=list)
    error_code: str | None = None

    @property
    def best(self) -> dict[str, float]:
        """Valori della coppia π-stack migliore (la prima nell'elenco), o {} se non ce ne sono."""
        stacks = [p for p in self.pairs if p.kind == "stack"]
        return dict(stacks[0].values) if stacks else {}


BatchResult = StructureEntry  # nome usato dalla prima versione


def analyse_one(path: str, limits: Limits) -> StructureEntry:
    """Analizza un file e ne riassume il risultato (eseguita nei processi di lavoro)."""
    from pistak.api import analyse_file
    from pistak.errors import describe_error, error_code

    name = Path(path).name
    try:
        a = analyse_file(path, limits)
    except Exception as e:  # noqa: BLE001
        return StructureEntry(path, name, error=describe_error(e), error_code=error_code(e))
    c = a.cif.cell
    counts = {k: 0 for k in ("stack", "parallel", "tilted")}
    for p in a.pairs:
        counts[p.kind] += p.count
    rows = [PairRow(name, p.index, p.kind, p.count, {k: p.get(k) for k in a.param_keys}) for p in a.pairs]
    cols = []
    for col in a.stacks:
        cols.append({"kind_key": col.kind_key, "rank": col.rank, "molecules": len(col.molecules), "repeat": col.repeat,
                     "interplanar": col.mean("interplanar"), "long_slip": col.mean("long_slip"), "lat_slip": col.mean("lat_slip")})
    return StructureEntry(
        path, name, None, a.cif.space_group, (c.a, c.b, c.c, c.alpha, c.beta, c.gamma),
        "; ".join(f"{g.count}×{g.text}" for g in a.molecule_groups()), len(a.pairs), counts, rows, cols,
    )


def run_batch(files: list[Path], limits: Limits, workers: int = 0,
              progress: Callable[[int, int, str], None] | None = None) -> list[StructureEntry]:
    """Analizza i file in parallelo. `progress(fatti, totale, nome)` viene chiamata a ogni file finito.
    L'ordine del risultato è quello di `files`."""
    workers = workers or max(1, min(len(files), (os.cpu_count() or 2) - 1))
    paths = [str(f) for f in files]
    out: dict[str, StructureEntry] = {}
    if workers == 1:
        for i, p in enumerate(paths, 1):
            out[p] = analyse_one(p, limits)
            if progress:
                progress(i, len(paths), out[p].name)
    else:
        try:
            with ProcessPoolExecutor(max_workers=workers) as ex:
                futures = {ex.submit(analyse_one, p, limits): p for p in paths}
                for i, fut in enumerate(as_completed(futures), 1):
                    p = futures[fut]
                    try:
                        out[p] = fut.result()
                    except Exception as e:  # noqa: BLE001  (processo morto, ecc.)
                        out[p] = StructureEntry(p, Path(p).name, error=f"{type(e).__name__}: {e}", error_code=type(e).__name__)
                    if progress:
                        progress(i, len(paths), out[p].name)
        except (RuntimeError, OSError, BrokenExecutor):
            # i processi non partono (ambiente senza multiprocessing, modulo principale senza guardia __main__…):
            # si prosegue in serie sui file non ancora fatti, invece di lasciare l'elaborazione a metà
            for p in paths:
                if p not in out:
                    out[p] = analyse_one(p, limits)
                    if progress:
                        progress(len(out), len(paths), out[p].name)
    return [out[p] for p in paths]


def find_cifs(folder: str | Path) -> list[Path]:
    """I .cif di una cartella (e sottocartelle), ignorando i metadati '._*' dei volumi exFAT."""
    folder = Path(folder)
    if folder.is_file():
        return [folder]
    return sorted(p for p in folder.rglob("*") if p.suffix.lower() == ".cif" and not p.name.startswith("."))


def write_summary_csv(results: list[StructureEntry], path: str | Path) -> Path:
    path = Path(path)
    head = ["file", "error", "space_group", "a", "b", "c", "alpha", "beta", "gamma", "molecules", "pairs",
            "n_stack", "n_parallel", "n_tilted", "columns"] + [f"stack1_{k}" for k in SUMMARY_KEYS]
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(head)
        for r in results:
            cell = [f"{x:.4f}" for x in r.cell] if r.cell else [""] * 6
            best = r.best
            vals = [("" if not math.isfinite(best.get(k, math.nan)) else f"{best[k]:.4g}") for k in SUMMARY_KEYS]
            w.writerow([r.name, r.error or "", r.space_group, *cell, r.formulas, r.n_pairs,
                        r.counts.get("stack", ""), r.counts.get("parallel", ""), r.counts.get("tilted", ""),
                        "; ".join(c["kind_key"] for c in r.columns), *vals])
    return path
