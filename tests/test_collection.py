"""Batch con righe per coppia, statistiche su più file e confronto tra strutture."""
import math
from pathlib import Path

import numpy as np
import pytest

from pistak.core.config import Limits
from pistak.io.batch import analyse_one, find_cifs, run_batch, write_summary_csv
from pistak.io.compare import compare, motif_distance, write_csv
from pistak.io.stats import PairTable, correlation, describe, histogram, summary_rows

EX = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture(scope="module")
def entries():
    files = find_cifs(EX)
    return run_batch(files, Limits(), workers=1)


def by_name(entries, name):
    return next(e for e in entries if e.name == name)


def test_find_cifs_ignora_i_metadati_di_macos(tmp_path):
    (tmp_path / "a.cif").write_text("x")
    (tmp_path / "._a.cif").write_text("x")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "B.CIF").write_text("x")
    assert [p.name for p in find_cifs(tmp_path)] == ["B.CIF", "a.cif"] or len(find_cifs(tmp_path)) == 2
    assert all(not p.name.startswith(".") for p in find_cifs(tmp_path))


def test_ogni_file_ha_coppie_e_colonne(entries):
    assert len(entries) == 6 and all(e.error is None for e in entries)
    demo = by_name(entries, "demo_naphthalene.cif")
    assert [p.kind for p in demo.pairs] == ["stack"]
    assert demo.pairs[0].values["interplanar"] == pytest.approx(3.4, abs=1e-3)
    assert demo.columns and demo.columns[0]["rank"] == 1
    assert by_name(entries, "infinite_chain.cif").pairs == []


def test_file_illeggibile_non_ferma_il_resto(tmp_path):
    bad = tmp_path / "bad.cif"
    bad.write_text("data_x\n_cell_length_a 5\n")
    res = run_batch([bad, EX / "demo_naphthalene.cif"], Limits(), workers=1)
    assert res[0].error and "cella" in res[0].error
    assert res[1].error is None


def test_riepilogo_csv(entries, tmp_path):
    out = write_summary_csv(entries, tmp_path / "s.csv")
    lines = out.read_text().splitlines()
    assert len(lines) == 7 and lines[0].startswith("file,errore")


def test_tabella_e_descrizione(entries):
    t = PairTable.from_entries(entries)
    assert len(t) == sum(len(e.pairs) for e in entries)
    stack = describe(t, ["interplanar"], "stack", weighted=False)[0]
    vals = t.column("interplanar")[t.kinds == "stack"]
    assert stack["pairs"] == len(vals) and stack["mean"] == pytest.approx(vals.mean())
    assert stack["min"] <= stack["q25"] <= stack["median"] <= stack["q75"] <= stack["max"]
    only_demo = describe(t, ["interplanar"], "stack", files={"demo_naphthalene.cif"})[0]
    assert only_demo["mean"] == pytest.approx(3.4, abs=1e-3) and only_demo["std"] == pytest.approx(0, abs=1e-6)


def test_istogramma_somma_i_pesi(entries):
    t = PairTable.from_entries(entries)
    edges, by_kind = histogram(t, "interplanar", bins=10, weighted=True)
    assert len(edges) == 11
    total = sum(h.sum() for h in by_kind.values())
    assert total == t.counts[np.isfinite(t.column("interplanar"))].sum()


def test_correlazione(entries):
    t = PairTable.from_entries(entries)
    c = correlation(t, ["interplanar", "centroid_distance", "long_slip"], None)
    assert c.shape == (3, 3) and np.allclose(np.diag(c)[np.isfinite(np.diag(c))], 1)
    assert np.allclose(c, c.T, equal_nan=True)
    # un parametro costante non ha correlazione definita
    small = PairTable.from_entries([by_name(entries, "demo_naphthalene.cif")])
    assert np.isnan(correlation(small, ["interplanar", "long_slip"], None)).all()


def test_confronto_con_riferimento(entries):
    ref = entries.index(by_name(entries, "demo_naphthalene.cif"))
    rows = compare(entries, ref)
    r_ref = rows[ref]
    assert r_ref.motif_distance == 0 and all(abs(v) < 1e-12 for v in r_ref.delta.values())
    other = next(r for r in rows if r.name == "herringbone_anthracene.cif")
    assert other.delta["interplanar"] == pytest.approx(other.best["interplanar"] - 3.4, abs=1e-3)
    expected = motif_distance(by_name(entries, "herringbone_anthracene.cif"), entries[ref])
    assert other.motif_distance == pytest.approx(expected) and math.isfinite(expected)
    nostack = next(r for r in rows if r.name == "infinite_chain.cif")
    assert math.isnan(nostack.motif_distance) and nostack.best == {}


def test_csv_di_confronto(entries, tmp_path):
    out = write_csv(compare(entries, 0), tmp_path / "c.csv")
    assert len(out.read_text().splitlines()) == 7


def test_riepilogo_per_file(entries):
    rows = summary_rows(entries)
    d = next(r for r in rows if r["file"] == "demo_naphthalene.cif")
    assert d["n_stack"] == 2 and d["mean_stack_interplanar"] == pytest.approx(3.4, abs=1e-3)


def test_analyse_one_diretto():
    e = analyse_one(str(EX / "herringbone_anthracene.cif"), Limits())
    assert e.counts["stack"] == 4 and e.counts["tilted"] == 8
