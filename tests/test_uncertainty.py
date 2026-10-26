"""Incertezze dal Monte Carlo sulle s.u. del CIF: lettura delle s.u., riproduzione della geometria nominale, confronto con
casi analitici, scala lineare, coerenza con la simmetria, export."""
import json
import math
import re
from pathlib import Path

import numpy as np
import pytest

from pistak import api
from pistak.core.analysis import analyse_text
from pistak.core.cif import cif_su
from pistak.core.config import Limits
from pistak.core.pairs import remeasure
from pistak.core.uncertainty import GEOMETRY_KEYS, monte_carlo, perturbed

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
ATOM = re.compile(r"^(\S+ [A-Z][a-z]?) (-?\d+\.\d+) (-?\d+\.\d+) (-?\d+\.\d+)$", re.M)


def with_su(text: str, coord_su: int = 0, cell_su: dict[str, str] | None = None) -> str:
    """Il CIF con una s.u. (in unità dell'ultima cifra) su tutte le coordinate e, a scelta, sui parametri di cella."""
    if coord_su:
        text = ATOM.sub(lambda m: f"{m.group(1)} {m.group(2)}({coord_su}) {m.group(3)}({coord_su}) {m.group(4)}({coord_su})", text)
    for tag, su in (cell_su or {}).items():
        text = re.sub(rf"^(_cell_{tag}\s+)(\S+)$", rf"\g<1>\g<2>({su})", text, flags=re.M)
    return text


def test_lettura_delle_su():
    assert cif_su("10.3528(16)") == pytest.approx(0.0016)
    assert cif_su("-0.15128(16)") == pytest.approx(0.00016)
    assert cif_su("96.30(2)") == pytest.approx(0.02) and cif_su("123(4)") == 4
    assert cif_su("0.5") == cif_su("?") == cif_su(None) == 0
    a = analyse_text(with_su((EX / "demo_naphthalene.cif").read_text(), 15, {"length_a": "3"}), "d", Limits())
    assert np.allclose(a.cif.atoms[0].su, 15e-5) and a.cif.cell_su[0] == pytest.approx(3e-5)


def test_senza_su_tutte_le_sigma_sono_zero():
    a = api.analyse_file(EX / "herringbone_anthracene.cif", Limits())
    u = api.uncertainty(a, samples=8)
    assert not u.has_su and a.uncertainty is u
    assert all(v == 0 for s in u.sigma.values() for v in s.values() if math.isfinite(v))


def test_senza_perturbazione_si_riottiene_la_geometria_nominale():
    a = analyse_text(with_su((EX / "herringbone_anthracene.cif").read_text(), 20), "h", Limits())
    a2 = perturbed(a, np.random.default_rng(1), scale=0.0)
    for p in a.pairs:
        q = remeasure(a2.structure, a2.fragments, p, a2.limits)
        for k in GEOMETRY_KEYS:
            assert getattr(q, k) == pytest.approx(getattr(p, k), abs=1e-9), k


def test_caso_analitico_distanza_tra_i_centroidi_dello_stack_lungo_a():
    # stack per traslazione lungo a: la distanza tra i centroidi è |a|, quindi σ = s.u.(a), se si perturba solo la cella
    text = with_su((EX / "demo_naphthalene.cif").read_text(), 0, {"length_a": "100"})  # s.u.(a) = 0.00100 Å
    a = analyse_text(text, "d", Limits())
    p = a.pairs[0]
    assert p.kind == "stack" and p.centroid_distance == pytest.approx(a.cif.cell.a, abs=1e-4)
    u = monte_carlo(a, samples=400, seed=3)
    assert u.get(p, "centroid_distance") == pytest.approx(0.001, rel=0.12)
    assert monte_carlo(a, samples=50, include_cell=False).get(p, "centroid_distance") == 0


def test_le_sigma_crescono_linearmente_con_le_su():
    base = (EX / "herringbone_anthracene.cif").read_text()
    s1 = monte_carlo(analyse_text(with_su(base, 20), "h", Limits()), samples=300, seed=5)
    s2 = monte_carlo(analyse_text(with_su(base, 40), "h", Limits()), samples=300, seed=5)
    for key in ("interplanar", "long_slip", "contact"):
        assert s2.get(0, key) / s1.get(0, key) == pytest.approx(2.0, rel=0.05)


def test_le_copie_per_simmetria_si_muovono_insieme():
    # lo stack dell'herringbone è tra una molecola e la sua traslata: piani sempre paralleli in ogni campione
    a = analyse_text(with_su((EX / "herringbone_anthracene.cif").read_text(), 30), "h", Limits())
    u = monte_carlo(a, samples=60, seed=2)
    p = a.pairs[0]
    assert p.kind == "stack" and p.shift != (0, 0, 0) and p.a == p.b
    assert u.get(p, "dihedral") < 1e-6 and u.get(p, "interplanar") > 1e-4
    tilted = next(q for q in a.pairs if q.kind == "tilted")
    assert u.get(tilted, "dihedral") > 1e-4  # tra molecole diverse l'angolo invece oscilla
    assert all(0 <= v <= 1 for v in u.kind_agreement.values()) and u.kind_agreement[p.index] == 1.0


def test_stesso_seme_stessi_risultati():
    a = analyse_text(with_su((EX / "dimers_naphthalene.cif").read_text(), 25), "d", Limits())
    u1, u2, u3 = (monte_carlo(a, samples=30, seed=s) for s in (7, 7, 8))
    assert u1.sigma == u2.sigma and u1.sigma != u3.sigma


def test_registro_e_export(tmp_path):
    from pistak.io import export

    a = api.analyse_text(with_su((EX / "herringbone_anthracene.cif").read_text(), 20), "h", Limits())
    u = api.uncertainty(a, samples=20, keys=["interplanar", "long_slip_signed", "beta_a"])
    assert set(u.keys) == {"interplanar", "long_slip_signed", "beta_a"} and u.has_su
    head = export.write_csv(a, tmp_path / "p.csv").read_text().splitlines()[0].split(",")
    assert {"interplanar_su", "long_slip_signed_su", "beta_a_su", "kind_agreement"} <= set(head)
    data = json.loads(export.write_json(a, tmp_path / "a.json").read_text())
    assert data["summary"]["uncertainty"]["samples"] == 20
    assert set(data["pairs"][0]["su"]) == {"interplanar", "long_slip_signed", "beta_a"}


def test_terminale(capsys):
    from pistak.cli import main

    assert main(["analyze", str(EX / "demo_naphthalene.cif"), "--mc", "5"]) == 0
    captured = capsys.readouterr()
    assert "Monte Carlo su 5 campioni" in captured.out and "non ha incertezze standard" in captured.err
