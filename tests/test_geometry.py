"""Fase 4a: geometria delle molecole (legami, angoli, torsioni) con le incertezze propagate dalle s.u. del CIF."""
import csv
import math
import re
from pathlib import Path

import numpy as np
import pytest

from pistak.core.analysis import analyse_file
from pistak.core.cif import _blocks, cif_number, cif_su, tokenize
from pistak.core.intramolecular import PRIME, angle, molecule_geometry, torsion, unique_molecules

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
REAL = ROOT / "DATI_REALI"


def test_misure_elementari():
    p = np.array([[1.0, 0, 0], [0, 0, 0], [0, 1.0, 0], [0, 1.0, 1.0]])
    assert angle(p[:3]) == pytest.approx(90)
    assert torsion(p) == pytest.approx(-90) and torsion(p[::-1]) == pytest.approx(-90)  # il verso non dipende dalla lettura
    p[3] = [1.0, 1.0, 0]
    assert torsion(p) == pytest.approx(0, abs=1e-9)


def test_molecola_sul_centro_di_inversione():
    a = analyse_file(EX / "demo_naphthalene.cif")
    assert unique_molecules(a.structure) == [0]
    g = molecule_geometry(a.structure, a.cif)
    bonds = [x for x in g if x.kind == "bond"]
    assert len(bonds) == 10  # 11 C–C + 8 C–H, metà per la simmetria interna (C3–C3′ una volta)
    assert any(PRIME in lab for x in bonds for lab in x.atoms)
    assert all(math.isnan(x.su) for x in g)  # il CIF dimostrativo non ha s.u.
    cc = [x.value for x in bonds if all(lab.startswith("C") for lab in x.atoms)]
    assert min(cc) > 1.3 and max(cc) < 1.5
    assert {x.kind for x in g} == {"bond", "angle", "torsion"}
    no_h = molecule_geometry(a.structure, a.cif, include_h=False)
    assert all(not lab.startswith("H") for x in no_h for lab in x.atoms) and len(no_h) < len(g)


def test_molecole_equivalenti_misurate_una_volta():
    a = analyse_file(EX / "herringbone_anthracene.cif")
    assert len(unique_molecules(a.structure)) < len(a.structure.molecules)


def test_csv(tmp_path):
    from pistak import api
    from pistak.io import export

    a = api.analyse_file(EX / "demo_naphthalene.cif")
    rows = list(csv.DictReader(export.write_geometry_csv(a, tmp_path / "g.csv").open(encoding="utf-8")))
    assert {r["kind"] for r in rows} == {"bond", "angle", "torsion"} and rows[0]["su"] == ""
    assert all(r["atom_4"] for r in rows if r["kind"] == "torsion") and not any(r["atom_3"] for r in rows if r["kind"] == "bond")


def _cif_geometry(path: Path) -> dict[tuple, str]:
    """Le tabelle _geom_bond / _geom_angle del CIF (solo le misure senza codici di simmetria)."""
    b = next(x for x in _blocks(tokenize(path.read_text(encoding="utf-8"))) if x.loops)
    out = {}
    for lp in b.loops:
        col = {t: i for i, t in enumerate(lp.tags)}
        if "_geom_bond_distance" in col:
            for r in lp.rows:
                out[("bond", *sorted((r[col["_geom_bond_atom_site_label_1"]], r[col["_geom_bond_atom_site_label_2"]])))] = \
                    r[col["_geom_bond_distance"]]
        if "_geom_angle" in col:
            for r in lp.rows:
                ends = sorted((r[col["_geom_angle_atom_site_label_1"]], r[col["_geom_angle_atom_site_label_3"]]))
                out[("angle", r[col["_geom_angle_atom_site_label_2"]], *ends)] = r[col["_geom_angle"]]
    return out


@pytest.mark.parametrize("name", ["4500022.cif", "UNITCELL.cif"])
def test_come_le_tabelle_del_cif_reale(name):
    path = REAL / name
    if not path.is_file():
        pytest.skip(f"{path} non c'è (CIF reali non distribuiti)")
    a = analyse_file(path)
    ref = _cif_geometry(path)
    matched = 0
    for g in molecule_geometry(a.structure, a.cif, kinds=("bond", "angle")):
        labs = [x.rstrip(PRIME) for x in g.atoms]
        key = ("bond", *sorted(labs)) if g.kind == "bond" else ("angle", labs[1], *sorted((labs[0], labs[2])))
        if key not in ref:
            continue
        matched += 1
        value, su = cif_number(ref[key]), cif_su(ref[key])
        unit = 10.0 ** -len(re.match(r"[-\d]*\.?(\d*)", ref[key]).group(1))  # l'ultima cifra scritta nel CIF
        # il programma di affinamento misura sulle coordinate non arrotondate: con un asse di 50 Å l'arrotondamento
        # delle coordinate frazionarie del CIF sposta le distanze fino a ~0,002 Å; gli idrogeni in posizione calcolata
        # (senza s.u.) hanno coordinate scritte con meno cifre
        tol = (0.003, 0.15) if math.isfinite(g.su) else (0.01, 0.6)
        assert g.value == pytest.approx(value, abs=max(tol[g.kind != "bond"], unit / 2 + 1e-9)), key
        if su and math.isfinite(g.su):
            # un'unità sull'ultima cifra (arrotondamento), o il 35% (il CIF include anche le s.u. di cella)
            assert abs(g.su - su) <= max(unit, 0.35 * su) + 1e-9, (key, g.su, su)
    assert matched > 100
