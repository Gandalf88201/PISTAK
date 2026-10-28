"""Fase 4a: geometria delle molecole (legami, angoli, torsioni, contatti di van der Waals) con le incertezze propagate
dalle s.u. del CIF e i codici di simmetria."""
import csv
import math
import re
from pathlib import Path

import numpy as np
import pytest

from pistak.core.analysis import analyse_file
from pistak.core.cif import _blocks, cif_number, cif_su, tokenize
from pistak.core.elements import vdw_radius
from pistak.core.intramolecular import angle, molecule_geometry, symmetry_codes, symop_text, torsion, unique_molecules, vdw_contacts

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
COD = ROOT / "tests" / "data" / "cod"
REAL = ROOT / "DATI_REALI"  # CIF non distribuibili (CSD): se mancano, i test che li usano si saltano


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
    assert len(bonds) == 10  # 11 C–C + 8 C–H, metà per la simmetria interna (C3–C3#1 una volta)
    labels, legend = symmetry_codes(bonds)
    assert legend == [("#1", "-x,-y,-z")] and ("C3", "C3#1") in labels
    assert sum(any(x.codes) for x in bonds) == 2  # tra le copie uguali vince quella con meno atomi trasformati
    assert all(math.isnan(x.su) for x in g)  # il CIF dimostrativo non ha s.u.
    cc = [x.value for x in bonds if all(lab.startswith("C") for lab in x.atoms)]
    assert min(cc) > 1.3 and max(cc) < 1.5
    assert {x.kind for x in g} == {"bond", "angle", "torsion"}  # una molecola sola in una cella larga: nessun contatto
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
    assert any(r["symop_2"] == "-x,-y,-z" for r in rows) and not any(r["vdw_sum"] for r in rows)
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


@pytest.mark.parametrize("path", [COD / "4500022.cif", COD / "2311088.cif", COD / "2100417.cif", REAL / "UNITCELL.cif"],
                         ids=lambda p: p.name)
def test_come_le_tabelle_del_cif_reale(path):
    if not path.is_file():
        pytest.skip(f"{path} non c'è (CIF reali non distribuiti)")
    a = analyse_file(path)
    ref = _cif_geometry(path)
    matched = 0
    for g in molecule_geometry(a.structure, a.cif, kinds=("bond", "angle")):
        labs = list(g.atoms)
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
            # un'unità sull'ultima cifra (arrotondamento), o il 35% (il CIF include anche le s.u. di cella); il 45% per le
            # misure tra un atomo e la sua immagine per simmetria, dove pesano le correlazioni tra coordinate che il CIF
            # non riporta (SHELXL usa la matrice di covarianza completa)
            rel = 0.45 if any(g.codes) else 0.35
            assert abs(g.su - su) <= max(unit, rel * su) + 1e-9, (key, g.su, su)
    assert matched >= 10


def test_testo_delle_operazioni():
    assert symop_text(-np.eye(3), np.array([1.0, 0.0, 0.5])) == "-x+1,-y,-z+1/2"
    assert symop_text(np.array([[1, -1, 0], [1, 0, 0], [0, 0, 1.0]]), np.array([0, 0, 1 / 3])) == "x-y,x,z+1/3"
    assert symop_text(np.eye(3), np.zeros(3)) == "x,y,z"


def _hbonds(path: Path) -> list[tuple[str, str, float]]:
    b = next(x for x in _blocks(tokenize(path.read_text(encoding="utf-8"))) if x.loops)
    lp = next(lp for lp in b.loops if "_geom_hbond_distance_ha" in lp.tags)
    col = {t: i for i, t in enumerate(lp.tags)}
    return [(r[col["_geom_hbond_atom_site_label_h"]], r[col["_geom_hbond_atom_site_label_a"]],
             cif_number(r[col["_geom_hbond_distance_ha"]])) for r in lp.rows]


def test_contatti_come_i_legami_h_del_cif():
    # ogni H···A della tabella _geom_hbond più corto della somma dei raggi di van der Waals è tra i contatti
    a = analyse_file(COD / "4500022.cif")
    contacts = vdw_contacts(a.structure, a.cif)
    assert contacts and all(c.delta < 0 for c in contacts)
    assert [c.delta for c in contacts] == sorted(c.delta for c in contacts)
    found = {(tuple(sorted(c.atoms)), round(c.value, 2)) for c in contacts}
    checked = 0
    for h, acc, d in _hbonds(COD / "4500022.cif"):
        if d >= vdw_radius("H") + vdw_radius(acc[0]) - 0.01:
            assert not any(k == tuple(sorted((h, acc))) and abs(v - d) <= 0.011 for k, v in found), (h, acc, d)
            continue
        checked += 1
        assert any(k == tuple(sorted((h, acc))) and abs(v - d) <= 0.011 for k, v in found), (h, acc, d)
    assert checked >= 10
    unique = set(unique_molecules(a.structure))
    assert all(c.molecule in unique for c in contacts)
    assert len({(tuple(sorted(c.atoms)), round(c.value, 3)) for c in contacts}) == len(contacts)  # niente doppioni
    _labels, legend = symmetry_codes(contacts)
    assert legend and all(op.count(",") == 2 for _n, op in legend)
    assert all(math.isfinite(c.su) for c in contacts if not any(x.startswith("H") for x in c.atoms))


def test_contatti_senza_idrogeni_e_nessun_contatto_nella_stessa_molecola():
    a = analyse_file(EX / "herringbone_anthracene.cif")
    all_c = vdw_contacts(a.structure, a.cif)
    no_h = vdw_contacts(a.structure, a.cif, include_h=False)
    assert no_h and len(no_h) < len(all_c) and all("H" not in x[0] for c in no_h for x in c.atoms)
    bonded = {tuple(sorted((m.labels[i], m.labels[j]))) for m in a.structure.molecules for i, j in m.bonds}
    assert all(c.value > 1.0 for c in all_c) and not any(tuple(sorted(c.atoms)) in bonded and c.value < 1.7 for c in all_c)


def test_json_completo(tmp_path):
    import json

    from pistak import api
    from pistak.io import export

    a = api.analyse_file(COD / "4500022.cif")
    data = json.loads(export.write_json(a, tmp_path / "a.json").read_text(encoding="utf-8"))
    assert data["version"] == export.JSON_VERSION == 4
    assert {"symops", "molecule_list", "stacks", "ring_pairs", "geometry"} <= set(data)
    g = data["geometry"]
    items = molecule_geometry(a.structure, a.cif)
    assert [len(g[k]) for k in ("bonds", "angles", "torsions", "contacts")] == \
        [sum(x.kind == k for x in items) for k in ("bond", "angle", "torsion", "contact")]
    hb = g["contacts"][0]
    assert hb["delta"] == pytest.approx(hb["value"] - hb["vdw_sum"]) and len(hb["symops"]) == 2 and hb["molecule_b"]
    assert sum(m["independent"] for m in data["molecule_list"]) == len(unique_molecules(a.structure))
    assert data["ring_pairs"] and all(r["pair"] >= 1 for r in data["ring_pairs"]) and len(data["symops"]) == 4
