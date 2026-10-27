"""Fase 3: frammenti π multipli per molecola, coppie frammento–frammento, riepilogo per molecola, asse lungo definito
dall'utente, slittamenti con segno."""
import json
import math
from pathlib import Path

import numpy as np
import pytest

import cifgen
from pistak.core.analysis import analyse_file, analyse_text
from pistak.core.cif import Cell
from pistak.core.config import Limits
from pistak.core.cores import pi_core, pi_fragments
from pistak.core.pairs import class_spread
from pistak.core.structure import cell_vectors
from pistak.labels import relation
from pistak.params import compute_all

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"


def naphthyl_benzene(twist: float = 50.0) -> list[tuple[str, np.ndarray]]:
    """2-fenilnaftalene: naftalene nel piano xy più un benzene legato al C2 e ruotato di `twist` gradi intorno al
    legame. Due sistemi di anelli separati (non condensati): due frammenti π, da 10 e da 6 atomi."""
    naph = cifgen.ring_system(2)
    carbons = [p for el, p in naph if el == "C"]
    ipso_n = max(carbons, key=lambda p: (round(p[0], 6), p[1]))  # C2 (x massimo, y positiva)
    u = np.array([math.cos(math.pi / 6), math.sin(math.pi / 6), 0.0])
    # via l'idrogeno del C2: al suo posto il benzene
    atoms = [(el, p) for el, p in naph if not (el == "H" and np.linalg.norm(p - ipso_n - 1.08 * u) < 1e-6)]
    v = np.array([-u[1], u[0], 0.0])
    w = np.cross(u, v)
    t = math.radians(twist)
    v2 = math.cos(t) * v + math.sin(t) * w
    center = ipso_n + (1.48 + 1.4) * u
    for k in range(6):
        ang = math.pi + k * math.pi / 3  # k = 0: il carbonio ipso, verso il naftalene
        d = math.cos(ang) * u + math.sin(ang) * v2
        atoms.append(("C", center + 1.4 * d))
        if k:
            atoms.append(("H", center + 2.48 * d))
    return atoms


def stacked_crystal(twist: float = 50.0) -> str:
    """P-1 con la molecola impilata per traslazione lungo a (≈3,75 Å, quasi perpendicolare al naftalene)."""
    mine = np.array([[1.2, 0.6, 3.5], [0.0, 16.0, 0.0], [-30.0, 0.0, 0.0]])
    ln = np.linalg.norm(mine, axis=1)
    ang = lambda i, j: math.degrees(math.acos(mine[i] @ mine[j] / (ln[i] * ln[j])))  # noqa: E731
    cell = Cell(ln[0], ln[1], ln[2], ang(1, 2), ang(0, 2), ang(0, 1))
    to_cell = np.linalg.inv(mine) @ cell_vectors(cell)  # stesse coordinate frazionarie nel sistema di cell_vectors
    if np.linalg.det(to_cell) < 0:
        raise AssertionError("terna sinistrorsa")
    mol = [(el, p @ to_cell) for el, p in naphthyl_benzene(twist)]
    return cifgen.make_cif(mol, cell, ["x, y, z", "-x, -y, -z"], "fenilnaftalene", "P -1", origin_frac=(0.25, 0.25, 0.25))


@pytest.fixture(scope="module")
def crystal():
    return stacked_crystal()


# ---- frammenti -------------------------------------------------------------------------------------------
def test_frammenti_di_una_molecola(crystal):
    a = analyse_text(crystal, "fn", Limits())
    m = a.structure.molecules[0]
    assert [len(f.atoms) for f in a.fragments[0]] == [10, 6]
    assert pi_core(m, a.limits).atoms.tolist() == a.fragments[0][0].atoms.tolist()  # il principale è il nucleo di sempre
    assert len(pi_fragments(m, a.limits.with_values(fragments="largest"))) == 1
    assert len(pi_fragments(m, a.limits.with_values(min_fragment_atoms=7))) == 1
    g = a.molecule_groups()[0]
    assert g.fragments == 2 and g.with_core == 2


def test_coppie_frammento_frammento(crystal):
    full = compute_all(analyse_text(crystal, "fn", Limits()))
    largest = analyse_text(crystal, "fn", Limits(fragments="largest"))
    assert all((p.fa, p.fb) == (0, 0) for p in largest.pairs)
    combos = {(p.fa, p.fb) for p in full.pairs if not p.intramolecular}
    assert {(0, 0), (1, 1)} <= combos
    # lo stack naftalene–naftalene lungo a c'è in entrambi i modi, con gli stessi valori
    n_full = next(p for p in full.pairs if (p.fa, p.fb) == (0, 0) and p.kind == "stack")
    n_big = next(p for p in largest.pairs if p.kind == "stack")
    assert n_full.interplanar == pytest.approx(n_big.interplanar) == pytest.approx(3.5, abs=1e-3)
    # fenile–fenile lungo a: piani paralleli (stessa molecola traslata)
    ph = next(p for p in full.pairs if (p.fa, p.fb) == (1, 1) and p.dihedral < 1)
    assert ph.get("fragment_a") == 2 and ph.get("fragment_b") == 2 and "π2–π2" in relation(ph)
    for p in full.pairs:
        assert class_spread(full.structure, full.fragments, full.limits, p) < 1e-6


def test_coppie_intramolecolari(crystal):
    a = analyse_text(crystal, "fn", Limits())
    intra = [p for p in a.pairs if p.intramolecular]
    assert len(intra) == 1  # π1–π2 nella stessa molecola; le due molecole della cella sono equivalenti
    p = intra[0]
    assert (p.fa, p.fb) == (0, 1) and p.count == 2 and p.kind == "tilted"
    assert p.dihedral == pytest.approx(50.0, abs=0.5) and relation(p).startswith("intramolecolare · π1–π2")


def test_riepilogo_per_molecola(crystal):
    a = analyse_text(crystal, "fn", Limits())
    contacts = a.molecule_contacts()
    assert sum(len(c.pairs) for c in contacts) == len(a.pairs)
    along_a = next(c for c in contacts if not c.intramolecular and c.a == c.b and tuple(abs(x) for x in c.shift) == (1, 0, 0))
    assert len(along_a.pairs) >= 2 and along_a.count == 2  # una per molecola della cella (Z = 2)
    assert along_a.best.kind == "stack" and along_a.best.fa == along_a.best.fb == 0
    intra = [c for c in contacts if c.intramolecular]
    assert len(intra) == 1 and intra[0].count == 2


def test_json_con_frammenti_e_contatti(crystal, tmp_path):
    from pistak.io import export

    a = compute_all(analyse_text(crystal, "fn", Limits()))
    data = json.loads(export.write_json(a, tmp_path / "a.json").read_text())
    assert data["version"] == 3
    assert [(f["molecule"], f["fragment"], f["atoms"]) for f in data["fragments"]] == [(1, 1, 10), (1, 2, 6), (2, 1, 10), (2, 2, 6)]
    assert {p["fragment_a"] for p in data["pairs"]} == {1, 2}
    assert sum(len(c["pairs"]) for c in data["molecule_contacts"]) == len(data["pairs"])


def test_esempi_con_un_solo_frammento_invariati():
    for cif in sorted(p for p in EX.glob("*.cif") if not p.name.startswith(".")):
        full, largest = analyse_file(cif, Limits()), analyse_file(cif, Limits(fragments="largest"))
        assert [(p.kind, round(p.interplanar, 9), p.count) for p in full.pairs] == \
               [(p.kind, round(p.interplanar, 9), p.count) for p in largest.pairs]


# ---- asse lungo definito dall'utente e slittamenti con segno --------------------------------------------
def _short_axis_labels(a) -> tuple[str, str]:
    """Due atomi del nucleo allineati lungo l'asse corto (stessa coordinata L, W opposta)."""
    m, core = a.structure.molecules[0], a.cores[0]
    W = np.cross(core.normal, core.long)
    best = None
    for i in core.atoms:
        for j in core.atoms:
            d = m.pos[j] - m.pos[i]
            score = abs(d @ W) - 10 * abs(d @ core.long)
            if best is None or score > best[0]:
                best = (score, i, j)
    return m.labels[best[1]], m.labels[best[2]]


def test_asse_lungo_dagli_atomi():
    cif = EX / "herringbone_anthracene.cif"
    auto = compute_all(analyse_file(cif, Limits()))
    i, j = _short_axis_labels(auto)
    turned = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{i},{j}")))
    assert all(c.axis_source == "atoms" for c in turned.cores)
    p0, p1 = auto.pairs[0], turned.pairs[0]
    assert p0.kind == p1.kind == "stack" and p1.interplanar == pytest.approx(p0.interplanar)
    # l'asse «lungo» ora è quello corto: i due slittamenti si scambiano
    assert p1.long_slip == pytest.approx(p0.lat_slip, abs=1e-3) and p1.lat_slip == pytest.approx(p0.long_slip, abs=1e-3)
    # etichette che non esistono: si resta sull'asse automatico
    assert all(c.axis_source == "pca" for c in analyse_file(cif, Limits(long_axis_atoms="X98,X99")).cores)


def test_invertire_gli_atomi_dell_asse_cambia_il_segno():
    # molecola fuori da centri di inversione (ogni etichetta una volta per molecola): il verso dell'asse è definito
    cif = EX / "dimers_naphthalene.cif"
    auto = analyse_file(cif, Limits())
    i, j = _short_axis_labels(auto)
    fwd = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{i},{j}")))
    rev = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{j},{i}")))
    for p, q in zip(fwd.pairs, rev.pairs):
        assert q.get("long_slip_signed") == pytest.approx(-p.get("long_slip_signed"), abs=1e-9)
        assert q.long_slip == pytest.approx(p.long_slip, abs=1e-9)


def test_slittamenti_con_segno():
    a = compute_all(analyse_file(EX / "dimers_naphthalene.cif", Limits()))
    for p in a.pairs:
        sl, sw = p.get("long_slip_signed"), p.get("lat_slip_signed")
        assert abs(sl) == pytest.approx(p.long_slip, abs=1e-9) and abs(sw) == pytest.approx(p.lat_slip, abs=1e-9)
        assert p.get("slip_direction") == pytest.approx(math.degrees(math.atan2(sw, sl)), abs=1e-9)
        assert 0 <= p.get("axis_aspect_a") <= 1


def test_isotropia_del_nucleo():
    naph = analyse_file(EX / "demo_naphthalene.cif", Limits()).cores[0].plane.aspect
    anth = analyse_file(EX / "herringbone_anthracene.cif", Limits()).cores[0].plane.aspect
    assert 0 < anth < naph < 1  # l'antracene è più allungato del naftalene


def test_asse_dalle_soglie():
    assert Limits.from_dict({"long_axis_atoms": " C1,C9 "}).long_axis_atoms == "C1,C9"
    assert Limits.from_dict({"long_axis_atoms": ["C1"]}).long_axis_atoms == ""
