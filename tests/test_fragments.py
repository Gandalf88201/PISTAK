"""Phase 3: several π fragments per molecule, fragment–fragment pairs, per-molecule summary, user-defined long axis,
signed slips."""
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
    """2-phenylnaphthalene: naphthalene in the xy plane plus a benzene bonded to C2 and turned by `twist` degrees about
    the bond. Two separate (not fused) ring systems: two π fragments, of 10 and 6 atoms."""
    naph = cifgen.ring_system(2)
    carbons = [p for el, p in naph if el == "C"]
    ipso_n = max(carbons, key=lambda p: (round(p[0], 6), p[1]))  # C2 (largest x, positive y)
    u = np.array([math.cos(math.pi / 6), math.sin(math.pi / 6), 0.0])
    # remove the C2 hydrogen: the benzene goes in its place
    atoms = [(el, p) for el, p in naph if not (el == "H" and np.linalg.norm(p - ipso_n - 1.08 * u) < 1e-6)]
    v = np.array([-u[1], u[0], 0.0])
    w = np.cross(u, v)
    t = math.radians(twist)
    v2 = math.cos(t) * v + math.sin(t) * w
    center = ipso_n + (1.48 + 1.4) * u
    for k in range(6):
        ang = math.pi + k * math.pi / 3  # k = 0: the ipso carbon, towards the naphthalene
        d = math.cos(ang) * u + math.sin(ang) * v2
        atoms.append(("C", center + 1.4 * d))
        if k:
            atoms.append(("H", center + 2.48 * d))
    return atoms


def stacked_crystal(twist: float = 50.0) -> str:
    """P-1 with the molecule stacked by translation along a (≈3.75 Å, nearly perpendicular to the naphthalene)."""
    mine = np.array([[1.2, 0.6, 3.5], [0.0, 16.0, 0.0], [-30.0, 0.0, 0.0]])
    ln = np.linalg.norm(mine, axis=1)
    ang = lambda i, j: math.degrees(math.acos(mine[i] @ mine[j] / (ln[i] * ln[j])))  # noqa: E731
    cell = Cell(ln[0], ln[1], ln[2], ang(1, 2), ang(0, 2), ang(0, 1))
    to_cell = np.linalg.inv(mine) @ cell_vectors(cell)  # same fractional coordinates in the cell_vectors frame
    if np.linalg.det(to_cell) < 0:
        raise AssertionError("left-handed axes")
    mol = [(el, p @ to_cell) for el, p in naphthyl_benzene(twist)]
    return cifgen.make_cif(mol, cell, ["x, y, z", "-x, -y, -z"], "phenylnaphthalene", "P -1", origin_frac=(0.25, 0.25, 0.25))


@pytest.fixture(scope="module")
def crystal():
    return stacked_crystal()


# ---- fragments -------------------------------------------------------------------------------------------
def test_fragments_of_a_molecule(crystal):
    a = analyse_text(crystal, "fn", Limits())
    m = a.structure.molecules[0]
    assert [len(f.atoms) for f in a.fragments[0]] == [10, 6]
    assert pi_core(m, a.limits).atoms.tolist() == a.fragments[0][0].atoms.tolist()  # the main one is the usual core
    assert len(pi_fragments(m, a.limits.with_values(fragments="largest"))) == 1
    assert len(pi_fragments(m, a.limits.with_values(min_fragment_atoms=7))) == 1
    g = a.molecule_groups()[0]
    assert g.fragments == 2 and g.with_core == 2


def test_fragment_fragment_pairs(crystal):
    full = compute_all(analyse_text(crystal, "fn", Limits()))
    largest = analyse_text(crystal, "fn", Limits(fragments="largest"))
    assert all((p.fa, p.fb) == (0, 0) for p in largest.pairs)
    combos = {(p.fa, p.fb) for p in full.pairs if not p.intramolecular}
    assert {(0, 0), (1, 1)} <= combos
    # the naphthalene–naphthalene stack along a is there in both modes, with the same values
    n_full = next(p for p in full.pairs if (p.fa, p.fb) == (0, 0) and p.kind == "stack")
    n_big = next(p for p in largest.pairs if p.kind == "stack")
    assert n_full.interplanar == pytest.approx(n_big.interplanar) == pytest.approx(3.5, abs=1e-4)
    # phenyl–phenyl along a: parallel planes (same molecule translated)
    ph = next(p for p in full.pairs if (p.fa, p.fb) == (1, 1) and p.dihedral < 1)
    assert ph.get("fragment_a") == 2 and ph.get("fragment_b") == 2 and "π2–π2" in relation(ph)
    for p in full.pairs:
        assert class_spread(full.structure, full.fragments, full.limits, p) < 1e-6


def test_intramolecular_pairs(crystal):
    a = analyse_text(crystal, "fn", Limits())
    intra = [p for p in a.pairs if p.intramolecular]
    assert len(intra) == 1  # π1–π2 in the same molecule; the two molecules of the cell are equivalent
    p = intra[0]
    assert (p.fa, p.fb) == (0, 1) and p.count == 2 and p.kind == "tilted"
    assert p.dihedral == pytest.approx(50.0, abs=1e-3) and relation(p).startswith("intramolecolare · π1–π2")


def test_summary_per_molecule(crystal):
    a = analyse_text(crystal, "fn", Limits())
    contacts = a.molecule_contacts()
    assert sum(len(c.pairs) for c in contacts) == len(a.pairs)
    along_a = next(c for c in contacts if not c.intramolecular and c.a == c.b and tuple(abs(x) for x in c.shift) == (1, 0, 0))
    assert len(along_a.pairs) >= 2 and along_a.count == 2  # one per molecule of the cell (Z = 2)
    assert along_a.best.kind == "stack" and along_a.best.fa == along_a.best.fb == 0
    intra = [c for c in contacts if c.intramolecular]
    assert len(intra) == 1 and intra[0].count == 2


def test_json_with_fragments_and_contacts(crystal, tmp_path):
    from pistak.io import export

    a = compute_all(analyse_text(crystal, "fn", Limits()))
    data = json.loads(export.write_json(a, tmp_path / "a.json").read_text())
    assert data["version"] == 5
    assert [(f["molecule"], f["fragment"], f["atoms"]) for f in data["fragments"]] == [(1, 1, 10), (1, 2, 6), (2, 1, 10), (2, 2, 6)]
    assert {p["fragment_a"] for p in data["pairs"]} == {1, 2}
    assert sum(len(c["pairs"]) for c in data["molecule_contacts"]) == len(data["pairs"])


def test_single_fragment_examples_unchanged():
    for cif in sorted(p for p in EX.glob("*.cif") if not p.name.startswith(".")):
        full, largest = analyse_file(cif, Limits()), analyse_file(cif, Limits(fragments="largest"))
        assert [(p.kind, round(p.interplanar, 9), p.count) for p in full.pairs] == \
               [(p.kind, round(p.interplanar, 9), p.count) for p in largest.pairs]


# ---- user-defined long axis and signed slips --------------------------------------------------------------
def _short_axis_labels(a) -> tuple[str, str]:
    """Two core atoms aligned along the short axis (same L coordinate, opposite W)."""
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


def test_long_axis_from_atoms():
    cif = EX / "herringbone_anthracene.cif"
    auto = compute_all(analyse_file(cif, Limits()))
    i, j = _short_axis_labels(auto)
    turned = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{i},{j}")))
    assert all(c.axis_source == "atoms" for c in turned.cores)
    p0, p1 = auto.pairs[0], turned.pairs[0]
    assert p0.kind == p1.kind == "stack" and p1.interplanar == pytest.approx(p0.interplanar)
    # the «long» axis is now the short one: the two slips swap
    assert p1.long_slip == pytest.approx(p0.lat_slip, abs=1e-4) and p1.lat_slip == pytest.approx(p0.long_slip, abs=1e-4)
    # labels that do not exist: the automatic axis stays
    assert all(c.axis_source == "pca" for c in analyse_file(cif, Limits(long_axis_atoms="X98,X99")).cores)


def test_swapping_axis_atoms_flips_sign():
    # molecule off inversion centres (each label once per molecule): the axis direction is defined
    cif = EX / "dimers_naphthalene.cif"
    auto = analyse_file(cif, Limits())
    i, j = _short_axis_labels(auto)
    fwd = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{i},{j}")))
    rev = compute_all(analyse_file(cif, Limits(long_axis_atoms=f"{j},{i}")))
    for p, q in zip(fwd.pairs, rev.pairs):
        assert q.get("long_slip_signed") == pytest.approx(-p.get("long_slip_signed"), abs=1e-9)
        assert q.long_slip == pytest.approx(p.long_slip, abs=1e-9)


def test_signed_slips():
    a = compute_all(analyse_file(EX / "dimers_naphthalene.cif", Limits()))
    for p in a.pairs:
        sl, sw = p.get("long_slip_signed"), p.get("lat_slip_signed")
        assert abs(sl) == pytest.approx(p.long_slip, abs=1e-9) and abs(sw) == pytest.approx(p.lat_slip, abs=1e-9)
        assert p.get("slip_direction") == pytest.approx(math.degrees(math.atan2(sw, sl)), abs=1e-9)
        assert 0 <= p.get("axis_aspect_a") <= 1


def test_core_isotropy():
    naph = analyse_file(EX / "demo_naphthalene.cif", Limits()).cores[0].plane.aspect
    anth = analyse_file(EX / "herringbone_anthracene.cif", Limits()).cores[0].plane.aspect
    assert 0 < anth < naph < 1  # anthracene is more elongated than naphthalene


def test_axis_from_thresholds():
    assert Limits.from_dict({"long_axis_atoms": " C1,C9 "}).long_axis_atoms == "C1,C9"
    assert Limits.from_dict({"long_axis_atoms": ["C1"]}).long_axis_atoms == ""


def test_dimers_xyz_collect_the_fragment_pairs_of_a_dimer(crystal, tmp_path):
    """Con due frammenti per molecola, le coppie π1–π1, π2–π2, π1–π2 tra le stesse due molecole sono un dimero solo:
    un file XYZ per dimero, non uno per coppia di frammenti. La coppia intramolecolare π1–π2 non è un dimero: resta
    fuori dall'esportazione dei dimeri e nell'XYZ della singola coppia la molecola c'è una volta sola."""
    from pistak.io import export

    a = analyse_text(crystal, "fn", Limits())
    inter = [p for p in a.pairs if not export.same_molecule(p)]
    ds = export.dimers(a)
    assert len(ds) < len(inter)  # più coppie di frammenti per lo stesso dimero
    assert any(len({(p.fa, p.fb) for p in d.pairs}) > 1 for d in ds)
    assert sorted(p.index for d in ds for p in d.pairs) == sorted(p.index for p in inter)
    # tutte le distanze tra gli atomi del dimero, ordinate: uguali a meno di traslazioni, rotazioni, scambio A/B
    def shape(x):
        return np.sort(np.linalg.norm(x[:, None] - x[None], axis=2).ravel())

    def same(u, v):
        return len(u) == len(v) and np.abs(u - v).max() < 1e-4

    def from_xyz(text):
        return shape(np.array([[float(v) for v in line.split()[1:4]] for line in text.splitlines()[2:]]))

    shapes = [from_xyz(export.dimer_xyz(a, d)) for d in ds]
    assert all(not same(shapes[i], shapes[j]) for i in range(len(ds)) for j in range(i))  # dimeri distinti
    for d, s in zip(ds, shapes):  # ogni coppia del dimero dà le stesse due molecole
        for p in d.pairs:
            ctx = a.pair_context(p)
            assert same(shape(np.vstack([ctx.pos_a, ctx.pos_b])), s)
    intra = next(p for p in a.pairs if export.same_molecule(p))
    lines = export.pair_xyz(a, intra).splitlines()
    assert int(lines[0]) == len(a.structure.molecules[intra.a].el) and "same_molecule" in lines[1]
    paths = export.write_dimers_xyz(a, tmp_path)
    assert len(paths) == len(ds) and all(p.name.startswith("fn-") for p in paths)


def test_dimer_column_and_scatter_rings(crystal):
    """La colonna «Dimero» della tabella e gli anelli dello scatter vengono da `export.dimers`, come i file XYZ:
    stessa etichetta per le coppie dello stesso dimero, trattino per due frammenti della stessa molecola; ordinando
    per la colonna, le coppie di un dimero stanno una dopo l'altra."""
    pytest.importorskip("PySide6")
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QApplication

    from pistak.gui.session import Session
    from pistak.gui.views.scatter_view import ScatterView
    from pistak.gui.views.table_view import FIXED, TableView
    from pistak.io import export

    _app = QApplication.instance() or QApplication([])
    s = Session()
    table, scatter = TableView(s, "A"), ScatterView(s, "A")
    a = compute_all(analyse_text(crystal, "fn", Limits()))  # come nella finestra: i parametri, per gli assi dello scatter
    s._apply(a)
    col = FIXED.index("Dimero")
    m = table.model
    labels = {m.pairs[r].index: m.data(m.index(r, col)) for r in range(m.rowCount())}
    ds = export.dimers(a)
    for d in ds:
        assert {labels[p.index] for p in d.pairs} == {d.label}
    assert all(labels[p.index] == "–" for p in a.pairs if export.same_molecule(p))
    assert m.headerData(col, Qt.Orientation.Horizontal) in ("Dimero", "Dimer")

    table.table.sortByColumn(col, Qt.SortOrder.AscendingOrder)
    shown = [table.proxy.mapToSource(table.proxy.index(r, 0)).row() for r in range(table.proxy.rowCount())]
    shown_labels = [labels[m.pairs[r].index] for r in shown]
    assert shown_labels == [d.label for d in ds for _p in d.pairs] + ["–"] * shown_labels.count("–")

    big = max(ds, key=lambda d: len(d.pairs))
    assert len(big.pairs) > 1
    scatter.choose(big.pairs[0].index)
    assert len(scatter.siblings.data) == len(big.pairs) - 1
    assert big.label in scatter._tip(0.0, 0.0, big.pairs[1].index)


def test_xyz_export_choice_warns_about_the_per_pair_files(crystal, tmp_path):
    """Un file per coppia della tabella: tanti file quante coppie, `dimer=` nella riga di commento, la molecola da sola
    per due frammenti della stessa molecola; i conteggi dell'avviso tornano con i file, e la finestra di scelta li
    mostra (avvisi solo dove servono)."""
    from pistak.io import export
    from pistak.labels import xyz_pairs_warnings

    a = analyse_text(crystal, "fn", Limits())
    c = export.xyz_counts(a)
    ds = export.dimers(a)
    assert (c.pairs, c.dimers) == (len(a.pairs), len(ds))
    assert c.same_molecule == sum(export.same_molecule(p) for p in a.pairs) > 0
    assert c.pairs == c.dimers + c.repeated + c.same_molecule and c.repeated > 0
    assert c.example is not None and len(c.example.pairs) == max(len(d.pairs) for d in ds)
    assert len(xyz_pairs_warnings(c)) == 2

    paths = export.write_pairs_xyz(a, tmp_path / "pairs")
    assert len(paths) == len(a.pairs) == len({p.name for p in paths})
    dimer_of = export.dimer_of_pairs(a)
    for p, path in zip(a.pairs, paths):
        head = path.read_text(encoding="utf-8").splitlines()[1]
        assert f"pair={p.index + 1} " in head
        if export.same_molecule(p):
            assert "same_molecule" in head and "dimer=" not in head
        else:
            assert head.endswith(f" dimer={dimer_of[p.index][1].label}")

    # una struttura senza frammenti multipli: nessun avviso, i due insiemi coincidono
    simple = analyse_file(EX / "herringbone_anthracene.cif")
    assert xyz_pairs_warnings(export.xyz_counts(simple)) == []

    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QLabel

    from pistak.gui.xyz_dialog import DIMERS, PAIRS, XyzExportDialog

    _app = QApplication.instance() or QApplication([])
    dlg = XyzExportDialog("fn.cif", c)
    assert dlg.choice() == DIMERS
    texts = [w.text() for w in dlg.findChildren(QLabel)]
    assert sum(t.startswith("⚠") for t in texts) == 2
    dlg.pairs.setChecked(True)
    assert dlg.choice() == PAIRS
    assert XyzExportDialog("fn.cif", c, PAIRS).choice() == PAIRS


def test_xyz_pairs_cli_writes_one_file_per_pair_and_warns(crystal, tmp_path, capsys):
    from pistak.cli import main

    cif = tmp_path / "fn.cif"
    cif.write_text(crystal, encoding="utf-8")
    assert main(["analyze", str(cif), "--xyz-pairs", str(tmp_path / "out")]) == 0
    a = analyse_file(cif)
    assert len(list((tmp_path / "out").glob("*.xyz"))) == len(a.pairs)
    err = capsys.readouterr().err
    assert err.count("\n") == 2  # i due avvisi, sullo standard error
