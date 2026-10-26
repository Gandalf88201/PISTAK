"""Fase 2: coppie uniche ed equivalenza per simmetria, celle adattive, più blocchi data_, disordine, operazioni dal
simbolo del gruppo, criterio di stack configurabile."""
import math
from pathlib import Path

import numpy as np
import pytest

import cifgen
from pistak.core.analysis import analyse_file, analyse_text
from pistak.core.cif import Cell, block_names, parse_cif
from pistak.core.config import Limits
from pistak.core.errors import CifError
from pistak.core.pairs import class_spread, measure_pair, neighbour_shells
from pistak.core.spacegroups import HM_TO_HALL, hall_operations, operations, triplet

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
EXAMPLES = sorted(p for p in EX.glob("*.cif") if not p.name.startswith("."))


# ---- coppie uniche ed equivalenza per simmetria -------------------------------------------------------
@pytest.mark.parametrize("cif", EXAMPLES, ids=[p.stem for p in EXAMPLES])
def test_ogni_contatto_una_volta(cif):
    site, canon = analyse_file(cif, Limits.site()), analyse_file(cif, Limits())
    # il sito conta ogni contatto due volte (a→b e b→a)
    assert sum(p.count for p in site.pairs) == 2 * sum(p.count for p in canon.pairs)
    seen = set()
    for p in canon.pairs:
        for a, b, sh in p.members:
            key = min((a, b, sh), (b, a, tuple(-x for x in sh)))
            assert key not in seen
            seen.add(key)


@pytest.mark.parametrize("cif", EXAMPLES, ids=[p.stem for p in EXAMPLES])
def test_le_coppie_equivalenti_hanno_gli_stessi_valori(cif):
    a = analyse_file(cif, Limits())
    for p in a.pairs:
        assert class_spread(a.structure, a.cores, a.limits, p) < 1e-6


@pytest.mark.parametrize("cif", EXAMPLES, ids=[p.stem for p in EXAMPLES])
def test_scambiare_a_e_b_non_cambia_i_valori(cif):
    a = analyse_file(cif, Limits())
    for p in a.pairs:
        q = measure_pair(a.structure, a.cores, p.b, p.a, tuple(-x for x in p.shift), a.limits)
        for attr in ("interplanar", "long_slip", "lat_slip", "slip", "dihedral", "contact", "overlap"):
            assert getattr(q, attr) == pytest.approx(getattr(p, attr), abs=1e-9), attr
        assert q.kind == p.kind


def test_herringbone_le_due_coppie_inclinate_del_sito_sono_un_solo_contatto():
    site = analyse_file(EX / "herringbone_anthracene.cif", Limits.site())
    canon = analyse_file(EX / "herringbone_anthracene.cif", Limits())
    assert [p.kind for p in site.pairs] == ["stack", "tilted", "tilted"]
    assert [p.kind for p in canon.pairs] == ["stack", "tilted"]
    # lo stack era già simmetrico: stessi valori nei due modi
    assert canon.pairs[0].long_slip == pytest.approx(site.pairs[0].long_slip)
    t = canon.pairs[1]
    assert t.count == 4 and t.interplanar == pytest.approx(site.pairs[1].interplanar)


def test_operazioni_che_non_mappano_le_molecole_non_uniscono_niente():
    a = analyse_file(EX / "dimers_naphthalene.cif", Limits())
    a.structure._images.clear()
    a.structure.symops = a.structure.symops[:1]  # solo l'identità: nessuna equivalenza oltre a quella banale
    from pistak.core.pairs import collect_pairs

    pairs, _ = collect_pairs(a.structure, a.cores, a.limits)
    assert len(pairs) >= len(a.pairs) and all(p.count >= 1 for p in pairs)


def test_colonne_uguali_nei_due_modi():
    for cif in EXAMPLES:
        site, canon = analyse_file(cif, Limits.site()), analyse_file(cif, Limits())
        assert [(c.rank, sorted(c.molecules)) for c in site.stacks] == [(c.rank, sorted(c.molecules)) for c in canon.stacks]


# ---- celle esplorate adattive ---------------------------------------------------------------------------
def oblique_stack() -> str:
    """Cella molto obliqua: la molecola impilata sta a 3a + b (3,8 Å), oltre le ±2 celle del sito."""
    R = cifgen.rotation("z", 90) @ cifgen.rotation("y", -90)  # piano xz (normale y), asse lungo z
    mol = [(el, R @ p) for el, p in cifgen.ring_system(2)]
    a1, a2 = np.array([7.0, 0, 0]), np.array([-20.5, 3.8, 0])
    gamma = math.degrees(math.acos(a1 @ a2 / (np.linalg.norm(a1) * np.linalg.norm(a2))))
    cell = Cell(7.0, float(np.linalg.norm(a2)), 12.0, 90.0, 90.0, gamma)
    return cifgen.make_cif(mol, cell, ["x, y, z"], "obliqua", "P 1", origin_frac=(0.3, 0.4, 0.5))


def test_celle_adattive_trovano_lo_stack_lontano_in_celle():
    text = oblique_stack()
    site = analyse_text(text, "obliqua", Limits.site())
    canon = analyse_text(text, "obliqua", Limits())
    assert neighbour_shells(canon.structure, canon.cores, canon.limits)[0] >= 3
    assert not any(p.kind == "stack" for p in site.pairs)  # ±2 celle non bastano
    stack = [p for p in canon.pairs if p.kind == "stack"]
    assert len(stack) == 1 and stack[0].interplanar == pytest.approx(3.8, abs=1e-3)
    assert stack[0].slip == pytest.approx(0.5, abs=1e-3) and tuple(abs(x) for x in stack[0].shift) == (3, 1, 0)


# ---- più blocchi data_ ----------------------------------------------------------------------------------
def two_blocks() -> str:
    first = (EX / "demo_naphthalene.cif").read_text().replace("data_demo_naphthalene_stack", "data_uno")
    second = (EX / "herringbone_anthracene.cif").read_text().replace("data_herringbone", "data_due")
    return "data_global\n_journal_name_full 'Prova'\n\n" + first + "\n" + second


def test_piu_blocchi(tmp_path):
    text = two_blocks()
    assert block_names(text) == ["uno", "due"]
    assert parse_cif(text).name == "uno"  # il primo con atomi, non data_global
    assert parse_cif(text, block="due").name == "due" and parse_cif(text, block=2).name == "due"
    with pytest.raises(CifError) as e:
        parse_cif(text, block="tre")
    assert e.value.code == "noBlock"
    f = tmp_path / "due_strutture.cif"
    f.write_text(text)
    from pistak.io.batch import run_batch

    entries = run_batch([f, EX / "demo_naphthalene.cif"], Limits(), workers=1)
    assert [e.name for e in entries] == ["due_strutture.cif [uno]", "due_strutture.cif [due]", "demo_naphthalene.cif"]
    assert [e.block for e in entries] == ["uno", "due", None] and all(e.error is None for e in entries)
    assert entries[1].space_group == "P 21/c"


def test_terminale_con_blocco(tmp_path, capsys):
    from pistak.cli import main

    f = tmp_path / "due.cif"
    f.write_text(two_blocks())
    assert main(["analyze", str(f)]) == 0
    assert "2 strutture" in capsys.readouterr().err
    assert main(["analyze", str(f), "--block", "due"]) == 0
    assert "P 21/c" in capsys.readouterr().out


# ---- disordine -------------------------------------------------------------------------------------------
DISORDERED = """data_dis
_cell_length_a 10 _cell_length_b 10 _cell_length_c 10
_cell_angle_alpha 90 _cell_angle_beta 90 _cell_angle_gamma 90
loop_
_atom_site_label _atom_site_type_symbol _atom_site_fract_x _atom_site_fract_y _atom_site_fract_z
_atom_site_occupancy _atom_site_disorder_assembly _atom_site_disorder_group
C1 C 0.1 0.1 0.1 1 . .
C2A C 0.2 0.1 0.1 0.3 A 1
C2B C 0.2 0.2 0.1 0.7 A 2
O3A O 0.5 0.5 0.5 0.6 B 1
O3B O 0.5 0.6 0.5 0.4 B 2
"""


@pytest.mark.parametrize("policy,kept", [
    ("first", ["C1", "C2A", "O3A"]),  # il sito: sempre il gruppo 1
    ("major", ["C1", "C2B", "O3A"]),  # per ogni assemblaggio il gruppo più occupato
    ("all", ["C1", "C2A", "C2B", "O3A", "O3B"]),
    ("2", ["C1", "C2B", "O3B"]),
])
def test_regola_di_disordine(policy, kept):
    c = parse_cif(DISORDERED, disorder=policy)
    assert [a.label for a in c.atoms] == kept and c.skipped_disorder == 5 - len(kept) and c.disorder == policy


def test_disordine_dalle_soglie():
    a = analyse_text(DISORDERED, "d", Limits(disorder="all"))
    assert a.cif.disorder == "all" and a.structure.n_cell_atoms == 5


# ---- operazioni dal simbolo del gruppo ---------------------------------------------------------------
EXPECTED_ORDER = {"P -1": 2, "P 21/c": 4, "C 2/c": 8, "P 21 21 21": 4, "P b c a": 8, "P n a 21": 4, "F d d d": 32,
                  "I 41/a": 16, "R -3": 18, "R -3 c": 36, "P 61 2 2": 12, "P 63/m m c": 24, "P a -3": 24,
                  "I a -3 d": 96, "F d -3 m": 192}


@pytest.mark.parametrize("hm,order", EXPECTED_ORDER.items())
def test_ordine_dei_gruppi(hm, order):
    assert len(hall_operations(HM_TO_HALL[hm])) == order


@pytest.mark.parametrize("hall,ops", [
    ("-P 2ybc", {"x,y,z", "-x,y+1/2,-z+1/2", "-x,-y,-z", "x,-y+1/2,z+1/2"}),
    ("-P 2yn", {"x,y,z", "-x+1/2,y+1/2,-z+1/2", "-x,-y,-z", "x+1/2,-y+1/2,z+1/2"}),
    ("P 2ac 2ab", {"x,y,z", "-x+1/2,-y,z+1/2", "-x,y+1/2,-z+1/2", "x+1/2,-y+1/2,-z"}),
    ("P 2c -2n", {"x,y,z", "-x,-y,z+1/2", "x+1/2,-y+1/2,z", "-x+1/2,y+1/2,z+1/2"}),
    # origine spostata e assi binari diagonali (International Tables, gruppo 178)
    ("P 61 2 (0 0 -1)", {"x,y,z", "-y,x-y,z+1/3", "-x+y,-x,z+2/3", "-x,-y,z+1/2", "y,-x+y,z+5/6", "x-y,x,z+1/6",
                         "y,x,-z+1/3", "x-y,-y,-z", "-x,-x+y,-z+2/3", "-y,-x,-z+5/6", "-x+y,y,-z+1/2", "x,x-y,-z+1/6"}),
])
def test_operazioni_note(hall, ops):
    assert {triplet(o) for o in hall_operations(hall)} == ops


def test_simboli_hm_in_varie_forme():
    ref = set(operations(hm="P 21/n")[0])
    for name in ("P21/n", "P 1 21/n 1", "p 21/n"):
        assert set(operations(hm=name)[0]) == ref
    assert operations(hm="R -3 :H")[1] == "table" and operations(hm="P 42/n c m") is None


def test_cif_con_solo_il_simbolo_di_hall():
    text = (EX / "herringbone_anthracene.cif").read_text()
    start, end = text.index("loop_\n_symmetry_equiv_pos_as_xyz"), text.index("loop_\n_atom_site_label")
    no_ops = text[:start] + text[end:]
    hall = no_ops.replace("_symmetry_space_group_name_H-M 'P 21/c'", "_space_group_name_Hall '-P 2ybc'")
    for t, source in ((no_ops, "table"), (hall, "hall")):
        c = parse_cif(t)
        assert c.symops_source == source and len(c.symops) == 4
    ref = analyse_file(EX / "herringbone_anthracene.cif", Limits())
    got = analyse_text(hall, "h", Limits())
    assert [(p.kind, round(p.interplanar, 6)) for p in got.pairs] == [(p.kind, round(p.interplanar, 6)) for p in ref.pairs]


# ---- criterio di stack ---------------------------------------------------------------------------------
def test_criterio_di_stack_configurabile():
    cif = EX / "herringbone_anthracene.cif"
    base = analyse_file(cif, Limits()).pairs[0]
    assert base.kind == "stack" and 0 < base.overlap < 1
    assert analyse_file(cif, Limits(min_overlap=base.overlap + 0.01)).pairs[0].kind == "parallel"
    assert analyse_file(cif, Limits(max_slip=base.slip - 0.1)).pairs[0].kind == "parallel"
    assert analyse_file(cif, Limits(max_slip=base.slip + 0.1)).pairs[0].kind == "stack"
    area = analyse_file(cif, Limits(overlap_measure="area"))
    ratio = area.pairs[0]
    assert ratio.kind == "stack"
    from pistak.params import compute_all

    r = compute_all(area).pairs[0].get("overlap_ratio")
    assert analyse_file(cif, Limits(overlap_measure="area", min_overlap=r + 0.01)).pairs[0].kind == "parallel"


# ---- soglie ----------------------------------------------------------------------------------------------
def test_soglie_da_dizionario():
    lim = Limits.from_dict({"pair_mode": "SITE", "disorder": "3", "max_dihedral": "15", "overlap_measure": "volume",
                            "max_slip": True, "sconosciuta": 1})
    assert lim.pair_mode == "site" and lim.disorder == "3" and lim.max_dihedral == 15.0
    assert lim.overlap_measure == "atoms" and lim.max_slip == 0.0  # valori non validi ignorati
    assert Limits.site() == Limits(pair_mode="site", disorder="first")
