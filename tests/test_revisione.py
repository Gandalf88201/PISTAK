"""Difetti trovati dalla revisione del codice (validation/REGISTRO.md), ciascuno con un oracolo che non usa pistak.core:
topologia nota dei policicli, normali ricalcolate con numpy, limite teorico delle celle da esplorare."""
import itertools
import math
import random
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import cifgen
from pistak import api
from pistak.api import analyse_text, finish, reanalyse
from pistak.core.cif import Cell
from pistak.core.config import Limits
from pistak.core.cores import small_rings
from pistak.core.pairs import MAX_SHELLS

EX = Path(__file__).resolve().parents[1] / "examples"

ROOT = Path(__file__).resolve().parents[1]


# ---- tokenize: spaziatura diversa da spazio e tab ----------------------------------------------------------------
@pytest.mark.parametrize("space", ["\xa0", "\x0c", "\x0b", " ", "\x1c"])
def test_tokenize_non_si_blocca_sugli_spazi_unicode(space):
    # in un processo a parte: il difetto era un ciclo infinito che riempie la memoria
    code = f"from pistak.core.cif import tokenize; print(tokenize('_cell_length_a{space}5.0 # nota'))"
    out = subprocess.run([sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True, timeout=30)
    assert out.stdout.strip() == "['_cell_length_a', '5.0']"


# ---- small_rings: anelli con tutti i legami condivisi -------------------------------------------------------------
def honeycomb(centres):
    """Atomi, legami ed esagoni attesi di un benzenoide dato dai centri dei suoi anelli (coordinate assiali)."""
    cc, h = 1.4, 1.4 * math.sqrt(3) / 2
    atoms, rings = [], []
    for q, r in centres:
        c = np.array([q * 2 * h + r * h, r * 1.5 * cc])
        ring = []
        for k in range(6):
            p = c + cc * np.array([math.cos(math.pi / 6 + k * math.pi / 3), math.sin(math.pi / 6 + k * math.pi / 3)])
            i = next((i for i, a in enumerate(atoms) if np.linalg.norm(a - p) < 1e-6), None)
            if i is None:
                atoms.append(p)
                i = len(atoms) - 1
            ring.append(i)
        rings.append(frozenset(ring))
    bonds = [(i, j) for i, j in itertools.combinations(range(len(atoms)), 2) if np.linalg.norm(atoms[i] - atoms[j]) < 1.45]
    return len(atoms), bonds, rings


CORONENE = [(0, 0), (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
HEXABENZOCORONENE = CORONENE + [(2, -1), (1, 1), (-1, 2), (-2, 1), (-1, -1), (1, -2)]


@pytest.mark.parametrize("centres", [CORONENE, HEXABENZOCORONENE], ids=["coronene", "esabenzocoronene"])
def test_small_rings_trova_gli_anelli_interni_in_ogni_ordine(centres):
    n, bonds, expected = honeycomb(centres)
    assert len(bonds) - n + 1 == len(expected)  # rango ciclico = numero di esagoni
    rng = random.Random(0)
    for _ in range(300):
        perm = list(range(n))
        rng.shuffle(perm)
        inv = {old: new for new, old in enumerate(perm)}
        b = np.array(sorted(tuple(sorted((inv[i], inv[j]))) for i, j in bonds))
        found = {frozenset(perm[i] for i in r) for r in small_rings(n, b, lambda i: True, 8)}
        assert found == set(expected)


# ---- mean_atom_height: distanza media dal piano di A --------------------------------------------------------------
def test_quota_media_dal_piano_di_a():
    an = analyse_text((ROOT / "examples/tilted_anthracene.cif").read_text(encoding="utf-8"), "tilted", Limits())
    finish(an)
    tilted = 0
    for p in an.pairs:
        ctx = an.pair_context(p)
        pa, pb = ctx.pos_a[ctx.core_a.atoms], ctx.pos_b[ctx.core_b.atoms]
        ca = pa.mean(axis=0)
        na = np.linalg.svd(pa - ca)[2][-1]
        assert p.get("mean_atom_height") == pytest.approx(np.abs((pb - ca) @ na).mean(), abs=1e-9)
        if p.dihedral > 30:  # non è più la distanza interplanare sotto altro nome
            tilted += 1
            assert abs(p.get("mean_atom_height") - p.interplanar) > 0.1
    assert tilted


# ---- MAX_SHELLS: il taglio della ricerca si segnala ---------------------------------------------------------------
def long_acene_column(a: float = 3.5) -> str:
    """Undecacene (estensione ~14 Å) impilato lungo un asse a corto: servono più di MAX_SHELLS celle."""
    mol = [(el, cifgen.rotation("y", 90) @ p) for el, p in cifgen.ring_system(11)]
    return cifgen.make_cif(mol, Cell(a, 8.0, 34.0, 90, 90, 90), ["x, y, z"], "undecacene", "P 1", origin_frac=(0.5, 0.5, 0.5))


def _truncated(an):
    return [x for x in an.alerts if x.code == "shellsTruncated"]


def test_celle_tagliate_danno_un_avviso():
    an = finish(analyse_text(long_acene_column(), "undecacene", Limits()))
    extent = an.structure.molecules[0].extent
    need = lambda cutoff: math.floor((2 * extent + cutoff) / 3.5 + 1)  # limite teorico lungo a (cella ortogonale)
    assert need(4.0) > MAX_SHELLS
    alerts = _truncated(an)
    assert len(alerts) == 1 and alerts[0].level == "B" and alerts[0].values["needed"].startswith(f"±{need(4.0)},")
    assert any(p.kind == "stack" for p in an.pairs)
    # con un altro contatto massimo servono altre celle: l'avviso riusato da reanalyse si ricalcola
    smaller = reanalyse(an, an.limits.with_values(contact_cutoff=0.5))
    assert need(0.5) < need(4.0) and _truncated(smaller)[0].values["needed"].startswith(f"±{need(0.5)},")


def test_nessun_avviso_sulle_celle_normali():
    an = finish(analyse_text(long_acene_column(a=6.0), "undecacene", Limits()))
    assert not _truncated(an)


# ---- contact_cutoff: ogni stack entro max_interplanar viene misurato ---------------------------------------------
def naphthalene_column(interplanar: float, slip: float) -> str:
    """Colonna infinita di naftalene a piani esattamente paralleli: quota e slittamento imposti dalla cella monoclina."""
    mol = cifgen.ring_system(2)
    c = math.hypot(interplanar, slip)
    beta = math.degrees(math.atan2(interplanar, slip))
    return cifgen.make_cif(mol, Cell(20.0, 20.0, c, 90, beta, 90), ["x, y, z"], "colonna", "P 1", origin_frac=(0.5, 0.5, 0.0))


@pytest.mark.parametrize("slip", [0.0, 1.4])
@pytest.mark.parametrize("d", [3.5, 4.05, 4.1, 4.19])
def test_stack_tra_4_0_e_4_2_angstrom(d, slip):
    an = analyse_text(naphthalene_column(d, slip), "colonna", Limits())
    stacks = [p for p in an.pairs if p.kind == "stack"]
    assert len(stacks) == 1 and stacks[0].interplanar == pytest.approx(d, abs=1e-4)


def test_oltre_max_interplanar_la_coppia_e_parallela():
    an = analyse_text(naphthalene_column(4.3, 0.0), "colonna", Limits())
    assert [p.kind for p in an.pairs] == ["parallel"]


def test_il_modo_del_sito_tiene_4_angstrom():
    assert Limits.site().contact_cutoff == 4.0 and Limits.site(contact_cutoff=5.0).contact_cutoff == 5.0
    an = analyse_text(naphthalene_column(4.1, 0.0), "colonna", Limits.site())
    assert not an.pairs


# ---- beta_a / beta_b → beta / gamma -------------------------------------------------------------------------------
def test_nomi_vecchi_di_beta_e_gamma():
    from pistak.params import REGISTRY

    an = finish(analyse_text((ROOT / "examples/tilted_anthracene.cif").read_text(encoding="utf-8"), "tilted", Limits()))
    p = an.pairs[0]
    assert "beta" in p.values and "gamma" in p.values and "beta_a" not in p.values
    assert p.get("beta_a") == p.get("beta") and p.get("beta_b") == p.get("gamma")
    assert REGISTRY.get("beta_b").key == "gamma" and REGISTRY.get("gamma").label.startswith("γ")


# ---- distanze X–H secondo l'ibridazione ---------------------------------------------------------------------------
def test_idrogeni_secondo_ibridazione():
    from pistak.io.qc import normalize_hydrogens

    # C0 sp³ (CH3–), C4 sp² (=CH–), N7 con quattro vicini (N⁺H3), P11 con un H: X–H tutti a 0,95 Å lungo assi diversi
    el = ["C", "H", "H", "H", "C", "H", "C", "N", "H", "H", "H", "P", "H"]
    bonds = np.array([(0, 1), (0, 2), (0, 3), (0, 4), (4, 5), (4, 6), (6, 7), (7, 8), (7, 9), (7, 10), (11, 12)])
    pos = np.zeros((len(el), 3))
    rng = np.random.default_rng(3)
    for x, h in bonds:
        if el[h] == "H":
            v = rng.normal(size=3)
            pos[h] = pos[x] + 0.95 * v / np.linalg.norm(v)
    out, moved = normalize_hydrogens(el, pos, bonds)
    length = lambda x, h: float(np.linalg.norm(out[h] - out[x]))
    # valori di Allen & Bruno, Acta Cryst. B66, 380 (2010), medie a T ≤ 60 K
    assert [length(0, h) for h in (1, 2, 3)] == pytest.approx([1.095] * 3, abs=1e-9)
    assert length(4, 5) == pytest.approx(1.085, abs=1e-9)
    assert [length(7, h) for h in (8, 9, 10)] == pytest.approx([1.040] * 3, abs=1e-9)
    assert length(11, 12) == pytest.approx(0.95, abs=1e-9) and moved == 7  # P–H senza fonte: resta com'è


# ---- simboli H-M con l'impostazione esplicita (:1/:2, :H/:R) -------------------------------------------------------
def test_suffissi_dell_impostazione():
    from pistak.core.spacegroups import operations

    rhombo, source = operations(hm="R -3 :R")
    assert source == "table" and len(rhombo) == 6 and "z,x,y" in rhombo  # assi romboedrici: niente centratura R
    assert len(operations(hm="R -3")[0]) == len(operations(hm="R -3 :H")[0]) == 18
    # origine 2 (predefinita): centro di inversione nell'origine; origine 1: no
    assert "-x,-y,-z" in operations(hm="F d d d")[0] and "-x,-y,-z" in operations(hm="F d d d :2")[0]
    assert "-x,-y,-z" not in operations(hm="F d d d :1")[0]
    assert "-x,-y,-z" not in operations(hm="I 41/a :1")[0] and "-x,-y,-z" in operations(hm="I 41/a")[0]


# ---- element_of: etichette ambigue e siti fittizi -----------------------------------------------------------------
@pytest.mark.parametrize("symbol, label, expected", [
    (None, "Cl1", "Cl"), (None, "CL1", "Cl"), (None, "BR2", "Br"), (None, "CA1", "C"), (None, "CD1", "C"),
    (None, "Co1", "Co"), (None, "D1", "H"), (None, "Cg1", ""), (None, "Q1", ""), ("Q", "Q1", ""),
    ("Fe2+", "Fe1", "Fe"), ("C", "Co1", "C"), ("Cl", "CL1", "Cl"),
])
def test_elemento_dall_etichetta(symbol, label, expected):
    from pistak.core.cif import element_of

    assert element_of(symbol, label) == expected


def test_siti_fittizi_scartati():
    from pistak.core.cif import parse_cif

    text = "\n".join([
        "data_t", "_cell_length_a 10", "_cell_length_b 10", "_cell_length_c 10", "_cell_angle_alpha 90",
        "_cell_angle_beta 90", "_cell_angle_gamma 90", "loop_", "_symmetry_equiv_pos_as_xyz", "'x, y, z'",
        "loop_", "_atom_site_label", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z",
        "_atom_site_calc_flag", "C1 0.1 0.1 0.1 d", "CL1 0.2 0.1 0.1 .", "X1 0.3 0.1 0.1 dum", "Cg1 0.4 0.1 0.1 .",
        "Q1 0.5 0.1 0.1 .", ""])
    assert [(a.label, a.el) for a in parse_cif(text).atoms] == [("C1", "C"), ("CL1", "Cl")]


# ---- operazioni elencate contro il simbolo del gruppo --------------------------------------------------------------
def test_operazioni_incoerenti_con_il_simbolo():
    ok = cifgen.make_cif(cifgen.ring_system(2), Cell(8, 6, 9, 90, 100, 90), cifgen.P21C, "ok", "P 21/c", origin_frac=(0.3, 0.3, 0.3))
    wrong = cifgen.make_cif(cifgen.ring_system(2), Cell(8, 6, 9, 90, 100, 90), ["x, y, z", "-x, -y, -z"], "no", "P 21/c",
                            origin_frac=(0.3, 0.3, 0.3))
    assert not [x for x in finish(analyse_text(ok, "ok", Limits())).alerts if x.code == "symopsMismatch"]
    (alert,) = [x for x in finish(analyse_text(wrong, "no", Limits())).alerts if x.code == "symopsMismatch"]
    assert alert.level == "C" and alert.values == {"declared": "P 21/c", "n_cif": 2, "n_symbol": 4}


# ---- contenuto della cella: atomi ripetuti dal CIF e posizioni speciali --------------------------------------------
def test_contenuto_della_cella_con_atomi_ripetuti():
    from pistak.core.cif import parse_cif
    from pistak.core.crystal import cell_contents

    text = "\n".join([
        "data_t", "_cell_length_a 10", "_cell_length_b 10", "_cell_length_c 10", "_cell_angle_alpha 90",
        "_cell_angle_beta 90", "_cell_angle_gamma 90", "loop_", "_symmetry_equiv_pos_as_xyz", "'x, y, z'", "'-x, -y, -z'",
        "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
        "_atom_site_fract_z", "_atom_site_occupancy",
        "C1 C 0.1 0.2 0.3 0.5", "C1i C -0.1 -0.2 -0.3 0.5",  # immagine scritta di nuovo, metà occupanza: in tutto 2 C
        "N1 N 0.3 0.3 0.1 1", "N1i N -0.3 -0.3 -0.1 1",       # immagine scritta di nuovo, occupanza piena: 2 N, non 4
        "O1 O 0 0 0 0.5",                                      # sul centro di inversione, mezzo occupato: 0,5 O
        ""])
    assert cell_contents(parse_cif(text, disorder="all")) == pytest.approx({"C": 2.0, "N": 2.0, "O": 0.5})


def _cif_p1(atoms: list[tuple[str, str, tuple[float, float, float]]], a: float = 15.0) -> str:
    rows = "\n".join(f"{lab} {el} {x / a:.6f} {y / a:.6f} {z / a:.6f}" for lab, el, (x, y, z) in atoms)
    return (f"data_x\n_cell_length_a {a}\n_cell_length_b {a}\n_cell_length_c {a}\n_cell_angle_alpha 90\n"
            f"_cell_angle_beta 90\n_cell_angle_gamma 90\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\nloop_\n"
            f"_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n{rows}\n")


@pytest.mark.parametrize("name, atoms, expected", [
    # C terminale con triplo legame a N o O: nessun idrogeno mancante
    ("isonitrile", [("C1", "C", (7, 7, 7)), ("N1", "N", (8.42, 7, 7)), ("C2", "C", (9.58, 7, 7)),
                    ("H1", "H", (6.64, 8.03, 7)), ("H2", "H", (6.64, 6.49, 7.89)), ("H3", "H", (6.64, 6.49, 6.11))], set()),
    ("cianuro", [("C1", "C", (7, 7, 7)), ("N1", "N", (8.15, 7, 7))], set()),
    ("carbonile", [("Cr1", "Cr", (7, 7, 7)), ("C1", "C", (8.91, 7, 7)), ("O1", "O", (10.05, 7, 7))], set()),
    # una catena di C senza idrogeni: l'avviso resta
    ("senza H", [("C1", "C", (7, 7, 7)), ("C2", "C", (8.5, 7, 7)), ("C3", "C", (9.25, 8.3, 7))], {"noHydrogens"}),
])
def test_idrogeni_mancanti_senza_falsi_positivi(name, atoms, expected):
    from pistak.core.analysis import analyse_text
    from pistak.core.validation import validate

    a = analyse_text(_cif_p1(atoms), name)
    assert {x.code for x in validate(a)} & {"noHydrogens", "missingH"} == expected


def _benzene_p1(occ: float) -> str:
    """Un benzene e la sua immagine per inversione, scritti in P1: la simmetria P-1 è mancata."""
    import math

    rows = []
    for k in range(6):
        t = 2 * math.pi * k / 6
        x, y, z = 3 + 1.39 * math.cos(t), 3 + 1.39 * math.sin(t), 3.0
        rows.append((f"C{k + 1}", (x / 12, y / 12, z / 12)))
        rows.append((f"C{k + 7}", (1 - x / 12, 1 - y / 12, 1 - z / 12)))
    atoms = "\n".join(f"{lab} C {f[0]:.6f} {f[1]:.6f} {f[2]:.6f} {occ}" for lab, f in rows)
    return ("data_x\n_cell_length_a 12\n_cell_length_b 12\n_cell_length_c 12\n_cell_angle_alpha 90\n_cell_angle_beta 90\n"
            "_cell_angle_gamma 90\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\nloop_\n_atom_site_label\n"
            "_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n"
            f"{atoms}\n")


def test_addsym_lascia_fuori_gli_atomi_disordinati():
    pytest.importorskip("spglib")
    from pistak.core.analysis import analyse_text
    from pistak.core.validation import validate

    full = {x.code: x for x in validate(analyse_text(_benzene_p1(1.0), "ordinato"))}
    assert full["missedSymmetry"].level == "A"  # spglib trova Immm: il benzene scritto così ha anche specchi
    half = {x.code for x in validate(analyse_text(_benzene_p1(0.5), "disordinato"))}
    assert "addsymDisorder" in half and "missedSymmetry" not in half


def test_valori_mancanti_nel_cif():
    # densità e Z sconosciuti («?»): nessun errore, la formula non si controlla
    from pistak.core.analysis import analyse_text
    from pistak.core.validation import validate

    text = _benzene_p1(1.0).replace("_cell_angle_gamma 90\n", "_cell_angle_gamma 90\n_exptl_crystal_density_diffrn ?\n"
                                    "_cell_formula_units_Z ?\n_chemical_formula_sum 'C6 H6'\n")
    codes = {x.code for x in validate(analyse_text(text, "mancanti"))}
    assert "formulaUnchecked" in codes and "density" not in codes and "checkFailed" not in codes


# ---- sezione 9: comportamenti silenziosi -------------------------------------------------------------------------
def test_spglib_che_fallisce_non_si_scambia_per_spglib_mancante(monkeypatch):
    """Se spglib c'è ma non riesce, l'avviso è checkFailed con il messaggio, non «serve il pacchetto spglib»."""
    import spglib

    from pistak.core.validation import validate

    def broken(*_a, **_k):
        return None  # come spglib < 2.5 quando fallisce

    monkeypatch.setattr(spglib, "get_symmetry_dataset", broken)
    a = api.analyse_file(EX / "demo_naphthalene.cif")
    codes = {x.code: x for x in validate(a)}
    assert "addsymUnavailable" not in codes
    assert codes["checkFailed"].values["check"] == "addsym" and "spglib" in codes["checkFailed"].values["error"]


def test_monte_carlo_conta_i_campioni_mancanti():
    """Un parametro che manca in metà dei campioni: la σ viene dagli altri e `incomplete` lo dice."""
    import math

    from pistak.core.uncertainty import geometry_values, monte_carlo

    a = api.analyse_file(EX / "demo_naphthalene.cif")
    calls = {"n": 0}

    def flaky(analysis, pair):
        v = geometry_values(analysis, pair)
        calls["n"] += 1
        if calls["n"] > len(a.pairs) and calls["n"] % 2:  # i nominali sempre, poi un campione su due manca
            v["interplanar"] = math.nan
        return v

    u = monte_carlo(a, samples=20, seed=1, evaluate=flaky)
    i = a.pairs[0].index
    assert 0 < u.incomplete[i]["interplanar"] < 20 == u.valid[i]["slip"]
    assert "slip" not in u.incomplete.get(i, {})


def test_raggio_vdw_non_tabulato_si_segnala():
    """Fe non ha un raggio di van der Waals nelle tabelle (Bondi, Mantina): l'avviso G dice che si usa il valore medio."""
    from pistak.core.elements import DEFAULT_VDW, VDW_RADII
    from pistak.core.validation import validate

    assert "Fe" not in VDW_RADII and VDW_RADII["Ra"] == 2.83  # Ra da Mantina 2009, tab. 12
    cif = cifgen.make_cif([("Fe", np.array([1.0, 1.0, 1.0])), ("O", np.array([1.0, 1.0, 2.8]))],
                          Cell(10, 10, 10, 90, 90, 90), ["x, y, z"], "ferro")
    alerts = {x.code: x for x in validate(analyse_text(cif, "ferro"))}
    assert alerts["vdwDefault"].values == {"elements": "Fe", "r": DEFAULT_VDW}
    assert "vdwDefault" not in {x.code for x in validate(api.analyse_file(EX / "demo_naphthalene.cif"))}
