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
from pistak.api import analyse_text, finish, reanalyse
from pistak.core.cif import Cell
from pistak.core.config import Limits
from pistak.core.cores import small_rings
from pistak.core.pairs import MAX_SHELLS

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
