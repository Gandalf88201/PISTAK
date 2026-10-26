"""Analisi a stadi con cache: ogni soglia sta in uno stadio, rifare da uno stadio dà lo stesso risultato di
un'analisi da zero, il criterio di stack non rimisura nulla e i parametri geometrici si riusano."""
import math
from dataclasses import fields
from pathlib import Path

import numpy as np
import pytest

from pistak import api
from pistak.core import analysis as core
from pistak.core.config import Limits
from pistak.core.errors import InputError
from pistak.core.stages import STAGE_FIELDS, STAGES, first_stage, stage_of
from pistak.params import Registry

EX = Path(__file__).resolve().parents[1] / "examples"
#: Una modifica per stadio (e, per `pairs`, anche il cambio di criterio di sovrapposizione).
CHANGES = {
    "cif": {"disorder": "all"},
    "structure": {"bond_tolerance": 0.3},
    "fragments": {"ring_rms": 0.05, "long_axis_atoms": "C1,C5"},
    "contacts": {"contact_cutoff": 3.6},
    "pairs": {"max_interplanar": 3.6, "max_dihedral": 8.0},
    "pairs-area": {"overlap_measure": "area", "min_overlap": 0.3},
}


def summary(a) -> list:
    """Tutto ciò che conta delle coppie, con i NaN confrontabili."""
    def clean(values):
        return {k: (None if isinstance(v, float) and math.isnan(v) else round(v, 12)) for k, v in values.items()}

    return [(p.key, p.kind, p.count, p.members, clean(p.values)) for p in a.pairs]


def test_ogni_soglia_sta_in_uno_e_un_solo_stadio():
    names = [f.name for f in fields(Limits)]
    placed = [n for s in STAGES for n in STAGE_FIELDS[s]]
    assert sorted(placed) == sorted(names), "una soglia nuova va assegnata a uno stadio in core/stages.py"
    assert stage_of("max_slip") == "pairs" and stage_of("disorder") == "cif"


def test_primo_stadio_da_rifare():
    base = Limits()
    assert first_stage(base, base) is None
    assert first_stage(base, base.with_values(min_overlap=0.2)) == "pairs"
    assert first_stage(base, base.with_values(min_overlap=0.2, contact_cutoff=3.5)) == "contacts"
    assert first_stage(base, base.with_values(pair_mode="site")) == "contacts"
    assert first_stage(base, base.with_values(fragments="largest", max_slip=2)) == "fragments"


@pytest.mark.parametrize("mode", ["canonical", "site"])
@pytest.mark.parametrize("change", list(CHANGES))
@pytest.mark.parametrize("cif", sorted(p.name for p in EX.glob("*.cif") if not p.name.startswith(".")))
def test_rifare_da_uno_stadio_come_da_zero(cif, change, mode):
    base = Limits() if mode == "canonical" else Limits.site()
    new = base.with_values(**CHANGES[change])
    text = (EX / cif).read_text()
    a = api.analyse_text(text, cif, base)
    b = api.reanalyse(a, new)
    ref = api.analyse_text(text, cif, new)
    assert summary(b) == summary(ref)
    assert [(c.molecules, c.rank) for c in b.stacks] == [(c.molecules, c.rank) for c in ref.stacks]
    assert [(c.a, c.b, c.shift, c.count) for c in b.molecule_contacts()] == \
        [(c.a, c.b, c.shift, c.count) for c in ref.molecule_contacts()]
    # e tornando alle soglie di partenza si ritrova l'analisi originale, che intanto non è cambiata
    assert summary(api.reanalyse(b, base)) == summary(a) == summary(api.analyse_text(text, cif, base))


def test_il_criterio_di_stack_non_rimisura(monkeypatch):
    a = api.analyse_file(EX / "herringbone_anthracene.cif", Limits())

    def forbidden(*_a, **_k):
        raise AssertionError("lo stadio pairs non deve rimisurare le coppie")

    monkeypatch.setattr(core, "iter_neighbour_pairs", forbidden)
    monkeypatch.setattr(core, "build_structure", forbidden)
    b = api.reanalyse(a, a.limits.with_values(max_interplanar=3.0, overlap_measure="area"))
    assert all(p.kind != "stack" for p in b.pairs) and any(p.kind == "stack" for p in a.pairs)
    assert b.structure is a.structure and b.cache.contacts is a.cache.contacts
    assert b.cache.values is not a.cache.values  # cache dei valori copiata: le due analisi restano indipendenti
    with pytest.raises(AssertionError):
        api.reanalyse(a, a.limits.with_values(contact_cutoff=3.5))


def test_stesse_soglie_stessa_analisi():
    a = api.analyse_file(EX / "demo_naphthalene.cif", Limits())
    assert api.reanalyse(a, Limits()) is a


def test_i_parametri_geometrici_si_riusano():
    calls = {"geo": 0, "plain": 0}
    reg = Registry()

    @reg.parameter("geo", "geo", geometric=True)
    def geo(ctx):
        calls["geo"] += 1
        return float(ctx.pair.interplanar)

    @reg.parameter("plain", "plain")
    def plain(ctx):
        calls["plain"] += 1
        return float(ctx.pair.kind == "stack")

    a = reg.compute_all(core.analyse_file(EX / "herringbone_anthracene.cif", Limits()))
    n = len(a.pairs)
    assert calls == {"geo": n, "plain": n}
    b = reg.compute_all(core.reanalyse(a, a.limits.with_values(max_interplanar=3.0)))
    assert calls == {"geo": n, "plain": 2 * n}  # il geometrico dalla cache, l'altro di nuovo
    assert [p.get("geo") for p in b.pairs] == [p.interplanar for p in b.pairs]
    assert all(p.get("plain") == 0 for p in b.pairs)
    c = reg.compute_all(core.reanalyse(b, b.limits.with_values(contact_cutoff=3.9)))
    assert c.pairs and calls["geo"] == n + len(c.pairs)  # nuove misure: la cache riparte


def test_i_parametri_extra_sono_geometrici_quelli_dei_plugin_no():
    from pistak.params import REGISTRY

    assert all(s.geometric for s in REGISTRY.specs() if s.source == "extra")
    assert not any(s.geometric for s in REGISTRY.specs() if s.source not in ("extra",))


def test_disordine_senza_testo():
    from pistak.core.cif import parse_cif

    cif = parse_cif((EX / "demo_naphthalene.cif").read_text())
    a = core.analyse_cif(cif, "d", Limits())
    assert core.reanalyse(a, Limits().with_values(max_dihedral=10)).pairs  # gli altri stadi bastano
    with pytest.raises(InputError) as e:
        core.reanalyse(a, Limits().with_values(disorder="all"))
    assert e.value.code == "noSource"


def test_le_geometrie_perturbate_non_portano_la_cache():
    from pistak.core.uncertainty import perturbed

    a = api.analyse_file(EX / "herringbone_anthracene.cif", Limits())
    assert perturbed(a, np.random.default_rng(0)).cache is None
