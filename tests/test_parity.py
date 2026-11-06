"""Parity with the website: the Python port must give the same numbers as the original TypeScript algorithm.

`tests/fixtures/site_expected.json` is produced by `tests/parity/dump_site.js`, which runs the website sources
(src/lib/crystal/*.ts) on the CIFs in `examples/`. To regenerate it:
    SITE_DIR=/path/to/SITE node tests/parity/dump_site.js
"""
import json
from pathlib import Path

import pytest

from pistak.core.analysis import analyse_file
from pistak.core.config import Limits

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = json.loads((ROOT / "tests" / "fixtures" / "site_expected.json").read_text())
NAMES = sorted(k for k in EXPECTED if not k.startswith("."))

# name in the JSON (website) → Python attribute
FIELDS = {"interplanar": "interplanar", "longSlip": "long_slip", "latSlip": "lat_slip", "slip": "slip",
          "dihedral": "dihedral", "centroidDistance": "centroid_distance", "contact": "contact", "overlap": "overlap"}


@pytest.mark.parametrize("name", NAMES)
def test_same_numbers_as_website(name):
    exp = EXPECTED[name]
    a = analyse_file(ROOT / "examples" / name, Limits.site())  # website thresholds and conventions
    assert a.structure.n_cell_atoms == exp["cellAtoms"]
    assert len(a.structure.molecules) == len(exp["molecules"])
    for i, (m, em) in enumerate(zip(a.structure.molecules, exp["molecules"])):
        assert len(m) == em["atoms"] and len(m.bonds) == em["bonds"] and m.extended == em["extended"]
        core = a.cores[i]
        if em["core"] is None:
            assert core is None
        else:
            assert core is not None and core.kind == em["core"]["kind"] and len(core.atoms) == em["core"]["atoms"]
            assert core.rms == pytest.approx(em["core"]["rms"], abs=1e-6)
    assert len(a.pairs) == len(exp["pairs"])
    # Pairs with the same interplanar distance (a→b and b→a) can come out swapped because of ~1e-15 numerical
    # noise in the sort: they are compared after sorting them with a rounded key.
    pairs = sorted(a.pairs, key=lambda p: (p.kind, round(p.interplanar, 4), round(p.contact, 4), round(p.long_slip, 4)))
    expected = sorted(exp["pairs"], key=lambda e: (e["kind"], round(e["interplanar"], 4), round(e["contact"], 4), round(e["longSlip"], 4)))
    for p, ep in zip(pairs, expected):
        assert p.kind == ep["kind"] and p.count == ep["count"]
        for js, py in FIELDS.items():
            # an angle near 0° amplifies numerical noise (acos): 1e-6° between TypeScript and numpy is the same value
            tol = 1e-4 if js == "dihedral" else 1e-6
            assert getattr(p, py) == pytest.approx(ep[js], abs=tol), f"{name} pair {p.index + 1}: {js}"


def test_comparison_is_not_trivial():
    """At least one example with pairs of all three kinds, otherwise parity would prove little."""
    kinds = {p["kind"] for e in EXPECTED.values() if "pairs" in e for p in e["pairs"]}
    assert kinds == {"stack", "parallel", "tilted"}
