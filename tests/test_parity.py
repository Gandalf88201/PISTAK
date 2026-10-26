"""Parità con il sito: il porting Python deve dare gli stessi numeri dell'algoritmo TypeScript originale.

`tests/fixtures/site_expected.json` è prodotto da `tests/parity/dump_site.js`, che esegue i sorgenti del
sito (src/lib/crystal/*.ts) sui CIF di `examples/`. Per rigenerarlo:
    SITE_DIR=/percorso/SITOMIO node tests/parity/dump_site.js
"""
import json
from pathlib import Path

import pytest

from pistak.core.analysis import analyse_file
from pistak.core.config import Limits

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = json.loads((ROOT / "tests" / "fixtures" / "site_expected.json").read_text())
NAMES = sorted(k for k in EXPECTED if not k.startswith("."))

# nome in JSON (sito) → attributo Python
FIELDS = {"interplanar": "interplanar", "longSlip": "long_slip", "latSlip": "lat_slip", "slip": "slip",
          "dihedral": "dihedral", "centroidDistance": "centroid_distance", "contact": "contact", "overlap": "overlap"}


@pytest.mark.parametrize("name", NAMES)
def test_stessi_numeri_del_sito(name):
    exp = EXPECTED[name]
    a = analyse_file(ROOT / "examples" / name, Limits.site())  # soglie e convenzioni del sito
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
    # Coppie con lo stesso interplanare (a→b e b→a) possono uscire in ordine scambiato per rumore numerico di
    # ~1e-15 nell'ordinamento: si confrontano dopo averle ordinate con una chiave arrotondata.
    pairs = sorted(a.pairs, key=lambda p: (p.kind, round(p.interplanar, 4), round(p.contact, 4), round(p.long_slip, 4)))
    expected = sorted(exp["pairs"], key=lambda e: (e["kind"], round(e["interplanar"], 4), round(e["contact"], 4), round(e["longSlip"], 4)))
    for p, ep in zip(pairs, expected):
        assert p.kind == ep["kind"] and p.count == ep["count"]
        for js, py in FIELDS.items():
            # l'angolo vicino a 0° amplifica il rumore numerico (acos): 1e-6° tra TypeScript e numpy è lo stesso valore
            tol = 1e-4 if js == "dihedral" else 1e-6
            assert getattr(p, py) == pytest.approx(ep[js], abs=tol), f"{name} coppia {p.index + 1}: {js}"


def test_il_confronto_non_e_banale():
    """Almeno un esempio con coppie di tutti e tre i tipi, altrimenti la parità proverebbe poco."""
    kinds = {p["kind"] for e in EXPECTED.values() if "pairs" in e for p in e["pairs"]}
    assert kinds == {"stack", "parallel", "tilted"}
