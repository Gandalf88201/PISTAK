import math
from pathlib import Path

import numpy as np
import pytest

from pistak.core.analysis import analyse_file
from pistak.core.config import Limits
from pistak.core.demo import DEMO_STACK
from pistak.labels import rank_name

EX = Path(__file__).resolve().parents[1] / "examples"


def test_demo_e_una_colonna_lungo_a():
    a = analyse_file(EX / "demo_naphthalene.cif", Limits())
    assert len(a.stacks) == 1
    col = a.stacks[0]
    assert col.rank == 1 and rank_name(col.rank) == "colonna" and col.kind_key == "column" and col.molecules == [0]
    d = DEMO_STACK
    assert col.repeat == pytest.approx(math.sqrt(d["interplanar"] ** 2 + d["long_slip"] ** 2 + d["lat_slip"] ** 2), abs=1e-3)
    assert col.mean("interplanar") == pytest.approx(d["interplanar"], abs=1e-3)
    normal = a.cores[0].normal
    assert col.axis_normal_angle(normal) == pytest.approx(math.degrees(math.acos(d["interplanar"] / col.repeat)), abs=0.01)


def test_herringbone_ha_colonne_di_rango_1():
    a = analyse_file(EX / "herringbone_anthracene.cif", Limits())
    assert a.stacks and all(c.rank == 1 for c in a.stacks)
    assert all(np.linalg.norm(c.axis) > 0 for c in a.stacks)


def test_senza_stack_nessuna_colonna():
    assert analyse_file(EX / "infinite_chain.cif", Limits()).stacks == []
