"""Round-5 review of the paper (PISTAK_review_round5.md, N1): results for extended networks must not depend on how
the repeat unit is cut. Lattice of periodicities, same bonded network, bridging atoms in other copies of the unit,
stacking graph on the fragments of networks, symmetry images of fragments that share an atom."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from pistak.core.analysis import _assemble, _contacts, _fragments, read_cif_file
from pistak.core.cif import parse_cif
from pistak.core.config import Limits
from pistak.core.cores import small_rings
from pistak.core.elements import NONMETALS
from pistak.core.structure import (
    _dedupe,
    build_structure,
    cell_bonds,
    in_lattice,
    lattice_basis,
    network_links,
    same_network,
    to_cart,
)

COD = Path(__file__).resolve().parent / "data" / "cod"


def test_lattice_basis_and_membership():
    b = lattice_basis([(0, 0, 2), (0, 0, 3), (0, 0, 0)])
    assert b.tolist() == [[0, 0, 1]]
    b = lattice_basis([(2, 0, 0), (1, 1, 0)])
    assert in_lattice(b, (3, 1, 0)) and in_lattice(b, (0, 2, 0)) and not in_lattice(b, (1, 0, 0))
    assert not in_lattice(b, (0, 0, 1))
    assert lattice_basis([]).shape == (0, 3)


def _recut(s, limits, rng):
    """The same crystal with every network repeat unit cut differently: rings and non-metal ligands kept whole, each
    piece moved by a random period of its network (as paper/scripts/recut_invariance.py)."""
    edges = cell_bonds(s.cell_el, to_cart(s.vectors, s.cell_frac), s.vectors, limits)
    inv = np.linalg.inv(s.vectors)
    mols = []
    for m in s.molecules:
        if not m.extended or not len(m.periods):
            mols.append(m)
            continue
        n = len(m)
        parent = list(range(n))

        def find(x, parent=parent):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i, j in m.bonds:
            if m.el[int(i)] in NONMETALS and m.el[int(j)] in NONMETALS:
                parent[find(int(i))] = find(int(j))
        for r in small_rings(n, m.bonds, lambda _x: True):
            for x in r[1:]:
                parent[find(int(x))] = find(int(r[0]))
        cell = [tuple(int(x) for x in cs) for cs in m.cell_shift]
        groups: dict[int, list[int]] = {}
        for k in range(n):
            groups.setdefault(find(k), []).append(k)
        for g in groups.values():
            if rng.random() < 0.5:
                continue
            t = m.periods[0] * (1 if rng.random() < 0.5 else -1)
            for k in g:
                cell[k] = tuple(cell[k][d] + int(t[d]) for d in range(3))
        where = {int(ci): k for k, ci in enumerate(m.cell_index)}
        bonds = [(k, where[int(j)]) for k, ci in enumerate(m.cell_index) for j, es in edges[int(ci)]
                 if int(j) in where and where[int(j)] > k and cell[where[int(j)]] == tuple(cell[k][d] + es[d] for d in range(3))]
        cs = np.array(cell, dtype=int)
        pos = to_cart(s.vectors, s.cell_frac[m.cell_index] + cs)
        c = pos.mean(axis=0)
        fc = np.floor(c @ inv + 1e-6)
        mols.append(replace(m, pos=pos - fc @ s.vectors, centroid=c - fc @ s.vectors, bonds=_dedupe(bonds),
                            cell_shift=cs - fc.astype(int), links=network_links(m.cell_index, cs, edges)))
    return replace(s, molecules=mols, _images={}, _atom_maps={})


def _signature(a):
    rows = []
    for p in a.pairs:
        vals = [-1.0 if v != v else round(float(v), 6) for v in (p.interplanar, p.centroid_distance, p.dihedral,
                                                                 p.overlap_area, p.slip, p.contact)]
        own = tuple(sorted([tuple(-1.0 if v != v else round(float(v), 6) for v in (p.long_slip_a, p.lat_slip_a)),
                            tuple(-1.0 if v != v else round(float(v), 6) for v in (p.long_slip_b, p.lat_slip_b))]))
        cat = "network" if (p.same_network or p.intramolecular) else ("between" if p.polymer else "molecules")
        rows.append((p.kind, cat, p.bonded, p.bridged, p.clamped, p.count, tuple(vals), own))
    return sorted(rows), sorted((c.rank, len(c.nodes)) for c in a.stacks)


def _analyse(path, rng=None):
    lim = Limits()
    cif = parse_cif(read_cif_file(path), disorder=lim.disorder)
    s = build_structure(cif, lim)
    if rng is not None:
        s = _recut(s, lim, rng)
    f = _fragments(s, lim)
    return _assemble(path.name, cif, s, f, lim, _contacts(s, f, lim, None, None))


def test_same_network_follows_the_periods_of_the_network():
    a = _analyse(COD / "7024607.cif")
    s = a.structure
    nets = [i for i, m in enumerate(s.molecules) if m.extended]
    assert nets and all(len(s.molecules[i].periods) for i in nets)
    for p in a.pairs:
        expected = p.a == p.b and (not any(p.shift) or in_lattice(s.molecules[p.a].periods, p.shift))
        assert p.same_network == expected == same_network(s, p.a, p.b, p.shift)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_results_do_not_depend_on_the_cut_of_the_repeat_unit(seed):
    path = COD / "7024607.cif"
    assert _signature(_analyse(path)) == _signature(_analyse(path, np.random.default_rng(seed)))
