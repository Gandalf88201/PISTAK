"""Controlli rapidi della struttura: una parte di ciò che controlla checkCIF/PLATON, quella che si può fare con il solo
CIF (senza il file dei fattori di struttura .fcf). Non sostituisce checkCIF: serve a trovare in fretta i problemi che
falsano l'analisi delle coppie (atomi sovrapposti, idrogeni mancanti, disordine, simmetria mancata).

Ogni controllo restituisce avvisi (`Alert`) con un codice stabile, un livello e i valori che lo descrivono; il testo
per l'utente lo costruisce `pistak.labels.alert_text` fuori dal nucleo. Livelli, come in checkCIF:

    A   problema serio, da risolvere o spiegare      B   problema potenzialmente serio
    C   da controllare                               G   informazione

Le soglie sono indicative (dello stesso ordine di quelle di checkCIF, non identiche). La ricerca di simmetria mancata
(ADDSYM) usa il pacchetto facoltativo `spglib`; senza, il controllo viene saltato con un avviso G.
"""
from __future__ import annotations

import math
import warnings
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from pistak.core.crystal import cell_contents, cell_volume, cif_value, density, formula_text, reported_formula
from pistak.core.elements import atomic_number

if TYPE_CHECKING:
    from pistak.core.analysis import Analysis

LEVELS = ("A", "B", "C", "G")

#: Intervalli plausibili delle lunghezze di legame (Å), per coppia di elementi in ordine alfabetico; fuori → avviso.
BOND_RANGES: dict[tuple[str, str], tuple[float, float, str]] = {
    ("C", "C"): (1.15, 1.65, "B"), ("C", "N"): (1.10, 1.56, "B"), ("C", "O"): (1.15, 1.50, "B"),
    ("C", "H"): (0.75, 1.15, "C"), ("H", "N"): (0.70, 1.10, "C"), ("H", "O"): (0.65, 1.10, "C"),
}
#: Contatti intermolecolari corti (Å), con `Limits.contact_check = "fixed"`: coppia di elementi → (soglia, livello).
SHORT_CONTACTS: dict[tuple[str, str], tuple[float, str]] = {("C", "C"): (3.0, "C"), ("H", "H"): (1.8, "C")}
#: Con `contact_check = "vdw"` non si segnalano i legami a idrogeno plausibili: H legato a uno di questi donatori…
HBOND_DONORS = frozenset({"N", "O"})
#: …verso uno di questi accettori.
HBOND_ACCEPTORS = frozenset({"N", "O", "F"})
#: Due atomi distinti più vicini di così sono quasi sovrapposti (Å).
OVERLAP = 0.7
#: Tolleranze di ADDSYM (Å): simmetria trovata alla prima → A, solo alla seconda → C (pseudosimmetria).
ADDSYM_TOLERANCES = (0.1, 0.3)


@dataclass(frozen=True)
class Alert:
    """Un avviso: `code` stabile, `level` (A, B, C, G) e i valori del caso (numeri o testi brevi, chiavi fisse)."""

    code: str
    level: str
    values: dict = field(default_factory=dict)

    @property
    def rank(self) -> int:
        return LEVELS.index(self.level)


Check = Callable[["Analysis"], list[Alert]]
CHECKS: dict[str, Check] = {}


def check(name: str):
    """Decoratore: registra un controllo. I controlli girano nell'ordine di registrazione."""

    def deco(fn: Check) -> Check:
        CHECKS[name] = fn
        return fn

    return deco


def validate(analysis: Analysis) -> list[Alert]:
    """Tutti i controlli; un controllo che fallisce diventa un avviso G `checkFailed` invece di fermare gli altri."""
    out: list[Alert] = []
    for name, fn in CHECKS.items():
        try:
            out += fn(analysis)
        except Exception as e:  # noqa: BLE001
            out.append(Alert("checkFailed", "G", {"check": name, "error": f"{type(e).__name__}: {e}"}))
    return sorted(out, key=lambda a: a.rank)


def counts(alerts: list[Alert]) -> dict[str, int]:
    return {lev: sum(a.level == lev for a in alerts) for lev in LEVELS}


def _labels(items: list[str], n: int = 6) -> str:
    return ", ".join(items[:n]) + (f" (+{len(items) - n})" if len(items) > n else "")


# ---- contenuto della cella ------------------------------------------------------------------------------------
@check("formula")
def _formula(a: Analysis) -> list[Alert]:
    rep = reported_formula(a.cif)
    z = cif_value(a.cif, "z")
    if not rep or not math.isfinite(z) or z <= 0:
        return [Alert("formulaUnchecked", "G")]
    model = {e: n / z for e, n in cell_contents(a.cif).items()}
    out = []
    for heavy in (True, False):
        diffs = []
        for e in sorted(set(rep) | set(model)):
            if (e != "H") != heavy:
                continue
            m, r = model.get(e, 0.0), rep.get(e, 0.0)
            if abs(m - r) > max(0.02 * r, 0.1):
                diffs.append(f"{e} {m:.2f}/{r:g}")
        if diffs:
            out.append(Alert("formulaHeavy" if heavy else "formulaH", "B" if heavy else "C",
                             {"diff": "; ".join(diffs), "model": formula_text(model), "cif": formula_text(rep)}))
    return out


@check("density")
def _density(a: Analysis) -> list[Alert]:
    rep = cif_value(a.cif, "density_cif")
    calc = density(cell_contents(a.cif), cell_volume(a.cif))
    if not (math.isfinite(rep) and math.isfinite(calc)):
        return []
    d = abs(calc - rep)
    if d <= 0.02:
        return []
    return [Alert("density", "B" if d > 0.1 else "C", {"calc": round(calc, 3), "cif": rep})]


# ---- idrogeni -------------------------------------------------------------------------------------------------
@check("hydrogens")
def _hydrogens(a: Analysis) -> list[Alert]:
    els = a.structure.cell_el
    if "C" in els and "H" not in els and not any(at.el == "H" for at in a.cif.omitted):
        return [Alert("noHydrogens", "B")]
    bare: list[str] = []
    for m in a.structure.molecules:
        nb = [[] for _ in m.el]
        for i, j in m.bonds:
            nb[i].append(j)
            nb[j].append(i)
        for k, e in enumerate(m.el):
            if e == "C" and len(nb[k]) <= 1 and not any(m.el[x] == "H" for x in nb[k]) and m.labels[k] not in bare:
                bare.append(m.labels[k])
    return [Alert("missingH", "C", {"n": len(bare), "atoms": _labels(bare)})] if bare else []


# ---- geometria del modello ------------------------------------------------------------------------------------
def _cell_neighbours(a: Analysis, cutoff: float):
    """Coppie di atomi della cella (i nella cella 0, j traslato di s) più vicine di `cutoff`: (i, j, s, d)."""
    s = a.structure
    cart = s.cell_frac @ s.vectors
    tree = cKDTree(cart)
    for shift in np.array([(x, y, z) for x in (-1, 0, 1) for y in (-1, 0, 1) for z in (-1, 0, 1)]):
        sm = tree.sparse_distance_matrix(cKDTree(cart + shift @ s.vectors), cutoff, output_type="ndarray")
        for i, j, d in zip(sm["i"], sm["j"], sm["v"]):
            if shift.any() or i < j:
                yield int(i), int(j), shift, float(d)


@check("overlap")
def _overlap(a: Analysis) -> list[Alert]:
    labels = a.structure.cell_labels
    found = sorted({(min(labels[i], labels[j]), max(labels[i], labels[j]), round(d, 2))
                    for i, j, _s, d in _cell_neighbours(a, OVERLAP) if d > 1e-3})
    if not found:
        return []
    return [Alert("overlap", "A", {"n": len(found), "atoms": _labels([f"{x}···{y} {d:.2f}" for x, y, d in found], 4)})]


@check("bonds")
def _bonds(a: Analysis) -> list[Alert]:
    bad: dict[str, list[str]] = {}
    seen: set[tuple[str, str]] = set()
    for m in a.structure.molecules:
        for i, j in m.bonds:
            pair = tuple(sorted((m.el[i], m.el[j])))
            rng = BOND_RANGES.get(pair)  # type: ignore[arg-type]
            key = tuple(sorted((m.labels[i], m.labels[j])))
            if rng is None or key in seen:
                continue
            seen.add(key)  # type: ignore[arg-type]
            d = float(np.linalg.norm(m.pos[i] - m.pos[j]))
            if not rng[0] <= d <= rng[1]:
                bad.setdefault(rng[2], []).append(f"{key[0]}–{key[1]} {d:.3f}")
    return [Alert("bondLength", lev, {"n": len(items), "bonds": _labels(items, 4)}) for lev, items in sorted(bad.items())]


@check("contacts")
def _contacts(a: Analysis) -> list[Alert]:
    if a.limits.contact_check == "vdw":
        return _vdw_contacts(a, a.limits.contact_check_margin)
    s = a.structure
    n = s.n_cell_atoms
    mol = np.full(n, -1)
    off = np.zeros((n, 3), dtype=int)
    for k, m in enumerate(s.molecules):
        if not m.extended:
            mol[m.cell_index] = k
            off[m.cell_index] = m.cell_shift
    found: dict[tuple[str, str], list[str]] = {}
    seen: set[tuple] = set()
    cutoff = max(t for t, _ in SHORT_CONTACTS.values())
    for i, j, shift, d in _cell_neighbours(a, cutoff):
        if mol[i] < 0 or mol[j] < 0 or (mol[i] == mol[j] and not (shift - off[j] + off[i]).any()):
            continue  # stessa molecola (o rete estesa)
        pair = tuple(sorted((s.cell_el[i], s.cell_el[j])))
        rule = SHORT_CONTACTS.get(pair)  # type: ignore[arg-type]
        key = (*sorted((s.cell_labels[i], s.cell_labels[j])), round(d, 2))
        if rule is None or d >= rule[0] or key in seen:
            continue
        seen.add(key)
        found.setdefault(pair, []).append(f"{key[0]}···{key[1]} {d:.2f}")  # type: ignore[arg-type]
    return [Alert("shortContact", SHORT_CONTACTS[p][1], {"pair": f"{p[0]}···{p[1]}", "limit": SHORT_CONTACTS[p][0],
                                                         "n": len(v), "contacts": _labels(v, 4)}) for p, v in found.items()]


def _vdw_contacts(a: Analysis, margin: float) -> list[Alert]:
    """Contatti intermolecolari più corti della somma dei raggi di van der Waals di oltre `margin` Å, per coppia di
    elementi (gli stessi dell'elenco della scheda Geometria), senza i legami a idrogeno plausibili."""
    from pistak.core.intramolecular import vdw_contacts

    element = {at.label: at.el for at in a.cif.atoms + a.cif.omitted}
    polar_h = {m.labels[i if m.el[i] == "H" else j] for m in a.structure.molecules for i, j in m.bonds
               if {m.el[i], m.el[j]} & {"H"} and ({m.el[i], m.el[j]} - {"H"}) & HBOND_DONORS}
    found: dict[tuple[str, str], list[str]] = {}
    for c in vdw_contacts(a.structure, a.cif):
        if c.delta >= -margin:
            break  # in ordine di differenza crescente
        els = tuple(element.get(x, "") for x in c.atoms)
        if any(lab in polar_h and other in HBOND_ACCEPTORS for lab, other in zip(c.atoms, els[::-1])):
            continue
        found.setdefault(tuple(sorted(els)), []).append(f"{c.atoms[0]}···{c.atoms[1]} {c.value:.2f} ({c.delta:+.2f})")
    return [Alert("vdwContact", "C", {"pair": f"{p[0]}···{p[1]}", "margin": f"{margin:.2f}", "n": len(v),
                                      "contacts": _labels(v, 4)}) for p, v in sorted(found.items())]


# ---- affinamento ----------------------------------------------------------------------------------------------
@check("refinement")
def _refinement(a: Analysis) -> list[Alert]:
    out = []
    r1, wr2, goof = (cif_value(a.cif, k) for k in ("r1", "wr2", "goof"))
    if math.isfinite(r1) and r1 > 0.075:
        out.append(Alert("highR1", "B" if r1 > 0.10 else "C", {"value": r1}))
    if math.isfinite(wr2) and wr2 > 0.20:
        out.append(Alert("highWR2", "B" if wr2 > 0.25 else "C", {"value": wr2}))
    if math.isfinite(goof) and not 0.8 <= goof <= 2.0:
        out.append(Alert("goof", "B" if not 0.6 <= goof <= 4.0 else "C", {"value": goof}))
    rho = max((abs(cif_value(a.cif, k)) for k in ("rho_max", "rho_min")), default=math.nan)
    heaviest = max((atomic_number(e) for e in a.structure.cell_el), default=0)
    if math.isfinite(rho) and heaviest <= 18 and rho > 1.0:  # solo strutture di atomi leggeri: con atomi pesanti è normale
        out.append(Alert("residualDensity", "B" if rho > 2.0 else "C", {"value": rho}))
    if not any(math.isfinite(v) for v in (r1, wr2, goof)):
        out.append(Alert("noRefinement", "G"))
    return out


# ---- disordine e simmetria ------------------------------------------------------------------------------------
@check("disorder")
def _disorder(a: Analysis) -> list[Alert]:
    out = []
    if a.cif.skipped_disorder:
        out.append(Alert("disorderSkipped", "G", {"n": a.cif.skipped_disorder, "policy": a.cif.disorder}))
    partial = sorted({at.label for at in a.cif.atoms if at.occ < 0.99})
    if partial:
        out.append(Alert("partialOccupancy", "C", {"n": len(partial), "atoms": _labels(partial)}))
    if a.cif.symops_source != "cif":
        out.append(Alert("symopsDerived", "G", {"source": a.cif.symops_source}))
    return out


def _spglib_get(ds, key: str):
    return getattr(ds, key) if hasattr(ds, key) else ds[key]


def find_symmetry(a: Analysis, tolerance: float) -> tuple[int, str, int] | None:
    """Gruppo spaziale del modello trovato da spglib: (operazioni nella cella data, simbolo, numero); None senza spglib
    o se spglib non trova nulla."""
    try:
        import spglib
    except ImportError:
        return None
    s = a.structure
    kinds = {e: k + 1 for k, e in enumerate(sorted(set(s.cell_el)))}
    cell = (s.vectors, s.cell_frac, [kinds[e] for e in s.cell_el])
    with warnings.catch_warnings():  # spglib ≥ 2.5 avvisa del vecchio modo di segnalare gli errori (None)
        warnings.simplefilter("ignore", DeprecationWarning)
        try:
            ds = spglib.get_symmetry_dataset(cell, symprec=tolerance)
        except Exception:  # noqa: BLE001  (SpglibError con la nuova gestione degli errori)
            return None
    if ds is None:
        return None
    return len(_spglib_get(ds, "rotations")), str(_spglib_get(ds, "international")), int(_spglib_get(ds, "number"))


@check("addsym")
def _addsym(a: Analysis) -> list[Alert]:
    declared = len(a.cif.symops)
    for k, tol in enumerate(ADDSYM_TOLERANCES):
        found = find_symmetry(a, tol)
        if found is None:
            return [Alert("addsymUnavailable", "G")] if k == 0 else []
        n, symbol, number = found
        if n > declared:
            return [Alert("missedSymmetry" if k == 0 else "pseudoSymmetry", "A" if k == 0 else "C",
                          {"found": f"{symbol} ({number})", "declared": a.cif.space_group or "?", "tolerance": tol,
                           "ops_found": n, "ops_declared": declared})]
    return []
