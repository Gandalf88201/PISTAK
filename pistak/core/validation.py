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
from pistak.core.elements import DEFAULT_VDW, VDW_RADII, atomic_number

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
#: Un C terminale legato solo a N o O più corto di così è un triplo legame (cianuro, isonitrile, carbonile metallico)
#: e non porta idrogeni; un doppio legame C=O è più lungo. F. H. Allen et al., J. Chem. Soc. Perkin Trans. 2 1987,
#: S1–S19, tabella: Csp≡N da 1,136 Å (nitrili C*–C≡N, n = 140) a 1,155 Å ((S–C≡N)⁻); C=O da 1,187 Å (anidridi, esteri
#: arilici; quartile inferiore 1,181) a 1,21–1,23 (chetoni, ammidi). Sui 200 CIF COD della revisione (registro,
#: voce 23; `check_struttura.py hthresh`) i 697 C=O di carbonili partono da 1,182 Å (mediana 1,221): nessuno sotto.
TRIPLE_CX = 1.18
#: Un C con due vicini e un angolo sotto questo valore è sp² (≈ 120°) e gli manca un H; un C sp (≈ 180°) no. Allen 1987
#: dà solo lunghezze di legame: la soglia si controlla sui CIF COD, dove i C sp² con un H arrivano a 134,5° e i C sp
#: (nitrili, alchini) stanno a 176° di mediana.
SP_ANGLE = 150.0


def _h_deficient(m, nb: list[list[int]], k: int) -> tuple[bool, bool]:
    """(il C `k` ha un solo vicino non idrogeno e non è un triplo legame terminale,
        il C `k` ha meno vicini di quanti ne vuole la sua geometria: il primo caso, o due vicini ad angolo < SP_ANGLE)."""
    if any(m.el[x] == "H" for x in nb[k]):
        return False, False
    if len(nb[k]) == 0:
        return True, True
    if len(nb[k]) == 1:
        x = nb[k][0]
        triple = m.el[x] in ("N", "O") and float(np.linalg.norm(m.pos[x] - m.pos[k])) <= TRIPLE_CX
        return not triple, not triple
    if len(nb[k]) == 2:
        u, v = m.pos[nb[k][0]] - m.pos[k], m.pos[nb[k][1]] - m.pos[k]
        ang = math.degrees(math.acos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))
        return False, ang < SP_ANGLE
    return False, False


@check("hydrogens")
def _hydrogens(a: Analysis) -> list[Alert]:
    """Idrogeni mancanti. Senza nessun H nel modello l'avviso B scatta solo se qualche C sembra volerne (terminale,
    o con due vicini ad angolo da sp²): un cianuro, un carbonile metallico, un fullerene non hanno idrogeni."""
    els = a.structure.cell_el
    bare: list[str] = []
    deficient = False
    for m in a.structure.molecules:
        nb = [[] for _ in m.el]
        for i, j in m.bonds:
            nb[i].append(j)
            nb[j].append(i)
        for k, e in enumerate(m.el):
            if e != "C":
                continue
            terminal, wants_h = _h_deficient(m, nb, k)
            deficient |= wants_h
            if terminal and m.labels[k] not in bare:
                bare.append(m.labels[k])
    if "C" in els and "H" not in els and not any(at.el == "H" for at in a.cif.omitted):
        return [Alert("noHydrogens", "B")] if deficient else []
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
        found.setdefault(pair, []).append(key)  # type: ignore[arg-type]
    # dal più corto, poi per etichette: l'ordine della ricerca dei vicini cambia da una macchina all'altra
    return [Alert("shortContact", SHORT_CONTACTS[p][1], {
        "pair": f"{p[0]}···{p[1]}", "limit": SHORT_CONTACTS[p][0], "n": len(v),
        "contacts": _labels([f"{x}···{y} {d:.2f}" for x, y, d in sorted(v, key=lambda k: (k[2], k[0], k[1]))], 4),
    }) for p, v in sorted(found.items())]


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
@check("vdwRadii")
def _vdw_radii(a: Analysis) -> list[Alert]:
    """Elementi senza raggio di van der Waals tabulato: contatti e parametri di vdW usano `DEFAULT_VDW`."""
    missing = sorted({e for e in a.structure.cell_el if e not in VDW_RADII}, key=atomic_number)
    return [Alert("vdwDefault", "G", {"elements": ", ".join(missing), "r": DEFAULT_VDW})] if missing else []


@check("oldDefaults")
def _old_defaults(a: Analysis) -> list[Alert]:
    """Soglie rimaste a un vecchio valore predefinito (di solito da un user.toml salvato prima del cambio)."""
    from pistak.core.config import superseded_values

    return [Alert("oldDefault", "C", {"key": k, "value": v, "default": d, "date": when})
            for k, v, d, when in superseded_values(a.limits)]


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
    if a.cif.atoms and not any(np.any(at.su > 0) for at in a.cif.atoms):
        out.append(Alert("noSU", "G"))  # revisione, voce 54: il Monte Carlo darebbe σ = 0
    else:
        out += _fragment_su(a)
    return out


def _fragment_su(a: Analysis) -> list[Alert]:
    """Atomi non H dei frammenti π senza s.u. in un CIF che le ha (gruppi rigidi, frammenti disordinati vincolati): nel
    Monte Carlo restano fermi e la σ dei descrittori di quei frammenti è sottostimata (revisione 2, voce 62)."""
    s, cif = a.structure, a.cif
    if not len(s.cell_origin) or not a.fragments:
        return []
    from pistak.core.uncertainty import site_projectors  # import locale: uncertainty importa l'analisi

    fixed = {i for i, p in site_projectors(cif).items() if np.abs(p).max() < 1e-6}  # tutte le coordinate fissate
    bare = set()
    for m, frags in zip(s.molecules, a.fragments):
        for f in frags:
            for i in f.atoms:
                if m.el[i] in ("H", "D"):
                    continue
                ia = int(s.cell_origin[m.cell_index[i], 0])
                at = cif.atoms[ia]
                if ia not in fixed and not np.any(at.su > 0):
                    bare.add(at.label)
    return [Alert("fragmentNoSU", "C", {"n": len(bare), "atoms": _labels(sorted(bare))})] if bare else []


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


def _op_set(ops) -> set:
    """Operazioni come insieme confrontabile: R intera e t in dodicesimi modulo il reticolo."""
    return {(tuple(np.rint(o.r).astype(int).ravel()), tuple(int(round(12 * x)) % 12 for x in o.t)) for o in ops}


@check("symops")
def _symops_vs_symbol(a: Analysis) -> list[Alert]:
    """Le operazioni elencate nel CIF (che PISTAK usa) contro quelle del simbolo del gruppo dichiarato."""
    from pistak.core.cif import parse_symop
    from pistak.core.spacegroups import operations

    if a.cif.symops_source != "cif" or not a.cif.space_group:
        return []
    hall = next((a.cif.items[t] for t in ("_symmetry_space_group_name_hall", "_space_group_name_hall") if a.cif.items.get(t)), None)
    found = operations(hall=hall, hm=a.cif.space_group)
    if found is None:
        return []
    from_symbol = _op_set(parse_symop(t) for t in found[0])
    if _op_set(a.cif.symops) == from_symbol:
        return []
    return [Alert("symopsMismatch", "C", {"declared": a.cif.space_group, "n_cif": len(a.cif.symops), "n_symbol": len(from_symbol)})]


# ---- ricerca delle coppie -------------------------------------------------------------------------------------
@check("shells")
def _shells(a: Analysis) -> list[Alert]:
    """Il limite `MAX_SHELLS` alle celle esplorate taglia la ricerca delle coppie: lo si dice invece di tacerlo."""
    from pistak.core.pairs import neighbour_shells, shells_needed

    cores = a.fragments or a.cores
    need, used = shells_needed(a.structure, cores, a.limits), neighbour_shells(a.structure, cores, a.limits)
    if need == used:
        return []
    return [Alert("shellsTruncated", "B", {"needed": ", ".join(f"±{n}" for n in need), "max": max(used)})]


def _spglib_get(ds, key: str):
    return getattr(ds, key) if hasattr(ds, key) else ds[key]


def find_symmetry(a: Analysis, tolerance: float) -> tuple[int, str, int] | None:
    """Gruppo spaziale del modello trovato da spglib: (operazioni nella cella data, simbolo, numero); None senza spglib.
    Se spglib non riesce (per esempio atomi sovrapposti) solleva un'eccezione, che `validate` riporta come avviso G."""
    try:
        import spglib
    except ImportError:
        return None
    s = a.structure
    keep = ordered_atoms(a)
    els = [e for e, k in zip(s.cell_el, keep) if k]
    kinds = {e: k + 1 for k, e in enumerate(sorted(set(els)))}
    cell = (s.vectors, s.cell_frac[keep], [kinds[e] for e in els])
    with warnings.catch_warnings():  # spglib ≥ 2.5 avvisa del vecchio modo di segnalare gli errori (None)
        warnings.simplefilter("ignore", DeprecationWarning)
        ds = spglib.get_symmetry_dataset(cell, symprec=tolerance)  # spglib ≥ 2.5: SpglibError se non riesce
    if ds is None:  # spglib < 2.5: None se non riesce. Non è "spglib manca": diventa l'avviso G checkFailed
        raise RuntimeError(f"spglib non trova la simmetria (symprec {tolerance} Å)")
    return len(_spglib_get(ds, "rotations")), str(_spglib_get(ds, "international")), int(_spglib_get(ds, "number"))


#: Come PLATON ADDSYM, gli atomi con occupanza parziale (componenti di disordine) restano fuori dalla ricerca: la
#: componente tenuta dalla regola del disordine può avere più simmetria del modello completo. Se resta fuori più di
#: questa frazione degli atomi non idrogeno, la ricerca non si fa.
ADDSYM_MAX_EXCLUDED = 0.5


def ordered_atoms(a: Analysis) -> np.ndarray:
    """Per ogni atomo della cella: vero se viene da un sito con occupanza piena."""
    s = a.structure
    if not len(s.cell_origin):
        return np.ones(len(s.cell_el), dtype=bool)
    occ = np.array([a.cif.atoms[int(i)].occ for i in s.cell_origin[:, 0]])
    return occ >= 0.99


@check("addsym")
def _addsym(a: Analysis) -> list[Alert]:
    heavy = np.array([e != "H" for e in a.structure.cell_el])
    excluded = int(np.count_nonzero(heavy & ~ordered_atoms(a)))
    if heavy.any() and excluded > ADDSYM_MAX_EXCLUDED * np.count_nonzero(heavy):
        return [Alert("addsymDisorder", "G", {"n": excluded})]
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
