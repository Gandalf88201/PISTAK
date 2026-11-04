"""Checklist, sezioni 4–5: parametri anello–anello di PISTAK (`ringpairs.py`) contro PLATON `CALC ALL` (oracolo esterno).

Gli output `.lis` stanno in `validation/platon/` (vedi `LEGGIMI.md`). Gli anelli si abbinano per insieme di etichette
degli atomi (PLATON li elenca in «n-Membered Ring (k)»), le coppie per anelli e distanza Cg–Cg. Per ogni coppia
abbinata si confrontano Cg–Cg, α, β, γ, le due distanze perpendicolari e lo slippage.

Corrispondenze: PLATON CgJ_Perp (Cg(J) sul piano I) = PISTAK perp_a; CgI_Perp = perp_b; slippage = proiezione di
Cg(J) sul piano I, come in PISTAK.

    .venv/bin/python validation/check_platon.py
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from pistak.api import analyse_file  # noqa: E402
from pistak.core.config import Limits  # noqa: E402
from pistak.core.ringpairs import ring_pairs  # noqa: E402

DIR = HERE / "platon"
NUM = re.compile(r"-?\d+\.\d+(?:\(\d+\))?")
ROW = re.compile(r"^Cg(\d+)\s+\[\s*\d+\]\s+->\s+Cg(\d+)\s+\[\s*(\d+)\.(\d+)\](.*)$")


def platon_rings(text: str) -> dict[int, frozenset[str]]:
    """Numero dell'anello → etichette degli atomi, dalle righe «6-Membered Ring ( 1)   N1 --> C23 --> …»."""
    rings = {}
    for m in re.finditer(r"(\d+)-Membered Ring\s*\(\s*(\d+)\)\s+(.*)", text):
        labels = [t for t in re.split(r"\s*-->\s*", m.group(3).strip()) if t]
        rings[int(m.group(2))] = frozenset(labels[: int(m.group(1))])
    return rings


def platon_pairs(text: str) -> list[dict]:
    """Le righe della tabella «Analysis of Short Ring-Interactions»."""
    start = text.find("Analysis of Short Ring-Interactions")
    if start < 0:
        return []
    block = text[start:text.find("Min or Max", start)]
    out = []
    for line in block.splitlines():
        m = ROW.match(line.strip())
        if not m:
            continue
        nums = [float(x.split("(")[0]) for x in NUM.findall(m.group(5))]
        # Cg-Cg, P, Q, R, S, alpha, beta, gamma, CgI_Perp, CgJ_Perp, [slippage]
        if len(nums) < 10:
            continue
        out.append({"i": int(m.group(1)), "j": int(m.group(2)), "cg": nums[0], "alpha": nums[5], "beta": nums[6],
                    "gamma": nums[7], "perp_i": abs(nums[8]), "perp_j": abs(nums[9]),
                    "slip": nums[10] if len(nums) > 10 else None, "intra": m.group(3) == "1" and m.group(4) == "555"})
    return out


def compare(stem: str) -> dict:
    text = (DIR / f"{stem}.lis").read_text(encoding="latin-1")
    rings = platon_rings(text)
    ref = [r for r in platon_pairs(text) if r["i"] in rings and r["j"] in rings]
    an = analyse_file(DIR / f"{stem}.cif", Limits())
    mine = []
    for p in an.pairs:
        if p.intramolecular:
            continue
        ctx = an.pair_context(p)
        for rp in ring_pairs(ctx):
            la = frozenset(ctx.A.labels[i] for i in rp.atoms_a)
            lb = frozenset(ctx.B.labels[i] for i in rp.atoms_b)
            mine.append((la, lb, rp))
    diffs: dict[str, list[float]] = {k: [] for k in ("cg", "alpha", "beta", "gamma", "perp_a", "perp_b", "slip")}
    matched_ref = set()
    unmatched_mine = 0
    rows: list[dict] = []  # una per coppia abbinata, nella notazione di PLATON (anello I → anello J)
    for la, lb, rp in mine:
        hits = []
        for k, r in enumerate(ref):
            ri, rj = rings[r["i"]], rings[r["j"]]
            if abs(r["cg"] - rp.cg_distance) > 0.006:
                continue
            if (ri, rj) == (la, lb):
                hits.append((k, False))
            if (ri, rj) == (lb, la):  # anche con la == lb: un anello e la sua immagine, PLATON può dare l'altro verso
                hits.append((k, True))
        if not hits:
            unmatched_mine += 1
            continue

        def score(h, rp=rp):  # tra le righe con gli stessi anelli e la stessa Cg–Cg, quella più vicina negli angoli
            r, sw = ref[h[0]], h[1]
            b, g = (rp.gamma, rp.beta) if sw else (rp.beta, rp.gamma)
            return abs(rp.alpha - r["alpha"]) + abs(b - r["beta"]) + abs(g - r["gamma"])

        k, swapped = min(hits, key=score)
        r = ref[k]
        matched_ref.update(h[0] for h in hits if score(h) < 1.0)
        beta, gamma = (rp.gamma, rp.beta) if swapped else (rp.beta, rp.gamma)
        pa, pb = (rp.perp_b, rp.perp_a) if swapped else (rp.perp_a, rp.perp_b)
        diffs["cg"].append(rp.cg_distance - r["cg"])
        diffs["alpha"].append(rp.alpha - r["alpha"])
        diffs["beta"].append(beta - r["beta"])
        diffs["gamma"].append(gamma - r["gamma"])
        diffs["perp_a"].append(pa - r["perp_j"])
        diffs["perp_b"].append(pb - r["perp_i"])
        slip = float(np.sqrt(max(0.0, rp.cg_distance ** 2 - pa ** 2)))
        if r["slip"] is not None:
            diffs["slip"].append(slip - r["slip"])
        rows.append({"stem": stem, "i": r["i"], "j": r["j"], "size_i": len(rings[r["i"]]), "size_j": len(rings[r["j"]]),
                     "pistak": {"cg": rp.cg_distance, "alpha": rp.alpha, "beta": beta, "gamma": gamma,
                                "perp_i": pb, "perp_j": pa, "slip": slip},
                     "platon": {k: r[k] for k in ("cg", "alpha", "beta", "gamma", "perp_i", "perp_j", "slip")}})
    inter = [k for k, r in enumerate(ref) if not r["intra"]]
    # classi: stessi due anelli (in qualunque ordine) e stessa Cg–Cg; PLATON elenca anche le copie per simmetria
    cls = lambda r: (frozenset((rings[r["i"]], rings[r["j"]])), round(r["cg"], 2))  # noqa: E731
    found = {cls(ref[k]) for k in matched_ref}
    missed = [r for r in {cls(ref[k]): ref[k] for k in inter}.values() if cls(r) not in found]
    return {"stem": stem, "mine": len(mine), "matched": len(mine) - unmatched_mine, "ref_inter": len(inter),
            "ref_missed": missed, "diffs": diffs, "n_rings": len(rings), "rows": rows, "rings": rings,
            "space_group": an.cif.space_group}


def main() -> None:
    stems = [s for s in (DIR / "elenco.txt").read_text().split() if (DIR / f"{s}.lis").exists()]
    total: dict[str, list[float]] = {}
    print("== coppie di anelli: PISTAK contro PLATON CALC ALL ==")
    for stem in stems:
        res = compare(stem)
        if not res["n_rings"]:
            print(f"  {stem}: PLATON senza anelli (CIF non letto?)")
            continue
        for k, v in res["diffs"].items():
            total.setdefault(k, []).extend(v)
        worst = {k: max((abs(x) for x in v), default=0) for k, v in res["diffs"].items()}
        print(f"  {stem}: PISTAK {res['mine']} coppie di anelli, abbinate {res['matched']}; PLATON {res['ref_inter']} "
              f"intermolecolari, classi non trovate da PISTAK {len(res['ref_missed'])}; scarto max "
              + ", ".join(f"{k} {v:.4f}" for k, v in worst.items()))
    print("  totale, scarto massimo e mediano (PISTAK − PLATON):")
    for k, v in total.items():
        a = np.abs(v)
        if len(a):
            print(f"    {k:7s} n = {len(a):4d}  max {a.max():.4f}  mediana {np.median(a):.4f}")


if __name__ == "__main__":
    main()
