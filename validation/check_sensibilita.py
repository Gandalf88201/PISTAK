"""Checklist, sezione 10 (passo 5): quanto i risultati dipendono dalle soglie.

Ogni soglia si sposta attorno al valore predefinito, una alla volta, sui CIF COD di `cod_cache/`; per ogni variante si
confrontano le coppie con quelle delle soglie predefinite. Una coppia è la stessa se ha le stesse due molecole, la
stessa traslazione e gli stessi frammenti π. Si contano:

    cambiano tipo   coppie presenti in entrambe le analisi che passano da stack a parallel/tilted o viceversa
    nuove / perse   coppie che compaiono o spariscono (nuclei π o ricerca delle coppie diversi)
    stack           numero totale di coppie stack, contro quello con le soglie predefinite
    strutture       strutture in cui il numero di stack cambia

Non è un oracolo esterno: dice quanto peso ha ciascuna scelta, da dichiarare nel paper.

    .venv/bin/python validation/check_sensibilita.py [N]     # N = numero di CIF (predefinito: tutti)
"""
from __future__ import annotations

import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from pistak.core.analysis import analyse_text, reanalyse  # noqa: E402
from pistak.core.config import Limits  # noqa: E402

#: (soglia, valori provati). I predefiniti sono in `Limits`.
VARIANTS: list[tuple[str, tuple]] = [
    ("max_dihedral", (10.0, 15.0, 25.0, 30.0)),
    ("max_interplanar", (3.8, 4.0, 4.4, 4.6)),
    ("contact_cutoff", (4.0, 5.0)),
    ("ring_rms", (0.05, 0.15)),
    ("molecule_rms", (0.10, 0.20)),
    ("bond_tolerance", (0.3, 0.5)),
    ("max_ring_size", (6, 10)),
    ("min_fragment_atoms", (4, 6)),
]


def _pairs(a) -> dict[tuple, str]:
    return {(p.a, p.b, tuple(p.shift), p.fa, p.fb): p.kind for p in a.pairs}


def one(path: Path) -> dict | None:
    """Coppie con le soglie predefinite e con ogni variante, per un CIF (None se non si analizza)."""
    try:
        base = analyse_text(path.read_text(encoding="utf-8", errors="replace"), path.stem, Limits())
    except Exception:  # noqa: BLE001  (CIF che PISTAK rifiuta: fuori dal confronto, contati a parte)
        return None
    out = {"base": _pairs(base)}
    for name, values in VARIANTS:
        for v in values:
            try:
                out[(name, v)] = _pairs(reanalyse(base, replace(base.limits, **{name: v})))
            except Exception as e:  # noqa: BLE001
                out[(name, v)] = e
    return out


def main(n: int | None) -> None:
    files = sorted((HERE / "cod_cache").glob("[0-9]*.cif"))[:n]  # [0-9]: non i file ._ di macOS
    with ProcessPoolExecutor() as ex:
        results = [r for r in ex.map(one, files) if r is not None]
    base_stacks = sum(Counter(r["base"].values())["stack"] for r in results)
    base_pairs = sum(len(r["base"]) for r in results)
    print(f"== sensibilità alle soglie: {len(results)} strutture su {len(files)} ({len(files) - len(results)} rifiutate), "
          f"{base_pairs} coppie, {base_stacks} stack con le soglie predefinite ==")
    d = Limits()
    print(f"  {'soglia':<20}{'valore':>8}  {'cambiano tipo':>13}  {'nuove':>6}  {'perse':>6}  {'stack':>12}  {'strutture':>9}")
    for name, values in VARIANTS:
        print(f"  {name:<20}{getattr(d, name)!s:>8}  (predefinito)")
        for v in values:
            flips = new = lost = stacks = structs = errors = 0
            for r in results:
                b, x = r["base"], r[(name, v)]
                if isinstance(x, Exception):
                    errors += 1
                    continue
                common = b.keys() & x.keys()
                flips += sum((b[k] == "stack") != (x[k] == "stack") for k in common)
                new += len(x.keys() - b.keys())
                lost += len(b.keys() - x.keys())
                s_b, s_x = Counter(b.values())["stack"], Counter(x.values())["stack"]
                stacks += s_x
                structs += s_b != s_x
            delta = f"{stacks} ({100 * (stacks - base_stacks) / base_stacks:+.0f}%)" if base_stacks else str(stacks)
            err = f"  [{errors} errori]" if errors else ""
            print(f"  {'':<20}{v!s:>8}  {flips:>13}  {new:>6}  {lost:>6}  {delta:>12}  {structs:>9}{err}")


if __name__ == "__main__":  # necessario: i processi reimportano questo modulo
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
