"""Fissa il campione COD del paper a N strutture (revisione 4: 500 al posto di 200, decisione dell'autore).

Il campione è costituito dai primi N ID trovati da `scarica_cod.py` nella sua estrazione a seme fisso (2026): i 200
del campione precedente ne sono i primi 200. Si leggono dal registro dello scaricamento (righe «k/2000 ID»), nell'ordine
in cui sono stati trovati.

    .venv/bin/python validation/congela_campione.py 500 paper/data/runs/rev4/download_2000.txt

Scrive `cod_sample.txt` (ordinato) e `cod_sample.sha256`, tiene i precedenti come `cod_sample_200.*`, e porta
`cod_cache/` a contenere esattamente quei CIF (gli script di verifica prendono il campione dall'elenco di quella
cartella).
"""
from __future__ import annotations

import hashlib
import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main(n: int, log: Path) -> None:
    found = [int(m.group(2)) for m in re.finditer(r"^(\d+)/\d+ (\d+)$", log.read_text(), re.M)]
    assert len(found) >= n, f"only {len(found)} structures in {log}"
    ids = found[:n]
    assert len(set(ids)) == n
    old = [int(x) for x in (HERE / "cod_sample.txt").read_text().split()]
    if len(old) != n:
        missing = sorted(set(old) - set(ids))
        assert not missing, f"structures of the previous sample not among the first {n}: {missing}"
        assert set(found[: len(old)]) == set(old), "the previous sample is not the first part of the new one"
        shutil.copy2(HERE / "cod_sample.txt", HERE / f"cod_sample_{len(old)}.txt")
        shutil.copy2(HERE / "cod_sample.sha256", HERE / f"cod_sample_{len(old)}.sha256")
    src, cache = HERE / "cod_cache_2000", HERE / "cod_cache"
    for i in ids:
        if not (cache / f"{i}.cif").exists():
            shutil.copy2(src / f"{i}.cif", cache / f"{i}.cif")
    extra = [f for f in cache.glob("[0-9]*.cif") if int(f.stem) not in set(ids)]
    assert not extra, f"cod_cache holds files outside the sample: {[f.name for f in extra][:5]}"
    assert len(list(cache.glob("[0-9]*.cif"))) == n
    (HERE / "cod_sample.txt").write_text("\n".join(map(str, sorted(ids))) + "\n")
    lines = [f"{hashlib.sha256((cache / f'{i}.cif').read_bytes()).hexdigest()}  {i}.cif" for i in sorted(ids)]
    (HERE / "cod_sample.sha256").write_text("\n".join(lines) + "\n")
    # the checksums of the previous sample must not have changed
    prev = HERE / f"cod_sample_{len(old)}.sha256"
    if prev.exists():
        now = dict(reversed(x.split("  ")) for x in lines)
        changed = [x.split("  ")[1] for x in prev.read_text().splitlines() if now.get(x.split("  ")[1]) != x.split("  ")[0]]
        assert not changed, f"files of the previous sample changed: {changed}"
    print(f"sample of {n} structures frozen: cod_sample.txt, cod_sample.sha256; previous {len(old)} kept; cod_cache has {n} CIFs")


if __name__ == "__main__":
    main(int(sys.argv[1]), Path(sys.argv[2]))
