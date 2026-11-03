"""Campione casuale di CIF della COD (dominio pubblico) per le verifiche della sezione 1 della checklist.

Gli ID si estraggono con un seme fisso dalle serie 1xxxxxx, 2xxxxxx, 4xxxxxx e 7xxxxxx (strutture organiche e
metallo-organiche); quelli che non esistono si saltano. L'elenco finale va in `cod_sample.txt` (nel repository), i
file in `cod_cache/` (fuori dal repository). Con `cod_sample.txt` già presente si scaricano solo quegli ID.

    .venv/bin/python validation/scarica_cod.py [quanti]
"""
from __future__ import annotations

import random
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cod_cache"
SAMPLE = HERE / "cod_sample.txt"
RANGES = [(1000000, 1561000), (2000000, 2242000), (4000000, 4346000), (7000000, 7246000)]
URL = "https://www.crystallography.net/cod/{}.cif"


def fetch(cod_id: int) -> bool:
    path = CACHE / f"{cod_id}.cif"
    if path.exists():
        return True
    try:
        with urllib.request.urlopen(URL.format(cod_id), timeout=30) as r:
            data = r.read()
    except urllib.error.HTTPError:
        return False
    if b"_atom_site_fract_x" not in data:
        return False
    path.write_bytes(data)
    time.sleep(0.2)  # una richiesta alla volta, con calma
    return True


def main(n: int) -> None:
    CACHE.mkdir(exist_ok=True)
    if SAMPLE.exists():
        ids = [int(x) for x in SAMPLE.read_text().split()]
        missing = [i for i in ids if not fetch(i)]
        print(f"{len(ids) - len(missing)}/{len(ids)} CIF in {CACHE}" + (f"; non trovati: {missing}" if missing else ""))
        return
    rng = random.Random(2026)
    got: list[int] = []
    tried = 0
    while len(got) < n and tried < 20 * n:
        lo, hi = RANGES[tried % len(RANGES)]
        cod_id = rng.randrange(lo, hi)
        tried += 1
        if cod_id not in got and fetch(cod_id):
            got.append(cod_id)
            print(f"{len(got)}/{n} {cod_id}", flush=True)
    SAMPLE.write_text("\n".join(map(str, sorted(got))) + "\n")
    print(f"{len(got)} CIF su {tried} tentativi; elenco in {SAMPLE}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200)
