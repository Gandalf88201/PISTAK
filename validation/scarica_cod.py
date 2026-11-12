"""Campione casuale di CIF della COD (dominio pubblico) per le verifiche della sezione 1 della checklist.

Gli ID si estraggono con un seme fisso dalle serie 1xxxxxx, 2xxxxxx, 4xxxxxx e 7xxxxxx (strutture organiche e
metallo-organiche); quelli che non esistono si saltano. L'elenco finale va in `cod_sample.txt` (nel repository), i
file in `cod_cache/` (fuori dal repository). Con `cod_sample.txt` già presente si scaricano solo quegli ID.

    .venv/bin/python validation/scarica_cod.py [quanti] [elenco]

Con un secondo argomento l'elenco va in quel file invece di `cod_sample.txt` (revisione 4, M7.1: il campione di 2000
strutture, `cod_sample_2000.txt`, estratto con lo stesso seme: gli ID di `cod_sample.txt` sono i suoi primi 200
trovati, purché esistano ancora nella COD).
"""
from __future__ import annotations

import http.client
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
    data = b""
    for attempt in range(6):
        try:
            with urllib.request.urlopen(URL.format(cod_id), timeout=30) as r:
                data = r.read()
            break
        except urllib.error.HTTPError:
            return False  # l'ID non esiste
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError, http.client.HTTPException):
            # errore di rete (connessione chiusa dal server): non è «l'ID non esiste», si riprova con calma; se si
            # rinunciasse, l'estrazione con seme fisso darebbe un campione diverso a ogni esecuzione
            time.sleep(5 * (attempt + 1))
    else:
        raise RuntimeError(f"COD {cod_id}: rete non raggiungibile dopo 6 tentativi; rilancia lo script (riparte dalla cache)")
    if b"_atom_site_fract_x" not in data:
        return False
    path.write_bytes(data)
    time.sleep(0.2)  # una richiesta alla volta, con calma
    return True


def main(n: int, sample: Path = SAMPLE) -> None:
    global CACHE
    if sample != SAMPLE:
        # un altro campione ha la sua cartella: gli script di verifica prendono «il campione» dall'elenco dei file di
        # cod_cache, e file in più lì dentro lo cambierebbero (revisione 4: è successo durante lo scaricamento)
        CACHE = HERE / f"cod_cache_{sample.stem.rsplit('_', 1)[-1]}"
    CACHE.mkdir(exist_ok=True)
    if sample.exists():
        ids = [int(x) for x in sample.read_text().split()]
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
    sample.write_text("\n".join(map(str, sorted(got))) + "\n")
    print(f"{len(got)} CIF su {tried} tentativi; elenco in {sample}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200, HERE / sys.argv[2] if len(sys.argv) > 2 else SAMPLE)
