"""Checklist, sezione 9: riproducibilità.

    batch   le stesse 40 strutture COD analizzate in un processo e in quattro, e due volte in quattro: risultati
            identici (l'ordine di completamento dei processi non entra nei risultati)
    mc      il Monte Carlo con lo stesso seme dà le stesse σ; con un seme diverso, σ diverse

L'ambiente in cui sono stati eseguiti i controlli è in `validation/ambiente.txt`.

    .venv/bin/python validation/check_riproducibilita.py [batch mc]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

from pistak.core.analysis import analyse_file  # noqa: E402
from pistak.core.config import Limits  # noqa: E402
from pistak.core.uncertainty import monte_carlo  # noqa: E402
from pistak.io.batch import run_batch  # noqa: E402


def _dump(entries) -> str:
    rows = [{k: round(v, 12) if isinstance(v, float) else v for k, v in vars(e).items() if k != "path"} for e in entries]
    return json.dumps(rows, sort_keys=True, default=str)


def check_batch() -> None:
    print("== cartella: 1 processo contro 4 ==")
    files = sorted((HERE / "cod_cache").glob("[0-9]*.cif"))[:40]  # [0-9]: non i file ._ di macOS sui dischi exFAT
    one, four, again = (_dump(run_batch(files, Limits(), w)) for w in (1, 4, 4))
    print(f"  {len(files)} file; 1 = 4 processi: {one == four}; due esecuzioni a 4 uguali: {four == again}")


def check_mc() -> None:
    print("\n== Monte Carlo: seme ==")
    a = analyse_file(ROOT / "tests" / "data" / "cod" / "2100417.cif", Limits())
    u1, u2, u3 = (monte_carlo(a, 30, seed=s) for s in (7, 7, 8))
    # NaN (s_L, s_W non definiti, voce 55) non è uguale a sé stesso: si confrontano le rappresentazioni
    same = lambda x, y: json.dumps(x, sort_keys=True) == json.dumps(y, sort_keys=True)  # noqa: E731
    print(f"  stesso seme, stesse σ: {same(u1.sigma, u2.sigma)}; seme diverso, σ diverse: {not same(u1.sigma, u3.sigma)}")


CHECKS = {"batch": check_batch, "mc": check_mc}

if __name__ == "__main__":  # necessario: i processi del batch reimportano questo modulo
    for name in sys.argv[1:] or CHECKS:
        CHECKS[name]()
