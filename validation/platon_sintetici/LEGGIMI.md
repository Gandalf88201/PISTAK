# PLATON sui tre cristalli sintetici del controllo Monte Carlo

I cristalli sono quelli di `../check_incertezze.py` (controllo `mc`) e di `../../paper/scripts/mc_data.py`, ricostruiti
con lo stesso generatore (seme 7, stesso ordine): il testo CIF è identico.

| File | Cristallo |
|---|---|
| `p1_benzeni.cif` | *P*1, due benzeni indipendenti a 6° |
| `pm1_benzene_inversione.cif` | *P*1̄, un benzene e la sua immagine per inversione |
| `pm1_naftalene.cif` | *P*1̄, naftalene sul centro di inversione (*Z*′ = ½), coppia con la copia traslata lungo *a* |

Per ognuno ci sono:

- `<nome>.cif`: il CIF del controllo Monte Carlo (senza idrogeni), letto da PISTAK;
- `<nome>_H.cif`: lo stesso CIF con gli H in posizione calcolata, senza s.u. Senza idrogeni PLATON non riconosce un
  anello di C con due soli legami: lo elenca come piano («Plan») e salta la tabella Cg–Cg. Gli H non entrano nei piani
  né nei centroidi degli anelli, quindi nessuna grandezza degli anelli cambia;
- `<nome>_H.lis`: il listato di PLATON (nel repository con il permesso di A. L. Spek; si rifà come descritto sotto);
- `<nome>_su<k>.cif`, `<nome>_su<k>_H.cif`, `<nome>_su<k>_H.lis`: lo stesso cristallo con tutte le s.u. delle coordinate
  moltiplicate per *k*. PLATON stampa una s.u. con una sola cifra significativa, salvo quando la cifra è 1: 0,3° può
  valere da 0,25° a 0,35°. Con *k* scelto in modo che la s.u. di α (o di Cg···Cg) cada fra 10 e 19 unità dell'ultima
  cifra, PLATON la stampa con due cifre (±3 %). A questo livello tutto è lineare nelle s.u., quindi i rapporti non
  dipendono da *k*.

## Come rifarli

Dalla radice del repository (PLATON compilato dal sorgente, V-200626):

```bash
.venv/bin/python paper/scripts/platon_synthetic.py --platon /percorso/di/platon 4000 \
    > paper/data/runs/rev3/platon_synthetic.txt
```

Lo script scrive i CIF, esegue `printf 'CALC ALL\nEXIT\n' | platon -o <nome>_H.cif` in una cartella temporanea, copia
qui il `.lis`, e confronta le s.u. della tabella *Analysis of Short Ring-Interactions* con il Monte Carlo di PISTAK
(*N* = 4000, seme 1, senza le s.u. della cella) e con una propagazione lineare scritta in numpy. Senza `--platon` legge i
listati già presenti. Il risultato è in `paper/data/platon_sigma_synthetic.csv`.

## Risultato

Eseguito il 9 novembre 2026 con PLATON V-200626 (`IGBL(4) = 200626` nel sorgente). La s.u. di α stampata da PLATON è
(σ_t,A² + σ_t,B²)^½, con σ_t l'inclinazione quadratica media della normale (entrambe le componenti), anche quando i due
anelli sono la stessa molecola e una sua immagine di simmetria; vedi `paper/revision/round3/platon.md`.
