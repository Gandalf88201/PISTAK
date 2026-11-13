# Riferimento PLATON per le sezioni 4–5 della checklist

I 20 CIF di `elenco.txt` vengono dal campione COD (`../cod_sample.txt`). Sono strutture organiche senza metalli,
senza avvisi A, con coppie di anelli entro 6 Å e almeno uno stack. Le copie dei CIF si rigenerano da `../cod_cache/`;
gli output di PLATON (`*.lis`) sono il riferimento esterno e sono nel repository con il permesso di A. L. Spek
(13 novembre 2026); si rifanno come descritto sotto.

## Cosa serve

Per ogni CIF, il file `.lis` di `CALC ALL`. Contiene la tabella *Analysis of Short Ring-Interactions with Cg-Cg
Distances < 6.0 Ångstrom* (Cg–Cg, α, β, γ, distanze perpendicolari, slippage) da confrontare con
`pistak/core/ringpairs.py`.

## Come ottenerlo

Con PLATON installato e nel PATH, da questa cartella:

```bash
for f in *.cif; do printf 'CALC ALL\nEXIT\n' | platon -o "$f"; done
```

`-o` evita la finestra grafica; i comandi arrivano da standard input. Se la tua versione non accetta l'input così,
basta aprire ogni CIF in PLATON e dare `CALC ALL` dal prompt (o dal pulsante nel menu principale): il file
`nome.lis` compare accanto al CIF.

Nessuna opzione da cambiare: i valori predefiniti di PLATON (Cg–Cg < 6,0 Å, β < 60°) sono quelli da confrontare.

## Risultato

Eseguito il 4 novembre 2026 con PLATON V-200626 (`printf 'CALC ALL\nEXIT\n' | platon -o nome.cif`); confronto in
`confronto.txt`, prodotto da `../check_platon.py`, e riassunto in `../REGISTRO.md`. PLATON non legge 2205356 (tabella
anisotropa prima degli atomi) e su 7242996 si ferma su ADDSYM.

Il rapporto PDF del confronto, con dati e tabelle, è `PISTAK_vs_PLATON.pdf` (in inglese `PISTAK_vs_PLATON_en.pdf`), prodotto da
`../rapporto_platon.py`; tutte le coppie abbinate sono in `coppie_abbinate.csv`.
