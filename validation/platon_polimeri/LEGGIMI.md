# Riferimento PLATON per le reti estese (polimeri di coordinazione)

I 26 CIF di `elenco.txt` sono le strutture del campione COD (`../cod_sample.txt`) in cui la rete di legami non si
chiude e che hanno anelli: polimeri di coordinazione e MOF. Servono a verificare le coppie π tra frammenti dell'unità
ripetuta (registro, voce 47). Come per `../platon/`, gli output di PLATON (`*.lis`) sono il riferimento esterno.
Sono nel repository con il permesso di A. L. Spek (13 novembre 2026). Le copie dei CIF e i file `.spf` e `.ckf` che
PLATON scrive accanto non sono nel repository: si rigenerano come descritto sotto.

## Come si rifà

Con i CIF del campione in `../cod_cache/` (`validation/scarica_cod.py`) e PLATON nel PATH, da questa cartella:

```bash
while read -r id; do cp "../cod_cache/$id.cif" .; done < elenco.txt
for f in *.cif; do printf 'CALC ALL\nEXIT\n' | platon -o "$f"; done
```

e poi, dalla cartella del progetto:

```bash
.venv/bin/python validation/check_platon.py platon_polimeri > validation/platon_polimeri/confronto.txt
```

## Risultato (5 novembre 2026, PLATON V-200626)

PLATON non legge 6 CIF (1521307, 1538535, 1540198, 2108013, 4326811, 7022074: nessun anello nel `.lis`). Sui 20
restanti, 217 coppie di anelli abbinate; tolta la riga di 7117874 descritta sotto, lo scarto massimo PISTAK − PLATON è
0,0013 Å su Cg···Cg e sulle distanze perpendicolari, 0,10° su α, 0,07° su β e γ, 0,002 Å sullo slippage: la precisione
con cui PLATON stampa.

Classi di coppie che PLATON elenca e PISTAK no: 114 con un anello non planare per PISTAK o con un metallo (chelati,
metallacicli: non sono frammenti π), più la riga di 7117874.

**7117874 (P6₃, Ni).** Per Cg2→Cg2 [−x+y, 1−x, z] PLATON stampa α 54,1°, β 34,2°, γ 0,0°: valori incompatibili tra
loro (con β 34,2° e γ 0° l'angolo tra i piani non può superare 34,2°), con una normale trasformata non unitaria
(−1,33, −1,72, 0,56) e CgI_Perp 4,825 Å, che con Cg···Cg 5,462 Å darebbe γ 28°. Ricalcolata con numpy sugli atomi del
CIF (gemmi per la cella), la coppia ha α 88,1°, β 34,2°, γ 61,6°: i valori di PISTAK. È un difetto di PLATON con
l'operazione esagonale.
