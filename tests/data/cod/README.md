# Strutture reali per i test

CIF della [Crystallography Open Database](https://www.crystallography.net/cod/), che li distribuisce come dominio
pubblico (CC0). Il nome del file è il codice COD; la pagina di origine è `https://www.crystallography.net/cod/<codice>.html`.

Rispetto ai file scaricati sono stati tolti solo i campi di testo `_shelx_res_file`, `_shelx_hkl_file` e
`_shelx_fab_file` e le loro checksum (il file .res e le riflessioni incorporati), che PISTAK non legge.

| Codice | Composto | Gruppo | T | Perché è qui |
|---|---|---|---|---|
| 2311088 | naftalene | P2₁/c, Z′ = ½ | 293 K | molecola su un centro di inversione; GoF dichiarato 0,25 (avviso B) |
| 2012157 | pentacene | P1̄, Z′ = 1 | 293 K | spina di pesce tipica dei semiconduttori organici |
| 2100417 | rubrene ortorombico | Cmca, Z′ = ¼ | 293 K | simbolo nella notazione vecchia (Cmca = Cmce); molecola su un sito di simmetria mm |
| 4500022 | C₁₄H₁₂N₆ · 2 H₂O | P2₁/n, Z′ = 3 | 100 K | tre molecole indipendenti, acqua, legami H; tabelle `_geom_` complete con le s.u. |
| 1550133 | C₂₄H₃₂O₈ | P2₁, Z′ = 1 | 100 K | anelli tutti saturi: nessun nucleo π, nessuna coppia |
| 7024607 | [Mn₂(C₁₄H₁₃N₅O₉)]ₙ | P1̄, rete 3D | 293 K | polimero di coordinazione: le coppie π si cercano tra frammenti dell'unità ripetuta; una piridina restava spezzata tra due celle (registro, voce 47) |
| 2210645 | (C₂₅H₂₁NO₂P)₂[Ni(C₄N₂S₂)₂] · 2 H₂O | P1̄, Z = 1 | 293 K | i due anelli chelati NiS₂C₂, complanari e con il Ni in comune, non sono un π-stack (frammenti legati tra loro) |
