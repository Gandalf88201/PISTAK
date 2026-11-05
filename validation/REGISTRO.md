# Registro delle verifiche

Registro della revisione del codice (checklist del 3 novembre 2026, sezione 10). Una riga per problema trovato o
verifica chiusa; sotto la tabella, per ogni voce, l'oracolo usato e il risultato. Le verifiche si rieseguono con

```bash
.venv/bin/python validation/check_registro.py            # voci del registro: cell rings beta height unit shells cutoff
.venv/bin/python -m pip install -r validation/requirements.txt   # gemmi, per la sezione 1
.venv/bin/python validation/check_simmetria.py           # sezione 1: simboli di Hall e H-M contro gemmi
.venv/bin/python validation/scarica_cod.py               # sezione 1: i 200 CIF COD di cod_sample.txt in cod_cache/
.venv/bin/python validation/check_cif.py                 # sezione 1: lettura, s.u., elementi, disordine
.venv/bin/python validation/check_struttura.py           # sezioni 1–3: cella, formule, cella obliqua, anelli, raggi
.venv/bin/python validation/check_platon.py              # sezioni 4–5: coppie di anelli contro PLATON CALC ALL
.venv/bin/python validation/check_incertezze.py          # sezione 6: Monte Carlo, geometria contro le tabelle _geom_
.venv/bin/python validation/check_export.py              # sezione 7: file per ORCA/Psi4/Gaussian e codici periodici
.venv/bin/python validation/check_snapshot.py            # sezione 8: istantanee COD rifatte con numpy
.venv/bin/python validation/check_struttura.py oblique origin   # sezione 8: invarianza per cella e origine
.venv/bin/python validation/check_riproducibilita.py    # sezione 9: batch 1 contro 4 processi, seme del Monte Carlo
```

Gli oracoli non usano `pistak.core`: geometria ricalcolata con numpy, topologia nota dei policicli, tensore metrico,
conti a mano.

| # | Modulo | Problema o verifica | Priorità | Stato | Azione |
| --- | --- | --- | --- | --- | --- |
| 1 | params/extra.py | `beta_b` etichettato come γ | P1 → P3 | **Corretto**: i numeri erano giusti, la chiave no | Chiavi `beta` e `gamma`; `beta_a`, `beta_b` valgono ancora in `Pair.get` e nel registro |
| 2 | params/extra.py | `mean_atom_height`: descrizione (piano di A) diversa dal codice (normale media) | P1 | **Corretto**: il vecchio codice dava la distanza interplanare sotto altro nome | Ora media di \|(x − Cg_A)·n_A\| sugli atomi di B; test in `tests/test_revisione.py` |
| 3 | core/cores.py | `small_rings` può perdere anelli centrali | P1 | **Corretto**: era confermato sul coronene; pirene, trifenilene, perilene erano già giusti | Per ogni legame tutti i cicli più corti; test su coronene ed esabenzocoronene in 300 ordinamenti |
| 4 | core/geometry.py | `unit()` restituisce il vettore nullo senza errore | P1 → P3 | Non raggiungibile: dimostrato e verificato | Nessuna sui numeri; facoltativo un `assert` |
| 5 | core/pairs.py | `MAX_SHELLS` = 8 taglia la ricerca senza avviso | P1 | **Corretto** (avviso) | Nuovo avviso B `shellsTruncated` nei controlli rapidi, ricalcolato anche da `reanalyse`; il limite resta |
| 6 | core/pairs.py | Coppie cercate entro 4,0 Å ma stack fino a 4,2 Å | P1 | **Corretto**: `contact_cutoff` = 4,5 Å (il modo del sito resta a 4,0) | Motivazione in `Limits`; test su colonne a 4,05–4,19 Å |
| 7 | io/qc.py | Distanze X–H per P e Si senza fonte; C–H unico per sp² e sp³ | P1 | **Corretto** sul testo di Allen & Bruno 2010 | Medie a T ≤ 60 K, C–H secondo l'ibridazione; P–H tolto (nessuna fonte) |
| 8 | repository | Manca LICENSE; manca `license` in pyproject.toml | P1 | **Corretto** | MIT: `LICENSE`, `license = "MIT"` in pyproject.toml, sezione nei due README |
| 9 | tests/cifgen.py | Oracolo circolare? usa `cell_vectors` di PISTAK | P1 | Verificato: `cell_vectors` è corretto per un oracolo indipendente, quindi `cifgen` non è circolare | Nessuna; la verifica resta in `check_registro.py` |
| 11 | core/spacegroups.py | **Nuovo.** Il suffisso dell'impostazione (`:1`, `:2`, `:R`) veniva scartato: `R -3 :R` riceveva le operazioni esagonali, `F d d d :1` e `I 41/a :1` l'origine 2; senza suffisso la tabella sceglieva l'origine 2 e gemmi la 1 | P1 (solo CIF senza operazioni elencate) | **Corretto** | Suffissi rispettati; senza suffisso sempre origine 2 (convenzione dichiarata); test |
| 12 | core/cif.py | **Nuovo.** `element_of`: `BR2` senza type_symbol diventava boro, `CL1` carbonio; centroidi `Cg1` letti come carbonio, picchi `Q1` tenuti come atomi; `_atom_site_calc_flag dum` ignorato | P2 | **Corretto** | CL, BR a due lettere; Cg, Q e dum scartati; altre etichette maiuscole (CA, CD, HG) a una lettera |
| 13 | core/validation.py | **Nuovo.** Operazioni elencate e simbolo del gruppo incoerenti: vincevano le operazioni, senza avviso | P2 | **Corretto** | Avviso C `symopsMismatch`; nessun falso positivo sui 200 CIF COD |
| 14 | core/cif.py | Disordine su un elemento di simmetria (gruppi negativi): le immagini simmetriche del gruppo si sovrappongono alle posizioni originali | P1 | **Aperto** | Oggi lo segnala l'avviso A `overlap` (COD 4113834: C32A···C32B 0,38 Å); serve togliere le immagini che coincidono con l'altra parte |
| 15 | core/crystal.py | **Nuovo.** `cell_contents` contava una volta sola una posizione occupata da due atomi diversi del CIF: un CIF che scrive anche le immagini per simmetria con occupanza ½ (COD 2100449) dava metà del contenuto della cella, e densità e formula sbagliate negli avvisi | P2 (avvisi, non coppie) | **Corretto** | Occupanze di atomi diversi sommate fino a 1; immagini dello stesso atomo (posizioni speciali) una volta; test |
| 16 | core/elements.py | Raggi diversi dalle tabelle di gemmi: covalente C 0,76 (sp³, dichiarato) contro 0,73; Lu 1,87 e Hf 1,75 scambiati rispetto a gemmi; vdW H 1,10 (Rowland e Taylor) mentre il commento cita Bondi (1,20) | P1 (H: contatti vdW) | **Corretto** (commento) | Valori giusti: i 96 covalenti uguali alla tabella di Cordero 2008 (le differenze erano di gemmi); i vdW uguali alla tab. 12 di Mantina 2009 per i gruppi principali (H 1,10 di Rowland e Taylor 1996), Br 1,85 e metalli di transizione di Bondi 1964. Corretti il commento in `elements.py` e il manuale, che citavano solo Bondi |
| 17 | core/structure.py | Formule ricostruite diverse da `_chemical_formula_moiety` non spiegate: 7208500 (due acetonitrili in un'unica "molecola"), 2201728 (un CH₃ in più) | P2 | **Aperto** | Da esaminare a mano |
| 10 | core/cif.py | **Nuovo.** `tokenize` va in ciclo infinito (memoria che cresce senza limite) su un carattere di spaziatura diverso da spazio e tab: NBSP, form feed, tab verticale, U+2003… | P2 | **Corretto** | `if ch.isspace()` al posto di `if ch in " \t"`; test con 5 caratteri in un processo separato |
| 18 | core/intramolecular.py | **Nuovo.** Torsioni riportate anche attraverso tre atomi allineati (X–M–X su un centro di inversione, angolo 180°): valore indefinito, che dipende dagli arrotondamenti | P2 | **Corretto** | Torsioni con un angolo di legame oltre 175° non riportate (`MAX_TORSION_ANGLE`); test, manuale |
| 19 | core/intramolecular.py | **Nuovo.** Misure con le stesse etichette considerate uguali per simmetria: L1–M–L2 e L1–M–L2′ (θ e 180° − θ attorno a un centro di inversione), torsioni cis e trans dello stesso frammento, ne restava una | P1 (geometria riportata) | **Corretto** | Uguali solo se legate da un'operazione del cristallo (posizioni relative, simmetria del sito del primo atomo); 2100417: una torsione in più; misure CIF senza corrispondente da 539 a 27 |
| 20 | core/validation.py | `_hydrogens`: falsi positivi su isonitrile (C terminale) e su composti senza H per natura (cianuri, carbonili metallici: avviso B) | P2 | **Corretto** | C terminale con triplo legame a N/O (≤ 1,18 Å) escluso; `noHydrogens` solo se qualche C sembra volerne; test |
| 21 | core/validation.py | ADDSYM su strutture disordinate: la componente tenuta da `major` aveva più simmetria del modello, avviso A falso (COD 1530955, 1531236; PLATON non trova nulla) | P2 | **Corretto** | Come PLATON, atomi con occupanza parziale esclusi; oltre metà esclusi, avviso G `addsymDisorder` |
| 22 | manuale, README | «pre-checkCIF» e livelli A/B/C/G fanno pensare agli alert ufficiali IUCr | P2 | **Corretto** | Detto nel manuale (it, en) e nei README che soglie e controlli sono di PISTAK |
| 23 | core/validation.py | `TRIPLE_CX` = 1,18 Å e `SP_ANGLE` = 150°: valori motivati (C≡N, C≡O 1,13–1,17 Å; C=O da 1,19 Å; sp 180°, sp² 120°) ma senza una tabella citata | P3 | **Aperto** | Controllare sulle tabelle di Allen et al. 1987 |
| 24 | io/qc.py | Molteplicità del dimero fissata all'accoppiamento ad alto spin: per due radicali l'energia dipende dalla scelta singoletto/tripletto | P1 (energie di radicali) | **Corretto** | Opzione `dimer_multiplicity` (API, `--dimer-mult`, finestra), controllata contro gli spin dei monomeri; il predefinito resta l'alto spin, scritto nel file |
| 25 | io/qc.py | `%maxcore` ORCA = memoria totale / processi, senza margine (ORCA a volte lo supera) | P2 | **Corretto** | 75 % della memoria per processo, come consigliano i centri di calcolo; MB per processo verificato |
| 26 | io/qc.py, io/periodic.py | Nei file non era scritto che i monomeri sono molecole intere, quanti H restano come nel CIF, che `ecutrho` = 8 × presuppone pseudopotenziali ultrasoft/PAW, che punti k, cutoff e Γ di CP2K sono valori di comodo | P2 | **Corretto** | Righe di commento in ogni file (CRYSTAL: nel `.notes.txt`); manuale |
| 27 | io/qc.py, io/periodic.py | Esecuzione reale con ORCA, Psi4, Gaussian, VASP, QE, CP2K, CRYSTAL | P1 | **Aperto** | I programmi non sono installati qui; servono una prova per programma (dimero di benzene, cella del naftalene) e il confronto CP ORCA contro Psi4 |
| 28 | tests/ | Tolleranze larghe (`abs=1e-3` Å, 0,5° e 0,01°) su cristalli sintetici costruiti a mano | P1 | **Corretto** | Scarti reali misurati: ≤ 5·10⁻⁶ Å e 5·10⁻⁵°; tolleranze portate a 10⁻⁴ Å e 10⁻³° (16 asserzioni in 5 file); quelle contro i CIF (densità, tabelle `_geom_`) restano legate alle cifre del CIF |
| 29 | tests/snapshots | Istantanee COD scritte dal programma stesso | P1 | **Verificato** | Le 95 coppie intermolecolari delle 5 istantanee rifatte con numpy: scarti ≤ 5·10⁻⁹ Å e 2·10⁻⁶° (precisione delle istantanee) |
| 30 | validation | Invarianza cristallografica | P1 | **Verificato** | Cella obliqua equivalente e origine spostata di (0,137; 0,291; 0,443): 40/40 strutture con le stesse coppie entro 10⁻⁶ Å (prima la verifica arrotondava a 3 decimali e segnalava differenze di arrotondamento) |
| 31 | configurazione utente | **Nuovo.** `analyse_file` senza soglie usa il file dell'utente (`user.toml`): su questa macchina `contact_cutoff` = 4,0 (il vecchio predefinito, salvato prima della voce 6), `max_dihedral` = 5, `overlap_measure` = area, `contact_check` = vdw. I test sono isolati (`PISTAK_CONFIG_DIR` temporanea), le analisi dell'utente no | P1 (risultati dell'utente) | **Aperto** | Decidere se aggiornare `user.toml` e se avvisare quando un valore salvato coincide con un vecchio predefinito |
| 32 | tests/ | Copertura (`pytest --cov`) e test basati su proprietà (hypothesis) | P2 | **Aperto** | Servono `pytest-cov` e `hypothesis`, non installati |
| 33 | tests/ | Nomi dei test in italiano | P3 | Aperto | Da decidere: docstring in inglese o traduzione |
| 34 | core/validation.py | **Nuovo.** Se spglib c'è ma non riesce, `find_symmetry` restituiva None come quando spglib manca: avviso G «serve il pacchetto spglib» (falso) a 0,1 Å, nessun avviso a 0,3 Å | P2 | **Corretto** | L'errore di spglib diventa l'avviso G `checkFailed` con il messaggio; test. Su 200 CIF COD nessun nuovo avviso |
| 35 | core/uncertainty.py | **Nuovo.** Un parametro NaN in parte dei campioni del Monte Carlo (un plugin che fallisce, una geometria degenere) dava una σ calcolata sui soli campioni rimasti, senza dirlo | P2 | **Corretto** | `Uncertainty.valid` e `incomplete`: avviso nel terminale, riga nel Dettaglio, campo `su_samples` nel JSON; test |
| 36 | core/elements.py | **Nuovo.** 48 elementi su 96 senza raggio di van der Waals (Fe, Co, Mn, Ru, Ti, lantanidi, …): contatti e parametri di vdW usavano 1,9 Å in silenzio | P2 | **Corretto** (avviso) | Avviso G `vdwDefault` con gli elementi; Po, At, Rn, Fr, Ra aggiunti dalla tab. 12 di Mantina 2009. Resta da decidere se adottare una tabella per i metalli (es. Alvarez 2013) |
| 37 | gui/flow_layout.py | **Nuovo.** «Adattato dall'esempio Flow Layout di Qt» senza nota di copyright né licenza (gli esempi di Qt for Python sono BSD-3-Clause, che chiede di conservarle) | P1 (licenza) | **Corretto** | Nota di copyright e testo BSD-3 nell'intestazione; README aggiornati. Da confermare sull'esempio originale (copyright Riverbank 2013, Qt Company 2022) |
| 38 | README | Licenze delle dipendenze: mancavano gemmi (MPL-2.0, facoltativo) e vtk (BSD, tramite pyvista) | P3 | **Corretto** | Aggiunte; le altre (numpy, scipy, spglib BSD; pyqtgraph, pyvista MIT; PySide6 LGPL) confermate dai metadati installati |
| 39 | requirements | Solo versioni minime (`>=`), nessun file di versioni bloccate: i numeri del paper non sono legati a un ambiente | P2 | **Corretto** (documentato) | `validation/ambiente.txt`: versioni esatte con cui girano test e controlli. Batch e Monte Carlo deterministici (`check_riproducibilita.py`) |
| 40 | pistak/ | Leggibilità: ruff pulito (nessun import o variabile inutilizzati; 11 argomenti non usati, tutti segnature di slot Qt). vulture non installato; il mypy di sistema (1.19) non legge le stub di numpy 2.5 | P3 | **Aperto** | Installare mypy e vulture nel venv (insieme alla voce 32) |
| 41 | paper | Dichiarazione sull'uso dell'IA | P1 (paper) | **Bozza** | `validation/dichiarazione_ia.md`: fatti dal repository (29 commit con Claude co-autore), testo IT/EN, parti da completare dall'autore |

## 1. `beta_b` etichettato γ

**Oracolo.** Definizioni di PLATON (`CALC ALL`, sezione Cg–Cg): β = angolo tra il vettore Cg(I)→Cg(J) e la normale
al piano I; γ = angolo tra lo stesso vettore e la normale al piano J. Normali ricalcolate con una SVD indipendente.

**Risultato.** Su 104 coppie (5 CIF COD e 6 esempi) lo scarto massimo è 2·10⁻¹³°. Il valore di `beta_b` è esattamente il γ di PLATON, e l'etichetta «γ (Cg→normale B)» è
giusta. Il difetto è solo nella chiave: chi legge il CSV trova una colonna `beta_b` che contiene γ. `ringpairs.py`
usa già `beta` e `gamma` in modo corretto, quindi i due moduli usano nomi diversi per la stessa grandezza.

**Proposta.** Chiave `gamma` (e `beta_a` → `beta`), con le vecchie chiavi accettate in lettura delle sessioni.
Cambia l'intestazione dei CSV: va rigenerato lo snapshot.

## 2. `mean_atom_height`

**Descrizione.** «Distanza media degli atomi del nucleo B dal piano medio di A».

**Codice.** `abs(mean(proj_b[atomi del nucleo B, 2]))`, cioè la media della coordinata lungo la normale **media** n
della coppia, con origine nel baricentro di A. Coincide con la descrizione solo se i due piani sono paralleli.

**Oracolo.** media di (x − Cg_A)·n_A sugli atomi del nucleo B, con n_A dalla SVD degli atomi del nucleo A.

**Risultato.** Su 104 coppie lo scarto massimo tra codice e descrizione è 0,67 Å per le coppie stack e parallel
(piani entro 20°) e 4,44 Å per le inclinate. Caso peggiore: COD 2100417, coppia 3, angolo tra i piani 80,8°: il codice
dà 0,055 Å, la distanza media dal piano di A è 4,498 Å. Non è un arrotondamento: sono due grandezze diverse.

Due differenze tra codice e descrizione:

1. normale media invece di normale di A (scarto che cresce con l'angolo tra i piani);
2. `abs(media)` invece di `media(abs)`: per un nucleo B che attraversa il piano di A (coppie inclinate) gli atomi
   sopra e sotto si compensano.

Insieme, le due cose rendevano il parametro identico alla distanza interplanare: la media con segno delle quote
degli atomi di B lungo n è la quota del baricentro di B, cioè d·n. **Corretto** con media(|(x − Cg_A)·n_A|), che è
la descrizione; la versione con abs(media) sarebbe stata `height_b_over_a`, già nel registro.

## 3. `small_rings` e l'anello centrale del coronene

**Oracolo.** Policicli benzenoidi costruiti sul reticolo a nido d'ape (`honeycomb` in `check_registro.py`): gli
anelli attesi sono gli esagoni, e il loro numero coincide con il rango ciclico E − V + 1. Si permuta l'ordine degli
atomi e dei legami, perché la ricerca in ampiezza dipende dall'ordine dei vicini.

**Risultato.**

| Molecola | Anelli attesi | Esito su 200 ordinamenti |
| --- | --- | --- |
| naftalene, fenantrene, pirene, trifenilene, perilene, benzo[ghi]perilene | 2–6 | sempre corretti |
| coronene | 7 | anello centrale perso in 8/200 (legami in ordine casuale), 121/2000 ≈ 6 % (legami ordinati come in `structure.py`) |

Il trifenilene e il pirene, citati nella checklist come casi critici, **non** sono affetti: il loro anello centrale
ha almeno un legame che non appartiene ad altri anelli. Il difetto riguarda gli anelli con **tutti** i legami
condivisi (coronene, e in generale anelli interni completamente circondati: esabenzocoronene, circumaceni).

**Conseguenze.** Il frammento π fuso non cambia (gli atomi ci sono comunque). Cambia la tabella anello–anello
(`ringpairs.py`): manca un centroide. Il risultato dipende dall'ordine degli atomi nel CIF, quindi due CIF della
stessa struttura con atomi in ordine diverso possono dare tabelle Cg–Cg diverse.

**Proposta.** Per ogni legame raccogliere *tutti* i cicli di lunghezza minima che lo contengono, non solo il primo
trovato (approssimazione dei «cicli rilevanti» di Vismara). Per i benzenoidi dà esattamente gli esagoni.

## 4. `unit()` e il vettore nullo

`unit()` si usa solo in `pair_frame`, due volte:

- n = unit(n_A + n_B'), con n_B' orientata in modo che n_A·n_B' ≥ 0: |n_A + n_B'|² = 2 + 2 n_A·n_B' ≥ 2, quindi
  |n_A + n_B'| ≥ √2;
- L = unit(l_A − n (l_A·n)): l_A ⟂ n_A, quindi l_A·n = l_A·n_B' / |n_A + n_B'| e |l_A·n| ≤ 1/√2; ne segue
  |l_A proiettato| ≥ 1/√2. Nel modo simmetrico l_B è orientato in modo che l_A·l_B ≥ 0, quindi |l_A + l_B| ≥ |l_A| e
  il ramo di riserva (`> 1e-6`) non scatta mai.

Verifica numerica su 20 000 coppie casuali, compreso il caso peggiore n_B = l_A: minimi 1,414 e 0,707, come previsto.
Il rischio descritto nella checklist (L = W = 0 con slittamenti nulli che sembrano validi) non si presenta nel codice
attuale. Resta un rischio solo per chi riusa `unit()` altrove.

## 5. `MAX_SHELLS`

**Oracolo.** Il limite teorico della docstring di `neighbour_shells`, |k| ≤ D·|a*| + 1 con D = 2·estensione +
`contact_cutoff`, calcolato a parte.

| Asse corto | Estensione molecola | Celle necessarie | Esplorate |
| --- | --- | --- | --- |
| 3,5 Å | 7 Å | ±6 | ±6 |
| 3,5 Å | 12 Å | ±9 | ±8, **senza avviso** |
| 3,4 Å | 15 Å | ±11 | ±8, **senza avviso** |

Riguarda molecole grandi (estensione > ~12 Å: discotici con catene, macrocicli) con un asse di stack corto. Il limite
è prudente (sovrastima le celle necessarie), quindi un taglio non significa per forza una coppia persa, ma oggi non
lo si può sapere.

## 6. `contact_cutoff` = 4,0 Å contro `max_interplanar` = 4,2 Å

Per due piani paralleli ogni distanza tra un atomo di A e uno di B è maggiore o uguale alla distanza interplanare.
Quindi, se la distanza interplanare supera 4,0 Å, il contatto più corto supera `contact_cutoff` e la coppia non
viene nemmeno misurata. **Per piani paralleli la soglia effettiva è 4,0 Å, non 4,2 Å**; la fascia 4,0–4,2 Å è
raggiungibile solo da coppie con angolo tra i piani diverso da zero.

**Verifica.** Colonne infinite di naftalene a piani esattamente paralleli, CIF scritto a mano (cella monoclina, quota
e slittamento imposti):

| Interplanare | Slittamento | Risultato |
| --- | --- | --- |
| 3,50 / 3,90 / 3,99 Å | 0 Å | stack, interplanare esatta |
| 3,99 Å | 1,5 Å | nessuna coppia (contatto più corto 4,0003 Å) |
| 4,05 / 4,10 / 4,19 Å | 0 e 1,5 Å | nessuna coppia |

Il caso 3,99 Å con slittamento 1,5 Å mostra che anche sotto 4,0 Å lo stack dipende dal contatto più corto, non dalla
distanza interplanare. Va deciso prima di produrre le tabelle: o `contact_cutoff` ≥ `max_interplanar` più un margine
per lo slittamento, o nel paper la soglia dichiarata è il contatto atomo–atomo di 4,0 Å.

## 7. Distanze X–H

Prima: `{C: 1.083, N: 1.009, O: 0.983, B: 1.19, S: 1.338, P: 1.42, Si: 1.48}`. Controllato sul testo di Allen & Bruno,
Acta Cryst. B66, 380 (2010):

| Legame | Prima | Articolo | Ora |
| --- | --- | --- | --- |
| C–H | 1,083 (un solo valore) | C(ar)–H 1,083 e C=Csp²–H 1,082 (tab. 2, tutte le T); Csp²–H 1,085 e Csp³–H 1,095 a T ≤ 60 K (conclusioni, tab. 4); CH₃ accorciato di ~0,01 Å dalla librazione; Csp–H 1,055 | sp² e aromatico 1,085; sp³ 1,095; sp 1,055 |
| N–H | 1,009 (valore CSD del 1987) | Z₂–N–H 1,015 (tab. 2); 1,020 a T ≤ 60 K; N⁺–H 1,040 | 1,020; N⁺ (4 vicini) 1,040 |
| O–H | 0,983 | Z–O–H 0,983 (tab. 2); ossidrili 0,989 a T ≤ 60 K; acidi carbossilici 1,018 | 0,989 (gli acidi non si distinguono) |
| B–H | 1,19 | terminale 1,185 (tab. 2) | 1,185 |
| S–H | 1,338 | 1,338(2), una sola struttura (NALCYS02) | 1,338 |
| Si–H | 1,48 | 1,506(2), una sola struttura (COQYUC01) | 1,506 |
| P–H | 1,42 | assente | tolto: quegli H restano come nel CIF |

Si usano i valori a bassa temperatura perché servono geometrie per calcoli quantistici, dove l'accorciamento dovuto
alla librazione è un artefatto. L'ibridazione del carbonio si ricava dal numero di vicini (4 = sp³, 3 = sp² o
aromatico, 2 = sp); i due casi richiesti, aromatico e sp³, danno 1,085 e 1,095 Å. Lo stesso vale per le celle
periodiche (`io/periodic.py` usa la stessa funzione), quindi dimeri e celle restano coerenti.

## 9. `tests/cifgen.py`

`make_cif` passa da coordinate cartesiane a frazionarie con `cell_vectors` di PISTAK, la stessa funzione che poi
rilegge il CIF. Se `cell_vectors` fosse sbagliata, scrittura e lettura si compenserebbero e i test sintetici
passerebbero comunque. **Oracolo.** Tensore metrico G_ij = a_i·a_j dai sei parametri di cella e volume dalla formula
chiusa, su 5000 celle casuali: max |M·Mᵀ − G| = 2·10⁻¹³ Å², errore relativo sul volume 2·10⁻¹⁵, a lungo x e b nel
piano xy. `cell_vectors` è corretta, quindi i CIF di `cifgen` hanno la geometria dichiarata. `pistak/core/demo.py`
(caso 3,4 / 1,5 / 0,8 Å) costruisce la cella con numpy senza codice di PISTAK: oracolo indipendente.

## 10. `tokenize` e la spaziatura non ASCII

Trovato durante la revisione: sul disco esterno (exFAT) macOS crea un file `._nome.cif` per ogni CIF, con byte come
0x0B, 0x1C, 0x1E. Leggendone uno l'analisi non finisce più. La causa è in `tokenize`:

```python
if ch in " \t":          # salta solo spazio e tab
    ...
while k < n and not line[k].isspace():   # ma la parola finisce su qualunque spazio Unicode
```

Con un carattere che `str.isspace()` considera spazio ma che non è spazio o tab, la parola è vuota, `j` non avanza e
il ciclo aggiunge `""` alla lista per sempre. Riproduzione: `tokenize("_cell_length_a\xa05.0")` (NBSP). Il batch
salta i file che iniziano con un punto, quindi i `._*.cif` non entrano in un'analisi di cartella; ma un NBSP o un form
feed in un CIF passato a mano o copiato da una pagina web blocca la GUI. Non cambia nessun numero pubblicato.

## Correzioni (3 novembre 2026)

| Voce | File | Cosa cambia nei risultati |
| --- | --- | --- |
| 2 | `pistak/params/extra.py` | `mean_atom_height` cambia per tutte le coppie. Prima valeva esattamente la distanza interplanare: la media con segno delle quote degli atomi di B lungo n è la quota del loro baricentro. Ora è la distanza media degli atomi del nucleo B dal piano medio di A |
| 3 | `pistak/core/cores.py` | Molecole con anelli interni completamente circondati (coronene e simili): la tabella anello–anello ha sempre tutti gli anelli, in qualunque ordine siano gli atomi del CIF |
| 5 | `pistak/core/pairs.py`, `validation.py`, `api.py` | Nessun numero; un avviso B quando le celle esplorate sono meno di quelle necessarie (`shells_needed`) |
| 10 | `pistak/core/cif.py` | Nessun numero; un CIF con NBSP, form feed o simili ora si legge invece di bloccare il programma |

Snapshot rigenerati e confrontati campo per campo con quelli precedenti (14 file). Sono cambiati solo
`mean_atom_height` (121 valori) e la tabella anello–anello. Confrontate come insiemi, le 513 coppie di anelli hanno
gli stessi valori: cambiano l'ordine e la numerazione degli anelli, più un rumore numerico di 10⁻⁶° su α (COD 2100417).
Suite completa: 376 test superati, 3 saltati (servono i CIF di `DATI_REALI/`).

## Decisioni dell'autore (3 novembre 2026)

| Voce | Decisione | Effetto sui risultati |
| --- | --- | --- |
| 6 | `contact_cutoff` = 4,5 Å | Con piani paralleli il contatto non scende sotto la distanza interplanare, e un atomo di B dentro il contorno di A dista in proiezione al più ~1,4 Å (raggio di un anello) dall'atomo di A più vicino: √(4,2² + 1,4²) ≈ 4,43 Å < 4,5. Negli snapshot canonici: nessuna coppia persa, nessun tipo cambiato, nessuno stack nuovo; 4 coppie parallel in più nel naftalene sintetico, 10 coppie in più (9 inclinate) tra frammenti secondari in COD 2100417 e 4500022, perché la loro portata `contact_cutoff` + 1 Å passa da 5,0 a 5,5 Å. Il modo del sito resta a 4,0 Å |
| 1 | Chiavi `beta`, `gamma` | Cambiano le intestazioni dei CSV; i valori no |
| 8 | Licenza MIT | Nessuno |
| 7 | X–H da Allen & Bruno 2010, C per ibridazione | Solo gli input per ORCA, Psi4, Gaussian e i file periodici |

## Sezione 1: lettura CIF e simmetria (3 novembre 2026)

Oracolo: gemmi 0.7.5 (`validation/requirements.txt`). Campione: 200 CIF della COD estratti a caso (seme fisso) dalle
serie 1, 2, 4, 7 (`validation/cod_sample.txt`; i file in `validation/cod_cache/`, fuori dal repository).

| Verifica | Risultato |
| --- | --- |
| `hall_operations` su tutte le impostazioni di gemmi | 564/564 uguali (insiemi di R e t modulo il reticolo) |
| `HM_TO_HALL` contro gemmi | 97/97 uguali (dopo la correzione 11; prima 3 diverse: F d d d, I 41/a, F d -3 m davano l'origine 2 a chi chiedeva l'1) |
| `operations(hm=…)` con 43 scritture (complete, compatte, :1/:2/:H/:R) | 43/43 uguali (prima: `R -3 :R`, `F d d d :1`, `I 41/a :1` sbagliati) |
| `parse_cif` contro gemmi sui 200 CIF | 200/200: cella, numero di siti, coordinate, occupanze, elemento di ogni sito, operazioni |
| `cif_su` su 12 scritture calcolate a mano (`1.2345(12)`, `123(4)`, `.5(1)`, `1.2e-3(2)`…) | 12/12 |
| `element_of` su 25 etichette ambigue | tutte come attese dopo la correzione 12 |
| Tokenizzatore contro gemmi: CRLF, tab, `#` e apici dentro le stringhe, campi `;`, `data_` dentro un testo, blocchi multipli, save frame | uguale; solo l'a capo iniziale dei campi `;` viene tolto (innocuo) |
| Regola `major` ricalcolata indipendentemente dalle colonne lette da gemmi | 41 CIF con disordine, 0 differenze; esempi stampati da `check_cif.py disorder` per il controllo a mano |
| `symopsMismatch` sui 200 CIF | 0 avvisi (tutti con operazioni elencate e coerenti) |

Note:

- 5 dei 200 CIF non hanno `_atom_site_type_symbol`: l'elemento viene dall'etichetta.
- COD 7236705 mette tutti i 65 atomi nel gruppo di disordine 4: con `disorder = "first"` (il modo del sito) non resta
  nessun atomo e la lettura si ferma con `noAtoms`; con `major` (predefinito) si leggono tutti. Comportamento del modo
  di compatibilità, lasciato com'è.
- I gruppi con un solo membro in un assemblaggio (solventi a occupanza parziale, 1513296: occupanze 0,1–0,5) si
  tengono tutti: non c'è un'alternativa da scegliere.

## Sezioni 1–3: cella, molecole, anelli (3 novembre 2026)

Script: `validation/check_struttura.py`; stesso campione di 200 CIF COD. Oracoli: gemmi (`get_all_unit_cell_sites`),
i valori dichiarati nei CIF, networkx (`minimum_cycle_basis`), una cella obliqua costruita a mano.

| Verifica | Risultato |
| --- | --- |
| Atomi della cella (`fill_cell`) contro gemmi, CIF senza occupanze parziali | 132/132 uguali, elemento per elemento |
| Contenuto della cella pesato con l'occupanza (`cell_contents`) contro gemmi | prima 192/200, dopo la correzione 15 197/200; nei 3 rimasti 2 volte PISTAK coincide con Z × formula e gemmi no, 1 (4117065, Br disordinati a meno di 0,5 Å) nessuno dei due |
| Densità contro `_exptl_crystal_density_diffrn` (182 CIF) | scarto mediano 0,021 %; 21 oltre 0,5 %, tutti con contenuto della cella uguale a gemmi: solvente tolto con SQUEEZE ma dichiarato, formula o Z del CIF diversi dal modello |
| Molecole ricostruite contro `_chemical_formula_moiety` (114 CIF) | 74 uguali, 4 diverse solo negli H; delle 36 diverse, quasi tutte per convenzione del CIF (formula per unità asimmetrica o per coppia di ioni, controione tenuto separato da PISTAK) o per reti estese (polimeri di coordinazione); da disordine non risolto 4064128 (gruppi negativi, voce 14), 7036933 e 4330346 (disordine senza gruppi nel CIF: solo l'avviso `partialOccupancy`); voce 17 aperta |
| Cella obliqua equivalente (a' = a, b' = b + 2a, c' = c + a + b), 40 strutture | 39 identiche; in 1514883 due coppie differiscono di 0,001 Å (arrotondamento al terzo decimale). Nessun legame perso: la ricerca a ±1 cella regge anche celle non ridotte di questo tipo |
| `small_rings` contro la base minima dei cicli (287 molecole distinte) | trovati tutti gli anelli della base tranne 2 anelli da 5 in un cluster Ag–Pt (anelli i cui legami stanno tutti in cicli da 3 col metallo); gli anelli in più sono cicli più corti di qualche legame, quasi tutti triangoli tra i Br disordinati di 4117065 |

Il raddoppio dei solventi in 7036933 e 4330346 non è un difetto del programma: il CIF descrive un disordine (atomi a
occupanza parziale vicini fra loro o alle loro immagini) senza dire quali atomi stanno insieme. PISTAK lo segnala con
`partialOccupancy`; per escludere queste molecole dalle coppie serve un criterio nuovo (voce 14).

## Sezioni 4–5: coppie di anelli contro PLATON (4 novembre 2026)

Oracolo: PLATON V-200626 (A. L. Spek), compilato dal sorgente, `CALC ALL` con le impostazioni predefinite su 20 CIF
COD organici con stack (`validation/platon/elenco.txt`; output `.lis` nella stessa cartella, confronto completo in
`confronto.txt`). PLATON non ha letto 2 CIF: 2205356 (tabella anisotropa prima di quella degli atomi) e 7242996 (si
ferma su ADDSYM). Restano 18 strutture.

| Grandezza | Coppie abbinate | Scarto massimo PISTAK − PLATON | Mediano |
| --- | --- | --- | --- |
| Cg···Cg | 667 | 0,0011 Å | 0,0001 Å |
| α | 667 | 0,052° | 0,004° |
| β, γ | 667 | 0,056° | 0,03° |
| distanze perpendicolari (CgJ_Perp = `perp_a`, CgI_Perp = `perp_b`) | 667 | 0,0009 Å | 0,0001 Å |
| slippage | 339 (PLATON lo dà solo con α < 20°) | 0,0011 Å | 0,0002 Å |

Gli scarti stanno nella precisione con cui PLATON stampa (0,1° per gli angoli, 0,001 Å per le distanze): **i parametri
anello–anello di PISTAK coincidono con quelli di PLATON**. Le definizioni di PLATON stampate nel `.lis` confermano
quelle di PISTAK: β rispetto alla normale del piano I, γ a quella del piano J (voce 1), slippage come proiezione di
Cg(J) sul piano I.

Coppie che PLATON elenca e PISTAK no (classi di simmetria, 279):

- 249: un anello che PISTAK non considera. Quasi tutti sono anelli da 8 atomi che sono il perimetro di due anelli
  condensati (5 + 5, 4 + 6): PLATON li elenca come anelli, PISTAK per costruzione no (solo gli anelli più corti di
  qualche legame). Gli altri sono anelli attraverso un elemento di simmetria (PLATON chiama `C1_a` gli atomi generati,
  quindi l'abbinamento per etichetta non li trova) o di una componente di disordine scartata da `major`.
- 30: molecole senza atomi pesanti entro `contact_cutoff` = 4,5 Å; PLATON guarda solo Cg···Cg < 6 Å. È una scelta di
  PISTAK, da dichiarare: la tabella anello–anello copre solo le coppie di molecole a contatto.

PLATON elenca soltanto le coppie con β < 60° (in almeno un verso); PISTAK le elenca tutte entro 6 Å, per questo ha
più righe (725 contro le 667 abbinate).

### Distanza interplanare: PISTAK contro le distanze perpendicolari di PLATON

PISTAK misura la distanza interplanare lungo la normale media dei due nuclei; PLATON dà le due distanze perpendicolari
(il baricentro di ciascun piano dall'altro piano). Sulle 1381 coppie stack o parallel dei 200 CIF COD:

| Angolo tra i piani | Coppie | |⟂A − ⟂B| mediano (max) | interplanare PISTAK − media(⟂A, ⟂B), mediano (max) |
| --- | --- | --- | --- |
| 0–5° | 1096 | 0,000 (1,10) Å | 0,000 (0,12) Å |
| 5–10° | 132 | 0,54 (1,94) Å | 0,007 (0,36) Å |
| 10–15° | 88 | 0,72 (3,79) Å | 0,020 (0,75) Å |
| 15–20° | 65 | 1,10 (4,60) Å | 0,033 (0,94) Å |

Per il paper: con piani non paralleli le due distanze «alla PLATON» possono differire di qualche ångström, mentre la
distanza lungo la normale media è un solo numero, vicino alla loro media. Va scritto che è questa la definizione, con
la formula |d · n̄|, n̄ = (n_A + n_B)/|n_A + n_B|.

### Rapporto PDF

Il confronto, con dati e tabelle, è documentato in `validation/platon/PISTAK_vs_PLATON.pdf` (in inglese
`validation/platon/PISTAK_vs_PLATON_en.pdf`), prodotto da
`validation/rapporto_platon.py` (tutte le coppie abbinate anche in `validation/platon/coppie_abbinate.csv`; i dati
della distanza interplanare in `validation/platon/interplanare.json`, ricalcolabili con `--rifai` dai CIF di
`cod_cache/`). Le 249 classi con un anello che PISTAK non considera si dividono così: 205 anelli perimetrali da 8 atomi,
40 anelli con atomi generati per simmetria (etichette `_a` di PLATON), 4 anelli di una componente di disordine
scartata (tutti in 2204537, atomi con suffisso A). Rispetto alla prima stesura della tabella qui sopra, due massimi
erano arrotondati per eccesso (0,7546 e 0,9446 Å, ora 0,75 e 0,94).

## Sezione 6: incertezze, controlli rapidi, geometria intramolecolare (5 novembre 2026)

Script: `validation/check_incertezze.py`. Oracoli: propagazione lineare e geometria della coppia ricalcolate con numpy su
cristalli scritti a mano come testo CIF; tabelle `_geom_` dei CIF COD lette con gemmi; PLATON ADDSYM (V-200626).

### Monte Carlo

| Verifica | Risultato |
| --- | --- |
| σ Monte Carlo (N = 4000) contro propagazione lineare, P1 con due benzeni | scarti −0,7…+1,3 % (errore statistico atteso ±1,1 %) |
| idem, benzene e sua immagine per inversione (P-1) | −0,2…+2,8 %; angolo tra i piani nullo in ogni campione (σ = 0), come deve |
| idem, naftalene sul centro di inversione, Z′ = ½ (P-1) | −0,5 %: gli spostamenti delle due metà sono coerenti |
| σ in funzione di N (50…1000, 5 semi, 3 coppie COD) | dispersione relativa della σ stimata come 1/√(2(N−1)): 10 % a N = 50, 5 % a 200, 3 % a 500, 2 % a 1000 |
| media − valore nominale (bias), N = 1000 | entro 0,02 σ per distanza interplanare, Cg···Cg, slittamento e angolo tra piani non nullo |
| σ degli slittamenti L/W, nucleo isotropo (benzene) contro allungato (naftalene) | 0,48 Å contro 0,003–0,005 Å con le stesse s.u.: σ alta su L/W indica davvero un asse non definito, come dice il manuale |
| `kind_agreement` su una coppia costruita a 20° = `max_dihedral` | 0,503 (atteso 0,5 ± 0,011) |

**Scelta di N per il paper (da decidere).** Con N = 500 la σ ha un'incertezza relativa del 3 %, abbondante per
riportarla con una o due cifre; costa circa 0,04 s per campione per una struttura con 50 coppie. N = 200 (5 %) basta
per una cifra significativa.

Il bias è riportato nei campi `mean` del Monte Carlo ma non nell'esportazione: per gli angoli tra piani quasi
paralleli la media è spostata in alto (l'angolo è un valore assoluto), e lì conviene riportare valore(σ) senza
correggere.

### Geometria contro le tabelle `_geom_` (182 CIF COD)

| Misure con s.u. nel CIF | n | \|Δ\| mediano | oltre 3 s.u. |
| --- | --- | --- | --- |
| legami | 9840 | 0,0004 Å | 1 (riga malformata in 7206203: il codice di simmetria è nella colonna sbagliata) |
| angoli | 18353 | 0,022° | 0 |
| torsioni | 7438 | 0,029° | 0; nessun segno opposto (IUPAC/Klyne–Prelog, come SHELXL) |

Il 63–66 % degli scarti sta entro mezza unità dell'ultima cifra scritta nel CIF; il resto è l'arrotondamento delle
coordinate del CIF rispetto a quelle del raffinamento. La checklist chiedeva 0,001 Å e 0,1°: si rispettano per l'83 %
dei legami e il 91–94 % di angoli e torsioni, e lo scarto è sempre entro 3 s.u.

**s.u.** PISTAK/CIF: mediana 1,06 (legami), 1,01 (angoli), 1,03 (torsioni), 5–95 % tra 0,84 e 1,56. PISTAK propaga
solo la diagonale delle s.u. delle coordinate; il programma di raffinamento usa la covarianza completa e le s.u. di
cella. Lo scarto atteso, quindi il limite da dichiarare, è di questo ordine: in genere entro il ±15 %, fino al 50 % in
qualche legame.

**Codici di simmetria.** 20954 misure con atomi trasformati rifatte dalle coordinate del CIF con le operazioni scritte
da PISTAK: tutte uguali. Il formato (`C2#1` con legenda `#1 = -x+1,-y,-z`) è quello delle tabelle di SHELXL. Delle
misure del CIF con atomi trasformati, 2192 hanno in PISTAK la stessa relazione di simmetria, 813 lo stesso valore su
un'altra copia equivalente, 27 nessun corrispondente: in tutte c'è un legame a un atomo trasformato che il programma
di raffinamento considera tale e PISTAK no (Cu–Cl a 3,0 Å, Na–O, Ag–O, un atomo di disordine).

### Controlli rapidi

* **Idrogeni.** Sui 200 CIF COD `missingH` scatta in 6. Su cristalli costruiti a mano davano falsi positivi
  l'isonitrile (C terminale) e i composti senza idrogeni per natura (cianuro, Cr(CO)₆: avviso B); corretto (voce 20).
  Nitrili e acetiluri erano già giusti.
* **ADDSYM.** spglib con tolleranze 0,1 e 0,3 Å. Sei strutture COD riscritte in P1 (P2₁2₁2₁, Pbca, P2₁/c, P-1, P2₁/n,
  C2/c): PISTAK e PLATON ritrovano lo stesso gruppo in tutti e sei i casi. Sui 200 CIF originali PISTAK dava due
  avvisi A in strutture disordinate dove PLATON non trova nulla; corretto (voce 21).
* **Valori mancanti** (`?` per densità e Z): nessun errore, la formula non si controlla (avviso G); test.

## Sezione 7: export per chimica quantistica e celle periodiche (5 novembre 2026)

Script: `validation/check_export.py`. I programmi non sono installati su questa macchina: si è verificato che i file
descrivano la struttura giusta (parser propri, confronto con la cella espansa da gemmi) e che le parole chiave
corrispondano ai manuali; resta da farli girare (voce 27).

### Geometria nei file

| Verifica | Risultato |
| --- | --- |
| POSCAR, pw.x, CP2K, CRYSTAL contro la cella del CIF espansa da gemmi (7 strutture COD, 68–592 atomi, H come nel CIF) | stessa composizione, ogni sito del CIF a meno di 5·10⁻¹⁵ Å da un atomo del file, volume uguale (10⁻¹⁰ relativo) |
| Con gli H normalizzati (POSCAR) | atomi pesanti fermi; C–H 1,085 (aromatico) e 1,095 Å (sp³), O–H 0,989 Å, come `XH_NEUTRON` |
| Dimero ORCA (8 strutture) | ogni atomo è un sito del cristallo (operazione + traslazione intera) entro 10⁻⁸ Å |

### Parole chiave contro i manuali

| Programma | Verificato |
| --- | --- |
| ORCA 6.1 | atomi fantasma `C :` (manuale, sezione Counterpoise Corrections: spazio e due punti, come scrive PISTAK); `%maxcore` in MB per processo |
| Psi4 | `energy('…-d3bj', bsse_type='cp')`, SAPT con `--` tra i monomeri (già nei test) |
| Gaussian 16 | `Counterpoise=2`, riga carica/molteplicità totale e per frammento, `(Fragment=n)`; D4 non disponibile: errore `qcUnsupported` esplicito |
| Quantum ESPRESSO | `vdw_corr = 'grimme-d3'`, `dftd3_version` 3 = D3 zero, 4 = D3(BJ) (INPUT_PW); `ecutrho` predefinito 4 × ecutwfc, 8–12 × per ultrasoft |
| VASP | `IVDW` 11 = D3 zero, 12 = D3(BJ), 13 = D4 |
| CRYSTAL23 | `BASISSET` dopo la fine della geometria, `DFTD3`/`VERSION`/`4`/`END` dopo `DFT … END`: stesso ordine dell'esempio del manuale (urea, cap. 5.1); `END` e `ENDG` equivalenti (contano i primi tre caratteri); `OPTGEOM` ultima parola chiave della geometria |

La normalizzazione degli H è la stessa per dimeri e celle periodiche (`normalize_hydrogens`), e il file lo dice.
Il dimero esporta sempre le molecole intere: per un macrociclo (COD 4032966, 260 atomi per molecola) il dimero ha 520
atomi, come deve.

## Sezione 8: i test (5 novembre 2026)

### Che cosa prova ogni file

| Tipo di oracolo | File | Basta per il paper? |
| --- | --- | --- |
| Valori impostati a mano su cristalli sintetici (`cifgen`, verificato nella voce 9) | test_core (parte nuova), test_symmetry, test_fragments, test_stacks, test_uncertainty, test_geometry (casi sintetici) | Sì: oracolo indipendente; tolleranze ora 10⁻⁴ Å (voce 28) |
| Riferimento esterno | test_geometry (tabelle `_geom_` dei CIF COD), test_revisione (topologia dei policicli, PLATON, gemmi, conti a mano), test_checks (densità e Z dichiarati nei CIF) | Sì |
| Coerenza interna, invarianti | test_stages (analisi a stadi = analisi da zero), test_collection, test_no_pairs, test_qc e test_periodic (file contro la struttura di PISTAK) | Parziale: un errore sistematico passa; per gli export c'è `check_export.py` (sezione 7) |
| Parità con il sito | test_parity, test_core (porting dei test del sito) | No da sola: stesso algoritmo |
| Istantanee | test_snapshots | Ora sì per le COD: valori rifatti a parte (voce 29) |
| Interfaccia, lingue, manuale | test_gui, test_i18n, test_layout, test_manual, test_measure3d, test_session, test_structure | Non riguardano i numeri |

### Tolleranze

Misurate dentro ogni test con un plugin di pytest che registra lo scarto reale delle asserzioni con tolleranza ≥ 10⁻³:
sui cristalli sintetici gli scarti sono 10⁻⁶–10⁻⁵ Å (arrotondamento delle coordinate a 6 decimali nei CIF di prova),
quindi le tolleranze erano 200 volte più larghe del necessario: un errore di 0,5 mÅ sarebbe passato. Strette a 10⁻⁴ Å
e 10⁻³°. Restano larghe, motivatamente: densità contro il CIF (scarto 3·10⁻⁴, il CIF ha 3 decimali), tabelle `_geom_`
(cifre del CIF), σ del Monte Carlo (rel 0,05–0,12: errore statistico con poche centinaia di campioni).

### CIF non distribuibili

`DATI_REALI/` (CIF della CSD) serve a test_geometry, test_checks e test_manual; senza, quei casi si saltano con un messaggio
esplicito («CIF reali non distribuiti»). Le stesse verifiche girano sui CIF COD in `tests/data/cod/` e, in grande,
negli script di `validation/`: la copertura dei numeri non dipende da `DATI_REALI/`.

## Sezione 9: rischi del codice generato (5 novembre 2026)

### Costanti e fonti

- **Raggi covalenti**: i 96 valori di `COVALENT_RADII` uguali, uno per uno, alla tabella di Cordero et al. 2008 (C sp³,
  Mn/Fe/Co a basso spin come dichiarato). Le differenze con gemmi della voce 16 sono di gemmi.
- **Raggi di van der Waals**: gruppi principali uguali alla tab. 12 di Mantina et al. 2009 (H 1,10 di Rowland e Taylor
  1996), Br 1,85 di Bondi (Mantina 1,83), metalli di transizione di Bondi 1964 come li riporta WebElements. Il commento
  citava solo Bondi: corretto, con il manuale.
- **X–H da neutroni**: già controllate sul testo di Allen e Bruno 2010 (voce 7). Riferimenti citati (Hall 1981, Acta
  Cryst. A37, 517; International Tables B §1.4; Cordero, Dalton Trans. 2008, 2832; Mantina, J. Phys. Chem. A 2009, 113,
  5806; Bondi, J. Phys. Chem. 1964, 68, 441; Allen e Bruno, Acta Cryst. B66, 380): volumi e pagine corretti.
- **Soglie che sono scelte scientifiche** (da dichiarare nel paper, non hanno una fonte unica): `max_dihedral` 20°,
  `max_interplanar` 4,2 Å, `ring_rms` 0,10 e `molecule_rms` 0,15 Å, `bond_tolerance` 0,4 Å, `same_site` 0,3 Å,
  `max_ring_size` 8, `ring_pair_max_distance` 6,0 Å (come PLATON), `min_fragment_atoms` 5, soglie fisse dei contatti
  (C···C 3,0, H···H 1,8 Å) e margine 0,5 Å, sovrapposizione 0,7 Å, tolleranze ADDSYM 0,1/0,3 Å, margine dei frammenti
  1,0 Å. Motivate nel codice: `contact_cutoff`, `MAX_TORSION_ANGLE`, `TRIPLE_CX`, `SP_ANGLE`, `ADDSYM_MAX_EXCLUDED`.
- **Limiti tecnici**: `MAX_SHELLS` 8 (segnalato, avviso `shellsTruncated`), `max_cell_atoms` 6000, `MAX_BYTES` 50 MB
  (errori espliciti), `neighbour_shells` 2.

### Blocchi `except Exception` (23)

| Esito | Dove | Giudizio |
| --- | --- | --- |
| Messaggio all'utente | cli.py (4), gui/session.py (3), gui/main_window.py (2), batch.py (3, la riga della struttura riporta l'errore), collection_window.py, view3d.py (vista 3D sostituita da un riquadro con il messaggio) | Corretti |
| Avviso o registro degli errori | validation.py:83 (G `checkFailed`), params/registry.py e params/structure.py (`param_errors`, scheda Avvisi), params/plugins.py (`load_errors`) | Corretti |
| Ripiego dichiarato | i18n.py (configurazione illeggibile → italiano), spacegroups.py (gemmi non accetta il simbolo → «gruppo sconosciuto» più avanti), main_window.py:314 (l'errore lo mostra poi l'analisi) | Corretti |
| Silenziosi, ora segnalati | validation.py:361 spglib (voce 34), params/uncertainty.py:29 NaN nei campioni (voce 35) | Corretti |

Altri ripieghi silenziosi: `DEFAULT_COVALENT` non si usa mai (il lettore scarta gli elementi sconosciuti);
`DEFAULT_VDW` sì, ora con l'avviso `vdwDefault` (voce 36).

### Licenze e proprietà intellettuale

Il parser dei simboli di Hall e quello dei CIF sono propri (nessuna traccia di sginfo, cctbx, PyCifRW); l'unico codice
derivato è `gui/flow_layout.py` (voce 37). Licenze delle dipendenze dai metadati installati (voce 38).

### Riproducibilità

`check_riproducibilita.py`: 40 strutture COD in 1 e in 4 processi, e due volte in 4, danno risultati identici; il Monte
Carlo con lo stesso seme dà le stesse σ. Versioni esatte in `ambiente.txt` (voce 39).

### Dichiarazione sull'uso dell'IA

Bozza in `dichiarazione_ia.md` (voce 41).

