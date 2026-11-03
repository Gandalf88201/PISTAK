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
| 16 | core/elements.py | Raggi diversi dalle tabelle di gemmi: covalente C 0,76 (sp³, dichiarato) contro 0,73; Lu 1,87 e Hf 1,75 scambiati rispetto a gemmi; vdW H 1,10 (Rowland e Taylor) mentre il commento cita Bondi (1,20) | P1 (H: contatti vdW) | **Aperto** | Controllare sul testo di Cordero 2008 e Bondi 1964; correggere il commento su H |
| 17 | core/structure.py | Formule ricostruite diverse da `_chemical_formula_moiety` non spiegate: 7208500 (due acetonitrili in un'unica "molecola"), 2201728 (un CH₃ in più) | P2 | **Aperto** | Da esaminare a mano |
| 10 | core/cif.py | **Nuovo.** `tokenize` va in ciclo infinito (memoria che cresce senza limite) su un carattere di spaziatura diverso da spazio e tab: NBSP, form feed, tab verticale, U+2003… | P2 | **Corretto** | `if ch.isspace()` al posto di `if ch in " \t"`; test con 5 caratteri in un processo separato |

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
