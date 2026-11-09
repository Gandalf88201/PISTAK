# Registro delle verifiche

*English summary and glossary of the Italian file names: `validation/REVIEW_LOG_EN.md`.*

Registro della revisione del codice (checklist del 3 novembre 2026, sezione 10). Una riga per problema trovato o
verifica chiusa; sotto la tabella, per ogni voce, l'oracolo usato e il risultato. Le verifiche si rieseguono con

```bash
.venv/bin/python validation/check_registro.py            # voci del registro: cell rings beta height unit shells cutoff
.venv/bin/python -m pip install -r validation/requirements.txt   # gemmi, networkx, RDKit, …: gli oracoli
.venv/bin/python validation/check_simmetria.py           # sezione 1: simboli di Hall e H-M contro gemmi
.venv/bin/python validation/scarica_cod.py               # sezione 1: i 200 CIF COD di cod_sample.txt in cod_cache/
.venv/bin/python validation/check_cif.py                 # sezione 1: lettura, s.u., elementi, disordine
.venv/bin/python validation/check_struttura.py           # sezioni 1–3: cella, formule, cella obliqua, anelli, raggi
.venv/bin/python validation/check_platon.py              # sezioni 4–5: coppie di anelli contro PLATON CALC ALL
.venv/bin/python validation/check_platon.py platon_polimeri   # voce 47: le stesse coppie nei polimeri
.venv/bin/python validation/check_incertezze.py          # sezione 6: Monte Carlo, geometria contro le tabelle _geom_
.venv/bin/python validation/check_export.py              # sezione 7: file per ORCA/Psi4/Gaussian e codici periodici
.venv/bin/python validation/check_snapshot.py            # sezione 8: istantanee COD rifatte con numpy
.venv/bin/python validation/check_struttura.py oblique origin   # sezione 8: invarianza per cella e origine
.venv/bin/python validation/check_struttura.py rdkit    # sezione 10: legami e anelli contro RDKit
.venv/bin/python validation/check_struttura.py hthresh  # voce 23: soglie del controllo sugli idrogeni
.venv/bin/python validation/check_riproducibilita.py    # sezione 9: batch 1 contro 4 processi, seme del Monte Carlo
.venv/bin/python validation/check_sensibilita.py        # sezione 10: sensibilità alle soglie (200 CIF, ~2 min)
```

Gli oracoli non usano `pistak.core`: geometria ricalcolata con numpy, topologia nota dei policicli, tensore metrico,
conti a mano.

| # | Modulo | Problema o verifica | Priorità | Stato | Azione |
| --- | --- | --- | --- | --- | --- |
| 1 | params/extra.py | `beta_b` etichettato come γ | P1 → P3 | **Corretto**: i numeri erano giusti, la chiave no | Chiavi `beta` e `gamma`; `beta_a`, `beta_b` valgono ancora in `Pair.get` e nel registro |
| 2 | params/extra.py | `mean_atom_height`: descrizione (piano di A) diversa dal codice (normale media) | P1 | **Corretto**: il vecchio codice dava la distanza interplanare sotto altro nome | Ora media di \|(x − Cg_A)·n_A\| sugli atomi di B; test in `tests/test_review.py` |
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
| 14 | core/cif.py, core/structure.py | Gruppi di disordine negativi: due casi diversi. (a) 1 e −1 usati come due alternative (COD 4113834: 0,61/0,39, C32A···C32B 0,38 Å): la regola `major` toglieva il segno e li teneva entrambi. (b) Gruppo −n davvero su un elemento di simmetria (4064128, 7042307): le sue immagini sono l'altra orientazione e venivano tenute tutte | P1 | **Corretto** | (a) −n separato da n se i loro atomi distano < 1 Å (`_signed_alternatives`); (b) per ogni gruppo negativo si tiene una copia per sito, scartando intera ogni copia che si sovrappone a una già tenuta (`special_disorder_copies`). Su 200 CIF cambiano solo i 3 interessati: spariscono gli avvisi `overlap` di 4113834 e 7042307, e in 4064128 il CH₂Cl₂ torna di 5 atomi (prima 10, con un falso anello e 4 coppie inclinate fittizie). Test |
| 15 | core/crystal.py | **Nuovo.** `cell_contents` contava una volta sola una posizione occupata da due atomi diversi del CIF: un CIF che scrive anche le immagini per simmetria con occupanza ½ (COD 2100449) dava metà del contenuto della cella, e densità e formula sbagliate negli avvisi | P2 (avvisi, non coppie) | **Corretto** | Occupanze di atomi diversi sommate fino a 1; immagini dello stesso atomo (posizioni speciali) una volta; test |
| 16 | core/elements.py | Raggi diversi dalle tabelle di gemmi: covalente C 0,76 (sp³, dichiarato) contro 0,73; Lu 1,87 e Hf 1,75 scambiati rispetto a gemmi; vdW H 1,10 (Rowland e Taylor) mentre il commento cita Bondi (1,20) | P1 (H: contatti vdW) | **Corretto** (commento) | Valori giusti: i 96 covalenti uguali alla tabella di Cordero 2008 (le differenze erano di gemmi); i vdW uguali alla tab. 12 di Mantina 2009 per i gruppi principali (H 1,10 di Rowland e Taylor 1996), Br 1,85 e metalli di transizione di Bondi 1964. Corretti il commento in `elements.py` e il manuale, che citavano solo Bondi |
| 17 | core/structure.py | Formule ricostruite diverse da `_chemical_formula_moiety` non spiegate: 7208500 (due acetonitrili in un'unica "molecola"), 2201728 (un CH₃ in più) | P2 | **Corretto** | 2201728: errore dell'oracolo, non di PISTAK (il lettore della formula ignorava il «2» di `(C H3)2`; corretto, ora coincide). 7208500: l'acetonitrile disordinato aveva le due orientazioni in assemblaggi diversi (A gruppo 1, 0,65; B gruppo 2, 0,35); `major` sceglieva per assemblaggio e le teneva entrambe, sovrapposte a 0,31 Å e legate in un falso C₄H₆N₂. Ora gli assemblaggi i cui gruppi si sovrappongono (< 1 Å) si uniscono prima della scelta (`_merged_assemblies`, decisione dell'autore). Stesso difetto, ora corretto, in 7025278 (BF₄⁻ su 3 orientazioni), 7034397, 7232723 (due molecole intere raddoppiate, 12 sovrapposizioni) e 7051759 (Yb su due siti a 0,79 Å). Formule contro `_chemical_formula_moiety`: 79 uguali su 114. Test |
| 10 | core/cif.py | **Nuovo.** `tokenize` va in ciclo infinito (memoria che cresce senza limite) su un carattere di spaziatura diverso da spazio e tab: NBSP, form feed, tab verticale, U+2003… | P2 | **Corretto** | `if ch.isspace()` al posto di `if ch in " \t"`; test con 5 caratteri in un processo separato |
| 18 | core/intramolecular.py | **Nuovo.** Torsioni riportate anche attraverso tre atomi allineati (X–M–X su un centro di inversione, angolo 180°): valore indefinito, che dipende dagli arrotondamenti | P2 | **Corretto** | Torsioni con un angolo di legame oltre 175° non riportate (`MAX_TORSION_ANGLE`); test, manuale |
| 19 | core/intramolecular.py | **Nuovo.** Misure con le stesse etichette considerate uguali per simmetria: L1–M–L2 e L1–M–L2′ (θ e 180° − θ attorno a un centro di inversione), torsioni cis e trans dello stesso frammento, ne restava una | P1 (geometria riportata) | **Corretto** | Uguali solo se legate da un'operazione del cristallo (posizioni relative, simmetria del sito del primo atomo); 2100417: una torsione in più; misure CIF senza corrispondente da 539 a 27 |
| 20 | core/validation.py | `_hydrogens`: falsi positivi su isonitrile (C terminale) e su composti senza H per natura (cianuri, carbonili metallici: avviso B) | P2 | **Corretto** | C terminale con triplo legame a N/O (≤ 1,18 Å) escluso; `noHydrogens` solo se qualche C sembra volerne; test |
| 21 | core/validation.py | ADDSYM su strutture disordinate: la componente tenuta da `major` aveva più simmetria del modello, avviso A falso (COD 1530955, 1531236; PLATON non trova nulla) | P2 | **Corretto** | Come PLATON, atomi con occupanza parziale esclusi; oltre metà esclusi, avviso G `addsymDisorder` |
| 22 | manuale, README | «pre-checkCIF» e livelli A/B/C/G fanno pensare agli alert ufficiali IUCr | P2 | **Corretto** | Detto nel manuale (it, en) e nei README che soglie e controlli sono di PISTAK |
| 23 | core/validation.py | `TRIPLE_CX` = 1,18 Å e `SP_ANGLE` = 150°: valori motivati (C≡N, C≡O 1,13–1,17 Å; C=O da 1,19 Å; sp 180°, sp² 120°) ma senza una tabella citata | P3 | **Verificato e citato** (6 novembre 2026) | Allen et al., J. Chem. Soc. Perkin Trans. 2 1987, S1–S19 (PDF fornito dall'autore): Csp≡N 1,136 Å nei nitrili C*–C≡N (σ 0,010, n = 140), 1,138 negli aril-nitrili, fino a 1,155 in (S–C≡N)⁻; C=O 1,192 nelle aldeidi, 1,210 nei chetoni (n = 474), 1,196 negli esteri, 1,187 in anidridi ed esteri arilici (le medie più corte; quartile inferiore minimo 1,181). 1,18 Å sta sopra ogni media C≡N e sotto ogni media C=O. Sui 200 CIF COD: 697 C=O di carbonili da 1,182 Å in su (mediana 1,221), nessuno sotto 1,18; 8980 C sp² con un H ad angoli ≤ 134,5°, nessuno sopra 150; i C con due vicini e nessun H stanno a 176° di mediana. Allen dà solo lunghezze: `SP_ANGLE` è verificato sui dati COD (`check_struttura.py hthresh`). Citazione nei commenti del codice |
| 24 | io/qc.py | Molteplicità del dimero fissata all'accoppiamento ad alto spin: per due radicali l'energia dipende dalla scelta singoletto/tripletto | P1 (energie di radicali) | **Corretto** | Opzione `dimer_multiplicity` (API, `--dimer-mult`, finestra), controllata contro gli spin dei monomeri; il predefinito resta l'alto spin, scritto nel file |
| 25 | io/qc.py | `%maxcore` ORCA = memoria totale / processi, senza margine (ORCA a volte lo supera) | P2 | **Corretto** | 75 % della memoria per processo, come consigliano i centri di calcolo; MB per processo verificato |
| 26 | io/qc.py, io/periodic.py | Nei file non era scritto che i monomeri sono molecole intere, quanti H restano come nel CIF, che `ecutrho` = 8 × presuppone pseudopotenziali ultrasoft/PAW, che punti k, cutoff e Γ di CP2K sono valori di comodo | P2 | **Corretto** | Righe di commento in ogni file (CRYSTAL: nel `.notes.txt`); manuale |
| 27 | io/qc.py, io/periodic.py | Esecuzione reale con ORCA, Psi4, Gaussian, VASP, QE, CP2K, CRYSTAL | P1 | **Rinviato** (decisione dell'autore, 5 novembre 2026: le prove le farà l'autore) | I programmi non sono installati qui; servono una prova per programma (dimero di benzene, cella del naftalene) e il confronto CP ORCA contro Psi4 |
| 28 | tests/ | Tolleranze larghe (`abs=1e-3` Å, 0,5° e 0,01°) su cristalli sintetici costruiti a mano | P1 | **Corretto** | Scarti reali misurati: ≤ 5·10⁻⁶ Å e 5·10⁻⁵°; tolleranze portate a 10⁻⁴ Å e 10⁻³° (16 asserzioni in 5 file); quelle contro i CIF (densità, tabelle `_geom_`) restano legate alle cifre del CIF |
| 29 | tests/snapshots | Istantanee COD scritte dal programma stesso | P1 | **Verificato** | Le 95 coppie intermolecolari delle 5 istantanee rifatte con numpy: scarti ≤ 5·10⁻⁹ Å e 2·10⁻⁶° (precisione delle istantanee) |
| 30 | validation | Invarianza cristallografica | P1 | **Verificato** | Cella obliqua equivalente e origine spostata di (0,137; 0,291; 0,443): 40/40 strutture con le stesse coppie entro 10⁻⁶ Å (prima la verifica arrotondava a 3 decimali e segnalava differenze di arrotondamento) |
| 31 | configurazione utente | **Nuovo.** `analyse_file` senza soglie usa il file dell'utente (`user.toml`): su questa macchina `contact_cutoff` = 4,0 (il vecchio predefinito, salvato prima della voce 6), `max_dihedral` = 5, `overlap_measure` = area, `contact_check` = vdw. I test sono isolati (`PISTAK_CONFIG_DIR` temporanea), le analisi dell'utente no | P1 (risultati dell'utente) | **Corretto** (decisione dell'autore, 6 novembre 2026) | `user.toml` aggiornato (copia in `user.toml.bak-2026-11-06`): tolto `contact_cutoff` = 4,0; restano le tre scelte personali. La finestra «Impostazioni» ora salva solo le soglie diverse dai predefiniti, così un predefinito cambiato in seguito vale anche per l'utente. Nuovo controllo rapido `oldDefault` (C): una soglia uguale a un vecchio predefinito (`config.SUPERSEDED_DEFAULTS`; non nel modo del sito). Test `test_old_default_reported_and_user_toml_only_differences` |
| 32 | tests/ | Copertura (`pytest --cov`) e test basati su proprietà (hypothesis) | P2 | **Copertura misurata** (revisione 2, 8 novembre 2026); hypothesis aperto | pytest-cov: 90 % delle righe in tutto, 97 % senza l'interfaccia grafica (core, io, params), 83 % l'interfaccia. Test basati su proprietà non scritti |
| 33 | tests/ | Nomi dei test in italiano | P3 | **Corretto** (decisione dell'autore: tradurre) | Nomi, docstring, commenti e messaggi di tests/ in inglese (251 test, 24 file; `test_revisione.py` → `test_review.py`). Restano in italiano solo i testi che i test confrontano con l'interfaccia italiana. Controllo: l'albero sintattico di ogni file, tolti nomi e testi, è identico a quello di prima; stessi 426 test superati |
| 34 | core/validation.py | **Nuovo.** Se spglib c'è ma non riesce, `find_symmetry` restituiva None come quando spglib manca: avviso G «serve il pacchetto spglib» (falso) a 0,1 Å, nessun avviso a 0,3 Å | P2 | **Corretto** | L'errore di spglib diventa l'avviso G `checkFailed` con il messaggio; test. Su 200 CIF COD nessun nuovo avviso |
| 35 | core/uncertainty.py | **Nuovo.** Un parametro NaN in parte dei campioni del Monte Carlo (un plugin che fallisce, una geometria degenere) dava una σ calcolata sui soli campioni rimasti, senza dirlo | P2 | **Corretto** | `Uncertainty.valid` e `incomplete`: avviso nel terminale, riga nel Dettaglio, campo `su_samples` nel JSON; test |
| 36 | core/elements.py | **Nuovo.** 48 elementi su 96 senza raggio di van der Waals (Fe, Co, Mn, Ru, Ti, lantanidi, …): contatti e parametri di vdW usavano 1,9 Å in silenzio | P2 | **Corretto** (decisione dell'autore: Alvarez 2013) | Blocchi d e f da S. Alvarez, Dalton Trans. 2013, 42, 8617, tabella 1 (letta dal PDF in accesso aperto e confrontata con la trascrizione di ASE, `ase.data.vdw_alvarez`: 93 elementi uguali); sostituiscono anche i valori provvisori di Bondi per Ni, Cu, Zn, Pd, Ag, Cd, Pt, Au, Hg. Gruppi principali invariati (Mantina 2009). Senza raggio restano Pm e Fm–Lr (avviso `vdwDefault`). Sulle istantanee cambia solo 7024607 (Mn): sparisce l'avviso G |
| 37 | gui/flow_layout.py | **Nuovo.** «Adattato dall'esempio Flow Layout di Qt» senza nota di copyright né licenza (gli esempi di Qt for Python sono BSD-3-Clause, che chiede di conservarle) | P1 (licenza) | **Corretto e verificato** | Nota di copyright e testo BSD-3 nell'intestazione; confrontata con l'originale (pyside-setup, `examples/widgets/layouts/flowlayout/flowlayout.py`, rami 6.5, 6.8, dev): stessi titolari (Riverbank Computing 2013, The Qt Company 2022); l'originale è «LicenseRef-Qt-Commercial OR BSD-3-Clause» e PISTAK usa l'opzione BSD-3-Clause, ora scritto nell'intestazione |
| 38 | README | Licenze delle dipendenze: mancavano gemmi (MPL-2.0, facoltativo) e vtk (BSD, tramite pyvista) | P3 | **Corretto** | Aggiunte; le altre (numpy, scipy, spglib BSD; pyqtgraph, pyvista MIT; PySide6 LGPL) confermate dai metadati installati |
| 39 | requirements | Solo versioni minime (`>=`), nessun file di versioni bloccate: i numeri del paper non sono legati a un ambiente | P2 | **Corretto** (documentato) | `validation/ambiente.txt`: versioni esatte con cui girano test e controlli. Batch e Monte Carlo deterministici (`check_riproducibilita.py`) |
| 40 | pistak/ | Leggibilità: ruff pulito (nessun import o variabile inutilizzati; 11 argomenti non usati, tutti segnature di slot Qt). vulture non installato; il mypy di sistema (1.19) non legge le stub di numpy 2.5 | P3 | **Aperto** | Installare mypy e vulture nel venv (insieme alla voce 32) |
| 41 | paper | Dichiarazione sull'uso dell'IA | P1 (paper) | **Approvata come bozza** (6 novembre 2026) | `validation/dichiarazione_ia.md`: fatti dal repository, testo IT/EN; l'autore ha approvato l'impostazione; le parti segnate da completare restano sue |
| 42 | repository | Sezione 10, passo 1: nessuna CI (`.github/workflows` assente); pytest gira solo in locale | P2 | **Corretto** (run 42055379324: 4 job verdi) | `.github/workflows/test.yml`: job `calcolo` (Python 3.11–3.13, solo numpy/scipy/spglib: ruff e `pytest -m "not gui"`) e job `gui` (Qt offscreen). Simulato in locale bloccando i pacchetti grafici: trovato `test_manual` che importava la finestra del manuale senza PySide6, ora si salta |
| 43 | core/cores.py | Sezione 10, passo 4: anelli contro RDKit | P2 | **Verificato** | `check_struttura.py rdkit` (RDKit 2025.09): su 204 molecole organiche distinte del campione COD, legami rifatti da RDKit dalle sole coordinate uguali in 193; anelli (SSSR simmetrizzata, ≤ 8) uguali in 200. Le 4 diverse: 3 polibromuri (triangoli Br₃ in più, non π) e 4112371, un anello a 6 in più in un sistema 4+6 fuso (rms 0,31 Å, mai un nucleo π). Nessun anello planare diverso |
| 44 | core/pairs.py | Sezione 10, passo 5: sensibilità dei risultati alle soglie | P1 | **Verificato** | `check_sensibilita.py` su 200 CIF COD: `max_dihedral` è la scelta che pesa (±5° → ±6 % di stack, ±10° → −18/+15 %); `max_interplanar` ±0,2 Å → ≤ 2 %; le altre soglie ≤ 2 %, tranne `min_fragment_atoms` = 6 (−5 %: frammenti di 5 atomi). Da riportare nel paper |
| 45 | repository | Sezione 10, passo 8: versione per il paper (tag attuali v0.0.1, v0.0.2); istantanee da rigenerare e congelare | P1 | **Chiuso** (6 novembre 2026) | `pistak.__version__` = 1.0.0, tag e release v1.0.0 su GitHub dopo le voci 23, 31, 33, 36, 37 e 41; istantanee rigenerate e confrontate a ogni cambio (14 + 7024607), tutte uguali al codice della release. Restano aperte solo 32 e 40 (pytest-cov, hypothesis, mypy, vulture nel venv) |
| 46 | core/ringpairs.py, core/intramolecular.py, core/validation.py, core/pairs.py, params/extra.py | **Nuovo** (dalla prima CI su Linux). L'ordine di elementi equivalenti per simmetria dipendeva dall'ultima cifra del calcolo, diversa tra macOS e Linux: coppie di anelli alla stessa Cg···Cg, contatti H1···H8/H2···H7, i 4 contatti mostrati nell'avviso `shortContact` (nell'ordine della ricerca dei vicini). Il diedro da `arccos` vicino a 0° (10⁻⁶° contro 0) e `slip_direction` −180 contro +180 (zero con segno) | P2 (riproducibilità) | **Corretto** | Ordinamenti con spareggio sugli indici/etichette; avviso con i 4 contatti più corti; diedro da atan2(|a×b|, |a·b|); `slip_direction` in (−180, 180]. Istantanee rigenerate: confrontate con le vecchie a meno dell'ordine, cambiano solo 10 diedri ~10⁻⁶° → 0 e 2 direzioni −180 → 180 |
| 47 | core/structure.py, core/pairs.py, core/analysis.py, core/diagnosis.py, io/qc.py, io/export.py, labels.py | **Nuovo** (dalla prova dell'autore su COD 7109544, un polimero di Cu₂ senza nessuna coppia). Le reti di legami estese (polimeri di coordinazione, MOF) erano escluse: nessuna coppia π. Nella ricostruzione dell'unità ripetuta, inoltre, la visita in ampiezza poteva lasciare un anello spezzato tra due celle (7024607: una piridina su due senza frammento π) | P1 | **Corretto** | L'unità ripetuta si ricompone unendo prima i legami degli anelli di non metalli, poi quelli dei leganti, poi quelli con i metalli: anelli e leganti restano interi. Coppie tra frammenti dell'unità e delle unità vicine entro `contact_cutoff` + 1 Å, segnate `polymer` (contatto tra i frammenti, escluse dall'export QC); simmetria per immagini di frammenti. Contro PLATON su 20 polimeri: 217 coppie di anelli, scarti entro la stampa di PLATON (sotto). Su 200 CIF COD nessuna struttura non polimerica cambia |
| 48 | params/extra.py | **Nuovo**. Nelle coppie intramolecolari (frammenti della stessa molecola) `cc_contact` e `n_vdw_contacts` si misuravano tra la molecola e sé stessa: 0 Å e centinaia di contatti | P2 | **Corretto** | Si misurano tra gli atomi dei due frammenti (anche per le coppie di polimero). Istantanee: in 2100417 e 4500022 cambiano solo queste due colonne delle coppie intramolecolari |
| 49 | core/pairs.py, io/export.py, manuale | **Nuovo** (trovato scrivendo il paper, 6 novembre 2026). Il criterio di stack accettava due frammenti legati tra loro: in COD 2210645 (NiS₂C₂) e 4339142 (CuO₂C₃) i due anelli chelati complanari attorno allo stesso metallo erano «π-stack» a 0,04–0,18 Å, in 1538535 (Ta₂VC₂) anelli Ta/V/C della rete a 0,00 e 2,14 Å: 9 stack su 144 nel campione COD. Il manuale diceva inoltre che gli anelli con un metallo non sono frammenti π nelle reti estese, ma il codice li accettava | P1 | **Corretto** | Due frammenti con un atomo in comune o un legame tra loro (raggi covalenti + `bond_tolerance`, come in `structure`) non sono mai uno stack: la coppia resta parallel/tilted con tutti i parametri, `Pair.bonded`, chiave `bonded` nel JSON. Gli anelli con metalli restano frammenti π (manuale corretto). Sui 200 CIF COD: stack 144 → 135, cambiano solo le 9 coppie (nessuna tra molecole distinte); sensibilità rifatta (135 stack). Istantanee esistenti invariate; nuovo CIF di prova 2210645 con le sue istantanee; 3 test |
| 50 | core/config.py, core/pairs.py, params/extra.py, default.toml | **Nuovo** (revisione del paper, 7 novembre 2026). Il criterio di stack predefinito (almeno un atomo di un frammento dentro l'involucro convesso dell'altro) accettava sovrapposizioni minime: 25 dei 100 stack intermolecolari sotto il 10 % dell'area, slittamenti fino a 9,5 Å | P1 | **Corretto** | Nel modo canonico la sovrapposizione si misura sull'unione dei poligoni dei singoli anelli (non sull'involucro: le baie e le tacche alle fusioni restano fuori) e il predefinito diventa `overlap_measure = "area"`, `min_overlap = 0,10`. Il modo del sito (`Limits.site()`) resta com'era: involucro, criterio sugli atomi, qualsiasi sovrapposizione; le sue istantanee non cambiano (solo colonne nuove). Test: `tests/test_revision.py::test_union_of_rings_versus_hull` |
| 51 | core/cores.py | **Nuovo** (revisione). Negli organometallici η⁵/η⁶ l'anello aromatico non si trovava: i cicli più corti attraverso un legame C–C del Cp o dell'arene sono i triangoli M–C–C | P1 | **Corretto** | Anelli in due passate: prima tra i soli non metalli, poi quelli con un metallo. Ferrocene (COD 7033930: due Cp), Cr(η⁶-C₆H₆)(CO)₃ (4103814: C₆), Ag(C₆H₆)ClO₄ (2101898: C₆, senza Ag). Nessuna struttura senza metalli cambia. Test: `test_metal_bound_rings_are_found` |
| 52 | core/cores.py, core/pairs.py, core/config.py, core/stages.py | **Nuovo** (revisione). Sistemi coniugati non condensati (biarili complanari, porfirine e ftalocianine senza metallo) erano frammenti separati; la stessa ftalocianina era un frammento con il Cu e quattro senza | P2 | **Opzione** | `fragment_merge = "conjugated"` (predefinito "fused", invariato): unisce sistemi complanari collegati da un legame o da atomi sp/sp² se l'insieme resta entro `molecule_rms`; i macrocicli (anello interno di 16 atomi) entrano nel poligono come triangoli dal centroide. H₂Pc β (2008099): 4 → 1 frammento di 40 atomi, colonna lungo b (4,73 Å, area 0,41). Un biarile ruotato di 40° resta in due frammenti. Test: `test_phthalocyanines_free_base_and_metal`, `test_conjugated_merge_coplanar_and_twisted_biaryl` |
| 53 | io/qc.py, gui/qc_dialog.py, cli.py | **Nuovo** (revisione). Gli input QC non sono mai stati eseguiti (voce 27) e non avvertivano dei limiti del modello | P2 | **Corretto in parte** | Avvertenze in testa all'input: modello disordinato (atomi a occupanza parziale, regola usata), solvente mascherato (SQUEEZE/masks), input controllato solo nella sintassi. Nuovo metodo r2SCAN-3c (solo ORCA; D4 e gCP inclusi, niente counterpoise). L'esecuzione reale resta da fare (voce 27); il paper non presenta più gli export come risultato |
| 54 | core/uncertainty.py, core/validation.py, labels.py | **Nuovo** (revisione). La σ Monte Carlo descrive male le grandezze limitate (θ vicino a 0°); un CIF senza s.u. dava σ = 0 senza avviso | P2 | **Corretto** | `Uncertainty.interval`: percentili 2,5 e 97,5 dei campioni per coppia e parametro. Avviso G `noSU` se nessuna coordinata ha una s.u. |
| 55 | core/pairs.py, params/extra.py, core/config.py | **Nuovo** (revisione). L'asse lungo (PCA) non è definito per un nucleo quasi isotropo (asse Cn, n ≥ 3: benzene, coronene, porfirine) e il riferimento delle coppie è discontinuo per assi quasi ortogonali | P1 | **Corretto** | `max_axis_isotropy = 0,9` (√(λ₂/λ₁), con λ₁ ≥ λ₂ gli autovalori della dispersione nel piano; corretto nella revisione 2: prima era scritto λ₂/λ₁): sopra, asse non definito; `max_axis_twist = 75°`: sopra, riferimento non definito. In entrambi i casi s_L e s_W sono NaN (lo slittamento totale resta). Nuovi parametri: `twist` (angolo tra gli assi lunghi), `planes_cross` (la retta d'intersezione dei piani passa tra i centroidi: d⊥ = |h_A + h_B|/(2 cos θ/2) è allora minore delle due distanze perpendicolari), `local_interplanar` (distanza media dal piano dell'altro frammento degli atomi che vi cadono sopra in proiezione). Il modo del sito non cambia |
| 56 | core/cores.py | **Nuovo** (revisione). Il ripiego «molecola intera» senza anelli accettava scheletri lineari (normale non definita) | P3 | **Corretto** | Rifiutato se √(λ₂/λ₁) della dispersione nel piano è sotto 0,05 (`MIN_PLANE_ASPECT`). Test: `test_linear_molecule_has_no_fallback_plane` |
| 57 | core/elements.py | **Nuovo** (revisione). Raggio di van der Waals di Br da Bondi (1,85 Å) mentre gli altri non metalli vengono da Mantina 2009 | P3 | **Corretto** | Br = 1,83 Å (Mantina): un solo insieme per i non metalli. Nessuna istantanea cambia |
| 58 | core/pairs.py, io/export.py | **Nuovo** (revisione 2, V1). I due Cp del ferrocene (legati entrambi al Fe, non tra loro) erano uno π-stack intramolecolare a 3,31 Å | P1 | **Corretto** | Regola `bridged`: due frammenti legati entrambi a uno stesso atomo (non H) non sono mai uno stack, come i `bonded`. Il ciclofano (nessun atomo comune) resta uno stack. Test: `test_metallocene_rings_are_bridged_not_stacked` |
| 59 | core/analysis.py, core/config.py, core/stages.py | **Nuovo** (revisione 2, C11). Le strutture senza C (reti inorganiche) non hanno stack ma sono le sole che rompono l'invarianza alla scelta della cella | P2 | **Corretto** | `pi_scope = "organic"` (predefinito): nessun frammento π se nessuna molecola ha un legame C–C o C–H (definizione che funziona anche per carburi, carbonati, cianuri). `"all"` ripristina la ricerca; `Limits.site()` usa "all". Test: `test_structures_without_organic_bonds_have_no_fragments` |
| 60 | core/structure.py | **Nuovo** (revisione 2, C7). Riempiendo la cella si fondevano tutte le posizioni dello stesso elemento entro 0,3 Å, anche atomi diversi (posizioni alternative di un atomo sdoppiato) | P2 | **Corretto** | Si fondono le immagini dello stesso sito e gli atomi ordinati ripetuti nel CIF; mai un atomo a occupanza parziale o in un gruppo di disordine con un altro atomo. Test: `test_split_positions_are_not_merged` |
| 61 | core/uncertainty.py | **Nuovo** (revisione 2, C8). Nel Monte Carlo gli spostamenti degli atomi su posizioni speciali rompevano i vincoli lineari (x = y, z = 0) | P2 | **Corretto** | Lo spostamento si proietta sulla simmetria del sito: media di R·δ sulle operazioni che mandano il sito in sé. Test: `test_monte_carlo_keeps_special_position_constraints` |
| 62 | core/validation.py, labels.py | **Nuovo** (revisione 2, C9). `noSU` scattava solo se nessuna coordinata aveva una s.u.: gruppi rigidi o frammenti vincolati senza s.u. davano σ = 0 senza avviso | P3 | **Corretto** | Avviso C `fragmentNoSU` con gli atomi non H dei frammenti π senza s.u.; esclusi gli atomi con tutte le coordinate fissate dalla simmetria |
| 63 | core/pairs.py | **Nuovo** (revisione 2, C2: trovato con l'oracolo indipendente `paper/scripts/geometry_oracle.py`, shapely). (a) La sovrapposizione somma le intersezioni anello per anello: se due anelli dello stesso frammento si sovrappongono in proiezione (gabbie, pile di anelli in una rete) l'area comune conta due volte e il rapporto supera 1 (trovato su COD 2100449, 66 anelli a 4 e rapporto 1,7; lì la causa erano atomi doppi, voce 65, ma la protezione resta per i casi veri). (b) I triangoli del macrociclo si aggiungevano anche dove coprivano già un anello. (c) `planes_cross` cambiava segno per frammenti complanari (altezze ≈ 0) | P2 | **Corretto** | (a) Se due anelli piccoli dello stesso frammento si sovrappongono per più del 2 % del più piccolo, il contorno è l'involucro convesso del frammento. (b) Del macrociclo entrano solo i triangoli che non coprono un anello piccolo. (c) `planes_cross` richiede (d·n_A)(d·n_B) < −10⁻⁶. Oracolo sulle coppie con θ ≤ 20° del campione COD: scarto massimo dell'area ≈ 10⁻¹⁴, nessun disaccordo su `planes_cross`. Test: `test_self_overlapping_rings_fall_back_to_the_hull` |
| 64 | core/cores.py | **Nuovo** (revisione 2). Con `fragment_merge = "conjugated"` un sistema già unito attraverso un atomo della frontiera si aggiungeva una seconda volta da un altro atomo: H₂Pc aveva ogni isoindolo due volte, e per la voce 63(a) il contorno diventava l'involucro (area 0,663 invece di 0,351) | P2 | **Corretto** | Un sistema già unito si salta. H₂Pc β (2008099): anelli 4 × 6, 4 × 5 e il 16, area 0,351 come l'unione esatta calcolata con shapely. Test: `test_conjugated_macrocycle_rings_are_not_duplicated` |
| 65 | core/structure.py | **Nuovo** (revisione 2, dal confronto degli anelli con RDKit dopo la voce 60). Con la regola della voce 60 due atomi del CIF a occupanza parziale non si fondevano mai: COD 2100449, che elenca anche le immagini per inversione come atomi distinti con occupanza ½, aveva ogni atomo due volte a 10⁻⁷ Å e 66 anelli a 4 inesistenti | P1 | **Corretto** | Si fondono anche due atomi del CIF che coincidono entro `DUPLICATE` = 0,05 Å, qualunque sia l'occupanza (lo stesso atomo scritto due volte); restano distinte le posizioni alternative di un atomo sdoppiato, più lontane. 2100449: due anelli a 6 per molecola; nel campione COD cambia solo questa struttura (contenuto della cella confrontato con e senza la regola). Test: `test_split_positions_are_not_merged` (esteso) |
| 66 | core/pairs.py | **Nuovo** (CI su Linux, Python 3.12, dopo la voce 46). Le coppie si ordinavano per tipo, distanza interplanare e contatto: due coppie con la stessa distanza (COD 2008099: 2,3655 Å, diversa solo nell'ultima cifra binaria) cambiavano posto da una macchina all'altra, e con loro gli indici | P2 (riproducibilità) | **Corretto** | Distanze arrotondate a 10⁻⁶ Å nella chiave e spareggio su distanza tra i centroidi, indici delle molecole e dei frammenti e traslazione. Cambia solo l'ordine di due coppie di COD 7033930 (stesse 11 coppie, resto dell'istantanea identico); istantanea rigenerata |

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
| 17 | Gruppi di disordine che si sovrappongono sono alternative anche se il CIF li mette in assemblaggi diversi (5 novembre 2026) | Su 200 CIF COD cambiano 8 strutture, tutte con alternative tenute insieme per errore: spariscono 19 sovrapposizioni e le molecole raddoppiate; nessuno stack cambia |
| Monte Carlo | N = 500 campioni per i risultati del paper (5 novembre 2026) | Errore relativo sulla σ ≈ 1/√(2·499) ≈ 3 % (sezione 6) |
| 27 | Esecuzione reale degli input QC e periodici rinviata: la farà l'autore (5 novembre 2026) | Nessuno sui numeri dell'analisi |
| 23 | Controllare `TRIPLE_CX` e `SP_ANGLE` su Allen et al. 1987 e citarlo (6 novembre 2026) | Nessuno: le soglie separano già le classi sui dati COD |
| 31 | Aggiornare `user.toml` e avvisare dei vecchi predefiniti (6 novembre 2026) | Le analisi dell'autore usano `contact_cutoff` = 4,5 Å |
| 33 | Tradurre in inglese nomi e testi dei test (6 novembre 2026) | Nessuno |
| 36 | Raggi di van der Waals di Alvarez 2013 per i metalli dei blocchi d e f (6 novembre 2026) | Contatti e parametri di vdW con i metalli: raggi più grandi (Cu 1,40 → 2,38 Å) |
| 45 | Tag v1.0.0 quando queste voci sono chiuse (6 novembre 2026) | — |

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
| Molecole ricostruite contro `_chemical_formula_moiety` (114 CIF) | 74 uguali, 4 diverse solo negli H; delle 36 diverse (rifatto il 5 novembre dopo le voci 14 e 17: 79 uguali, 3 solo negli H, 32 diverse), quasi tutte per convenzione del CIF (formula per unità asimmetrica o per coppia di ioni, controione tenuto separato da PISTAK) o per reti estese (polimeri di coordinazione); da disordine non risolto 4064128 (gruppi negativi, voce 14), 7036933 e 4330346 (disordine senza gruppi nel CIF: solo l'avviso `partialOccupancy`); voce 17 aperta |
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
| Riferimento esterno | test_geometry (tabelle `_geom_` dei CIF COD), test_review (topologia dei policicli, PLATON, gemmi, conti a mano), test_checks (densità e Z dichiarati nei CIF) | Sì |
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

`check_riproducibilita.py`: 40 strutture COD (361 coppie, nessun errore) in 1 e in 4 processi, e due volte in 4, danno
risultati identici; il Monte Carlo con lo stesso seme dà le stesse σ. (La prima versione dello script prendeva anche i
file `._*.cif` che macOS crea sui dischi exFAT, analizzando così soprattutto righe d'errore: corretta e rieseguita.) Versioni esatte in `ambiente.txt` (voce 39).

### Dichiarazione sull'uso dell'IA

Bozza in `dichiarazione_ia.md` (voce 41).

## Sezione 10: sensibilità alle soglie (5 novembre 2026)

`check_sensibilita.py` sposta una soglia alla volta attorno al predefinito su 200 CIF COD (3816 coppie, 116 stack con
le soglie predefinite) e confronta le coppie, identificate da molecole, traslazione e frammenti π.

| Soglia (predefinito) | Valore | Cambiano tipo | Nuove | Perse | Stack | Strutture con stack diversi |
| --- | --- | --- | --- | --- | --- | --- |
| `max_dihedral` (20°) | 10 / 15 / 25 / 30° | 21 / 7 / 7 / 17 | 0 | 0 | 95 / 109 / 123 / 133 (−18 / −6 / +6 / +15 %) | 15 / 5 / 7 / 13 |
| `max_interplanar` (4,2 Å) | 3,8 / 4,0 / 4,4 / 4,6 | 5 / 2 / 1 / 2 | 0 | 0 | 111 / 114 / 117 / 118 (−4 / −2 / +1 / +2 %) | 5 / 2 / 1 / 2 |
| `contact_cutoff` (4,5 Å) | 4,0 / 5,0 | 0 / 0 | 0 / 440 | 640 / 3 | 114 / 116 | 2 / 0 |
| `ring_rms` (0,10 Å) | 0,05 / 0,15 | 9 / 1 | 174 / 212 | 203 / 133 | 117 / 115 | 4 / 4 |
| `molecule_rms` (0,15 Å) | 0,10 / 0,20 | 1 / 0 | 25 / 45 | 71 / 10 | 114 / 116 | 1 / 0 |
| `bond_tolerance` (0,4 Å) | 0,3 / 0,5 | 0 / 2 | 294 / 93 | 279 / 144 | 116 / 116 | 0 / 2 |
| `max_ring_size` (8) | 6 / 10 | 0 | 0 | 0 | 116 | 0 |
| `min_fragment_atoms` (5) | 4 / 6 | 0 | 66 / 6 | 0 / 379 | 116 / 110 | 0 / 4 |

Lettura:

- **`max_dihedral` è la soglia che decide**: ±5° spostano il 6 % degli stack, ±10° il 15–18 %. Nel paper conviene
  dare i risultati anche a 15° e 25°, o la distribuzione degli angoli, così che il lettore veda dove cade il taglio.
- **`max_interplanar`** pesa poco: ±0,2 Å cambiano al più il 2 % degli stack (le coppie quasi parallele sono quasi
  tutte ben sotto 4,2 Å o ben sopra).
- **`contact_cutoff` = 4,0 Å perde 2 stack** (coppie con distanza interplanare ≤ 4,2 Å ma contatto tra atomi pesanti
  > 4,0 Å): conferma la voce 6 (4,5 Å). A 5,0 Å arrivano 440 coppie, nessuna stack. Le 3 «perse» a 5,0 Å non cambiano
  il numero di stack: con più candidate la coppia rappresentativa di una classe di equivalenza può essere un'altra.
- **Nuclei e molecole** (`ring_rms`, `molecule_rms`, `bond_tolerance`): molte coppie nuove o perse, ma il numero di stack
  cambia al più dell'1–2 %. Una parte delle «nuove/perse» è solo rinumerazione: se cambia un nucleo cambiano gli indici
  dei frammenti, quindi il confronto per chiave le conta due volte.
- **`min_fragment_atoms` = 6** toglie 6 stack in 4 strutture: frammenti di 5 atomi (anelli a cinque, eteroaromatici).
  Il predefinito 5 li tiene; va detto nel paper.
- `max_ring_size` 6–10 non cambia nulla in questo campione.

## Sezione 10: legami e anelli contro RDKit (5 novembre 2026)

RDKit è indipendente da PISTAK in entrambi i passi: rifà i legami dalle sole coordinate (`DetermineConnectivity`, raggi
propri) e trova gli anelli con la SSSR simmetrizzata. Il confronto si fa due volte, sui legami di PISTAK (si prova solo
l'algoritmo degli anelli) e su quelli di RDKit (si prova la catena intera).

| Gruppo | Molecole | Legami uguali | Anelli uguali (legami di PISTAK) | Anelli uguali (legami di RDKit) |
| --- | --- | --- | --- | --- |
| organiche (solo non metalli) | 204 | 193 | 200 | 198 |
| con metalli | 83 | 71 | 76 | 72 |

- **Algoritmo**: dove le molecole differiscono, `small_rings` dà anelli *in più* (tutti i cicli più corti per legame,
  non una base): migliaia di triangoli nei cluster di metalli e nei polibromuri, e in 4112371 un anello a 6 di un
  sistema 4+6 fuso, rms 0,31 Å contro la soglia di 0,10 Å. Nessuno di questi è planare, quindi nessuno diventa un nucleo π.
- **Legami**: le 11 molecole organiche con legami diversi hanno contatti lunghi che RDKit unisce e PISTAK no (Si···N,
  B···F a ponte). Cambiano gli anelli solo in 2 (4063970 Si, 7025278 B₃F₁₁), nessuno aromatico.
- Conclusione: per le molecole organiche gli anelli planari, quelli che contano per l'analisi, coincidono con RDKit.

## Reti estese: coppie π nei polimeri contro PLATON (5 novembre 2026)

Voce 47. Oracolo: PLATON V-200626, `CALC ALL` predefinito, sui 26 CIF del campione COD con una rete di legami estesa e
anelli (`validation/platon_polimeri/`, con `LEGGIMI.md` e il confronto completo in `confronto.txt`). PLATON non ne legge
6; restano 20 strutture.

| Grandezza | Coppie abbinate | Scarto massimo PISTAK − PLATON (senza 7117874) | Mediano |
| --- | --- | --- | --- |
| Cg···Cg | 217 | 0,0013 Å | 0,0002 Å |
| α | 217 | 0,10° | 0,005° |
| β, γ | 217 | 0,07° | 0,02–0,03° |
| distanze perpendicolari | 217 | 0,0013 Å | 0,0001 Å |
| slippage | 103 | 0,0021 Å | 0,0003 Å |

- **7117874**: per una coppia PLATON stampa valori incompatibili tra loro (α 54,1°, β 34,2°, γ 0,0°, normale trasformata
  non unitaria). Ricalcolata con numpy dalle coordinate del CIF: α 88,1°, β 34,2°, γ 61,6°, come PISTAK. Difetto di
  PLATON con l'operazione esagonale −x+y, 1−x, z.
- **Coppie di PLATON che PISTAK non elenca**: 114 classi, tutte con un anello non planare per PISTAK o con un metallo
  (chelati, metallacicli), che non sono frammenti π.
- **Difetto trovato e corretto** (7024607): la visita in ampiezza che ricompone l'unità ripetuta arrivava a una piridina
  da due lati in celle diverse e la lasciava spezzata; la coppia N5···N5′ (−x, −y, −z, Cg···Cg 5,247 Å, β 58,0°)
  mancava. Ora l'unità tiene interi anelli e leganti; il test `test_polymer_rings_whole_in_repeat_unit` (CIF in
  `tests/data/cod/`) ritrova le cinque coppie di PLATON e fallisce con la vecchia ricostruzione.
- **Campione COD (200 CIF)**: nessuna struttura non polimerica cambia (celle, coppie, stack). Delle 39 reti estese, 27
  hanno ora coppie π (32 stack in queste strutture); le altre 12 non hanno frammenti π abbastanza vicini.
- Test sintetici: catena di poli(p-fenilene) con stack noto a 3,6 Å (trovato, segnato polimero, non esportato) e asse
  2₁ lungo la catena (le due coppie equivalenti diventano una ×2 solo con le immagini dei frammenti).
