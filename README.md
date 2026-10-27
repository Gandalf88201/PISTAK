# PISTAK — analisi dei π-stack da file CIF

Applicazione indipendente, in Python, nata dall'algoritmo dello strumento «Analizzatore di π-stack da file CIF» di
moleculesfriends.eu. Fa le stesse analisi del sito (distanza interplanare, slittamento longitudinale e laterale,
angolo tra i piani, contatti, sovrapposizione) e aggiunge: parametri in più, colonne di stack, tabella di
anelli (stile Janiak), elaborazione di cartelle intere, viste 3D e di proiezione dei parametri sincronizzabili,
e un sistema di plugin.

**Manuale d'uso**: [pistak/manual/it.md](pistak/manual/it.md) (in inglese [pistak/manual/en.md](pistak/manual/en.md)),
con un esempio guidato, la spiegazione di ogni vista e di ogni grandezza, le soglie, il terminale e la libreria. È anche
nella finestra: *Aiuto ▸ Manuale* (F1). Questo README è la guida per chi sviluppa.

## Avvio

```bash
./start.py                       # finestra grafica (con il cristallo dimostrativo)
./start.py examples/herringbone_anthracene.cif
./start.py analyze a.cif --csv out.csv --json out.json     # solo terminale (--block NOME per i file con più strutture)
./start.py batch cartella/ --csv riepilogo.csv              # molti file in parallelo
./start.py compare a.cif b.cif c.cif --ref 1 --csv confronto.csv   # confronto del motivo di stack
./start.py stats cartella/ --kind stack --corr --csv coppie.csv --summary descrizione.csv
./start.py qc a.cif --program psi4 --method sapt -o input/  # dimeri degli stack per ORCA / Psi4 / Gaussian
./start.py check a.cif                                       # controlli rapidi (formula, densità, H, simmetria mancata)
./start.py geom a.cif --csv geom.csv                         # legami, angoli, torsioni con le s.u.
./start.py periodic a.cif --program qe -o calcoli/           # la cella per VASP / Quantum ESPRESSO / CP2K / CRYSTAL
./start.py sessione.pistak                                   # riapre una sessione salvata
./start.py params                # parametri disponibili (con quelli dei plugin)
./start.py config                # dove stanno impostazioni personali, ultima sessione e plugin
```

La prima volta `start.py` crea `.venv/` con un Python 3.11–3.13 (se il tuo predefinito è più nuovo ne cerca uno adatto,
oppure `./start.py --python /percorso/python3.12`) e vi installa `requirements.txt`. Rilanciando è immediato; se
modifichi i requisiti reinstalla da solo. Altri pacchetti tuoi: scrivili in `requirements-user.txt`.
`./start.py --dev` aggiunge pytest, `./start.py --reinstall` rifà l'installazione.

## Come è fatto

```
pistak/
  core/      calcolo puro (numpy/scipy), nessuna GUI: usabile come libreria
    cif.py          lettura del CIF (blocchi data_, cella, simmetria, atomi, regola di disordine)
    spacegroups.py  operazioni di simmetria dal simbolo di Hall o H-M (quando il CIF non le elenca)
    structure.py    cella → atomi simmetrizzati → legami → molecole intere
    cores.py        anelli e frammenti π di ogni molecola (sistemi di anelli planari), asse lungo dagli atomi
    pairs.py        coppie di frammenti π vicini, equivalenza per simmetria, classificazione (stack / parallela / inclinata)
    stacks.py       colonne di stack (rango 0–3) dal grafo delle coppie
    ringpairs.py    parametri anello–anello (Cg···Cg, α, β, γ)
    config.py       soglie (Limits) e preferenze: pistak/default.toml + user.toml dell'utente
    errors.py       errori con codice stabile (CifError, StructureError, InputError)
    stages.py       stadi dell'analisi e soglie da cui dipende ciascuno (rianalisi parziale con cache)
    analysis.py     la catena completa, `reanalyse`, il PairContext e il riepilogo per molecola (MoleculeContact)
    uncertainty.py  incertezze dei parametri con un Monte Carlo sulle s.u. del CIF
    crystal.py      contenuto della cella (con tutte le componenti disordinate), densità, Z′, valori dell'affinamento
    validation.py   controlli rapidi (pre-checkCIF): avvisi A/B/C/G con codici stabili; ADDSYM con spglib facoltativo
    intramolecular.py  legami, angoli e torsioni delle molecole indipendenti, s.u. propagate dalle coordinate
    diagnosis.py    perché un'analisi non ha coppie (anelli saturi, nessun anello, nuclei lontani…), per le viste vuote
    elements.py     numeri atomici, masse, raggi, colori
  params/    registro dei parametri misurabili (Registry; builtin = sito, extra = fase 2, plugins) e registro delle
             grandezze dell'intera struttura (structure.py: StructureRegistry, `@structure_parameter`)
  io/        export CSV/JSON/XYZ con chiavi fisse; batch (parallelo), stats (statistiche), compare (confronto);
             qc = input del dimero per ORCA, Psi4, Gaussian; periodic = cella per VASP, Quantum ESPRESSO, CP2K, CRYSTAL
  gui/       finestra (PySide6): session.py = stato condiviso e sincronizzazione; views/ = le viste;
             limits_panel.py = cursori delle soglie in tempo reale;
             collection_window.py = confronto e statistiche; session_file.py = formato .pistak
  errors.py  messaggi degli errori nella lingua corrente
  labels.py  testi per l'utente costruiti dai risultati (relazione della coppia, tipo di colonna)
  i18n.py    traduzioni (i18n_en.py: il dizionario inglese)
  cli.py     comandi da terminale
  default.toml  soglie e preferenze predefinite
  manual/    il manuale (it.md, en.md) e le sue figure (img/<lingua>/), mostrato anche in Aiuto ▸ Manuale
plugins/     i tuoi parametri (caricati all'avvio)
tools/       manual_figures.py: rigenera le figure del manuale dalla finestra vera (senza schermo)
examples/    CIF sintetici di prova (non sperimentali)
tests/       pytest, compresa la parità con l'algoritmo originale del sito
```

Regola di dipendenza (controllata da `tests/test_structure.py`): `core` importa solo sé stesso, niente traduzioni;
`params` usa `core` (e `i18n` per le etichette dei parametri); `io` usa `core` e `params` e scrive file con
intestazioni fisse in inglese, uguali in ogni lingua; `gui` e `cli` usano tutto e traducono con `labels`, `errors`, `i18n`.

**Impostazioni personali.** Stanno fuori dal progetto, nella cartella dell'utente (`./start.py config` la mostra):
macOS `~/Library/Application Support/PISTAK/`, Windows `%APPDATA%\PISTAK\`, Linux `~/.config/pistak/`, oppure quella
indicata da `PISTAK_CONFIG_DIR`. Contiene `user.toml` (soglie e preferenze che sovrascrivono `pistak/default.toml`),
`last.pistak` (l'ultima sessione) e `plugins/`. Le versioni precedenti usavano `config/` dentro il progetto: alla prima
apertura `user.toml` e `last.pistak` vengono copiati da lì nella cartella nuova.

## Finestre e sincronizzazione

Ogni vista (tabella, scatter di due parametri, istogramma, 3D, proiezione, dettaglio) è un pannello che puoi
trascinare, impilare o staccare in una finestra propria; da *Vista ▸ Nuova vista* ne apri altre copie. In alto a
sinistra di ogni vista c'è il selettore 🔗:

* **Gruppo A/B/C**: le viste dello stesso gruppo seguono la stessa coppia selezionata e lo stesso filtro. Clic su un punto
  dello scatter → la tabella evidenzia la riga, il 3D e la proiezione mostrano quella coppia.
* **Non sincronizzata**: la vista tiene una selezione propria (utile per confrontare due coppie affiancate).

Nello scatter «Selezione a riquadro» filtra le coppie dentro il rettangolo; nell'istogramma un clic su una barra filtra
quell'intervallo; il filtro vale per il gruppo (la tabella mostra solo quelle coppie).

## Confronto e statistiche su più strutture

*Analisi ▸ Confronto e statistiche…* (Ctrl+Maiusc+C) apre una finestra dove si trascinano file e cartelle di CIF
(analisi in parallelo). Tutte le coppie di tutti i file finiscono nella stessa tabella:

* **Strutture**: per ciascuna la coppia π-stack più vicina, le differenze (Δ) dal file di riferimento (tasto destro ▸
  «Imposta come riferimento») e la *distanza del motivo*: distanza in Å, nello spazio (interplanare, slittamento L,
  slittamento W), tra i due motivi di stack più simili. Doppio clic: apre quel file nella finestra principale.
* **Distribuzioni**: istogramma di qualunque parametro, impilato per tipo di coppia o per file.
* **Dispersione**: due parametri a scelta, colore per tipo, file o un terzo parametro; un clic su un punto e «Apri nella
  finestra principale» porta direttamente a quella coppia.
* **Correlazioni**: matrice di Pearson tra i parametri spuntati (solo coppie con valori finiti; NaN se un parametro è costante).
* **Riepilogo**: media, deviazione standard, minimo, quartili e massimo; più una riga per file.

Il filtro «Coppie» (π-stack / parallele / inclinate / tutte) e «Pesa per coppie equivalenti» valgono per tutte le schede.
Da terminale: `compare` e `stats` (vedi sopra). Tutto si esporta in CSV.

## Sessioni

*File ▸ Salva la sessione…* scrive un file `.pistak` (JSON leggibile) con: il CIF (incorporato, quindi la sessione si apre
anche se sposti il file), le soglie, la coppia selezionata e il filtro di ogni gruppo, le viste con le loro impostazioni
(assi scelti, colonne della tabella, opzioni e camera del 3D…), la disposizione e la geometria dei pannelli, e l'elenco dei
file della finestra di confronto. Si riapre da *File ▸ Apri una sessione…*, trascinando il file sulla finestra o con
`./start.py sessione.pistak`. Alla chiusura l'ultima sessione si salva da sola in `last.pistak`, nella cartella delle impostazioni personali
(*File ▸ Riapri l'ultima sessione*). Una vista nuova che scrivi tu salva le sue impostazioni implementando `extra_state()` e
`apply_state()`; le chiavi sconosciute di un file più nuovo si ignorano.

## Lingua e dimensioni delle finestre

* **Italiano / English:** menu *Lingua* (la finestra si ricostruisce da sola conservando sessione, selezione e pannelli), oppure
  `./start.py --lang en`, `PISTAK_LANG=en`, o `language = "en"` in `user.toml` (sezione `[gui]`; `auto` = lingua del sistema).
  I testi sono scritti in italiano dentro `tr("...")` e tradotti da `pistak/i18n_en.py`; `python -m pistak.i18n_check` elenca
  quelli senza traduzione. I plugin traducono le proprie etichette con `pistak.i18n.add_translations("en", {...})`.
* **Tema chiaro / scuro:** menu *Vista ▸ Tema* (automatico come il sistema, chiaro, scuro), oppure `background` in
  `user.toml` (sezione `[gui]`).
* **Finestre:** la finestra principale non supera mai lo schermo (94% dell'area utile), si può rimpicciolire fino a circa
  640×420 e ingrandire a piacere; le barre degli strumenti delle viste vanno a capo invece di allargare la finestra; le
  sessioni salvate su uno schermo più grande vengono riportate dentro lo schermo.

## Estendere

* **Un parametro nuovo**: copia `plugins/example_param.py` (in `plugins/` del progetto o in quella della cartella delle
  impostazioni personali), cambia la funzione, *Analisi ▸ Ricarica i plugin*. Compare
  ovunque (tabella, assi, export). Contesto disponibile: vedi la docstring del file e `core/analysis.py::PairContext`.
  Da codice si può anche calcolare un insieme diverso di parametri con un registro proprio: `reg = Registry()`,
  `@reg.parameter(...)`, `compute_all(analisi, reg)`; gli errori dei parametri finiscono in `analisi.param_errors`.
* **Una grandezza dell'intera struttura**: `@structure_parameter("chiave", "Etichetta", unit=…, group=…)` su una
  funzione che riceve l'`Analysis` (anche nei plugin). Compare nel Dettaglio ▸ Struttura, in `analyze`, nel JSON
  (`summary.structure`) e, se la aggiungi a `STRUCTURE_KEYS` in `io/batch.py`, nel riepilogo di `batch`.
* **Un controllo rapido nuovo**: `@check("nome")` in `core/validation.py` restituisce una lista di `Alert(codice,
  livello, valori)`; il testo per l'utente va in `ALERT_TEXT` di `labels.py` (e la traduzione in `i18n_en.py`).
* **Una vista nuova**: eredita da `gui/views/base.py::PairView` (la sincronizzazione è già fatta) e aggiungila a
  `VIEW_TYPES` in `gui/views/__init__.py`.
* **Soglie**: *Analisi ▸ Soglie…* o `user.toml` nella cartella delle impostazioni personali (`[limits]`).
* **Un altro renderer 3D**: la scena è descritta da `gui/scene.py` (solo numpy); `views/view3d.py` la disegna con PyVista.

## Convenzioni dell'analisi

Tutte si scelgono in *Analisi ▸ Soglie…*, in `user.toml` (`[limits]`) o da terminale (`--pair-mode`, `--disorder`,
`--overlap-measure`, `--min-overlap`, `--max-slip`). *Valori del sito* le riporta a quelle del sito (`Limits.site()`).

* **Coppie** (`pair_mode`). `canonical` (predefinito): ogni contatto tra due molecole compare una volta; due contatti sono
  la stessa coppia se un'operazione del gruppo spaziale porta l'uno nell'altro (non perché i valori arrotondati
  coincidono); l'asse lungo L è la bisettrice degli assi lunghi delle due molecole, quindi slittamento longitudinale e
  laterale non dipendono da quale si chiama A; «Coppie equivalenti» conta i contatti per cella. Le celle vicine esplorate
  si allargano da sole quanto serve (celle piccole o molto oblique), con `neighbour_shells` come minimo. `site`: come il
  sito — ogni contatto due volte (a→b e b→a, con slittamenti che possono differire), equivalenza per valori arrotondati,
  sempre ±`neighbour_shells` celle.
* **Frammenti π** (`fragments`, `min_fragment_atoms`). `all` (predefinito): ogni sistema di anelli planari condensati
  di una molecola è un frammento (π1 = il più grande, il nucleo π di sempre; gli altri con almeno `min_fragment_atoms`
  atomi). Le coppie sono tra frammenti: per due molecole a contatto si misurano π1–π1 e le altre combinazioni di
  frammenti a meno di `contact_cutoff` + 1 Å; con più frammenti nella stessa molecola anche le coppie intramolecolari.
  La relazione in tabella dice quali («molecola 2 [1 0 0] · π2–π1», «intramolecolare · π1–π3»); i parametri
  `fragment_a`/`fragment_b` li danno come numeri. Equivalenza per simmetria anche tra frammenti. `largest`: solo π1
  (il sito). Il **riepilogo per molecola** (*Dettaglio ▸ Contatti*, `molecule_contacts` nel JSON) raccoglie per ogni
  contatto tra molecole le sue coppie di frammenti, quanti stack, il contatto più corto e la coppia migliore (doppio
  clic la seleziona).
* **Asse lungo** (`long_axis_atoms`). Vuoto: la direzione di massima estensione del frammento; per nuclei quasi
  isotropi (parametro *Isotropia del nucleo* vicino a 1: coronene, porfirine) è mal definita e così la divisione dello
  slittamento in L e W. Con due etichette del CIF («C1,C9») l'asse va dal primo atomo al secondo, proiettato sul piano,
  per i frammenti che li contengono (o per π1). Se un'etichetta compare più volte nella molecola (molecola su un
  elemento di simmetria) si prende la coppia più lontana: «C1,C1» va bene.
* **Slittamenti con segno**: `long_slip_signed`, `lat_slip_signed` (componenti di s da A verso B lungo L e W) e
  `slip_direction` (−180…180°). Il segno ha senso fisico se l'asse è definito dagli atomi e la molecola non sta su un
  centro di inversione (che rovescia l'asse); scambiando A e B cambia segno.
* **Incertezze (σ)**: *Analisi ▸ Incertezze (Monte Carlo)…*, `analyze --mc 200` (`--mc-seed`, `--mc-no-cell`) o
  `pistak.api.uncertainty(analisi, 200)`. Dalle s.u. del CIF (le cifre tra parentesi delle coordinate e della cella)
  si generano N geometrie perturbate: lo spostamento casuale di ogni atomo dell'unità asimmetrica passa alle sue copie
  per simmetria attraverso l'operazione che le genera, si ricalcolano piani e assi dei frammenti e si rimisurano tutte
  le coppie. La σ di ogni parametro del registro (plugin compresi) è la deviazione standard sui campioni; in più, per
  ogni coppia, la frazione di campioni in cui il tipo (stack/parallela/inclinata) resta lo stesso. Il Dettaglio mostra
  «valore ± σ», il CSV aggiunge le colonne `<chiave>_su` e `kind_agreement`, il JSON `su` per coppia. Approssimazioni:
  coordinate indipendenti (nessuna covarianza), posizioni speciali non vincolate, atomi senza s.u. fermi. Una σ
  grande sugli slittamenti L/W di un nucleo quasi isotropo (fenile, coronene) dice che l'asse lungo non è definito.
* **Dimero per calcoli quantistici**: *Esporta ▸ Coppia selezionata per calcoli quantistici…* o
  `./start.py qc file.cif` (`--pair 1,3-5` o `--kind stack|parallel|tilted|all`, `-o` file o cartella; senza `-o` e
  con una sola coppia l'input va sul terminale). I monomeri sono le **due molecole intere** della coppia nella geometria
  del cristallo (B già traslata nella cella giusta), non i soli frammenti π: per un contatto tra frammenti l'energia è
  comunque quella dell'intero dimero. Gli idrogeni da raggi X sono troppo corti: con la normalizzazione (predefinita,
  `--raw-h` per toglierla) ogni X–H viene allungato lungo il legame alle distanze medie da neutroni (Allen 1987/2010:
  C–H 1,083, N–H 1,009, O–H 0,983, B–H 1,19, S–H 1,338 Å); gli atomi pesanti non si muovono.
  - **DFT-D** (`--method dft-d`, predefinito B3LYP-D3(BJ)/def2-TZVP; `--functional`, `--dispersion d3bj|d3zero|d4|none`,
    `--basis`): energia di interazione con counterpoise (Boys–Bernardi). ORCA: tre lavori in un file (`$new_job`):
    dimero, A con B fantasma (`El :`), B con A fantasma; E_int = E(AB) − E(A^AB) − E(B^AB). Psi4:
    `energy('b3lyp-d3bj', bsse_type='cp')`, che fa tutto da sé. Gaussian: `Counterpoise=2` con `(Fragment=n)` sugli
    atomi (D4 non c'è: `EmpiricalDispersion=GD3BJ|GD3`). `--no-cp` dà il solo calcolo del dimero.
  - **SAPT** (`--method sapt`, solo Psi4; `--sapt-level sapt0|ssapt0|sapt2|sapt2+|sapt2+(3)|sapt2+3`, base predefinita
    jun-cc-pVDZ, la scelta usuale per SAPT0): energia separata in elettrostatica, scambio, induzione e dispersione.
  - Cariche e molteplicità dei monomeri: `--charges 0,0`, `--mult 1,1` (il dimero ha la somma delle cariche e, per
    accoppiamento ad alto spin, 2S = 2S_A + 2S_B). `--nprocs`, `--mem` (GB) vanno nelle intestazioni dei programmi.
  - Il file inizia con un commento che dice da dove viene (CIF, blocco, coppia, simmetria, distanza, idrogeni) ed è
    tutto ASCII. Le coppie intramolecolari (tra frammenti della stessa molecola) non si esportano: non sono un dimero.
    Le parole chiave seguono le versioni correnti (ORCA 5–6, Psi4 1.7+, Gaussian 16): controllale con la tua.
* **Analisi a stadi e soglie in tempo reale**: la catena è divisa in stadi (`core/stages.py`) — lettura del CIF
  (`disorder`), struttura (legami, cella), frammenti π (anelli, asse lungo), contatti (`contact_cutoff`,
  `neighbour_shells`, `pair_mode`: la misura di tutte le coppie, la parte costosa) e coppie (il criterio di stack:
  `max_dihedral`, `max_interplanar`, `min_overlap`, `max_slip`, `overlap_measure`). Cambiando le soglie si rifà solo
  dal primo stadio toccato (`pistak.api.reanalyse(analisi, soglie)`; la finestra lo fa da sola): con quelle del
  criterio di stack le coppie già misurate si riclassificano senza rifare la geometria (pochi ms anche con decine
  di coppie), e i parametri del registro dichiarati `geometric=True` (tutti quelli di `params/extra.py`) si riusano.
  *Analisi ▸ Soglie in tempo reale* (Ctrl+Shift+L) apre i cursori: tabella, scatter (con le soglie come linee
  tratteggiate), istogramma e colonne seguono il cursore, selezione e filtri restano sulle stesse coppie; il contatto
  massimo rimisura e si applica al rilascio, in un thread a parte. *Salva come predefinite* scrive `user.toml`. Le
  incertezze del Monte Carlo vanno ricalcolate dopo un cambio di soglie. Un plugin che dipende solo dalla
  geometria della coppia può dichiararlo con `@parameter(..., geometric=True)`.
* **Atomi disordinati** (`disorder`): `major` (predefinito) tiene, per ogni assemblaggio (`_atom_site_disorder_assembly`),
  il gruppo con l'occupanza media più alta; `first` il gruppo 1 (il sito); `all` tutti; un numero, quel gruppo.
* **Criterio di stack**: oltre ad angolo tra i piani e distanza interplanare, la sovrapposizione dei nuclei proiettati
  deve superare `min_overlap`, misurata come frazione di atomi dentro il contorno dell'altro nucleo (`atoms`, il sito) o
  come area comune / area del nucleo minore (`area`); con `max_slip` > 0 si scartano anche le coppie troppo slittate.
* **File con più strutture**: ogni blocco `data_` con atomi è una struttura. La finestra chiede quale aprire, il confronto
  e `batch`/`stats`/`compare` le analizzano tutte (nome «file.cif [blocco]»), `analyze` prende la prima o quella di `--block`.
* **Operazioni di simmetria mancanti**: se il CIF non le elenca si ricavano dal simbolo di Hall (`_space_group_name_Hall`)
  o dal simbolo H-M (tabella dei gruppi più comuni in `core/spacegroups.py`; con il pacchetto facoltativo `gemmi`
  installato, tutti). Il riepilogo JSON dice da dove vengono (`symops_source`).

## Parità con il sito

`tests/test_parity.py` confronta il porting, con `Limits.site()`, con i numeri prodotti dai sorgenti TypeScript originali
del sito su 6 cristalli sintetici (tutti e tre i tipi di coppia). Per rigenerare i valori di riferimento:
`SITE_DIR=/percorso/SITOMIO node tests/parity/dump_site.js`.

Differenze volute anche in modalità `site`: il limite di atomi in cella è 6000 (2500 nel browser); i versi degli assi
(L, W) sono fissati in modo riproducibile; in più ci sono i parametri aggiuntivi e le colonne.

## Sviluppo

`./start.py --dev --setup-only` installa pytest e ruff in `.venv/`. Il progetto è anche un pacchetto installabile
(`pyproject.toml`, comando `pistak`): `pip install -e ".[gui,dev]"`.

```bash
.venv/bin/python -m pytest -m "not gui"     # solo il calcolo, pochi secondi
.venv/bin/python -m pytest                  # tutto, comprese le finestre (offscreen)
.venv/bin/ruff check .                      # controllo statico (regole in pyproject.toml)
.venv/bin/python tools/manual_figures.py    # rigenera le figure del manuale (dopo cambiamenti dell'interfaccia)
```

**Manuale.** `pistak/manual/it.md` ed `en.md` vanno tenuti allineati: `tests/test_manual.py` controlla che abbiano gli
stessi capitoli e le stesse figure, che i link interni portino a un titolo e che l'esempio della libreria giri davvero.
Le figure si rigenerano con lo script sopra (solo i CIF di `examples/`, PNG a 256 colori).

Marcatori: `gui` (finestre), `parity` (confronto con il sito), `snapshot` (istantanee). **Istantanee**:
`tests/snapshots/site/*.json` e `tests/snapshots/canonical/*.json` contengono tutti i risultati sugli esempi nei due modi
(molecole, nuclei, coppie con tutti i parametri, coppie di anelli, colonne); ogni cambiamento di un numero fa fallire
`pytest -m snapshot`. Dopo un cambiamento voluto: `pytest -m snapshot --update-snapshots` e si controlla il diff. Per
strutture reali non distribuibili (CSD): `PISTAK_SNAPSHOT_DIR=/cartella/cif pytest -m snapshot --update-snapshots`
scrive in `tests/snapshots/local/<modo>/`, che resta fuori dal repository.

## Note per volumi exFAT/FAT (macOS)

macOS crea file `._nome` accanto ai file copiati. `start.py` li rimuove dalla `.venv` dopo ogni installazione (rompono
matplotlib) e i caricatori di plugin/CIF li ignorano. Se preferisci, sposta il progetto su un disco APFS.
