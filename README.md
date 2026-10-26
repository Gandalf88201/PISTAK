# PISTAK — analisi dei π-stack da file CIF

Applicazione indipendente, in Python, nata dall'algoritmo dello strumento «Analizzatore di π-stack da file CIF» di
moleculesfriends.eu. Fa le stesse analisi del sito (distanza interplanare, slittamento longitudinale e laterale,
angolo tra i piani, contatti, sovrapposizione) e aggiunge: parametri in più, colonne di stack, tabella di
anelli (stile Janiak), elaborazione di cartelle intere, viste 3D e di proiezione dei parametri sincronizzabili,
e un sistema di plugin.

## Avvio

```bash
./start.py                       # finestra grafica (con il cristallo dimostrativo)
./start.py examples/herringbone_anthracene.cif
./start.py analyze a.cif --csv out.csv --json out.json     # solo terminale
./start.py batch cartella/ --csv riepilogo.csv              # molti file in parallelo
./start.py compare a.cif b.cif c.cif --ref 1 --csv confronto.csv   # confronto del motivo di stack
./start.py stats cartella/ --kind stack --corr --csv coppie.csv --summary descrizione.csv
./start.py sessione.pistak                                   # riapre una sessione salvata
./start.py params                # parametri disponibili (con quelli dei plugin)
```

La prima volta `start.py` crea `.venv/` con un Python 3.11–3.13 (se il tuo predefinito è più nuovo ne cerca uno adatto,
oppure `./start.py --python /percorso/python3.12`) e vi installa `requirements.txt`. Rilanciando è immediato; se
modifichi i requisiti reinstalla da solo. Altri pacchetti tuoi: scrivili in `requirements-user.txt`.
`./start.py --dev` aggiunge pytest, `./start.py --reinstall` rifà l'installazione.

## Come è fatto

```
pistak/
  core/      calcolo puro (numpy/scipy), nessuna GUI: usabile come libreria
    cif.py          lettura del CIF (cella, simmetria, atomi, disordine)
    structure.py    cella → atomi simmetrizzati → legami → molecole intere
    cores.py        anelli e nucleo π di ogni molecola
    pairs.py        coppie di molecole vicine e loro classificazione (stack / parallela / inclinata)
    stacks.py       colonne di stack (rango 0–3) dal grafo delle coppie
    ringpairs.py    parametri anello–anello (Cg···Cg, α, β, γ)
    config.py       soglie (Limits) lette da config/default.toml e config/user.toml
    analysis.py     la catena completa e il PairContext
  params/    registro dei parametri misurabili (builtin = sito, extra = fase 2, plugins)
  io/        export CSV/JSON/XYZ; batch (parallelo), stats (statistiche), compare (confronto)
  gui/       finestra (PySide6): session.py = stato condiviso e sincronizzazione; views/ = le viste;
             collection_window.py = confronto e statistiche; session_file.py = formato .pistak
  cli.py     comandi da terminale
plugins/     i tuoi parametri (caricati all'avvio)
config/      soglie e preferenze
examples/    CIF sintetici di prova (non sperimentali)
tests/       pytest, compresa la parità con l'algoritmo originale del sito
```

Regola di dipendenza: `core` non importa nient'altro; `params` usa `core`; `io` e `gui` usano entrambi.

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
`./start.py sessione.pistak`. Alla chiusura l'ultima sessione si salva da sola in `config/last.pistak`
(*File ▸ Riapri l'ultima sessione*). Una vista nuova che scrivi tu salva le sue impostazioni implementando `extra_state()` e
`apply_state()`; le chiavi sconosciute di un file più nuovo si ignorano.

## Lingua e dimensioni delle finestre

* **Italiano / English:** menu *Lingua* (la finestra si ricostruisce da sola conservando sessione, selezione e pannelli), oppure
  `./start.py --lang en`, `PISTAK_LANG=en`, o `language = "en"` in `config/user.toml` (sezione `[gui]`; `auto` = lingua del sistema).
  I testi sono scritti in italiano dentro `tr("...")` e tradotti da `pistak/i18n_en.py`; `python -m pistak.i18n_check` elenca
  quelli senza traduzione. I plugin traducono le proprie etichette con `pistak.i18n.add_translations("en", {...})`.
* **Finestre:** la finestra principale non supera mai lo schermo (94% dell'area utile), si può rimpicciolire fino a circa
  640×420 e ingrandire a piacere; le barre degli strumenti delle viste vanno a capo invece di allargare la finestra; le
  sessioni salvate su uno schermo più grande vengono riportate dentro lo schermo.

## Estendere

* **Un parametro nuovo**: copia `plugins/example_param.py`, cambia la funzione, *Analisi ▸ Ricarica i plugin*. Compare
  ovunque (tabella, assi, export). Contesto disponibile: vedi la docstring del file e `core/analysis.py::PairContext`.
* **Una vista nuova**: eredita da `gui/views/base.py::PairView` (la sincronizzazione è già fatta) e aggiungila a
  `VIEW_TYPES` in `gui/views/__init__.py`.
* **Soglie**: *Analisi ▸ Soglie…* o `config/user.toml` (`[limits]`).
* **Un altro renderer 3D**: la scena è descritta da `gui/scene.py` (solo numpy); `views/view3d.py` la disegna con PyVista.

## Parità con il sito

`tests/test_parity.py` confronta il porting con i numeri prodotti dai sorgenti TypeScript originali del sito su 6 cristalli
sintetici (tutti e tre i tipi di coppia). Per rigenerare i valori di riferimento:
`SITE_DIR=/percorso/SITOMIO node tests/parity/dump_site.js`.

Differenze volute: il limite di atomi in cella è 6000 (2500 nel browser); i versi degli assi (L, W) sono fissati in
modo riproducibile; in più ci sono i parametri aggiuntivi e le colonne.

## Sviluppo

`./start.py --dev --setup-only` installa pytest e ruff in `.venv/`. Il progetto è anche un pacchetto installabile
(`pyproject.toml`, comando `pistak`): `pip install -e ".[gui,dev]"`.

```bash
.venv/bin/python -m pytest -m "not gui"     # solo il calcolo, pochi secondi
.venv/bin/python -m pytest                  # tutto, comprese le finestre (offscreen)
.venv/bin/ruff check .                      # controllo statico (regole in pyproject.toml)
```

Marcatori: `gui` (finestre), `parity` (confronto con il sito), `snapshot` (istantanee). **Istantanee**:
`tests/snapshots/examples/*.json` contengono tutti i risultati attuali sugli esempi (molecole, nuclei, coppie con tutti i
parametri, coppie di anelli, colonne); ogni cambiamento di un numero fa fallire `pytest -m snapshot`. Dopo un
cambiamento voluto: `pytest -m snapshot --update-snapshots` e si controlla il diff. Per strutture reali non
distribuibili (CSD): `PISTAK_SNAPSHOT_DIR=/cartella/cif pytest -m snapshot --update-snapshots` scrive in
`tests/snapshots/local/`, che resta fuori dal repository.

## Note per volumi exFAT/FAT (macOS)

macOS crea file `._nome` accanto ai file copiati. `start.py` li rimuove dalla `.venv` dopo ogni installazione (rompono
matplotlib) e i caricatori di plugin/CIF li ignorano. Se preferisci, sposta il progetto su un disco APFS.
