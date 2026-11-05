# Dichiarazione sull'uso dell'intelligenza artificiale (bozza per il paper)

Bozza da adattare alle regole della rivista. Le parti tra [quadre] le deve completare o confermare l'autore: qui c'è
solo ciò che risulta dal repository e dalla revisione.

## Fatti verificabili nel repository

- Tutti i 29 commit della storia (dal 26 ottobre 2026) hanno la riga `Co-Authored-By` di un modello Claude di
  Anthropic: 28 Claude Opus 5.5, 1 Claude Sonnet 5.5, usati tramite Claude Code.
- Ogni commit e ogni push sono stati approvati esplicitamente dall'autore.
- La revisione del codice (`validation/REGISTRO.md`, sezioni 1–9 della checklist) ha controllato i risultati numerici
  con oracoli indipendenti dal codice di PISTAK: PLATON (coppie di anelli, ADDSYM), gemmi (simmetria, celle espanse),
  le tabelle `_geom_` dei CIF COD (geometria intramolecolare), calcoli indipendenti in numpy, la propagazione lineare
  delle incertezze, i testi delle fonti per ogni costante (Cordero 2008, Mantina 2009, Bondi 1964, Allen e Bruno 2010).
- Gli script di verifica sono in `validation/` e si possono rieseguire; l'ambiente usato è in `validation/ambiente.txt`.

## Testo proposto (italiano)

Il codice di PISTAK è stato scritto con l'assistenza di modelli linguistici di grandi dimensioni (Claude Opus 5.5 e
Claude Sonnet 5.5, Anthropic, tramite Claude Code), sotto la direzione dell'autore, che ha definito i requisiti
scientifici, [ha scelto i metodi e le soglie / altro] e ha approvato ogni modifica. Poiché il codice generato da un
modello può contenere errori plausibili ma non evidenti, la correttezza dei risultati non è stata affidata ai test
scritti insieme al codice: ogni quantità riportata è stata confrontata con un riferimento indipendente (PLATON, gemmi,
le tabelle geometriche dei CIF depositati, calcoli separati e i valori delle fonti bibliografiche per ogni costante),
come documentato nel registro di revisione distribuito con il codice. [Il testo dell'articolo è stato scritto
dall'autore / con l'assistenza di …]. L'autore è responsabile del contenuto del software e dell'articolo.

## Proposed text (English)

The PISTAK code was written with the assistance of large language models (Claude Opus 5.5 and Claude Sonnet 5.5,
Anthropic, via Claude Code) under the direction of the author, who defined the scientific requirements, [chose the
methods and thresholds / other] and approved every change. Because model-generated code can contain plausible but
non-obvious errors, the correctness of the results was not left to the tests written together with the code: every
reported quantity was compared with an independent reference (PLATON, gemmi, the geometry tables of the deposited
CIFs, separate calculations, and the literature values for every constant), as documented in the review log
distributed with the code. [The text of the article was written by the author / with the assistance of …]. The author
takes full responsibility for the content of the software and of the article.

## Da decidere

- Se la rivista chiede la dichiarazione nei metodi, nei ringraziamenti o in una sezione apposita.
- Se citare la versione esatta degli strumenti (Claude Code 2.1, date) e il registro di revisione come materiale
  supplementare.
- Se dichiarare anche l'uso per il testo dell'articolo e per le figure.
