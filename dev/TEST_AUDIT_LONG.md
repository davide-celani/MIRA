# Audit dei test longitudinali — 26 ottobre 2026

## Estensione metodologica successiva — risultato finale

Questa sezione documenta l'ampliamento richiesto dopo l'audit originario,
partendo da HEAD `196a357ac47066a49e67470b8ecac0322edac90f` e da un working
tree pulito. I numeri e le dichiarazioni sull'ambiente nelle sezioni successive
descrivono l'audit precedente. Per lo stato attuale valgono i risultati qui
riportati, comprese le nuove dipendenze installate con autorizzazione.

**Aggiunte tre prospettive, senza rimuovere i metodi esistenti**, con la
catena analisi → nomenclatura/stampa/summary → report/RDS/CSV → testthat.
La suite ordinaria finale contiene **364 blocchi `test_that`, 11.223
aspettative passate, zero failure/error/warning e quattro skip**, sia dai
sorgenti sia dal namespace del pacchetto appena costruito e installato.
Sono 39 nuovi blocchi rispetto alla baseline di 325; il conteggio delle
aspettative include anche i controlli esistenti applicati alle nuove tabelle
e non rappresenta altrettante verifiche scientifiche indipendenti.

| Percorso del risultato | Aggiunta e utilità | Assunzioni e interpretazione |
|---|---|---|
| `advanced_tests$kenward_roger` | Test F globali Kenward–Roger del modello misto primario, con F scalata, df corretti, fattore di scala e campione effettivo | Inferenza gaussiana basata sul modello; TIME e ARM comprendono le rispettive interazioni. Non è una correzione robusta della covarianza. |
| `advanced_models$mmrm` | Modello marginale gaussiano con covarianza non strutturata, varianze specifiche per visita, REML e Satterthwaite; coefficienti, test congiunti, medie marginali e contrasti | Covarianza comune fra bracci, baseline selezionata inclusa come risposta, missingness delle risposte ignorabile MAR e specificazione adeguata; non è una sensitivity MNAR. |
| `advanced_tests$rank_longitudinal` | Effetti relativi marginali non pesati `H0p`, SE/IC logit puntuali e test longitudinali ATS/WTS tramite nparLD | Analisi non aggiustata per covariate, risposte mancanti MCAR, approssimazioni asintotiche. Un effetto > 0,5 indica risposte maggiori, indipendentemente dalla direzione clinica. |

Gli `emmeans` preesistenti potevano già utilizzare Kenward–Roger: la nuova
aggiunta è il modulo esplicito di **test globali**. MMRM distingue la struttura
marginale non strutturata dai modelli con intercetta/pendenza casuale, dai
precedenti `nlme` CS/AR1 e dai GEE. Il modulo rank aggiunge un'analisi
fattoriale marginale con profili incompleti, distinta dal Friedman su profili
completi. Le tre prospettive hanno obiettivi e, talvolta, campioni differenti;
non devono essere utilizzate come voto di maggioranza o scelta del p minimo.

### Adattamento, campioni e risultati non stimabili

- I tre moduli seguono lo switch `model` già risolto dalla funzione; nessun
  nuovo argomento pubblico. `analyses = "none"` li disabilita in assenza di
  override espliciti e `inspect_only` non esegue fit. Ogni adattatore è protetto
  separatamente: un errore non elimina i risultati degli altri metodi.
- KR riutilizza il campione del modello primario. Richiede convergenza,
  nonsingolarità e varianza residua positiva; conserva test nativi e matrici
  di restrizione. Le ipotesi ridotte per colonne eliminate sono etichettate
  con coefficienti richiesti/testati/esclusi, senza presentarle come complete.
- Rank usa il calendario originale soggetto × visita, con righe NA esplicite,
  senza escludere soggetti per covariate aggiuntive mancanti. Esclude i profili
  interamente non osservati. L'adattatore richiede almeno due osservazioni
  finite per cella e soggetti sufficienti per braccio; rifiuta duplicati
  soggetto/visita e ipotesi degeneri. L'assenza completa di variazione entro
  soggetto non produce un test temporale apparentemente valido.
- MMRM usa le righe complete per i predittori selezionati e mantiene tutti i
  livelli delle visite richieste. Una visita interamente assente, anche dopo
  il filtraggio delle covariate, rende il modulo non disponibile: non viene
  reinterpretato un follow-up come nuova baseline. Sono verificati disegno a
  rango pieno, unicità soggetto/visita, sovrapposizione di ogni coppia di visite,
  convergenza e rango numerico delle covarianze. Il fallback riguarda soltanto
  l'inizializzazione; non cambia il modello o il metodo inferenziale.
- Per nparLD la covarianza nativa della stima scalata per sqrt(N) viene divisa
  per N, numero di soggetti indipendenti osservati. Le SE pubblicate sono
  `sqrt(diag(Sigma)/N)`, non le SE native scalate. Il fit nativo resta intatto.
  ATS usa il denominatore convenzionale del motore `1e10`; WTS non ha df al
  denominatore. Non è stato attivato un test di permutazione.
- I criteri MMRM restano quelli nativi: AIC penalizza i parametri di covarianza
  e BIC usa il numero di soggetti. Metadati e colonna `criterion_convention`
  chiariscono che non sono armonizzati con lme4/nlme per un ranking comune.

### Catena del report e correzioni dimostrate

I nuovi campi statistici seguono il prefisso esatto dell'outcome dichiarato e
le mappe canoniche; identità, visite, campioni e metodi restano strutturali.
Stampa e summary comprendono le nuove tabelle. `sensitivity` aggiunge righe
finite con ipotesi, campione e percorso; le famiglie dei contrasti MMRM hanno
riferimenti alle tabelle canoniche e aggiustamenti separati. Nessuna nuova
famiglia aggrega tutti i p-value globali o metodi diversi.

Il report riconosce KR, nparLD e MMRM come oggetti nativi terminali, preserva
le informazioni del campione e le assunzioni, esporta le tabelle standardizzate
e mostra il motivo dei moduli richiesti ma non disponibili. Non attraversa
ricorsivamente l'infrastruttura TMB o le strutture native nei CSV.

Le prove hanno rilevato e corretto il filtraggio indebito del campione rank
per covariate/visite mancanti e il possibile cambio implicito della baseline
MMRM. Hanno inoltre esposto un difetto precedente della stampa GEE:
`geepack::print.summary.geeglm()` modificava `options(digits)`, alterando le
tabelle successive. Il reporter ora ripristina la precisione del chiamante
anche in caso di errore del printer. Prima della guardia il test dedicato
aveva nove failure; dopo passa tutte le 20 aspettative. Due precedenti test
di formattazione fallivano nella suite completa e ora passano.

Un confronto con il risultato salvato prima dell'estensione conferma valori
identici per tutti i 51 percorsi metrici preesistenti, campione e coefficienti
del modello, configurazione, dati longitudinali e vecchie righe sensitivity.
Si escludono dal confronto soltanto la nuova riga MMRM e la colonna esplicativa
del confronto modelli, le righe sensitivity aggiunte e l'identità degli
ambienti formula/terms ricreati dalla lettura RDS; i dati e gli altri attributi
sono confrontati esattamente.

### Forza dei nuovi test

Le verifiche combinano aritmetica indipendente sui soggetti e confronti con
le API reali, con ipotesi/campioni indicizzati per nome:

- KR: somme dei quadrati bilanciate, df, restrizioni TIME/ARM/interazione,
  F scalata rispetto a F non scalata, dati incompleti, ipotesi ridotte e
  fallimenti parziali. Nel test report la differenza numerica dal p analitico
  è circa `6,97e-12`; il limite assoluto esplicito è `1e-10`, con controllo
  separato e più stretto contro il motore.
- MMRM: medie di cella, covarianza campionaria e sua trasformazione, bracci
  sbilanciati e campioni incompleti. Nel caso manuale a otto soggetti il
  contrasto finale−baseline è `6,125`, SE circa `1,35538`, df `7`: usare le
  osservazioni ripetute come unità indipendenti darebbe df `21` ed è respinto.
- Rank: oracolo manuale a tre soggetti con effetti `5/12` e `7/12`, SE `1/12`,
  covarianza nota, IC logit e ATS/WTS; identità/associazione dei soggetti,
  ties, profili mancanti, statica completa e colonne nuisance.
- API pubblica: separazione dei campioni rank/MMRM, visite iniziale e finale
  interamente NA, disable/inspect senza fit, stampa/summary senza mutazioni
  o consumo RNG e conservazione dei risultati precedenti.
- Report: numeri e identità attraverso RDS, CSV/manifest e runtime esportato
  con tabelle HTML realmente generate da kable; nomi arbitrari con spazi,
  cronologia non lessicografica, mappe, famiglie e diagnostica reason-only.
  Un processo R separato elabora il QMD/runtime/RDS salvato senza caricare
  MIRA, mantenendo tutti i nuovi moduli e `digits`.

Mutazioni controllate, mai lasciate nei sorgenti finali: F KR non scalata
(tre failure); df MMRM calcolati sulle righe con p ricalcolati (quattro);
covarianza MMRM attribuita alle visite sbagliate (due); SE rank native scalate
(una); scambio TIME/ARM rank (quattro); associazione errata dei soggetti rank
con margini/N invariati (undici). Sono prove mirate di sensibilità dei test,
non una percentuale di mutation coverage.

### Esecuzioni finali e limiti pratici

| Verifica | Risultato |
|---|---|
| Baseline prima dell'ampliamento | 325 blocchi, 10.022 pass, 4 skip; zero failure/error/warning |
| Suite completa sorgenti | 364 blocchi, 11.223 pass, 4 skip; zero failure/error/warning; exit 0; 458,47 s |
| Suite completa namespace installato | Stessi conteggi; exit 0; 321,38 s; percorso del namespace verificato nella libreria isolata |
| Stato prima/dopo le suite | Seed/RNG, directory, dispositivi, ambiente e opzioni preesistenti preservati; il caricamento iniziale di namespace installati aggiunge soltanto default dei pacchetti |
| Build e installazione del tarball | Passati; tarball SHA256 `bb2740d3b046935a9f28af9b098ad16889eacea23b12aa92a56ef56df03c5094` |
| `R CMD check --no-manual --no-build-vignettes --no-tests` | `Status: OK`; esempi passati; test eseguiti integralmente nelle due esecuzioni separate sopra |
| Runtime salvato in processo R nuovo | Passato, senza caricare MIRA; nessun rendering Quarto implicato |
| Profilo reale HTML/DOCX/PDF installato | Tentato con Quarto 1.9.38 e TinyTeX disponibili; bloccato da `processx__create_pipe`, errore Windows 5 «Accesso negato» |

I quattro skip ordinari sono il sampling CmdStan e i tre profili Quarto
espliciti. Nessuno dipende dall'assenza di mmrm, nparLD o pbkrtest. Il profilo
Quarto separato riporta 93 pass, sei failure, tre errori e quattro warning
(12 blocchi), conseguenze del mancato avvio del processo e dell'assenza dei
documenti attesi. Non sono nascosti né assimilati a rendering riusciti.
Non è stata eseguita una nuova CI remota e non si dichiara validazione su
altri sistemi operativi, R 4.2 o tutte le versioni future delle dipendenze.

Ambiente eseguito: Windows 10, R 4.4.3, testthat 3.3.2 edition 3,
pbkrtest 0.5.4, mmrm 0.3.17, TMB 1.9.21, nparLD 2.3.3, multcomp 1.4-32,
TH.data 1.1-5. `DESCRIPTION` dichiara direttamente pbkrtest, mmrm >= 0.3.17 e
nparLD >= 2.3.3. mmrm/TMB usano binari CRAN ufficiali R 4.4: il tentativo
sorgente del TMB più recente incontrava il vincolo Windows sulle pipe.
Non viene rivendicata esecuzione locale di mmrm 0.3.19.

Le dipendenze nuove sono state prima installate e provate in una libreria
isolata e poi rese disponibili nella libreria utente R 4.4 standard con
permesso di scrittura. Lo script conserva pacchetti esistenti sufficienti e
non sovrascrive versioni già presenti. Nessun percorso personale è stato
aggiunto ai file del pacchetto o al profilo R. MIRA stesso è stato installato
soltanto nella libreria isolata di validazione.

Codice in quattro nuovi helper `R/mira_info_long_{extensions,kenward_roger,
mmrm,rank}.R`, hook nei due file principali, sette nuovi file testthat e
aggiornamenti minimi ai contratti/schema precedenti. Documentazione roxygen,
Rd rigenerati, README, NEWS, contratto degli output e disclosure AI aggiornati.
Artefatti riproducibili e log nel workspace `work/methodology-expansion-2026-10-26/`:
`full-source-final*`, `full-installed-final*`, `preservation-final.log`,
`cold-runtime-final.log`, `rendering-installed-explicit*`, `build/`,
`packages/*provenance.json` e staging degli oracoli/mutazioni.

Fonti primarie consultate: [pbkrtest e KRmodcomp](https://hojsgaard.github.io/pbkrtest/reference/kr-modcomp.html),
[metodologia MMRM](https://openpharma.github.io/mmrm/v0.3.17/articles/methodological_introduction.html),
[criteri nativi mmrm](https://openpharma.github.io/mmrm/v0.3.17/reference/mmrm_tmb_methods.html),
[manuale nparLD 2.3.3](https://cran.r-project.org/web/packages/nparLD/refman/nparLD.html)
e la vignetta del pacchetto sorgente per l'ipotesi MCAR. La moderna API nparLD
2.3.3 non viene confusa con quella del lavoro originale del 2012.

**Nessun commit, staging, push, PR o pubblicazione effettuato.** La revisione
scientifica del manutentore resta necessaria; questa espansione non dimostra
completezza per ogni disegno, correttezza universale o accettazione editoriale.

---

## Audit precedente — esito e perimetro

Audit di `mira_info_long()`, `mira_report_freq_long()`, helper effettivamente
coinvolti, metodi di stampa e contratti di output. Lettura integrale distribuita
fra quattro agenti di entrambe le implementazioni e di tutti i file della suite.
Nessuna modifica a Stan, inferenza bayesiana, prior, QTE o API estranee.
Nessun commit, staging, push, PR, release, installazione di dipendenze o modifica
globale dell'ambiente. Le modifiche sono locali e revisionabili.

Risultato ordinario: **325 blocchi `test_that`, 10.022 aspettative passate,
zero failure/error/warning testthat, quattro skip**, sia dai sorgenti sia dal
namespace installato in una libreria isolata. La baseline aveva 314 blocchi,
9.665 aspettative passate e gli stessi quattro skip, ma il runner terminava
con exit 1 per una perdita preesistente del seed globale. Il controllo RNG
finale passa. I numeri delle aspettative non sono numeri di test indipendenti.

Correzioni produttive dimostrate prima del fix:

1. Nei report frazionati per colonne, i confronti tra bracci a una singola
   visita potevano perdere `arm_a`/`arm_b` accanto a stime, intervalli e p-value.
2. Una RM-ANOVA interamente non stimabile poteva avere `performed = TRUE`,
   nessun motivo e righe non finite nella tabella di sensibilità. Il fit grezzo
   ora resta consultabile, ma non è dichiarato inferenza disponibile. Il motivo
   arriva al report anche quando non accompagnato da errori o warning.
3. Con braccio interamente mancante/vuoto, un `reference_arm` esplicito veniva
   ignorato prima della validazione, pur risultando specificato manualmente.
4. Impostazioni vuote di direzione/soglia producevano un errore interno di
   indicizzazione: ora l'errore identifica l'argomento pubblico non valido.

Rafforzati inoltre il verificatore PDF, l'isolamento RNG delle fixture e due
script di validazione che cercavano nomi generici dopo una dichiarazione
esplicita dell'outcome, ottenendo vettori vuoti invece delle metriche.

**Rendering pubblico Quarto ancora bloccato dall'ambiente**: HTML, DOCX e PDF
sono stati realmente tentati prima e dopo, ma nessun documento Quarto è stato
prodotto con successo. Non sono assimilati ai controlli ordinari con mock,
kable o knitr. Non si dichiara correttezza universale né accettazione JOSS.

## Stato iniziale e ambiente

- Repository: `C:/Users/david/OneDrive/Desktop/MIRA/MIRA`.
- HEAD iniziale: `f9a26b65954bf4a54f9a781b24ba12fb98df8cb3`.
- `git status --short` inizialmente vuoto: nessuna modifica o file non tracciato
  preesistente da incorporare. Cercati gli `AGENTS.md` nella catena degli
  antenati e nel repository; nessuno applicabile trovato.
- Windows 10 x64, build 19045; R 4.4.3 ucrt, processi `--vanilla`,
  `LC_ALL=C`, `LANG=C`, data della sessione riferita a Europe/Rome.
- Librerie delle dipendenze già presenti: libreria utente R 4.4 e libreria
  di R 4.4.3. Installato soltanto MIRA nella libreria isolata dell'audit.
- Versioni: testthat 3.3.2 (edizione 3), withr 3.0.2, pkgload 1.5.3,
  lme4 2.0-1, lmerTest 3.2-1, emmeans 2.0.3, afex 1.5-1,
  geepack 1.3.13, nlme 3.1-168, clubSandwich 0.7.0,
  performance 0.15.3, ggplot2 3.5.2, quarto R 1.5.1, knitr 1.50,
  pdftools 3.8.0; roxygen2 8.1.0 disponibile.
- Quarto CLI RStudio 1.9.38; TinyTeX/pdfTeX 1.40.29, TeX Live 2026;
  xelatex/lualatex presenti. La loro presenza non dimostra che possano essere
  avviati dal processo di rendering in questo ambiente.
- `covr` non installato: nessuna coverage misurata, percentuale o copertura
  dei rami rivendicata. Non sono state installate dipendenze per ottenerla.
- Le tre variabili `MIRA_RUN_CMDSTAN_TESTS`, `MIRA_RUN_QUARTO_TESTS`,
  `MIRA_RUN_QUARTO_PDF_TESTS` erano inizialmente non impostate. I runner
  ordinari le impostano a `false`; il profilo reale imposta le ultime due a
  `true`. Installazione/aggiornamento automatico TeX disabilitati in tutti i run.

Intermedi, log, risultati RDS/CSV dei reporter, build e copie di lavoro sono
fuori dal repository, in
`C:/Users/david/Documents/Codex/2026-10-26/files-mentioned-by-the-user-mira/work/audit-2026-10-26`.
Questa directory è indicata come `$audit` nei comandi sotto. I rapporti storici
del 12/16 ottobre sono stati letti; i loro risultati non sono risultati attuali.

## Inventario completo e ricostruzione dei contratti

Inventario iniziale: 42 file in `tests/testthat`, di cui 39 file di test e tre
helper; 151 blocchi info, 123 report, 40 altri. Nessun file separato di dati,
fixture o snapshot, setup o teardown. Le fixture sono funzioni/literal nei
file; l'eventuale directory `_snaps` vuota non contiene evidenze di snapshot.
Inventario finale: gli stessi file, più tre file di test; nessun test precedente
rimosso. Le classificazioni sotto riguardano il rischio realmente verificato,
non una certificazione universale del rispettivo modulo.

Tutti i nomi della tabella sono relativi a `tests/testthat/` e sono stati letti
integralmente dal gruppo di audit.

| File | Relazione e garanzia principale |
| --- | --- |
| `helper-mira-info-metrics.R` | Supporto diretto: vista canonica dei nomi per oracoli numerici; i test di naming restano separati. |
| `helper-mira-report-freq_long.R` | Supporto diretto: fixture manuali, parser HTML/DOCX, gate dei prerequisiti Quarto. |
| `helper-mira-targeted-cells.R` | Supporto diretto: chiavi dei contrasti e verificatori scientifici DOCX/PDF. |
| `test-mira-info_long-advanced.R` | Contratti dei moduli avanzati, formule, diagnostica e famiglie di p-value. |
| `test-mira-info_long-affine-properties.R` | Proprietà affini, riflessione e soggetto interamente mancante con configurazione fissata. |
| `test-mira-info_long-audit-2026-10-16.R` | Welch/Hedges/Holm manuali, campioni ID×visita esatti, quoting delle covariate. |
| `test-mira-info_long-boundaries-pathological.R` | Soglie, costanti, dati mancanti, ties e non stimabilità Friedman. |
| `test-mira-info_long-canonical-contract.R` | Nomi, ordine, riferimenti canonici, strutture native; isolamento RNG rafforzato. |
| `test-mira-info_long-configuration-contracts.R` | Detection, cronologia, etichette, provenienza, selettori e routing multi-outcome. |
| `test-mira-info_long-core.R` | Descrittivi, pairing, soglie di risposta, molteplicità, correlazioni, ICC e bracci. |
| `test-mira-info_long-detection.R` | Parser, override e validazione di ID/bracci/riferimento/covariate. |
| `test-mira-info_long-generated-properties.R` | Dodici casi generati limitati, aritmetica/BH indipendenti e pairing; seed locale. |
| `test-mira-info_long-model-oracles.R` | Oracoli del disegno e confronti diretti con motori per mixed/ANOVA/GEE/emmeans/nlme/CR2. |
| `test-mira-info_long-name-provenance.R` | Nome dichiarato/inferito, chiave, label, colonne sorgenti, prefissi e collisioni. |
| `test-mira-info_long-optional-dependencies.R` | Indisponibilità dei motori simulata e isolamento; non è un'installazione priva degli Imports. |
| `test-mira-info_long-properties.R` | Invarianze per righe, colonne, ID, label e configurazioni equivalenti. |
| `test-mira-info_long-public.R` | Argomenti pubblici, dispatch, visibilità e isolamento dei fallimenti fra outcome. |
| `test-mira-info_long-sampling-inference.R` | Campioni pairwise/complete-profile, intervalli e inferenza indisponibile. |
| `test-mira-info_long-schema-print.R` | Classi, tipi, ordine e semantica di print/summary/plot. |
| `test-mira-info_long-scientific-regressions.R` | Wald/rango/scala, termini del modello e ipotesi CR2 ridotte. |
| `test-mira-info_long-statistical-oracles.R` | Aritmetica, ranghi, conteggi, intervalli, ICC e bersagli Wald. |
| `test-mira-report-freq_long-bundle-contracts.R` | RDS/CSV/manifest completi, visualizzazione distinta da traversata e trasferimento. |
| `test-mira-report-freq_long-canonical.R` | Naming, ordine, selezione e strutture native nella filiera report. |
| `test-mira-report-freq_long-document-parser.R` | Selezione univoca di tabelle DOCX, blocchi e XML avversari. |
| `test-mira-report-freq_long-figure-oracles.R` | Dati/coordinate/gruppi/caption di tre figure prioritarie, con attesi manuali. |
| `test-mira-report-freq_long-files.R` | Filesystem, overwrite, sentinelle, export, payload e runtime. |
| `test-mira-report-freq_long-helpers.R` | Formattazione, escaping e helper di presentazione. |
| `test-mira-report-freq_long-html-parser.R` | Caption univoche, entità, righe irregolari, header e celle. |
| `test-mira-report-freq_long-integration.R` | Analisi reale→bundle e integrazioni HTML/DOCX Quarto opt-in. |
| `test-mira-report-freq_long-publication.R` | Etichette, diagnostica, suddivisione tabelle e runtime isolato. |
| `test-mira-report-freq_long-rendering-contracts.R` | Fault injection con renderer simulato e PDF reale opt-in. |
| `test-mira-report-freq_long-targeted-cells.R` | Corruzioni leggibili di chiavi e celle, inclusa la nuova cardinalità PDF. |
| `test-mira-report-freq_long-tree.R` | Sezioni eleggibili, vuoto/NULL/skipped, pruning e profondità. |
| `test-mira-report-freq_long-validation-regressions.R` | Output vecchi/nuovi, rendering parziale/vuoto e diagnostica. |
| `test-mira-report-freq_long-validation.R` | Input, metadata, percorsi, formati e forwarding simulato. |
| `test-mira-report-freq_long.R` | Omissione dei modelli, stampa GEE e conservazione analisi→report. |
| `test-package.R` | Indiretto: caricamento namespace, non oracolo scientifico. |
| `test-outcome-registry_long.R` | Estraneo al callgraph frequentista; preservato ed eseguito per regressioni globali. |
| `test-data-prior-long-oracles.R` | Fuori perimetro: mapping e prior bayesiani; letto e preservato. |
| `test-dynamic-covariates_long.R` | Fuori perimetro: design/prior/summary e mock; sampling CmdStan opt-in saltato. |
| `test-mira-fit-long-compilation.R` | Fuori perimetro: mock di compilatore/cache, nessuna compilazione reale rivendicata. |
| `test-qte-long.R` | Fuori perimetro: mock posteriori e oracoli QTE; nessun nuovo intervento. |

Nuovi file letti e verificati: `test-mira-info_long-empty-reference.R`,
`test-mira-info_long-estimability-contrasts.R`,
`test-mira-report-freq_long-scientific-keys.R`.

Letti inoltre: `tests/testthat.R`, `DESCRIPTION`, `NAMESPACE`, `README.md`,
`CONTRIBUTING.md`, roxygen nelle due implementazioni,
`man/mira_info_long.Rd`, `man/mira_detect_long.Rd`,
`man/mira_report_freq_long.Rd`, `man/print.mira_report_freq_long.Rd`,
`dev/mira-info-output-contract.md`, `dev/VALIDATION.md`, `dev/AI_USAGE.md`,
`dev/JOSS_READINESS.md`, `dev/RELEASE.md`, entrambi i rapporti datati,
`dev/test_mira_info_long.R`, `dev/test_mira_report_freq_long.R` e tutti i sei
`dev/validate-*.R`. Ispezionati i workflow `R-CMD-check.yaml`,
`cmdstan-smoke.yaml`, `citation.yaml`. Non si dichiara lettura integrale dei
restanti script dimostrativi bayesiani, non necessari a questo callgraph.

Flusso info ricostruito: validazione globale → parser/gruppi longitudinali →
ID/braccio/riferimento/covariate → adattamento per outcome → analisi singola
su dati finiti wide/long → moduli statistici → naming canonico e metodi di
presentazione. L'unità dell'input wide è normalmente il soggetto; `strict_id`
e la sua eccezione documentata non sono stati reinterpretati.

Flusso report ricostruito: dispatch `x`/`data` → vista transitoria dei nomi →
albero eleggibile e pruning → traversata completa prevista dei risultati →
CSV/manifest, payload RDS, runtime e QMD → rendering temporaneo/verifica/
promozione → descrittore. `sections` e limiti di visualizzazione non cambiano
arbitrariamente l'archivio; `max_depth` e reportabilità conservano il proprio
ruolo. Un modulo non eseguito può essere omesso dai CSV mentre il suo oggetto
grezzo resta completo nell'RDS: non è una promessa che ogni elemento interno
sia pubblicato in una tabella.

## Matrice delle garanzie e decisioni

La prima matrice è stata scritta fuori repository prima dei fix
(`initial-matrix.md`). Questa versione integra i risultati effettivi.
A = formule/dati teorici indipendenti; B = motore di riferimento chiamato
direttamente; C = coerenza interna o della filiera. C da solo non dimostra A.

| Modulo / rischio | Contratto e test esistenti | Oracolo / qualità iniziale | Dipendenze / profilo | Lacuna verificata, priorità e intervento | Verifica finale |
| --- | --- | --- | --- | --- | --- |
| Descrittivi, CI, pairing, ICC, correlazioni | `core`, `statistical-oracles`, `sampling-inference`, `generated-properties`: finite per visita, coppie per contrasto, profili completi solo dove richiesti; `to-from`, t CI e N pairwise | A/B, già discriminanti; identità campioni e valori asimmetrici presenti | Ordinario/base, modelli condizionali | Nessuna nuova lacuna prioritaria dimostrata; preservati senza duplicazione | Suite source/installed completa |
| Bracci e molteplicità | `audit-2026-10-16`, core/oracles: Welch/Hedges/Holm con gruppi e missingness asimmetrici, famiglie separate | A/B, forte | Ordinario/stats | Conservare segni, campioni e associazioni; nuova catena report con tre bracci | Welch 90% keyed→RDS/CSV/kable |
| emmeans marginali, semplici e interazioni | `advanced`, `model-oracles`, `canonical` eseguivano i rami; naming forte, attesi numerici dei simple effects insufficienti | B/C esistenti → nuovo A | emmeans, lm; ordinario | Alta: fixture satura 3×3 con celle diverse e N disuguali, chiavi complete prima dei numeri | Stima/SE/df/t/p grezzo/CI; due corruzioni rilevate |
| RM-ANOVA non stimabile | Fit e tabelle del motore disponibili; `performed`/sensibilità non distinguevano totale non stimabilità | B reale + A per caso parzialmente valido | afex; ordinario | Alta: gate minimo se nessun F/p/df valido; raw e warning conservati | Tre failure prima; 5 blocchi nuovi verdi inclusa filiera report |
| Report e chiavi dei confronti | `publication`, `rendering-contracts` forti per `from/to`; fallback primi due campi insufficiente per `time+arm_a+arm_b` | C → nuovo B/C | knitr; ordinario | Alta: ripetere visita e tutti i bracci se il limite lascia spazio a una statistica; fallback stretto documentato preservato | Otto failure prima; 187 aspettative dopo, limiti 7/4/3, HTML/LaTeX/pipe |
| Verificatore PDF | `targeted-cells` riconosceva valori corretti ma poteva ignorare un duplicato errato | Fixture manuale, negativo incompleto | Test ordinario su testo; PDF reale separato | Alta: contare prima righe strutturali/chiavi, poi controllare valori | Cinque nuove corruzioni prima accettate ora rifiutate; 48 pass nel file |
| Riferimento/argomenti vuoti | `public`, `detection`, `configuration`: override e gruppi osservati già forti | Contratto documentato, input negativi | Ordinario/base | Media: riferimento esplicito bypassato con tutti NA/blank; direzione/soglia vuote con errore interno | 35 pass nei tre workflow, automatico NULL conservato |
| RNG e isolamento | Generatori/QTE già protetti; due fixture canonical non protette | Stato prima/dopo, forte dopo fix | withr; ordinario | Media: preservare seed e RNGkind, anche senza seed o su errore; non mascherare nel runner | Sei assertion di stato fallivano; file finale 1508 pass, stato ripristinato |
| Bundle/runtime/filesystem/render mock | `bundle-contracts`, `files`, `validation*`, `tree`, `canonical`: RDS/CSV/manifest, collisioni, stale/symlink e sentinelle distinti | C, già forte sulle condizioni enumerate | Ordinario, tempdir; mock renderer | Preservati; nessun nuovo requisito transazionale, nessuna apertura reale di app | Tutti i report: 4236 pass/3 skip; runtime trasferito in processo pulito |
| Mixed, GEE, slope, nlme, CR2 | `model-oracles`, `advanced`, `scientific-regressions`, audit storico: formule, sample keys, termini e rango | A/B/C; forte su casi mirati, non validazione universale dei motori | Imports presenti; ordinario | Evitare equiparazioni fra ipotesi diverse; frontiere numeriche diagnosticate, senza soglie nuove arbitrarie | Suite completa; limite concreto descritto sotto |
| Script di verifica | `validate-science`, `validate-bundle` leggevano colonne generiche su outcome dichiarato | A previsto ma alcuni controlli vacui | MIRA installato, stats/knitr | Alta per evidenza: nomi `score_*`, forma/finitezza/chiavi prima dell'aritmetica; p appaiato indipendente | Entrambi falliti prima; bundle e 1000 repliche riusciti dopo |
| Quarto reale | `integration`, `rendering-contracts`: gating esplicito e verificatori semantici | A/C, opt-in e separati dai mock | CLI + R quarto; PDF anche TeX | Blocco ambientale Windows, non skip trasformato in successo | Baseline/finale entrambi 6 failure/3 error/4 warning; nessun formato pubblico verificato |

## Nuovi oracoli e contratti scientifici

### Contrasti emmeans

Modello gaussiano saturo `value ~ time_factor * arm_factor`, tre visite e tre
bracci con livelli intenzionalmente non lessicografici e conteggi per cella
3/4/5 ruotati. L'unità della fixture è l'osservazione indipendente del modello
lm, non un fit misto appaiato. Non si trasferisce questa varianza a un altro
modello: gli altri test mantengono la verifica del disegno mixed pubblico.

Medie di cella esplicite; SSE = 51, df residui = 27, sigma² = 17/9.
La varianza della media di cella è sigma²/N; le medie marginali usano pesi
uguali fra visite, distinti dai pesi osservati in un disegno sbilanciato.
I contrasti semplici sono differenze delle celle indicate dalle chiavi; le
interazioni sono differenze di differenze. Ipotesi nulla: contrasto zero,
due code, t con df 27, intervalli al 90%. Stime, SE, df, statistica, p grezzo
e limiti sono confrontati dopo cardinalità, unicità e corrispondenza completa
delle chiavi. Nessun atteso numerico viene calcolato da emmeans o helper MIRA.

L'oracolo nuovo non impone una nuova politica di molteplicità fra strati:
conserva i test esistenti delle famiglie e verifica qui i p grezzi. Il pooling
attuale delle tabelle simple effects è un contratto implementato da distinguere
dai default `by` del motore, non un nuovo assunto scientifico dell'audit.
Fonti primarie consultate: [contrasti emmeans](https://rvlenth.github.io/emmeans/reference/contrast.html),
[guida emmeans](https://rvlenth.github.io/emmeans/articles/AQuickStart.html).

Tolleranza 1e-10 per algebra lineare/stime/SE/CI/p a questa scala, 1e-12 per
df numerici esterni; differenze fra celle ben maggiori delle tolleranze.
Chiavi e numerosità discrete sono esatte. Non sono margini di equivalenza
inferenziale né tolleranze universali per tutti i modelli.

### ANOVA interamente e parzialmente non stimabile

Otto soggetti con valori `(3,5,4,7,6,8,5,9)`, invariati nelle tre visite:
variazione fra soggetti, nessuna variazione entro soggetto. Il vero afex
restituisce un fit con F/p non finiti e warning catturati. Prima del fix,
`performed=TRUE`, motivo assente, tre righe RM non finite nella sensibilità.
Ora serve almeno un effetto con F e p finiti e df positivi finiti.
Nessun p-value o modello è ricalcolato o sostituito; oggetto afex, tabelle
native/standardizzate e warning restano disponibili nell'RDS.

Controesempio adiacente con due bracci: SS_arm = 81/8, SS_error = 75/4,
F = 81/25 su df 1/6. Il valido effetto fra bracci mantiene `performed=TRUE`
anche con effetti entro soggetto indisponibili. Evitato un gate che rifiuti
qualunque warning o qualunque riga non finita.

La filiera pubblica controlla anche cambio osservato zero distinto dal p
appaiato indisponibile nei CSV; il fit RM non eseguito è omesso dalle tabelle
reportabili, conservato nell'RDS e spiegato nella diagnostica. Una variante
con lo stesso fit afex e solo `reason_skipped` isola il collector: due failure
prima della correzione, motivo visibile dopo. Un modulo solo disabilitato e
senza fit continua a essere omesso come documentato.

### Tre bracci, campioni e report frazionato

Dodici soggetti, gruppi control/zeta/alpha di dimensioni 3/4/5; due visite,
missingness differente e valori asimmetrici. N a baseline 3/4/5, alla finale
2/4/4. Le sei chiavi sono visita×arm_a×arm_b; baseline_balance ne ha tre.
Per ogni riga sono verificati campione, `arm_b-arm_a`, CI Welch al 90% e p
bilaterale con chiamata diretta `stats::t.test(b,a,var.equal=FALSE)`.
È un oracolo B del wrapper, distinto dalle formule Welch manuali già presenti.

Confrontati risultato analitico, payload identico e CSV selezionato attraverso
manifest univoco. Ogni blocco realmente passato a kable conserva chiavi e
quantità; l'HTML emesso è anche parsato indipendentemente per l'identità
scientifica delle righe. La selezione HTML/LaTeX/pipe è controllata nel test,
la serializzazione kable è reale; non è rendering completo Quarto/DOCX/PDF.
Il limite sette e il bordo quattro lasciano spazio alle tre chiavi; il limite
tre esercita il fallback documentato, senza imporre chiavi complete impossibili
in quel budget. Tutte le colonne rimangono distribuite nei blocchi.

Numeri raw: tolleranza 1e-12 per il confronto deterministico con stats.
I p pubblicati sono controllati nella loro forma a tre decimali; la verifica
non scambia un p arrotondato con il p raw. Gli altri verificatori documentali
mantengono il limite di arrotondamento 0.0005000001 e confronti esatti per N.

### Verificatori e isolamento

Il parser PDF è intenzionalmente limitato alla fixture t0/t1 con quattro
celle n/mean/sd/p. Ora rileva prima le righe numeriche plausibili, incluse
NA/Inf/NaN, poi pretende una riga per chiave e verifica i valori. Non pretende
di parsare qualsiasi PDF. Restano ammessi ordine inverso, pagine separate e
menzioni discorsive delle visite. I precedenti negativi restano intatti.

Le due fixture canonical conservano esattamente i seed originali ma chiamano
prima `withr::local_preserve_seed()`. Il nuovo blocco controlla seed presente
e assente, RNGkind alternativo, riproducibilità e fallimento per ID duplicato.
Un primo tentativo con il solo `local_seed()` non ripristinava lo stato in
questa installazione; quel tentativo fallito è conservato nei log e non è la
soluzione finale. Non modificati withr, il runner o le aspettative storiche.

## Sensibilità dimostrata

Nessuna mutazione lasciata nei sorgenti. Corruzioni tramite binding locali
in processi R separati; ripristino automatico, controllato esplicitamente per
la partizione del report. Le cifre sotto contano aspettative fallite, non bug
indipendenti. Rilevazioni dovute a chiavi/valori scientifici, senza errori
accidentali di parsing o caricamento nei tre esperimenti di mutazione.

| Esperimento | Evidenza e ragione della rilevazione |
| --- | --- |
| Scambio di tutti i valori numerici fra due righe simple_arm, chiavi intatte | Sei failure su stima/SE/limiti/t/p; 105 pass, zero error/warning/skip. Numeri ancora leggibili e plausibili, assegnati al confronto errato. |
| Duplicazione di una riga simple_time, dimensioni intatte | Dieci failure; la prima è unicità della chiave, poi chiave mancante e valori; 101 pass, zero error/warning/skip. |
| Omissione di arm_b dai blocchi | Dodici failure pertinenti all'identità dei confronti; nessun errore/warning/skip, binding originale ripristinato. |
| Riga PDF extra con media errata, SD errata, NA, Inf, o visita t2 | Cinque corruzioni prima accettate, tutte rifiutate dopo; artefatto testuale ancora leggibile. Non una prova di PDF reale. |
| Regressioni su implementazione originale | Otto failure chiavi report; tre failure RM; due failure collector; diciotto failure riferimento ignorato (con warning preesistenti); sei assertion RNG. Il test di impostazione vuota si interrompeva sul previsto errore interno non informativo prima del fix. |

I test numerici esistenti per pairing/denominatore/segno/BH e i loro runner
di mutazione sono stati letti e preservati. Non si attribuisce a questa sessione
la campagna di mutazione dei rapporti storici, né si inventa un mutation score.

## Baseline e verifiche finali

| Profilo effettivamente eseguito | Blocchi | Aspettative passate | Failure / error / warning test | Skip | Exit / esito |
| --- | ---: | ---: | --- | ---: | --- |
| Baseline completa source | 314 | 9665 | 0 / 0 / 0 | 4 | 1: seed globale non ripristinato |
| Finale completa source | 325 | 10022 | 0 / 0 / 0 | 4 | 0 |
| Finale namespace installato, libreria isolata | 325 | 10022 | 0 / 0 / 0 | 4 | 0 |
| Report adiacenti, prima dell'ultimo collector RM | 125 | 4236 | 0 / 0 / 0 | 3 | 0; collector coperto dai run finali |
| Nuovi contrasti/estimabilità/filiera RM | 5 | 111 | 0 / 0 / 0 | 0 | 0 |
| Nuove chiavi report | 1 | 187 | 0 / 0 / 0 | 0 | 0 |
| Riferimento e argomenti vuoti | 3 | 35 | 0 / 0 / 0 | 0 | 0 |
| Canonical e RNG | 8 | 1508 | 0 / 0 / 0 | 0 | 0 |
| Verificatori targeted-cells | 4 | 48 | 0 / 0 / 0 | 0 | 0 |
| Quarto reale baseline source, HTML/DOCX/PDF attivi | 12 | 93 | 6 / 3 / 4 | 0 | 1: blocco processo Windows |
| Quarto reale finale installed, stessi formati | 12 | 93 | 6 / 3 / 4 | 0 | 1: stesso blocco |

Incremento ordinario: 11 blocchi e 357 aspettative, associati alle garanzie
sopra. I confronti prima/dopo ordinari usano lo stesso runner e profilo.
Tempi elapsed: baseline 349.31 s; finale source 433.98 s; installed 312.14 s;
report adiacenti 182.30 s. Processi concorrenti: non sono benchmark.
Rendering reale: baseline source 25.06 s, finale installed 6.17 s; entrambi
si fermano prima di produrre i documenti.

I quattro skip sono: uno smoke CmdStan fuori perimetro; HTML reale;
DOCX reale; PDF reale. Il motivo identifica la rispettiva variabile opt-in.
Abilitando i profili reali si ottengono failure/error espliciti, zero skip:
nessun falso successo per mancanza di prerequisiti.

Il runner ordinario controlla seed, RNGkind, working directory, device,
opzioni preesistenti e TZ/LANG/LC_ALL/CMDSTAN/QUARTO_PATH. Tutti ripristinati
nei run finali. Source: nessuna nuova opzione. Installed: 18 opzioni nuove,
nessuna preesistente cambiata; un processo fresco che carica soltanto gli
Imports riproduce 18 nuove opzioni (afex e altri default dei pacchetti),
senza modificare opzioni esistenti. Non è un monitor di ogni risorsa OS.
I test nuovi non aprono connessioni persistenti o applicazioni, né richiedono
rete; cleanup e binding sono locali anche su errore.

Build e installazione completate con exit 0. Archivio
`MIRA_0.0.1.tar.gz`, SHA256:
`8E3D84533A7B5C8591CD82E71C82B901571507146712148C984C8892273E1812`.
Il runner installed stampa il path `$audit/library/MIRA`, i suoi export e
sessionInfo; non carica una vecchia installazione né sorgenti per aggirare il
namespace. Man/roxygen aggiornati in modo coerente; nessuna firma o export
cambiato, nessuna rigenerazione non pertinente di NAMESPACE.

Verificati anche i testi degli 11 file modificati inclusi nel tarball (R,
manuali, NEWS e test): corrispondono al working tree finale. Roxygen2 8.1.0
è stato eseguito su una copia isolata: il manuale info coincide esattamente;
per il report restano soltanto due spostamenti di delimitatori
`description/details` preesistenti. Ripetendo la generazione sul sorgente e
manuale originali di HEAD si riproducono gli stessi due spostamenti. I nuovi
paragrafi coincidono con quelli generati; conservata la collocazione storica
per evitare modifiche generate estranee al comportamento corretto. Il check
Rd e dei mismatch codice/documentazione passa. Log `docs-verification.log`
(confronto integrale inizialmente negativo) e `docs-baseline.log` conservati.

`R CMD check --no-manual --no-vignettes` completato con exit 0 e
**Status: OK**, zero ERROR/WARNING/NOTE formali. Include esempi e suite del
pacchetto installato dal tarball: 10.022 pass, zero failure/warning test,
quattro skip; processo dei test elapsed 262.64 s. Il console log conserva
cinque avvisi ambientali: quattro per indici CRAN/Bioconductor irraggiungibili
e uno per `du` (status 322), distinti dai warning formali e da quelli testthat;
nessuna dipendenza scaricata.

### Script di sviluppo, portabilità e integrazione reale

Le copie originali di `validate-science.R` e `validate-bundle.R`, eseguite con
caricamento source esplicito, falliscono rispettivamente con mismatch delle
righe nel data.frame e `is.numeric(csv$paired_t_p)` falso. Prima di fermarsi,
emettono warning `max(numeric(0))`: alcuni precedenti confronti erano vacui.
I nuovi guard controllano forma/finitezza/chiavi prima dell'aritmetica e usano
i nomi `score_*` giustificati dalla dichiarazione pubblica `outcomes="score"`.

Entrambi gli script corretti passano sul pacchetto installato:

- Bundle: risultato indipendente finale 33/8, varianza delle differenze 33/8,
  SE sqrt(33)/8, CI t su df 7 e p 2*pt(-sqrt(33),7); 26 CSV, manifest relativo,
  RDS, trasferimento in directory con spazi e runtime isolato; 15 figure.
- Studio già previsto dal progetto: n=24, 500 repliche per delta=0 e delta=2,
  seed 2026101201, cambiamenti normali indipendenti senza missingness. Oltre
  agli aggregati, ogni replica verifica stima/CI/p appaiato con formule t.
  Bias 0.01239/-0.01948; coverage 0.944/0.942; reiezioni 0.056/0.998;
  nessun segnale oltre i quattro MCSE definiti nello script. Durata 65.81 s.
  Non è una nuova simulazione generale né convalida di altri motori o di
  missingness informativa; i conteggi di replica non sono blocchi testthat.

Verifica aggiuntiva: il bundle trasferito è stato nuovamente knitted in un
**processo R pulito, senza caricare MIRA nel namespace**, con ambiente di
valutazione avente parent `baseenv()`. Risultati analitici (tutti i componenti
eccetto plots) identici al payload, dati/classi/nome delle figure conservati,
quantità scientifiche riconoscibili nel Markdown. Due verifiche preliminari
del driver sono fallite perché confrontavano l'identità di ambienti ggplot
ricreati da letture RDS separate, poi perché `all.equal` non gestiva un loro
ambiente; il controllo finale distingue quelle identità private dai dati
scientifici, senza normalizzare risultati analitici. Log intermedi conservati.
Questa prova stabilisce portabilità del runtime, non rendering Quarto.

Il profilo reale usa Quarto 1.9.38 con i tre formati attivati e TeX esistente.
La chiamata R fallisce in `processx_exec`, creazione della write pipe,
Windows error 5/access denied. Un tentativo diretto della CLI su QMD trasferito
fallisce con `Failed to spawn cmd.exe: Invalid handle`; indicando R_HOME
esplicitamente fallisce sullo spawn di Rscript con lo stesso problema di handle.
Nessuna installazione o modifica ACL tentata. Gli errori successivi sui file
assenti sono conseguenze di questo blocco. I successi storici ottenuti con
Pandoc alternativo non sono rivendicati come esecuzioni di questa sessione.

## Riproduzione

Comandi principali PowerShell. Usare una nuova directory esterna, con le
dipendenze già installate, e controllare ogni exit code; non usare un'installazione
MIRA preesistente come prova del working tree corrente.

```powershell
$repo = 'C:\Users\david\OneDrive\Desktop\MIRA\MIRA'
$audit = 'C:\path\outside-repository\new-long-audit'
$rscript = 'C:\Program Files\R\R-4.4.3\bin\Rscript.exe'
$r = 'C:\Program Files\R\R-4.4.3\bin\R.exe'
$env:LC_ALL='C'; $env:LANG='C'
New-Item -ItemType Directory -Path "$audit/build","$audit/library","$audit/bundle","$audit/science","$audit/rendering" | Out-Null
& $rscript --vanilla "$repo/dev/validate-tests.R" $repo "$audit/source" source
Push-Location "$audit/build"
& $r CMD build --no-build-vignettes --no-manual $repo
& $r CMD INSTALL "--library=$audit/library" MIRA_0.0.1.tar.gz
$env:R_LIBS="$audit/library"
& $rscript --vanilla "$repo/dev/validate-tests.R" $repo "$audit/installed" installed
$env:MIRA_RUN_CMDSTAN_TESTS='false'
$env:MIRA_RUN_QUARTO_TESTS='false'; $env:MIRA_RUN_QUARTO_PDF_TESTS='false'
$env:QUARTO_TEXLIVE_AUTO_INSTALL='false'; $env:QUARTO_TEXLIVE_AUTO_UPDATE='false'
$env:_R_CHECK_FORCE_SUGGESTS_='true'
$env:_R_CHECK_CRAN_INCOMING_REMOTE_='false'
& $r CMD check --no-manual --no-vignettes MIRA_0.0.1.tar.gz
Pop-Location
& $rscript --vanilla "$repo/dev/validate-bundle.R" "$audit/bundle"
& $rscript --vanilla "$repo/dev/validate-science.R" "$audit/science"
$env:QUARTO_PATH='C:\Program Files\RStudio\resources\app\bin\quarto\bin\quarto.exe'
$env:PATH='C:\Users\david\AppData\Roaming\TinyTeX\bin\windows;' + $env:PATH
$env:MIRA_RUN_QUARTO_PDF_TESTS='true'
& $rscript --vanilla "$repo/dev/validate-rendering.R" $repo "$audit/rendering"
```

`dev/validate-tests.R` è stato letto prima dell'uso e disattiva esplicitamente
i profili esterni. L'ultimo comando è separato e li attiva; per completare la
prova dei tre formati serve un ambiente in cui processx/Quarto possano avviare
i processi figli. Non risolvere il blocco installando automaticamente tool.

Comando mirato riproducibile per i nuovi file e verificatori:

```r
testthat::test_local(repo,
  filter = "empty-reference|estimability-contrasts|scientific-keys|canonical-contract|targeted-cells",
  stop_on_failure = TRUE)
```

Le prove prima/dopo sono nei log `reference-before/after`,
`empty-settings-before`, `configuration-final`, `estimability-before/after`,
`estimability-before-collector`, `report-keys-final-pre-fix/post-fix`,
`pdf-verifier-before/after`, `rng-before-error/rng-final`.
I runner locali di sensibilità e portabilità, con comandi esatti e RDS dei
reporter, sono conservati in `$audit`; nessun log grezzo grande nel repository.

Driver effettivamente utilizzati nella directory di questa sessione (qui
`$audit` è il percorso reale delle evidenze indicato sopra, non una directory
vuota; per una ripetizione conservare gli originali e cambiare le destinazioni
dei log/RDS nei driver):

```powershell
& $rscript --vanilla "$audit/run-rendering.R" $repo "$audit/baseline-rendering" source
& $rscript --vanilla "$audit/run-rendering.R" $repo "$audit/final-rendering" installed
& $rscript --vanilla "$audit/run-report-keys.R" "$audit/report-keys-post-fix"
& $rscript --vanilla "$audit/run-report-key-mutation.R" "$audit/report-key-mutation"
& $rscript --vanilla "$audit/run-statistical-mutation.R" 'swap-simple-arm-values'
& $rscript --vanilla "$audit/run-statistical-mutation.R" 'duplicate-simple-time-key'
& $rscript --vanilla "$audit/cold-runtime.R" "$audit/bundle/relocated bundle"
```

## CI, limiti residui e problemi separati

Il check ordinario CI raccoglie tutti i test, inclusi i tre nuovi file. Il job
Quarto separato filtra `integration|rendering-contracts` e abilita HTML/DOCX;
PDF è un profilo esplicito distinto. I nuovi test kable/verificatori restano
ordinari. Nessuna modifica CI necessaria; nessuna esecuzione CI remota osservata
o dispatch effettuato. Non verificati localmente macOS/Linux, R minimo 4.2,
altre versioni di testthat o altre combinazioni di dipendenze.

Limiti scientifici concreti: sugli stessi profili esattamente statici, alcuni
motori gaussiani possono restituire inferenze finite instabili alla frontiera
di varianza residua nulla. Il mixed principale segnala non convergenza e warning
di ottimizzazione ma un LRT può avere p circa 1.7e-32; GEE exchangeable può
restituire p molto piccolo con warning NaN catturato. Sono osservazioni su
questa toolchain, non costanti da cristallizzare nei test. Il fix RM non rende
affidabile ogni inferenza di frontiera degli altri motori. Non introdotta una
soglia residua arbitraria comune o un veto generale sui warning: cambierebbe
politica scientifica e bersagli diversi. Le diagnostiche restano visibili e
richiedono interpretazione; `performed` non certifica adeguatezza scientifica.

Restano limiti espliciti: oracoli dei wrapper non certificano ottimizzatori,
CR2/GEE/nlme non sono convalidati su ogni disegno, le figure controllate non
sono tutte le figure possibili, il verificatore PDF non è generale e nessun
PDF Quarto corrente è stato verificato. Non si inferisce correttezza dai soli
conteggi o da coverage non misurata. Nessuna nuova garanzia Bayes/Stan/QTE.
Le prove di errore del filesystem comprendono fault injection; non è stata
ottenuta una nuova prova sistematica di dinieghi ACL reali. Non si attribuisce
ai mock una verifica di tutti i permessi del filesystem Windows.

Discrepanze documentali minori preesistenti, separate e non usate come requisiti:
`dev/test_mira_info_long.R` descrive inspection-only nell'esempio 2 ma esegue
analisi complete; l'esempio roxygen corretto usa `inspect_only=TRUE`.
`mira_detect_long.Rd` eredita una descrizione di `outcomes=character(0)` che
nomina il contenitore info, mentre valore/codice restano detection. Non
introdotte decisioni scientifiche basate su queste descrizioni incoerenti.

La suite preserva tutte le garanzie precedenti; i nuovi controlli coprono
lacune dimostrate. Per la filiera pubblica Quarto rimane necessaria una prova
in ambiente funzionante: non dichiarare le due funzioni interamente chiuse.

## File modificati e revisione locale

- `R/mira_info_long.R`: estimabilità RM e validazione dei casi vuoti.
- `R/mira_report_freq_long.R`: chiavi dei blocchi e diagnostica RM.
- `man/mira_info_long.Rd`, `man/mira_report_freq_long.Rd`: contratti sincronizzati.
- Tre nuovi file test citati nell'inventario.
- `tests/testthat/helper-mira-targeted-cells.R`,
  `test-mira-report-freq_long-targeted-cells.R`: verificatore PDF e regressioni.
- `tests/testthat/test-mira-info_long-canonical-contract.R`: cleanup RNG e suo test.
- `dev/validate-bundle.R`, `dev/validate-science.R`: verifiche numeriche non vacue.
- `NEWS.md`, `dev/AI_USAGE.md`, questo rapporto: comportamento ed evidenza attuale.

Revisione finale: `git diff --check` superato; diff produttivo limitato ai
comportamenti descritti, tutti i test originali conservati. Nessun file
mutato, build, log o documento renderizzato aggiunto al repository. Le copie
di staging e gli esperimenti restano soltanto nella directory esterna
dell'audit. Lo stato Git finale contiene esclusivamente i 15 file locali
indicati sopra, inclusi tre nuovi test e questo rapporto; niente staged.
