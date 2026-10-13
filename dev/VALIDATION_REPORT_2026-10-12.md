# MIRA: validazione tecnica e scientifica del 12 ottobre 2026

## Esito e perimetro

Lavoro applicato direttamente al working tree di
`C:/Users/david/OneDrive/Desktop/MIRA/MIRA`, seguendo `test info.txt`.
Sono state rafforzate le verifiche scientifiche e di integrità dei report,
corretti bug riprodotti e riparati test che davano falsi positivi.
API pubblica, estimandi, metodi e famiglie di molteplicità sono preservati.
Le conclusioni valgono per le proprietà e gli scenari descritti, non costituiscono
una prova generale di correttezza né di ammissibilità a JOSS.

HEAD iniziale e finale: `74aa85293fad7c7e61d980a2e992d857ebc2fdfc`.
Il repository era inizialmente pulito, senza modifiche staged o untracked.
Nessun commit, staging, push, PR, tag, release, snapshot Git, reset, clean,
rebase, dispatch CI o submission è stato eseguito. Le modifiche sono lasciate
locali per la revisione del mantenitore. L'inventario finale e il diff sono
nell'archivio delle evidenze.

## Ambiente e riproduzione

Sistema locale Windows, PowerShell 7.6.5, R 4.4.3; processi principali con
`--vanilla`, `LC_ALL=C`, `LANG=C`. R 4.5.3 era presente ma mancavano Imports
obbligatori: non è stato presentato come profilo installato verificato.
Versioni principali: testthat 3.3.2, withr 3.0.2, pkgload 1.5.3, afex 1.5-1,
clubSandwich 0.7.0, cmdstanr 0.9.0, emmeans 2.0.3, geepack 1.3.13,
ggplot2 3.5.2, knitr 1.50, lme4 2.0-1, lmerTest 3.2-1, nlme 3.1-168,
performance 0.15.3, quarto R 1.5.1. Gli elenchi completi delle librerie e
sessionInfo sono conservati con ogni profilo.

Gli engine frequentisti, cmdstanr, knitr e quarto R sono Imports obbligatori;
le Suggests sono testthat, withr, yaml, loo e posterior. Quarto CLI, LaTeX e
CmdStan sono software esterno. Non sono stati rimossi Imports per simulare
assenze opzionali. Non sono state aggiunte dipendenze R, effettuate
installazioni automatiche o richieste risorse di rete dai test ordinari.
L'installazione del package MIRA costruito è stata eseguita esplicitamente in
una libreria temporanea; le dipendenze erano già disponibili.

La guida `dev/VALIDATION.md` contiene comandi PowerShell riproducibili per
source, installed, build/check, bundle, rendering, mutazioni e simulazione.
Gli script `dev/validate-*.R` salvano risultati fuori dal repository.
Il runner offline disabilita espressamente i tre flag di integrazione.
Le integrazioni abilitate falliscono quando mancano gli strumenti richiesti:
non sono trasformate in skip per ottenere un risultato verde.

## Risultati eseguiti

| Profilo | Blocchi | Aspettative passate | Failure / error / warning | Skip | Esito |
| --- | ---: | ---: | --- | ---: | --- |
| Baseline source, prima delle modifiche | 238 | 3757 | 0 / 0 / 0 | 3 | Suite passata |
| Source finale | 270 | 4556 | 0 / 0 / 0 | 4 | Passata |
| Namespace installato finale | 270 | 4556 | 0 / 0 / 0 | 4 | Passata |
| Test durante R CMD check finale | — | 4556 | 0 / 0 / 0 | 4 | Status: OK |

Un blocco è un `test_that`, distinto dalle singole aspettative.
Tempi wall-clock del runner finale: source 209.50 secondi, installed
200.67 secondi; includono caricamento e lavoro concorrente e non sono benchmark.
Il dettaglio per blocco è nei CSV. Gli skip offline finali sono CmdStan,
Quarto HTML, DOCX e PDF: sono esclusioni esplicite del profilo, non verifiche
riuscite. La baseline aveva tre skip perché il test DOCX esplicito è nuovo.

Build da sorgente con `R CMD build --no-build-vignettes --no-manual`,
installazione dell'archivio con `R CMD INSTALL --library=...` e
`R CMD check --no-manual --no-vignettes` sono stati realmente eseguiti.
Il package installato è stato caricato dalla libreria temporanea, con
registrazione del percorso namespace e delle API esportate; i test installed
usano `load_package="installed"`, non pkgload. Il check include esempi e
test, ma non il manuale PDF, il profilo as-cran/incoming remoto o una matrice di sistemi.
Il check finale termina con Status: OK: zero ERROR, WARNING o NOTE formali. Il console log registra quattro warning per indici CRAN/Bioconductor non raggiungibili durante il controllo delle dipendenze, e un warning del comando esterno du per la dimensione del package (status 322). Entrambi i controlli terminano comunque OK; tutte le dipendenze richieste erano già installate. Nessun package è stato scaricato o installato dal check. Queste cinque segnalazioni ambientali sono distinte dai warning dei test, che sono zero; la disponibilità degli indici remoti non è stata verificata.

RNG, RNGkind, directory, device grafici, variabili d'ambiente registrate e
opzioni già esistenti risultano ripristinati. Il primo caricamento degli
Imports nel processo installed aggiunge 18 opzioni predefinite: un processo
indipendente che carica solo gli Imports, senza test, riproduce esattamente
le stesse aggiunte. Non sono attribuite ai test. `show_praise=FALSE` evita
che il messaggio casuale del reporter testthat consumi RNG.

Trasparenza sui runner: il primo export CSV della baseline fallì dopo la
suite, tentando di scrivere una colonna list; il risultato RDS era già salvato
e il CSV è stato ricostruito. Un primo script di mutazione modificato mentre
R lo leggeva terminò con errore dopo le prove; il profilo canonico è stato
rieseguito da un file stabile in `mutations-final`, con riepilogo completo.
Questi errori del coordinamento non sono considerati difetti del package o
mutazioni rilevate. I primi processi agenti registrarono warning di startup
per locale Windows C.UTF-8: le esecuzioni finali usano C e ne sono prive.

## Bug riprodotti e correzioni minime

1. **Covariate sovrascritte internamente.** `patient_mean` e `time_ordinal`
   potevano essere accettate come covariate e poi sostituite dalla media
   dell'outcome o dall'indice di visita. Sono ora incluse nei nomi riservati:
   selezione manuale rifiutata, auto-detection esclusa. Aggiornata anche la
   documentazione ereditata da `mira_detect_long()`.
2. **Rango Wald dipendente dall'unità.** La soglia degli autovalori aveva un
   pavimento assoluto 1; outcome riscalati a 1e-6 diventavano falsamente non
   stimabili. La soglia è relativa alla massima grandezza degli autovalori.
   Oracoli esatti Q=2/rango2 e Q=1/rango1, matrice nulla non stimabile e GEE
   reale con independence/exchangeable/AR1 verificano il comportamento.
3. **Ipotesi GEE/CR2 contaminate da nomi delle covariate.** Substring come
   `age_time_factor`, `pre_arm_factor` o nomi di livelli potevano aggiungere
   coefficienti estranei. La selezione usa termini del modello e attributo
   `assign` della matrice del disegno. Entrambi i chiamanti produttivi passano
   il modello. TIME e ARM mantengono le rispettive interazioni già previste.
4. **Effetti ANOVA selezionati erroneamente.** `lifetime`/`farm` potevano
   essere letti come time/arm. Il matcher ora usa etichette canoniche esatte.
   La revisione finale ha inoltre riprodotto con afex reale la collisione
   case-sensitive `Time_factor`/`Arm_factor`: R distingue questi nomi dai
   fattori interni minuscoli. Rimossa la conversione indiscriminata in minuscolo.
5. **Vecchio report dichiarato nuovo successo.** Un renderer che non produce
   output poteva far riutilizzare il vecchio documento; un renderer parziale
   poteva distruggerlo. Ogni formato usa ora un target temporaneo fresco,
   verificato e promosso solo dopo successo. Fault injection copre no-op,
   scrittura parziale, percorso errato, conflitti e fallimento di promozione.
6. **Etichette interpretate come markup.** Heading e metadati inline trattano
   letteralmente caratteri Markdown/HTML, conservando Unicode. Il helper è
   incluso nel runtime esportato. Sei casi semantici indipendenti CommonMark
   completano le asserzioni sul testo; non sono presentati come test Pandoc.
7. **Diagnostica persa per moduli avanzati falliti.** Errori e warning di
   moduli non eseguiti appaiono ora nei diagnostici, con motivo registrato,
   senza ripubblicare stime stale. I moduli disabilitati con solo motivo
   continuano a essere omessi dalle sezioni risultati.
8. **Heading missingness ambiguo.** Il nome condiviso per la missingness
   qualificava erroneamente come “by arm” anche grafici senza arm. La label
   generica preserva il contesto corretto e le caption specifiche.
9. **Tabelle DOCX perse nella conversione.** knitr restituisce il pipe-table
   come vettore di righe; il cat produttivo lo concatenava senza separatori.
   Pandoc produceva testo, senza celle XML di tabella. Ora le righe sono unite
   con newline prima dell'emissione. La regressione usa knitr reale e verifica
   header, separatore e dati; HTML e LaTeX scalari mantengono lo stesso testo.

Le nuove regressioni scientifiche iniziali mostrarono 26 failure sostanziali
prima dei fix; un ulteriore errore era lo scaffolding del nuovo argomento
interno non ancora implementato e non conta come prova scientifica. Il caso
ANOVA maiuscole/minuscole mostrò nove failure reali, senza errori o warning.
Il primo gruppo di report mostrò dieci failure di contratto prima dei fix.
Il titolo missingness e la struttura pipe-table hanno mostrato tre failure
pertinenti ciascuno prima delle rispettive correzioni.
I log before/after e i nuovi test permettono di verificare queste affermazioni.

## Affidabilità della suite e indipendenza degli oracoli

Ventuno verifiche QTE di argomenti invalidi passavano chiamando un simbolo
`summarize` inesistente e usando `expect_error` generico. Ora invocano il
wrapper corretto e richiedono messaggi specifici dell'argomento. Il cleanup
QTE ripristina RNGkind e seed, anche quando il seed iniziale non esiste.
Due soppressioni generiche di warning sono sostituite da verifiche della causa.
Non sono state alzate tolleranze o ridotte verifiche numeriche.

Gli oracoli nuovi includono osservazioni Stan/ID/visita/arm elencate a mano,
tempi irregolari e riferimento, formule indipendenti di conversione dei prior
e mappatura di 16 campi con valori distinti; prior metadata e campi invalidi
sono rifiutati. Friedman con ties usa ranghi manuali: Q=9/7, W=9/56,
p=exp(-9/14). Per GEE independence il sandwich è formato da design OLS,
residui e somme per cluster, indipendentemente dal wrapper produttivo.
I confronti afex/CR2 dopo rinomina testano identità delle ipotesi e wrapper,
non certificano gli algoritmi degli engine.

Le trasformazioni affine verificano stime, SD, SE, CI, cambiamenti e contrasti
nelle proprie unità; inversione del segno con coerente direzione di
miglioramento conserva i responder. Un soggetto interamente missing cambia
i denominatori di disponibilità senza aggiungere coppie finite. I test
esistenti per cronologia custom, livelli inutilizzati, riferimenti, outcome
aggiuntivi e ordinamento sono conservati, evitando duplicazioni.

Il generatore rapido usa 12 fixture, seed 20261012, N=6..12 e 2..5 visite,
almeno tre profili completi e NA/NaN/Inf/-Inf nelle altre osservazioni.
Formule esplicite per t appaiato, CI e BH (ordinamento e cumulative-min)
sono indipendenti dagli helper produttivi. Verifica anche permutazioni di
righe/colonne e biiezioni degli ID con configurazione bloccata. Le diagnostiche
includono seed e caso; le regressioni deterministiche restano la protezione
dei bug trovati. Aritmetica moderata: tolleranze circa 1e-12; fit e trasformazioni
numeriche: 1e-10; display a tre decimali: mezzo dell'ultima cifra mostrata.

## Matrice contratto → rischio → test/oracolo → esito

| Contratto/proprietà | Rischio e copertura precedente | Lacuna corretta / nuova prova | Esito |
| --- | --- | --- | --- |
| Pairing ID/visita, N, segno, CI | Oracoli esistenti `statistical-oracles`, `sampling-inference` | Generatore con formule indipendenti e trasformazioni affine | Pass; mutanti pairing/N/segno rilevati |
| NA/NaN/Inf e disponibilità | `boundaries-pathological`, `configuration-contracts` | Coppie finite generate e soggetto interamente missing | Pass; mutante Inf rilevato |
| Covariate interne | Rifiuti generali già presenti | Nomi omessi patient_mean/time_ordinal | Riprodotto e corretto |
| GEE Wald alle diverse scale | Confronti engine già presenti | Q e rango analitici più GEE reale riscalato | Riprodotto e corretto |
| Ipotesi GEE/CR2 | `model-oracles` | Termini ingannevoli, livelli e rinomina; sandwich indipendente | Riprodotto e corretto |
| ANOVA sensibilità | Confronti di modello presenti | Canonicalità esatta e case sensitivity con afex reale | Riprodotto e corretto |
| Friedman con ties | Test globali generali presenti | Ranghi manuali Q/W/p | Pass indipendente |
| Molteplicità | BH già verificato in inferenza | Algoritmo BH indipendente su casi generati | Pass; omissione BH rilevata |
| Stan data e prior | Fixture core/shared presenti | Enumerazioni e formule indipendenti, metadata invalidi | Pass; pairing Stan mutato rilevato |
| QTE input e stato | 21 errori generici falsamente verdi | Funzione reale, causa specifica, cleanup completo | Pass affidabile; nessuna validazione sampling implicita |
| Payload e CSV completi | `bundle-contracts`, `files`, `integration` | Analisi vera, formule manuali, relocation e knitr isolato | Pass; mutanti truncation/section-loss rilevati |
| Nuovo rendering e overwrite | Mock di successo/fallimento presenti | Vecchio file+no-op, partial, wrong path, promotion | Riprodotto e corretto; Quarto pubblico bloccato |
| Label e diagnostica | `tree`, `publication`, helper LaTeX esistenti | Escaping semantico, advanced failure, missingness heading | Pass; conversione reale separata documentata |
| Stato/sessione | Cleanup locali parziali | Audit esterno source/installed e RNG QTE senza seed | Pass nei limiti registrati |
| Portabilità | Workflow release presente | Matrice oldrel/minima e job Quarto installato distinti | Configurato, non eseguito remotamente |

## Mutazioni mirate realmente eseguite

Copie isolate dell'intero working tree corrente, incluse modifiche locali,
senza `.git`, ciascuna in un processo R distinto. Nessuna modifica alla
produzione reale durante gli esperimenti. Rilevazione solo tramite failure
pertinenti, mai per parse/load/infrastruttura.

| Mutante | Failure pertinenti | Errori downstream | Esito |
| --- | ---: | ---: | --- |
| Segno del cambiamento invertito | 167 | 0 | Rilevato |
| Pairing soggetti errato | 269 | 0 | Rilevato |
| N appaiato errato | 246 | 0 | Rilevato |
| Aggiustamento BH omesso | 11 | 0 | Rilevato |
| Cronologia interna ignorata | 0 | 0 | Sopravvissuto |
| Inf conservati come osservazioni | 108 | 0 | Rilevato |
| CSV troncato come display | 32 | 0 | Rilevato |
| Soglia Wald assoluta ripristinata | 15 | 0 | Rilevato |
| Ipotesi GEE contaminate | 11 | 0 | Rilevato |
| Pairing osservazioni Stan errato | 3 | 0 | Rilevato |
| Falso successo renderer | 17 | 0 | Rilevato |
| Archivio limitato alle sezioni mostrate | 10 | 6 | Rilevato dalle dieci failure |

Tutti i mutanti hanno zero warning; alcuni filtri includono i consueti skip
d'integrazione. Gli errori downstream dell'ultimo mutante derivano da file
scientifici mancanti e non sono usati per classificarlo come rilevato.
Il totale 11/12 descrive sensibilità a questi errori concreti, non coverage
o probabilità di correttezza. Il sopravvissuto modifica il fallback di
ordinamento interno, ma nel percorso pubblico `.outcome_name` è già fornito
e `clear_order` è falso: il ramo che differisce non viene eseguito. È un
limite di questa mutazione; i test pubblici di cronologia custom restano attivi.
Questa campagna precede gli ultimi fix di label/case/pipe-table, che non cambiano i
bersagli mutati; non è dichiarata rieseguita su tali ultimi dettagli.

## Simulazione eseguita, distinta dalla suite rapida

Prima dell'esecuzione sono fissati nel codice: N=24 soggetti indipendenti,
baseline N(10,3²), cambiamento indipendente N(delta,2²), delta=0 oppure 2,
nessun missing, 500 repliche per scenario, seed 2026101201. Target:
E[follow-up − baseline]. Prima di aggregare, tutte le 1000 repliche
confrontano stima e CI con formule indipendenti a 1e-10.

| Target | Bias (MCSE) | Copertura 95% (MCSE) | Rifiuto bilaterale (MCSE) |
| ---: | --- | --- | --- |
| 0 | 0.01239 (0.01858) | 0.944 (0.01028) | 0.056 (0.01028), tipo I |
| 2 | −0.01948 (0.01857) | 0.942 (0.01045) | 0.998 (0.00200), potenza |

Nessun segnale supera i criteri di investigazione di quattro MCSE specificati
prima dell'esecuzione. Tempo registrato 88.72 secondi. Il metodo t è appropriato
al DGP normale appaiato scelto. Non si inferisce adeguatezza sotto missingness
informativa, altri DGP, tutti gli engine o modelli Bayesiani. I CSV contengono
tutte le repliche e i riepiloghi. Studio eseguito sul package installato
prima degli ultimi fix ANOVA/label/pipe-table, che non sono utilizzati in questo profilo.

## Catena report e integrazioni effettive

Il profilo `validate-bundle.R` usa otto soggetti e tre visite. Totali manuali
64,78,97 → medie 8,9.75,12.125; cambiamento finale-baseline: somma33,
somma dei quadrati165, media33/8, varianza33/8, SE=sqrt(33)/8,
t=sqrt(33), df7. Payload e CSV sono confrontati con questi valori.
Il bundle esporta 27 CSV con riferimenti relativi, viene spostato e il runtime
eseguito in ambiente isolato; knitr cattura 15 figure. Questo profilo passa.

Quarto 1.9.38 è presente nella distribuzione RStudio, fuori da PATH; TinyTeX
ha pdflatex/xelatex/lualatex. Con CLI esplicitamente configurata, i test
pubblici reali sono stati eseguiti: 103 aspettative, 93 passate, sei failure,
tre errori downstream e quattro warning, zero skip. Causa primaria:
`processx` non può creare una write pipe (OS error5, accesso negato), prima
dell'avvio del renderer. Gli errori di lettura derivano dai documenti assenti.
Anche Quarto CLI diretto è stato tentato, ma non avvia Rscript per invalid
handle. Non è stato aggirato il problema cambiando la produzione, installando
software o dichiarando verde l'integrazione. Questi tentativi usano lo stesso
codice report salvo la successiva correzione della label missingness.

È riuscita una catena reale alternativa sul package installato finale, senza
override di sorgenti: analisi → payload/runtime esportati → knitr → Pandoc
3.8.3 → HTML/DOCX; per PDF, XeLaTeX 2026 in modalità no-pdf e xdvipdfmx
eseguiti separatamente. Non è un rendering Quarto e non verifica la promozione
dei documenti o i riferimenti generati da Quarto. Gli script e README
nell'archivio descrivono esattamente le opzioni, inclusi echo=FALSE per knitr
diretto e LC_CTYPE=English_United States.utf8 per conservare Unicode su questo
Windows. Nessuna installazione o download.

Risultati: HTML sei tabelle e otto figure embedded; DOCX 14 tabelle XML e
otto immagini PNG valide (prima della correzione: zero tabelle); PDF sette
pagine e otto figure vettoriali con caption. Per tutti e tre i documenti sono
controllate celle N, media, SD, cambiamenti e p-value con aritmetica indipendente
e formula chiusa t con cinque gradi di libertà. Fixture separata di sei
soggetti: medie 3.5/4.5/6.5, SD=sqrt(3.5); cambiamenti medi 1/3/2 e p
0.0755868184/0.0034781651/0.0027649603, mostrati .076/.003/.003. Anche
caratteri Markdown/HTML, entity letterale e μ risultano conservati.

Tutte le pagine PDF sono state renderizzate con Poppler e ispezionate: nessuna
sovrapposizione o troncamento osservata; una tabella lunga prosegue alla pagina
successiva con header e riga conservati. La tipografia resta compatta: sono
artefatti interni di verifica, non documenti editoriali rifiniti. Knitr emette
otto warning DOCX perché fig.align non è supportato: allineamento ignorato,
otto figure presenti. Il secondo passaggio XeLaTeX non registra warning,
overfull/underfull, caratteri mancanti o errori. Log, controlli numerici,
manifest e immagini delle pagine sono nell'archivio; questi warning di
formattazione sono distinti dai zero warning della suite offline.

La fault injection ordinaria è separata dagli esperimenti reali. Non promette
transazioni dell'intero bundle: un errore di serializzazione/esportazione può
lasciare file parziali. La promozione dei documenti usa file.copy, senza
garanzia di atomicità del filesystem. File estranei e documenti precedenti
sono protetti nei contratti testati; nessuna nuova garanzia generale di
overwrite è stata introdotta.

## File e documentazione

Produzione: `R/mira_info_long.R`, `R/mira_report_freq_long.R`.
Documentazione: `NEWS.md`, manuali di mira_info_long/mira_detect_long/
mira_report_freq_long, `dev/JOSS_READINESS.md`, `dev/AI_USAGE.md`, questa
relazione e `dev/VALIDATION.md`. CI: `.github/workflows/R-CMD-check.yaml`.
Test modificati: helper report, integration, publication, rendering-contracts,
QTE; test nuovi: data-prior-long-oracles, info affine-properties,
generated-properties, scientific-regressions e report validation-regressions.
I cinque script nuovi `dev/validate-{tests,rendering,bundle,science,mutations}.R`
rendono ripetibili i profili. L'inventario finale nell'archivio riporta tutti
i nomi, stato Git e SHA256; include i file nuovi che `git diff` ordinario omette.

## Rischi residui e JOSS

1. **Validazione scientifica esterna incompleta.** CmdStan non è stato
   eseguito: directory locale inaccessibile e configurazione non utilizzabile.
   Compilazione/sampling e recupero Bayesiano, modelli mal specificati,
   missingness informativa e giustificazione dei threshold clinici richiedono
   studi distinti. Gli smoke test corti predisposti non li sostituiscono.
2. **Rendering pubblico Quarto non riuscito su questo host.** Gli artefatti
   alternativi validano solo la catena esplicitamente descritta. Il profilo
   Quarto deve essere rieseguito in un ambiente in grado di avviare i processi.
3. **Non stimabilità e assunzioni.** Friedman interamente tied può conservare
   performed=TRUE con Q/p/W NaN senza reason_skipped: quel flag non garantisce
   inferenza stimabile. La ridondanza tra random intercept e compound symmetry
   in alcuni fit nlme richiede valutazione metodologica; non sono cambiati
   silenziosamente metodo o estimando. Il fallback interno di selezione per
   nomi resta substring, ma tutti i chiamanti produttivi passano il modello.
4. **Portabilità da verificare.** Release Linux/Windows/macOS, oldrel-1 e
   minima R4.2.3 con Matrix1.6-5 sono configurati, non eseguiti. La risoluzione
   delle dipendenze transitive sulla minima e i comportamenti degli altri OS
   rimangono da verificare. Nessun esito remoto viene inventato.
5. **Revisione umana e limiti del bundle.** Il mantenitore deve riesaminare
   modifiche, interpretazione scientifica e documentazione. Bundle non
   transazionali e limiti di formattazione sono dichiarati, non mascherati.

Consultate le fonti ufficiali JOSS aggiornate: [review criteria](https://joss.readthedocs.io/en/latest/review_criteria.html)
e [submission guidance](https://joss.readthedocs.io/en/latest/submitting.html).
Installazione, documentazione e verifica oggettiva sono pertinenti alla review;
JOSS non prescrive questi conteggi, un mutation score o gli studi avanzati
richiesti qui. Storia pubblica, sviluppo attivo oltre sei mesi, uso scientifico,
impatto e disclosure AI richiedono evidenza separata. Un timestamp locale di
commit non prova la data di pubblicazione. Le sezioni del paper e la disclosure
sono aggiornate nei promemoria, senza inventare adozione, contributi o revisione
umana. Quest'ultima resta pendente; nessuna comunicazione editoriale effettuata.

## Evidenze e controllo finale

L'archivio associato contiene guida, script, inventario/diff/hash finali,
baseline, log e CSV source/installed/check, regressioni before/after,
repliche della simulazione, riepiloghi e failure delle mutazioni, e prove
di report. Le grandi copie dei mutanti e le librerie R temporanee sono escluse.
L'archivio costruito del package verificato è incluso e identificato con hash.

Controllo finale: HEAD invariato, indice vuoto, git diff --check riuscito. Stato finale: 14 file tracked modificati e 12 file nuovi untracked, tutti intenzionali; nessun file staged. Il dettaglio completo è in git/status.txt e working-tree/ nell'archivio. Le modifiche alla documentazione di validazione sotto dev/ sono escluse dalla build come da .Rbuildignore; il tarball verificato include codice, manuali e test finali.

I file temporanei sono nel workspace di audit, fuori dal repository.
Nessuna evidenza di test simulato o saltato viene presentata come integrazione
riuscita. Restano disponibili i comandi per ripetere i controlli mancanti.
