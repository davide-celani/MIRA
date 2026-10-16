# MIRA: validazione longitudinale mirata — 16 ottobre 2026

## Perimetro e provenienza

Eseguita la richiesta in `test info.txt` su `mira_info_long()`,
`mira_report_freq_long()` e percorsi direttamente necessari. HEAD iniziale
`56f3777f652e4ba9cd87a9946963525683365bdd`; repository inizialmente pulito,
senza staging o untracked. I rapporti del 12 ottobre e le correzioni CI del
13 ottobre sono materiale preesistente: sono stati letti e non sovrascritti.
Le loro esecuzioni non sono attribuite a questa sessione.

Nessun commit, push, PR, tag, release, snapshot Git, dispatch CI, reset,
clean, rebase, comunicazione o submission. Nuove modifiche locali non staged.
Intermedi, libreria MIRA, build, copie mutate e log sono fuori dal repository,
in `work/validation-2026-10-16` e nelle tre directory di audit della sessione.
Non sono state installate dipendenze né modificati ACL/configurazioni globali.

## Esito sintetico

Due bug produttivi riprodotti e corretti con modifiche minime:

1. La rimozione dei backtick dai nomi delle covariate poteva confondere due
   colonne selezionate distinte (`age\`years` e `ageyears`), riducendo il
   disegno e alterando il modello senza errore. I backslash potevano produrre
   errori di parsing. Il quoting ora usa `encodeString` con delimitatore
   backtick, conservando
   i simboli letterali e il testo precedente per nomi ordinari. Correzione
   condivisa dal modello primario e dagli helper dei modelli avanzati.
2. Quarto mocked poteva creare un file di zero byte senza errore: il report
   lo promuoveva e dichiarava successo, distruggendo il documento precedente.
   Prima della copia ora si richiede una dimensione finita positiva. La
   regressione copre HTML, DOCX, PDF, diagnostica, mancata apertura del browser,
   conservazione dei vecchi documenti e cleanup. La dimensione positiva è una
   condizione necessaria, non una verifica completa della validità del formato.

Rafforzati gli oracoli Welch/Hedges/Holm e le identità ID×visita dei modelli;
verificati dati, coordinate e caption di tre figure prioritarie. Riparati i
verificatori HTML/DOCX contro recycling silenzioso e selezioni non univoche.
Conservate le protezioni valide esistenti senza riscrivere la suite.

## Ambiente e risultati eseguiti

Windows 10 x64 build 19045, R 4.4.3 ucrt; processi `--vanilla`,
`LANG=C`, `LC_ALL=C`. Tutti gli Imports richiesti sono già disponibili.
R 4.5.3 è presente ma manca di afex, clubSandwich, cmdstanr, geepack,
lmerTest e quarto: profilo **non verificato/bloccato** senza installazioni.
Quarto CLI RStudio 1.9.38 e TinyTeX/pdfTeX 2026 sono presenti fuori PATH.
Versioni, librerie e sessionInfo completi sono nei log e CSV allegati.

| Profilo | Blocchi | Aspettative passate | Failure/error/warning test | Skip | Stato |
| --- | ---: | ---: | --- | ---: | --- |
| Baseline source completa | 273 | 4575 | 0/0/0 | 4 | Verificato |
| Baseline pertinente, estratta dalla suite completa | 233 | 3448 | 0/0/0 | 3 | Verificato |
| Finale source completa | 285 | 4737 | 0/0/0 | 4 | Verificato |
| Finale source pertinente | 245 | 3610 | 0/0/0 | 3 | Verificato |
| Finale namespace installato completa | 285 | 4737 | 0/0/0 | 4 | Verificato |
| Finale namespace installato pertinente | 245 | 3610 | 0/0/0 | 3 | Verificato |
| R CMD check, esempi e test inclusi | — | 4737 nei test | 0/0/0 formali | 4 | Status: OK |
| Quarto pubblico installed, HTML/DOCX/PDF abilitati | — | 93 | 6 failure/3 error/4 warning | 0 | Non verificato/bloccato |

Il sottoinsieme pertinente comprende i file `test-mira-info_long*` e
`test-mira-report-freq_long*`; è estratto dagli stessi run, senza ripeterli.
Le tre esclusioni pertinenti sono HTML/DOCX/PDF reali; la quarta della suite
completa è il sampling CmdStan fuori perimetro. Skip non equivalgono a successo.
Target scientifico finale: 79/79; target report congelato: 378 pass, zero
failure/error/warning, due skip Quarto. Tempi elapsed: baseline 199.45 s,
source finale 258.84 s, installed finale 226.53 s; comprendono concorrenza e
non sono benchmark. Il processo test durante check registra 181.60 s.

Build e installazione del solo MIRA completate. SHA256 dell'archivio finale
`MIRA_0.0.1.tar.gz`:
`0E924180EB23C69EC0AD527AEA086E68101669E92AB093E138C61B5955CE91BA`.
Il namespace installed è caricato da `verified/library/MIRA`, non da pkgload.
Il check ha zero ERROR/WARNING/NOTE formali; il console log contiene quattro
warning ambientali per indici CRAN/Bioconductor irraggiungibili e un warning
`du` status 322. Questi cinque avvisi, già presenti nel contesto storico,
sono distinti dai warning dei test (zero). Nessuna dipendenza scaricata.

RNG, RNGkind, directory, device, opzioni esistenti e le cinque variabili
monitorate dal runner (TZ/LANG/LC_ALL/CMDSTAN/QUARTO_PATH) sono ripristinati.
Source: nessuna nuova opzione. Installed: 18 opzioni nuove; un processo fresco
che carica soltanto gli Imports riproduce esattamente le stesse 18, senza
modificare opzioni preesistenti. Il target scientifico controlla anche l'intero
ambiente del proprio processo. I messaggi informativi di fit singolare/R2
sono conservati, non trasformati in evidenza di adeguatezza del modello.

L'integrazione Quarto reale fallisce prima dell'avvio del renderer:
`processx_exec`, creazione write pipe, errore Windows 5 (accesso negato).
Gli errori downstream riguardano documenti assenti. Lo stesso blocco è stato
riprodotto in un piccolo driver pubblico source e nel profilo installed finale;
CLI e TeX presenti non bastano a superare questo vincolo. Non cambiati ACL,
toolchain o implementazione per aggirarlo. Nessun HTML/DOCX/PDF Quarto riuscito.

La catena **alternativa** installed → bundle pubblico render=FALSE → knitr →
Pandoc 3.8.3 ha invece prodotto HTML e DOCX reali. I verificatori validati
controllano celle N, media, SD, SE, CI, cambiamenti e p-value contro somme e
varianze manuali, riallineando le chiavi dei contrasti e selezionando i blocchi
univocamente. Medie 3.5/4.5/6.5, SD=sqrt(3.5); cambiamenti 1/3/2, varianze
1.2/2/0.8, p da t con df=5. Sei tabelle HTML e 14 tabelle Word effettive sono
conservati nell'archivio. Non è una prova del percorso pubblico Quarto.
PDF pubblico: **non verificato/bloccato**, nessuna conversione alternativa
PDF eseguita in questa sessione.

Il bundle installed separato passa: contrasto finale indipendente, 27 CSV
raw/manifest relativi, spostamento in directory isolata, esecuzione del runtime
in ambiente con parent baseenv e 15 figure catturate. Nessuna installazione.

## Matrice dei contratti

A = formule/enumerazioni indipendenti; B = confronto con engine, senza
certificarlo; C = coerenza della catena, che può propagare un errore comune.
I riferimenti seguenti sono sotto `tests/testthat/`; il CSV per blocco conserva
gli esiti effettivi. Una copertura adeguata riguarda il rischio descritto,
non tutti i possibili dati.

| Contratto/output → rischio | Campione ed estimando | Test e origine expected | Copertura iniziale → intervento e sensibilità | Esito e limite residuo |
| --- | --- | --- | --- | --- |
| N/media/SD/SE/CI; pairing, segno e cambiamenti | Visite disponibili e coppie finite per ID, baseline e consecutive | `statistical-oracles`, `sampling-inference`, `generated-properties`: A, differenze/varianze e distribuzioni standard | Già adeguata; preservata. Mutazione pairing mirata corrente | Verificato con limite; non prova validità sotto missingness informativa |
| Correlazioni e N effettive | Coppie disponibili; correlazioni perfette positive/negative | `statistical-oracles:141`, `boundaries-pathological:95`: A | Già adeguata per quei rischi; nessuna duplicazione | Verificato con limite del campione enumerato |
| Variabilità e ICC | Variabilità su profili incompleti distinta da ICC su completi | `sampling-inference:87`, `model-oracles:118`: A per componenti ridotte, B/C per fit | Già protetta; preservata | Verificato con limite; ICC del modello non certifica ottimizzatore |
| Welch, CI e Hedges; Holm finito | Tre visite, arm asimmetrici; unità soggetto; differenza active−control; Welch s²/n, df Satterthwaite; Hedges con pooled SD e J | `audit-2026-10-16`: A, formule centrate, pt/qt, Holm ordinato/cummax | Parziale (t.test come riferimento) → oracolo indipendente, ID espliciti, p non stimabile escluso dalla famiglia; mutazione arm-sign | Verificato con limite; copre stimandi dichiarati, non tutti i test di rango |
| Modelli: campione effettivo | 12 ID/36 righe long originarie; 24 righe/9 ID analizzati; covariate, arm e outcome missing; soggetto con sola visita conservato | `audit-2026-10-16`: A per 24 chiavi esplicite, model.frame misto e GEE independence | Parziale (solo conteggi) → identità esatte; controcampione con stesso N differisce | Verificato con limite; non esteso automaticamente a ogni engine/covariance |
| Covariate additive letterali | Nomi con backtick/backslash, due colonne distinte e disegno asimmetrico | `audit-2026-10-16`: simboli attesi e celle della matrice esplicite; fit primario e GEE | Lacuna → regressione fallita prima del fix; 14 failure contrattuali, zero error/warning | Verificato con limite; nomi ordinari compatibili |
| Friedman e ties | Profili completi; Q tie-corrected e Kendall W | `scientific-regressions:188`: A, Q=9/7, W=9/56, p=exp(−9/14) | Già adeguata per ties parziali; conservata | Verificato con limite: interamente tied produce NaN e performed TRUE come tentativo engine, non successo inferenziale |
| Wald, rango, termini TIME/ARM/interazione | Direzioni stimabili, covariance singolare, unità 1e−6…1e6 | `scientific-regressions:35/53/123/152/170`, `statistical-oracles:367`: A per Q/rango/indici; B/C per rinomina | Già adeguata per difetti storici; documentazione termini riallineata al codice | Verificato con limite; non tutti i fit reali con coefficienti eliminati |
| GEE independence e sandwich | 8 cluster×3 visite, Gaussian identity | `model-oracles:47`: A, OLS e score cluster indipendenti | Già adeguata sul caso analitico ridotto | Verificato con limite; non certifica exchangeable/AR1, QIC/CIC o CR2 |
| Mixed/ANOVA/emmeans | Disegno bilanciato, REML distinto da LRT ML; trend ordinale | `model-oracles:18/80/118`: A per estimandi/SS/trasformazioni; B per engine | Già protetti; preservati | Verificato con limite; GG/HF/Mauchly, df e ottimizzazione non convalidati universalmente |
| Slope/nlme/CR2 e failure | Strutture e inferenze specifiche, senza sostituire engine o covariance | `advanced`, `model-oracles:155/168`, `scientific-regressions:152`: B/C | Parziale; diagnostica e wrapper protetti, nessun nuovo studio inferenziale | Verificato con limite; identificabilità/singolarità restano limiti del modello |
| NA/NaN/Inf/−Inf e non stimabilità | Singleton, costanti, visite tutte missing, arm sparsi | `boundaries-pathological:217`, `sampling-inference:210/235/256`, nuovo Welch: A/B | Già protetta sui casi enumerati; non assimilati genericamente i valori speciali | Verificato con limite per scale estreme/quasi costanti non enumerate |
| Metamorfismi e cronologia | Configurazione fissata, permutazioni/ID/affine/direzione/livelli inutilizzati | `properties`, `affine-properties`, `generated-properties`, `configuration-contracts`: A/C | Già protetti, conservati; mutante fallback chronology dichiarato separatamente | Verificato con limite; un ramo non raggiunto non dimostra equivalenza |
| Payload/CSV/manifest/sezioni/relocation | Tutte le righe e colonne raw, chiavi relative; output parzialmente falliti | `bundle-contracts`, `files`, `integration`: C + A sui valori elementari | Già adeguati; mutazione CSV; relocation reale installed | Verificato con limite; no garanzia transazionale nuova |
| Celle HTML/DOCX e verificatori | Caption o header signature univoci; celle rettangolari, blocchi separati | `html-parser`, `document-parser`, `integration`: fixture note e negative; A/C | Parziale → errori espliciti prima degli indici; old HTML 12 failure+1 warning, nessun errore | Verificato con limite; parser intenzionalmente limitati a tabelle semplici, merge rifiutati |
| Figure prioritarie | 8 ID, arm asimmetrici/missing, cronologia risolta; mezzi e differenze nelle unità outcome, missing su n arm=4 | `figure-oracles`: A, n/somme/somme quadrati manuali, Welch, ggplot_build e disegno/caption | Non verificata numericamente → tre figure, dati e mapping effettivi | Verificato con limite; nessuna certificazione visiva di tutte le figure |
| Rendering/overwrite/cleanup | Output nuovo, vecchio documento, renderer no-op/parziale/wrong-path/promozione | `validation-regressions`, `rendering-contracts`, `files`: fault injection C | Già adeguato salvo zero byte → regressione 10 failure prima del fix; mutante false-success | Verificato con limite; mock distinto da Quarto reale; size>0 non valida contenuto |
| QTE/prior/Stan | API indipendenti escluse dall'intervento | Suite completa eseguita solo per regressioni indirette | Non applicabile motivato al lavoro nuovo; nessuna mutazione/sampling dedicata | Nessuna nuova convalida di queste famiglie |

## Expected, tolleranze e compatibilità

Nuovo Welch: control C1:C4/active A1:A4 a baseline; C1:C3/A1,A5 alla
visita intermedia; C1:C3/A1:A5 alla finale. Cambiamenti appaiati usano
C1:C3/A1 e C1:C3/A1:A4. X senza arm entra soltanto nei descrittivi core.
Le osservazioni non appaiabili hanno valori estremi che rendono discriminante
una selezione errata. Un controcampione dello stesso N dimostra che i soli
conteggi non proteggono il modello.

Welch/Hedges/Holm e figure usano 1e−12 per aritmetica moderata e quantili t;
conteggi, chiavi, segni e testi coincidono esattamente. Il confronto delle
matrici di disegno richiede dimensioni e ogni cella, escludendo soltanto
attributi engine estranei all'expected manuale. Per celle pubblicate a tre
decimali vale il limite di arrotondamento 0.0005000001, separato dalla
tolleranza dei raw. Non aumentate tolleranze esistenti né indebolite assertion.

Il confronto del caso deterministico ordinario (8 ID, due arm, tre visite,
missingness, correlazioni e arm tests) è identico prima/dopo per firme,
classi/schema, metadata, long_data, quantità scientifiche, payload, manifest
e tabelle CSV. Escluso soltanto `report$created_at`, motivato timestamp.

Un secondo confronto sul modello con covariata ordinaria non sintattica
`age years` è identico a tolleranza zero per formula, campioni primario/GEE,
metadata, coefficienti, covariance, componenti di varianza e quantità core.
Il riferimento originale è ricostruito fuori dal repository con `git show`
del HEAD iniziale pulito, senza checkout/commit; vengono confrontate le
componenti scientifiche e non l'identità degli ambienti R interni ai fit.

Differenze intenzionali: i nomi con backtick/backslash ora rappresentano le
colonne realmente selezionate; la regressione pubblica prima aveva quattro
coefficienti invece dei cinque attesi e passa dopo il fix. Per un renderer
vuoto `rendered` diventa FALSE, `files` è vuoto, compare una diagnosi e resta
il vecchio documento. Nessun estimando, metodo, famiglia p, default o schema
pubblico è stato cambiato.

## Mutazioni, simulazioni e affidabilità dei runner

Protocollo locale fissato prima del run: sette mutazioni univoche, piccole
fixture e test mirati, una baseline sulla copia non mutata prima delle sette
copie; budget sette processi test più baseline, nessuna ricerca di score.
Seed dei test/generatori preesistenti 20261012; nessuna nuova simulazione
stocastica nella campagna. Tutte le copie includono le modifiche locali,
senza .git. Baseline copiata: 95 blocchi, zero failure/error/warning, uno skip
PDF previsto, 56.28 s. Risultati della campagna corrente:

| Mutante | Ramo raggiunto | Failure pertinenti | Classificazione | Tipo di evidenza |
| --- | --- | ---: | --- | --- |
| arm-sign | Sì | 2 | detected | A, stime active−control indipendenti |
| wrong-subject-pairing | Sì | 243 | detected | A, differenze e campioni discriminanti |
| missing-BH-family | Sì | 10 | detected | A, famiglia/aggiustamento indipendenti |
| truncate-scientific-CSV | Sì | 6 | detected | C, perdita righe scientifiche nell'archivio |
| false-render-success | Sì | 4 | detected | Comportamento pubblico, vecchio file/no-op |
| parser-ragged | Sì | 1 | detected | Tecnica valida: recycling delle celle |
| ignore-custom-chronology | No | 0 | non-reached | Nessuna conclusione sulla sensibilità del ramo |

Nessun mutante invalido, blocked o survived nella campagna finale; il caso
non raggiunto **non è un successo né una prova di equivalenza**. Il mutante
parser emette anche un warning di recycling, ma è rilevato dall'assertion
negativa pertinente; non dal warning. CSV e falso rendering hanno ulteriori
failure coerenti: sono disponibili tutte, senza conteggiare errori accidentali
come rilevazioni. Nessuna mutazione dedicata a Stan/prior/QTE.

Nessuna nuova simulazione: le lacune affrontate hanno expected deterministici
costruibili. Lo studio normale appaiato storico non è stato rieseguito né
generalizzato a Welch, missingness informativa o modelli avanzati. I generatori
ordinari preesistenti sono eseguiti nella suite, con seed e ripristino RNG già
previsti. Nessuna nuova preregistrazione esterna.

Il runner mutazioni ora prova prima una copia non mutata del working tree
corrente, senza .git; richiede marker sull'espressione alterata e assertion
pertinente dichiarata prima del run. Errori di parsing/caricamento non sono
rilevazioni. Output, marker e log preesistenti impediscono sovrascritture.
La classificazione automatica è stata controllata leggendo i failure CSV.

Nel preflight sono emersi tre confronti di attributi accessori nelle nuove
matrici e un mock preesistente che creava un file vuoto come successo. Sono
stati corretti conservando dimensioni/celle/assertion. La revisione indipendente
ha inoltre corretto l'ipotesi sbagliata che n e mean DOCX fossero nello stesso
blocco: ora ogni blocco è identificato univocamente e ha una fixture propria.
Il primo baseline mutazioni si è fermato correttamente senza avviare mutanti.
Log preliminari conservati, non presentati come esiti finali.

## Riproduzione, limiti e JOSS

Comandi principali (PowerShell; scegliere una nuova directory esterna per
ogni esecuzione e verificare ogni exit code):

```powershell
$repo = 'C:\Users\david\OneDrive\Desktop\MIRA\MIRA'
$audit = 'C:\path\outside-repository\new-mira-audit'
$rscript = 'C:\Program Files\R\R-4.4.3\bin\Rscript.exe'
$r = 'C:\Program Files\R\R-4.4.3\bin\R.exe'
$env:LC_ALL = 'C'; $env:LANG = 'C'
New-Item -ItemType Directory -Path "$audit/build","$audit/library","$audit/mutations" | Out-Null
& $rscript --vanilla "$repo/dev/validate-tests.R" $repo "$audit/source" source
Push-Location "$audit/build"
& $r CMD build --no-build-vignettes --no-manual $repo
& $r CMD INSTALL "--library=$audit/library" MIRA_0.0.1.tar.gz
$env:R_LIBS = "$audit/library"
& $rscript --vanilla "$repo/dev/validate-tests.R" $repo "$audit/installed" installed
$env:MIRA_RUN_CMDSTAN_TESTS='false'
$env:MIRA_RUN_QUARTO_TESTS='false'; $env:MIRA_RUN_QUARTO_PDF_TESTS='false'
$env:QUARTO_TEXLIVE_AUTO_INSTALL='false'; $env:QUARTO_TEXLIVE_AUTO_UPDATE='false'
& $r CMD check --no-manual --no-vignettes MIRA_0.0.1.tar.gz
Pop-Location
Push-Location $repo
& $rscript --vanilla dev/validate-mutations.R "$audit/mutations" 'arm-sign,wrong-subject-pairing,missing-BH-family,truncate-scientific-CSV,false-render-success,parser-ragged,ignore-custom-chronology'
Pop-Location
```

Guida `dev/VALIDATION.md` per bundle e Quarto. Nel run attuale Quarto usa
QUARTO_PATH al CLI RStudio, QUARTO_R a Rscript 4.4.3, PATH di processo con
R 4.4.3 e TinyTeX, PDF esplicitamente abilitato e auto-install/update TeX
disabilitati. Il check non include manuale PDF, as-cran/incoming o CI remota.
Nell'archivio `validation/scripts` conserva i driver dei confronti e della
catena alternativa; `REPRODUCE.md` descrive conversione e ispezione dei documenti.

Rischi residui prioritari:

1. Inferenza sotto missingness informativa, covariance/df complessi, fit
   singolari/rank deficient e identificabilità nlme non sono certificati da
   queste fixture. Non cambiati metodi, correlazioni o p-value per convergere.
2. Quarto pubblico bloccato dall'ambiente; una catena alternativa riuscita non
   dimostra la promozione e il funzionamento del percorso pubblico Quarto.
3. Friedman interamente tied conserva NaN per Q/p/W con performed TRUE e
   reason NULL: risultato non stimabile, semantica storica del flag da chiarire
   prima di una futura modifica pubblica. Nessun p-value inventato.
4. R 4.2.3, Linux/macOS, matrice CI remota e CmdStan non eseguiti. Pin Matrix
   1.6-5 e Deriv 4.2.0 conservati; i workflow non sono stati dispatchati.
5. Parser stretti ma non generali; dimensione positiva non garantisce che un
   renderer produca un documento valido non vuoto. Nessuna promessa nuova
   di atomicità o rollback dell'intero bundle.

Consultati il 16 ottobre 2026 i [criteri ufficiali JOSS](https://joss.readthedocs.io/en/latest/review_criteria.html)
e le [istruzioni di submission](https://joss.readthedocs.io/en/latest/submitting.html).
Installazione da archivio, test oggettivi, esempi e istruzioni riproducibili
sono le evidenze pertinenti a installazione/funzionalità/test/documentazione.
Correttezza implementativa, validità statistica, qualità dei test e ammissibilità
editoriale sono distinte. Storia pubblica, uso scientifico, contributi umani e
revisione umana non sono attestati da questa sessione. Assistenza tramite Codex
per audit, codice/test e rapporto; nessuna dichiarazione inventata su modello
specifico, revisione umana o idoneità editoriale.

## File e stato finale

Modifiche di questa sessione soltanto:

- Produzione: `R/mira_info_long.R` (due punti quoting e roxygen termini),
  `R/mira_report_freq_long.R` (guard dimensione prima della promozione).
- Documentazione: `NEWS.md`, `man/mira_info_long.Rd`, `dev/VALIDATION.md`,
  nuovo `dev/VALIDATION_REPORT_2026-10-16.md`.
- Runner: `dev/validate-mutations.R`.
- Test modificati: helper report; `test-mira-report-freq_long-files.R`,
  `-html-parser.R`, `-integration.R`, `-validation-regressions.R`.
- Test nuovi: `test-mira-info_long-audit-2026-10-16.R`,
  `test-mira-report-freq_long-document-parser.R`,
  `test-mira-report-freq_long-figure-oracles.R`.

DESCRIPTION, NAMESPACE, CI/pin e materiale storico sono invariati.
HEAD finale uguale all'iniziale. Staging vuoto; nuovi file untracked
intenzionali. `git diff --check` passa. Gli avvisi Git LF→CRLF sono informativi
e non sono warning dei test; configurazione Git non modificata.
L'archivio contiene stato iniziale/finale, diff, inventario SHA256 dei file
testati e finali, hash del tarball, log/tempi/CSV/RDS e sorgenti finali interessati.
Il rapporto è aggiunto dopo i run; `dev` è esclusa dal build e nessun file R,
test, manuale o NEWS è cambiato dopo le verifiche canoniche.


