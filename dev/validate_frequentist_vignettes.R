# Validate the installed package and the complete frequentist collection.
# Rscript --vanilla dev/validate_frequentist_vignettes.R <outside-output-dir>
# Optional: MIRA_RUN_ADVANCED_VIGNETTES=true
# Optional, separate run: MIRA_RENDER_VIGNETTE_FORMATS=html,docx,pdf
# MIRA_VALIDATION_LIBRARY asserts the installed MIRA location.
# Every render, extracted-code run and reader-visible-code run has its own
# --vanilla R process. No pkgdown or external installation is performed.
args <- commandArgs(trailingOnly = TRUE)
compare_path <- function(x) if (.Platform$OS.type == "windows") tolower(x) else x
inside <- function(path, root) {
  path <- compare_path(sub("/+$", "", path))
  root <- compare_path(sub("/+$", "", root))
  identical(path, root) || startsWith(path, paste0(root, "/"))
}
visible_code <- function(file) {
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  result <- character()
  active <- FALSE
  keep <- FALSE
  for (line in lines) {
    if (!active && grepl("^```\\{r([ ,}]|$)", line)) {
      active <- TRUE
      keep <- !grepl("(include|echo|eval)[[:space:]]*=[[:space:]]*FALSE", line)
      if (keep) result <- c(result, paste0("# ", line))
    } else if (active && grepl("^```[[:space:]]*$", line)) {
      active <- FALSE
      keep <- FALSE
      result <- c(result, "")
    } else if (active && keep) result <- c(result, line)
  }
  if (active) stop("Unclosed R chunk in ", file)
  result
}
worker <- function(config) {
  cfg <- readRDS(config)
  .libPaths(cfg$libraries)
  setwd(cfg$stage)
  stopifnot(requireNamespace("MIRA", quietly = TRUE),
    identical(compare_path(normalizePath(find.package("MIRA"), winslash = "/")),
              compare_path(cfg$package_path)))
  writeLines(c(paste("Mode:", cfg$mode), paste("Input:", cfg$source),
    paste("Input MD5:", unname(tools::md5sum(cfg$source))),
    paste("Installed package:", find.package("MIRA")),
    paste("Advanced:", Sys.getenv("MIRA_RUN_ADVANCED_VIGNETTES", "false")),
    paste("External formats:", Sys.getenv("MIRA_RENDER_VIGNETTE_FORMATS", "")),
    capture.output(sessionInfo())), file.path(cfg$stage, "context.txt"))
  execution <- new.env(parent = globalenv())
  if (identical(cfg$mode, "render")) {
    target <- rmarkdown::render(cfg$source, output_dir = cfg$output,
      intermediates_dir = cfg$stage, envir = execution, clean = TRUE,
      quiet = FALSE)
    stopifnot(file.exists(target), file.info(target)$size > 1000L)
    # These compact, named analytical ledgers are useful for independent review.
    keys <- intersect(c("module_status", "sample_ledger", "formula_ledger",
      "sample_audit", "endpoint", "adjusted_endpoint", "adjusted_conclusion",
      "render_status", "transfer_solution", "transfer", "covariate_sample_audit",
      "primary_contrast", "sample_counts"),
      ls(execution))
    for (key in keys) {
      value <- get(key, execution)
      if (is.data.frame(value)) utils::write.csv(value,
        file.path(cfg$stage, paste0(key, ".csv")), row.names = FALSE)
      else if (is.character(value)) writeLines(value,
        file.path(cfg$stage, paste0(key, ".txt")))
    }
    # Inspect the rendered DOM, not the number of kable calls in the source.
    if (requireNamespace("xml2", quietly = TRUE)) {
      doc <- xml2::read_html(target)
      tables <- xml2::xml_find_all(doc, "//table")
      cells <- list()
      for (i in seq_along(tables)) {
        caption <- xml2::xml_text(xml2::xml_find_first(tables[[i]], "./caption"))
        rows <- xml2::xml_find_all(tables[[i]], ".//tr")
        for (j in seq_along(rows)) {
          entries <- xml2::xml_find_all(rows[[j]], "./th|./td")
          if (length(entries)) cells[[length(cells) + 1L]] <- data.frame(
            table = i, caption = caption, row = j,
            column = seq_along(entries), header = xml2::xml_name(entries) == "th",
            value = trimws(xml2::xml_text(entries)), stringsAsFactors = FALSE)
        }
      }
      if (length(cells)) utils::write.csv(do.call(rbind, cells),
        file.path(cfg$stage, "rendered-table-cells.csv"), row.names = FALSE)
      writeLines(paste("Rendered tables:", length(tables)),
                 file.path(cfg$stage, "dom-inspection.txt"))
    } else writeLines("DOM cell extraction unavailable: xml2 is not installed.",
                      file.path(cfg$stage, "dom-inspection.txt"))
  } else {
    code <- file.path(cfg$stage, paste0(cfg$name, ".R"))
    if (identical(cfg$mode, "extracted")) {
      knitr::purl(cfg$source, output = code, documentation = 0L, quiet = TRUE)
    } else if (identical(cfg$mode, "visible")) {
      writeLines(visible_code(cfg$source), code)
    } else stop("Unknown worker mode.")
    invisible(parse(file = code))
    # A script-only run may draw figures; keep its device outside the repository.
    grDevices::pdf(file.path(cfg$stage, "script-figures.pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    sys.source(code, envir = execution)
  }
  invisible(TRUE)
}
if (length(args) == 2L && identical(args[[1L]], "--worker")) {
  worker(args[[2L]])
  quit(save = "no", status = 0L)
}
if (length(args) != 1L) stop("Supply one output directory outside the package source tree.")
repo <- normalizePath(".", winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(repo, "DESCRIPTION"))) stop("Run from the MIRA source directory.")
# Resolve the existing ancestor and normalize .. components before mkdir.
ancestor <- path.expand(args[[1L]])
suffix <- character()
while (!dir.exists(ancestor)) {
  parent <- dirname(ancestor)
  if (identical(parent, ancestor)) stop("Cannot resolve the output directory.")
  suffix <- c(basename(ancestor), suffix)
  ancestor <- parent
}
out <- normalizePath(ancestor, winslash = "/", mustWork = TRUE)
for (component in suffix) {
  if (identical(component, "..")) out <- dirname(out)
  else if (!identical(component, ".")) out <- file.path(out, component)
}
out <- gsub("\\\\", "/", out)
if (inside(out, repo)) stop("The validation output directory must be outside the package tree.")
# The reading index is the architecture. Orphan frequentist articles are errors;
# unrelated future Bayesian articles are neither rendered nor rejected here.
index <- readLines("vignettes/MIRA.Rmd", warn = FALSE, encoding = "UTF-8")
links <- unlist(regmatches(index,
  gregexpr("frequentist-long-[A-Za-z0-9_-]+[.]html", index)))
expected <- c("MIRA", unique(sub("[.]html$", "", links)))
available <- tools::file_path_sans_ext(list.files("vignettes",
  pattern = "^(MIRA|frequentist-long-.+)[.]Rmd$"))
if (!setequal(expected, available)) stop("Index/source mismatch. Missing sources: ",
  paste(setdiff(expected, available), collapse = ", "), "; unlisted sources: ",
  paste(setdiff(available, expected), collapse = ", "))
files <- file.path("vignettes", paste0(expected, ".Rmd"))
stopifnot(requireNamespace("MIRA", quietly = TRUE),
          requireNamespace("rmarkdown", quietly = TRUE))
package_path <- normalizePath(find.package("MIRA"), winslash = "/", mustWork = TRUE)
requested_library <- Sys.getenv("MIRA_VALIDATION_LIBRARY", "")
if (nzchar(requested_library)) stopifnot(identical(
  compare_path(dirname(package_path)),
  compare_path(normalizePath(requested_library, winslash = "/", mustWork = TRUE))))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/", mustWork = TRUE)
if (inside(out, repo)) stop("Resolved output directory is inside the package tree.")
evidence <- file.path(out, "validation")
dir.create(evidence, showWarnings = FALSE)
utils::write.csv(data.frame(order = seq_along(expected), article = expected,
  source = files), file.path(evidence, "architecture.csv"), row.names = FALSE)
source_files <- c("DESCRIPTION", "NAMESPACE", list.files("R", full.names = TRUE),
  files, list.files("tests", recursive = TRUE, full.names = TRUE),
  "dev/validate_frequentist_vignettes.R")
source_md5 <- unname(tools::md5sum(source_files))
utils::write.csv(data.frame(file = source_files, md5 = source_md5),
  file.path(evidence, "validated-source-md5.csv"), row.names = FALSE)
writeLines(c(paste("MIRA version:", utils::packageVersion("MIRA")),
  paste("Installed path:", package_path), paste("R:", R.version.string),
  paste("Libraries:", paste(.libPaths(), collapse = ";")),
  paste("Advanced branch:", Sys.getenv("MIRA_RUN_ADVANCED_VIGNETTES", "false")),
  paste("External formats:", Sys.getenv("MIRA_RENDER_VIGNETTE_FORMATS", "")),
  "Isolation: one --vanilla process per article and execution mode",
  "HTML layout: flat collection; temporary inputs/intermediates under validation/"),
  file.path(evidence, "validation-context.txt"))
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
runner <- normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = TRUE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
status <- list()
for (i in seq_along(files)) {
  name <- expected[[i]]
  for (mode in c("render", "extracted", "visible")) {
    stage <- file.path(evidence, name, mode)
    dir.create(stage, recursive = TRUE, showWarnings = FALSE)
    stage <- normalizePath(stage, winslash = "/", mustWork = TRUE)
    copied <- file.path(stage, basename(files[[i]]))
    stopifnot(file.copy(files[[i]], copied, overwrite = TRUE),
      identical(unname(tools::md5sum(copied)), source_md5[match(files[[i]], source_files)]))
    config <- file.path(stage, "worker-config.rds")
    saveRDS(list(source = copied, name = name, stage = stage, mode = mode,
      output = out, package_path = package_path, libraries = .libPaths()), config)
    started <- proc.time()[["elapsed"]]
    log <- file.path(stage, "execution.log")
    result <- tryCatch(system2(rscript, args = c("--vanilla", shQuote(runner),
      "--worker", shQuote(config)), stdout = log, stderr = log),
      error = function(e) { writeLines(conditionMessage(e), log); 127L })
    status[[length(status) + 1L]] <- data.frame(vignette = name, mode = mode,
      success = identical(as.integer(result), 0L), exit_status = as.integer(result),
      seconds = proc.time()[["elapsed"]] - started,
      source_md5 = unname(tools::md5sum(files[[i]])),
      log = substring(log, nchar(out) + 2L), stringsAsFactors = FALSE)
    utils::write.csv(do.call(rbind, status), file.path(evidence, "vignette-status.csv"), row.names = FALSE)
    cat(name, mode, if (result == 0L) "PASS" else "FAIL", "\n")
  }
}
status <- do.call(rbind, status)
# Verify local targets from each actual flat HTML location, including anchors.
link_rows <- list()
if (requireNamespace("xml2", quietly = TRUE)) {
  for (name in expected) {
    html <- file.path(out, paste0(name, ".html"))
    if (!file.exists(html)) next
    doc <- xml2::read_html(html)
    hrefs <- unique(xml2::xml_attr(xml2::xml_find_all(doc, "//a[@href]"), "href"))
    hrefs <- hrefs[grepl("^(MIRA|frequentist-long-).*[.]html(#.*)?$", hrefs)]
    for (href in hrefs) {
      target <- sub("#.*$", "", href)
      path <- file.path(dirname(html), target)
      exists <- file.exists(path)
      anchor <- if (grepl("#", href, fixed = TRUE)) sub("^[^#]*#", "", href) else ""
      anchor_ok <- TRUE
      if (exists && nzchar(anchor)) {
        ids <- xml2::xml_attr(xml2::xml_find_all(xml2::read_html(path), "//*[@id]"), "id")
        anchor_ok <- utils::URLdecode(anchor) %in% ids
      }
      link_rows[[length(link_rows) + 1L]] <- data.frame(article = name,
        href = href, target_exists = exists, anchor_exists = anchor_ok)
    }
  }
  if (length(link_rows)) {
    link_rows <- do.call(rbind, link_rows)
    utils::write.csv(link_rows, file.path(evidence, "flat-html-links.csv"), row.names = FALSE)
    stopifnot(all(link_rows$target_exists), all(link_rows$anchor_exists))
  }
}
stopifnot(identical(unname(tools::md5sum(source_files)), source_md5))
writeLines(capture.output(sessionInfo()), file.path(evidence, "session-info.txt"))
print(status)
if (!all(status$success)) stop("Vignette validation failed; inspect per-mode logs.")
if (identical(Sys.getenv("MIRA_VALIDATE_INSTALLED_VIGNETTES"), "true")) {
  installed <- utils::vignette(package = "MIRA")$results
  stopifnot(all(expected %in% installed[, "Item"]))
  installed_doc <- file.path(package_path, "doc")
  stopifnot(all(file.exists(file.path(installed_doc, paste0(expected, ".html")))),
            all(file.exists(file.path(installed_doc, paste0(expected, ".Rmd")))))
  utils::write.csv(installed, file.path(evidence, "installed-vignettes.csv"), row.names = FALSE)
  # The generated index.html and future Bayesian articles are allowed.
}
