list_stan_models <- function(package = NULL) {

  # Se non viene specificato il package, usa il package corrente
  if (is.null(package)) {
    package <- utils::packageName()
  }

  # Percorso della cartella stan installata nel package
  stan_dir <- system.file("stan", package = package)

  if (stan_dir == "") {
    stop("La cartella 'inst/stan' non è stata trovata nel package.")
  }

  # Trova tutti i file .stan
  files <- list.files(
    stan_dir,
    pattern = "\\.stan$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    return(data.frame(
      name = character(),
      file = character(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    name = tools::file_path_sans_ext(basename(files)),
    file = files,
    stringsAsFactors = FALSE
  )
}
