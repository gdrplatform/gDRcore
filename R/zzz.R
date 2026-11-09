
.drugNameRegex <- NULL
.untreated_tag_patterns <- NULL
.untreatedDrugNameRegex <- NULL
patterns <- NULL

#' onload function
#'
#' @param libname library name
#' @param pkgname package name
#' @noRd
.onLoad <- function(libname, pkgname) {
  # scientific notation was disabled due to the problem with unit tests
  options(scipen = 999)

  cores <- Sys.getenv("NUM_CORES")
  # based on https://github.com/Bioconductor/BiocParallel/issues/98
  if (.Platform$OS.type != "windows" && cores != "") {
    BiocParallel::register(
      BiocParallel::MulticoreParam(workers = as.numeric(cores)),
      default = TRUE
    )
  } else {
    BiocParallel::register(
      BiocParallel::SerialParam(),
      default = TRUE
    )
  }
  # CONS
  drugs_id <- gDRutils::get_env_identifiers("drug_name")
  utils::assignInNamespace(".drugNameRegex", sprintf("^%s$|^%s_[[:digit:]]+$", drugs_id, drugs_id), ns = pkgname)

  untreated_tag_patterns <- vapply(
    gDRutils::get_env_identifiers("untreated_tag"),
    sprintf,
    fmt = "^%s$",
    character(1)
  )
  utils::assignInNamespace(".untreated_tag_patterns", untreated_tag_patterns, ns = pkgname)
  utils::assignInNamespace(".untreatedDrugNameRegex", paste(untreated_tag_patterns, collapse = "|"), ns = pkgname)

  # data.table compatible
  utils::assignInNamespace("patterns", data.table:::patterns, ns = pkgname)
  utils::globalVariables(
    c(
      ".",
      "..cotrt_var",
      "..present_ref_cols",
      "..y",
      "bliss_score",
      "column",
      "CorrectedReadout",
      "cotrt_value",
      "Duration",
      "hsa_score",
      "isDay0",
      "LogFoldChange",
      "normalization_type",
      "priority1",
      "priority2",
      "ratio",
      "ReadoutValue",
      "ReadoutValue_T0",
      "record_id",
      "rn",
      "smooth",
      "x",
      "x.N"
    ),
    pkgname)

  .register_fit_fns()
}


# The synergy fits this package exports, announced to the gDRutils registry so
# that a caller can list what is available instead of hardcoding it (GDR-3788).
# Registration is by name, so nothing is loaded here beyond the entry itself.
#
# MuSyC carries musyc_available as its predicate: it fits through the Python
# synergy package, so it is registered everywhere but only offered where that
# is provisioned.
#' @keywords internal
.register_fit_fns <- function() {
  if (!isTRUE(.fit_fn_registry_available())) {
    return(invisible(NULL))
  }
  gDRutils::register_fit_fn("bliss_fit_fn", package = "gDRcore",
                            data_type = "combination")
  gDRutils::register_fit_fn("hss_fit_fn", package = "gDRcore",
                            data_type = "combination")
  gDRutils::register_fit_fn("musyc_fit_fn", package = "gDRcore",
                            data_type = "combination",
                            available = musyc_available)
  invisible(NULL)
}


# gDRutils is a hard dependency, but the registry arrived in 1.11.24 and the
# DESCRIPTION floor is lower. Registering is an optional extra, so an older
# gDRutils should mean no catalogue rather than a package that will not load.
#' @keywords internal
.fit_fn_registry_available <- function() {
  exists("register_fit_fn", envir = asNamespace("gDRutils"), inherits = FALSE)
}
