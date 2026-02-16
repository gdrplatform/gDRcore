
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
  .drugNameRegex <<- sprintf("^%s$|^%s_[[:digit:]]+$", drugs_id, drugs_id)
  
  .untreated_tag_patterns <<- vapply(
    gDRutils::get_env_identifiers("untreated_tag"),
    sprintf,
    fmt = "^%s$",
    character(1)
  )
  .untreatedDrugNameRegex <<- paste(.untreated_tag_patterns, collapse = "|")
  
  # data.table compatible
  patterns <<- data.table:::patterns
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
}
