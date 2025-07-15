#' @import S4Vectors
#' @import checkmate
#' @useDynLib gDRcore, .registration = TRUE

# data.table awareness
.datatable.aware <- TRUE
patterns <- data.table:::patterns

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
}

# Prevent R CMD check from complaining about the use of pipe expressions
# standard data.table variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "normalization_type",
      ".",
      "rn",
      "column",
      "CorrectedReadout",
      "cotrt_value",
      "x",
      "Duration",
      "isDay0",
      "record_id",
      "ratio",
      "smooth",
      "priority1",
      "priority2",
      "x.N"
    ), 
    utils::packageName())
}

