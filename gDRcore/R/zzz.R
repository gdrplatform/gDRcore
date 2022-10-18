.onLoad <- function(libname, pkgname) {
  options(scipen = 999) # scientific notation was disabled due to the problem with unit tests
  cores <- Sys.getenv("NUM_CORES")
  if (cores != "") {
    BiocParallel::register(BiocParallel::MulticoreParam(workers = as.numeric(cores)), default = TRUE)
  }
}

# data.table awareness
.datatable.aware <- TRUE
