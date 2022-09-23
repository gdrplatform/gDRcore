.onLoad <- function(libname, pkgname) {
  options(scipen = 999) # scientific notation was disabled due to the problem with unit tests
  # set default options for BiocParallel
  options(MulticoreParam = BiocParallel::MulticoreParam(BiocParallel::multicoreWorkers() - 2))
}

# data.table awareness
.datatable.aware <- TRUE
