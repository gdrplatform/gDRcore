.onLoad <- function(libname, pkgname) {
  options(scipen = 999) # scientific notation was disabled due to the problem with unit tests
  # set default options for BiocParallel
  options(MulticoreParam = BiocParallel::MulticoreParam(
    workers = max(1,
                  BiocParallelmulticoreWorkers() - 5)))
}

# data.table awareness
.datatable.aware <- TRUE
