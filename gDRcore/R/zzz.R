.onLoad <- function(libname, pkgname) {
  options(scipen = 999) # scientific notation was disabled due to the problem with unit tests
}

# data.table awareness
.datatable.aware <- TRUE
cores <- as.numeric(Sys.getenv("NUM_CORES"))
if (is.na(cores)) {
  cores <- parallel::detectCores() - 1
}
