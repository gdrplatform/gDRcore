.onLoad <- function(libname, pkgname) {
  options(scipen = 999) # scientific notation was disabled due to the problem with unit tests
}

# data.table awareness
.datatable.aware <- TRUE
cores <- parallel::detectCores() - 1
