.onLoad <- function(libname, pkgname) {
  # scientific notation was disabled due to the problem with unit tests
  options(scipen = 999) 
  
  cores <- Sys.getenv("NUM_CORES")
  if (cores != "") {
    BiocParallel::register(
      BiocParallel::MulticoreParam(workers = as.numeric(cores)), 
      default = TRUE
    )
  }
}

# data.table awareness
.datatable.aware <- TRUE
