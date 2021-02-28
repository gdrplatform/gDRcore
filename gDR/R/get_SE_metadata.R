#' @export
#'
set_SE_keys <- function(se, value) {
  .set_SE_metadata(se, name = "Keys", value)
}


#' @export
#'
get_SE_experiment_metadata <- function(se) {
  .get_SE_metadata(se, name = "experiment_metadata")
}


#' @export
#'
get_SE_keys <- function(se, key_type = NULL) {
  .get_SE_metadata(se, name = "Keys", subname = key_type)
}


###############
# Internals
###############

#' The primary purpose of this function is to allow other functions to create getter functions off of this one.
#'
.get_SE_metadata <- function(se, name, subname = NULL) {
  v <- S4Vectors::metadata(se)[[name]]
  if (!is.null(subname)) {
    v <- v[[subname]]
  }
  v
}


#' The primary purpose of this function is to allow other functions to create getter functions off of this one.
#'
.set_SE_metadata <- function(se, name, value) {
  if (!is.null(.get_SE_metadata(se, name))) {
    warning(sprintf("overwriting existing metadata entry: '%s'", name))
  }
  S4Vectors::metadata(se)[[name]] <- value
}
