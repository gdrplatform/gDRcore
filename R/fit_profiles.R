####
# Fit profile registry
####
#
# Profiles define the default slicing behaviour for each experiment type.
# They are loaded from inst/extdata/fit_profiles.json at package load time
# and held in a package-local environment so users can extend them at runtime.
#
# Fields per profile:
#   slicing_cols   — character vector; columns used to split each BumpyMatrix
#                    cell into sub-experiments (e.g. "normalization_type")
#   slicing_values — character vector or null; values to iterate; null = all
#                    unique values found in each cell
#   input_assay    — string; default source assay name
#   description    — string; human-readable description (optional)
#
# Public API:
#   get_fit_profiles()         — list all registered profiles
#   get_fit_profile(name)      — retrieve one profile
#   register_fit_profile(name, ...) — add or update a profile at runtime

.fit_profile_env <- new.env(parent = emptyenv())

#' @keywords internal
.load_fit_profiles <- function() {
  json_path <- system.file(
    "extdata", "fit_profiles.json",
    package = "gDRcore", mustWork = TRUE
  )
  profiles <- jsonlite::read_json(json_path, simplifyVector = TRUE)
  for (nm in names(profiles)) {
    p <- profiles[[nm]]
    # slicing_values may be null in JSON → keep as NULL in R
    sv <- if (length(p$slicing_values) == 0L) NULL else p$slicing_values
    .fit_profile_env[[nm]] <- list(
      slicing_cols   = p$slicing_cols,
      slicing_values = sv,
      input_assay    = p$input_assay,
      description    = p$description %||% ""
    )
  }
  invisible(NULL)
}

# Tiny helper — R has no built-in %||%
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Get all registered fit profiles
#'
#' Returns the list of all experiment-type profiles used by
#' \code{\link{apply_fit}} and \code{\link{apply_fits}} to resolve default
#' slicing configuration.
#'
#' @return Named list of profile definitions. Each element contains
#'   \code{slicing_cols}, \code{slicing_values}, \code{input_assay}, and
#'   \code{description}.
#'
#' @examples
#' get_fit_profiles()
#'
#' @export
get_fit_profiles <- function() {
  as.list(.fit_profile_env)
}


#' Get a single fit profile by name
#'
#' @param name string; profile name (e.g. \code{"single-agent"}).
#'
#' @return Named list with \code{slicing_cols}, \code{slicing_values},
#'   \code{input_assay}, and \code{description}.  Stops if the profile
#'   does not exist.
#'
#' @examples
#' get_fit_profile("single-agent")
#'
#' @export
get_fit_profile <- function(name) {
  checkmate::assert_string(name)
  p <- .fit_profile_env[[name]]
  if (is.null(p)) {
    stop(sprintf(
      "Unknown fit profile '%s'. Available: %s",
      name, paste(ls(.fit_profile_env), collapse = ", ")
    ))
  }
  p
}


#' Register or update a fit profile
#'
#' Adds a new experiment-type profile or overwrites an existing one.
#' The profile is available immediately in the current R session and is
#' picked up by \code{\link{apply_fit}} and \code{\link{apply_fits}}.
#'
#' @param name string; unique profile identifier (e.g.
#'   \code{"biochemical"}).
#' @param slicing_cols character vector of column names used to split each
#'   BumpyMatrix cell into sub-experiments.
#' @param slicing_values character vector of values to iterate, or
#'   \code{NULL} to iterate all unique values found in each cell.
#' @param input_assay string; name of the default source assay.
#' @param description string; human-readable description (optional).
#'
#' @return Invisibly returns the registered profile list.
#'
#' @examples
#' register_fit_profile(
#'   "biochemical",
#'   slicing_cols   = "assay_type",
#'   slicing_values = c("Ki", "IC50"),
#'   input_assay    = "Biochemical",
#'   description    = "Biochemical activity assay (Ki / IC50)"
#' )
#' get_fit_profile("biochemical")
#'
#' @importFrom checkmate assert_string assert_character
#' @export
register_fit_profile <- function(name,
                                 slicing_cols,
                                 slicing_values = NULL,
                                 input_assay,
                                 description = "") {
  checkmate::assert_string(name)
  checkmate::assert_character(slicing_cols, min.len = 1L)
  checkmate::assert_character(slicing_values, null.ok = TRUE, min.len = 1L)
  checkmate::assert_string(input_assay)
  checkmate::assert_string(description)

  profile <- list(
    slicing_cols   = slicing_cols,
    slicing_values = slicing_values,
    input_assay    = input_assay,
    description    = description
  )
  .fit_profile_env[[name]] <- profile
  invisible(profile)
}
