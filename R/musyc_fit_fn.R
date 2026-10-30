.ref_musyc <- paste(
  "Meyer et al. 2019, Cell Systems 8(2):97-108, doi:10.1016/j.cels.2019.01.003;",
  "Wooten et al. 2021, Nature Communications 12:4607, doi:10.1038/s41467-021-24789-z"
)


#' musyc_fit_fn
#'
#' Fit function for the MuSyC (multi-dimensional synergy of combinations)
#' model on combination dose-response data, for use with
#' \code{\link[gDRutils]{apply_fit}}.
#'
#' The surface is fitted with the Python \code{synergy} package through
#' \code{reticulate}; install it with
#' \code{reticulate::py_install("synergy", pip = TRUE)}.  The fit provenance
#' records the Python, \code{synergy}, \code{numpy} and \code{scipy} versions.
#'
#' The whole grid is fitted, single-agent edges included.  The untreated
#' corner is set to 1 (added when absent), parameters are bounded
#' (\code{E} from -1 for GR and from 0 for RV, up to 1.5) and the fit is
#' restarted from random initial values until \code{r2} exceeds
#' \code{r2_threshold} or \code{patience} restarts in a row bring no
#' improvement.
#'
#' \code{alpha > 1} means synergistic potency, \code{beta > 0} synergistic
#' efficacy.
#'
#' @param avg_dt \code{data.table} for one (drug1 x drug2 x cell line x
#'   normalization_type) combination.  Must contain columns \code{Concentration},
#'   \code{Concentration_2}, \code{x}, and \code{normalization_type}.
#' @param fit_gamma logical; whether to fit the cooperativity parameters
#'   \code{gamma12} and \code{gamma21}.
#' @param n_restarts integer; maximum number of restarts from random initial
#'   values.
#' @param r2_threshold numeric; restarts stop once \code{r2} exceeds this value.
#' @param patience integer; restarts stop after this many restarts in a row
#'   without improvement.
#' @param seed integer; seed for the random initial values.  The global random
#'   number generator is left as it was.
#'
#' @return Named list with \code{normalization_type}, the parameters
#'   \code{E0}, \code{E1}, \code{E2}, \code{E3}, \code{h1}, \code{h2},
#'   \code{C1}, \code{C2}, \code{alpha12}, \code{alpha21}, the efficacy synergy
#'   \code{beta}, \code{r2} and \code{converged}.  Numeric fields are
#'   \code{NA_real_} when the fit does not converge.
#'
#' @examples
#' if (requireNamespace("reticulate", quietly = TRUE) &&
#'     reticulate::py_module_available("synergy")) {
#'   mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#'   combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#'   combo_se_out <- gDRutils::apply_fit(
#'     combo_se, musyc_fit_fn, "combination",
#'     output_assay = "musyc_params",
#'     summary_fn = musyc_summary_fn, summary_assay = "musyc_summary",
#'     fit_source = "musyc"
#'   )
#' }
#'
#' @seealso \code{\link{musyc_summary_fn}}
#' @export
musyc_fit_fn <- function(avg_dt,
                         fit_gamma = FALSE,
                         n_restarts = 100L,
                         r2_threshold = 0.5,
                         patience = 20L,
                         seed = 1L) {
  checkmate::assert_data_table(avg_dt)
  checkmate::assert_flag(fit_gamma)
  checkmate::assert_count(n_restarts)
  checkmate::assert_number(r2_threshold)
  checkmate::assert_count(patience)
  checkmate::assert_int(seed)

  norm_type <- as.character(avg_dt$normalization_type[1])
  conc1_col <- gDRutils::get_env_identifiers("concentration")
  conc2_col <- gDRutils::get_env_identifiers("concentration2")
  na_result <- .musyc_na_result(norm_type)

  obs <- .musyc_grid(avg_dt[[conc1_col]], avg_dt[[conc2_col]], avg_dt[["x"]])
  conc1 <- unique(obs$d1)
  conc2 <- unique(obs$d2)
  if (length(conc1) < 2L || length(conc2) < 2L || NROW(obs) < 4L) {
    warning(sprintf(
      "musyc_fit_fn: insufficient data (%d x %d concentrations, %d points) for %s",
      length(conc1), length(conc2), NROW(obs), norm_type
    ))
    return(na_result)
  }

  .check_synergy_available()
  e_lower <- if (norm_type == "RV") 0 else -1
  bounds <- .musyc_bounds(e_lower)
  corners <- .musyc_corners(obs, e_lower)
  med_d1 <- stats::median(conc1[conc1 > 0])
  med_d2 <- stats::median(conc2[conc2 > 0])
  p0_default <- c(corners, 1, 1, med_d1, med_d2, 1, 1)
  log_low <- log10(c(0.01, 0.01, med_d1 / 100, med_d2 / 100, 0.001, 0.001))
  log_high <- log10(c(100, 100, med_d1 * 100, med_d2 * 100, 1000, 1000))

  synergy <- reticulate::import("synergy.combination", convert = FALSE)
  warnings_mod <- reticulate::import("warnings", convert = FALSE)
  warnings_mod$filterwarnings("ignore")

  best <- .musyc_fit_once(synergy, obs, bounds, p0_default, fit_gamma)
  random_p0 <- .with_local_seed(seed, lapply(seq_len(n_restarts), function(i) {
    c(corners, 10^stats::runif(length(log_low), log_low, log_high))
  }))
  no_improve <- 0L
  for (p0 in random_p0) {
    if (isTRUE(best$r2 > r2_threshold) || no_improve >= patience) {
      break
    }
    candidate <- .musyc_fit_once(synergy, obs, bounds, p0, fit_gamma)
    if (is.na(best$r2) || isTRUE(candidate$r2 > best$r2)) {
      best <- candidate
      no_improve <- 0L
    } else {
      no_improve <- no_improve + 1L
    }
  }

  if (is.null(best$params)) {
    return(na_result)
  }
  c(
    list(normalization_type = norm_type),
    best$params,
    list(beta = best$beta, r2 = best$r2, converged = TRUE)
  )
}

musyc_fit_fn <- gDRutils::set_fit_fn_info(
  musyc_fit_fn,
  model = "MuSyC (multi-dimensional synergy of combinations)",
  description = paste(
    "Fits the bounded MuSyC dose-response surface to the whole grid with the Python",
    "synergy package, restarting from random initial values, and returns its parameters,",
    "the efficacy synergy beta and the R-squared of the fit."
  ),
  reference = .ref_musyc,
  env = function() .musyc_runtime()
)


#' musyc_summary_fn
#'
#' Summary function for \code{\link{musyc_fit_fn}}: called once per
#' (drug1 x drug2 x cell line) cell on the fits of all normalization types.
#'
#' @param fit_dt \code{data.table} of the \code{musyc_fit_fn} results for one
#'   cell.
#'
#' @return Named list with \code{mean_beta}, \code{mean_log2_alpha} (mean of
#'   \code{log2(alpha12)} and \code{log2(alpha21)}), \code{mean_r2},
#'   \code{any_converged} and \code{synergy_call} (\code{"synergistic"},
#'   \code{"antagonistic"}, \code{"additive"} or \code{"unknown"}, from
#'   \code{mean_beta} with a dead zone of 0.1).
#'
#' @examples
#' fit_dt <- data.table::data.table(
#'   normalization_type = c("GR", "RV"), alpha12 = c(2, 1.5),
#'   alpha21 = c(1.2, 1), beta = c(0.3, 0.2), r2 = c(0.9, 0.95),
#'   converged = TRUE
#' )
#' musyc_summary_fn(fit_dt)
#'
#' @seealso \code{\link{musyc_fit_fn}}
#' @export
musyc_summary_fn <- function(fit_dt) {
  checkmate::assert_data_table(fit_dt)

  converged <- fit_dt[["converged"]] %in% TRUE
  mean_beta <- mean(fit_dt[["beta"]][converged], na.rm = TRUE)
  mean_log2_alpha <- mean(
    log2(c(fit_dt[["alpha12"]][converged], fit_dt[["alpha21"]][converged])),
    na.rm = TRUE
  )
  synergy_call <- if (is.na(mean_beta)) {
    "unknown"
  } else if (mean_beta > 0.1) {
    "synergistic"
  } else if (mean_beta < -0.1) {
    "antagonistic"
  } else {
    "additive"
  }

  list(
    mean_beta = mean_beta,
    mean_log2_alpha = mean_log2_alpha,
    mean_r2 = mean(fit_dt[["r2"]][converged], na.rm = TRUE),
    any_converged = any(converged),
    synergy_call = synergy_call
  )
}

musyc_summary_fn <- gDRutils::set_fit_fn_info(
  musyc_summary_fn,
  model = "MuSyC synergy summary",
  description = paste(
    "Averages beta, log2 alpha and R-squared of the MuSyC fits over the normalization types",
    "and calls the pair synergistic, antagonistic or additive from the mean beta."
  )
)


#' @keywords internal
.musyc_na_result <- function(norm_type) {
  list(
    normalization_type = norm_type,
    E0 = NA_real_, E1 = NA_real_, E2 = NA_real_, E3 = NA_real_,
    h1 = NA_real_, h2 = NA_real_, C1 = NA_real_, C2 = NA_real_,
    alpha12 = NA_real_, alpha21 = NA_real_,
    beta = NA_real_, r2 = NA_real_, converged = FALSE
  )
}

# Observed points of the grid without NAs, with the untreated corner set to 1
#' @keywords internal
.musyc_grid <- function(d1, d2, x) {
  obs <- data.table::data.table(d1 = d1, d2 = d2, E = x)
  obs <- obs[!is.na(obs$E), ]
  obs <- obs[!(obs$d1 == 0 & obs$d2 == 0), ]
  obs <- rbind(data.table::data.table(d1 = 0, d2 = 0, E = 1), obs)
  obs[order(obs$d1, obs$d2), ]
}

#' @keywords internal
.musyc_bounds <- function(e_lower) {
  list(
    E0_bounds = c(0.99, 1.01),
    E1_bounds = c(e_lower, 1.5),
    E2_bounds = c(e_lower, 1.5),
    E3_bounds = c(e_lower, 1.5),
    h1_bounds = c(0.01, 100),
    h2_bounds = c(0.01, 100),
    C1_bounds = c(1e-4, 1e4),
    C2_bounds = c(1e-4, 1e4),
    alpha12_bounds = c(1e-3, 1e3),
    alpha21_bounds = c(1e-3, 1e3)
  )
}

# Initial E0..E3 from the grid corners, clipped to [e_lower, 1]; a missing
# corner falls back to the median response
#' @keywords internal
.musyc_corners <- function(obs, e_lower) {
  d1_max <- max(obs$d1)
  d2_max <- max(obs$d2)
  corner <- function(a, b) {
    value <- obs$E[obs$d1 == a & obs$d2 == b]
    if (length(value) == 0L) stats::median(obs$E) else value[1]
  }
  values <- c(corner(0, 0), corner(d1_max, 0), corner(0, d2_max), corner(d1_max, d2_max))
  pmin(pmax(values, e_lower), 1)
}

#' @keywords internal
.musyc_fit_once <- function(synergy, obs, bounds, p0, fit_gamma) {
  failed <- list(params = NULL, beta = NA_real_, r2 = NA_real_)
  bounds <- lapply(bounds, function(b) reticulate::tuple(b[1], b[2]))
  model <- do.call(synergy$MuSyC, c(list(fit_gamma = fit_gamma), bounds))
  if (fit_gamma) {
    p0 <- c(p0, 1, 1)
  }
  fitted <- tryCatch({
    model$fit(
      reticulate::np_array(obs$d1),
      reticulate::np_array(obs$d2),
      reticulate::np_array(obs$E),
      p0 = reticulate::r_to_py(as.list(p0))
    )
    reticulate::py_to_r(model$is_specified)
  }, error = function(e) FALSE)
  if (!isTRUE(fitted)) {
    return(failed)
  }

  param_names <- c("E0", "E1", "E2", "E3", "h1", "h2", "C1", "C2", "alpha12", "alpha21")
  params <- lapply(reticulate::py_to_r(model$get_parameters())[param_names], as.numeric)
  predicted <- as.numeric(reticulate::py_to_r(
    model$E(reticulate::np_array(obs$d1), reticulate::np_array(obs$d2))
  ))
  ss_res <- sum((obs$E - predicted)^2)
  ss_tot <- sum((obs$E - mean(obs$E))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  if (!is.finite(r2)) {
    return(failed)
  }
  list(params = params, beta = as.numeric(reticulate::py_to_r(model$beta)), r2 = r2)
}

#' @keywords internal
.check_synergy_available <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE) ||
      !reticulate::py_module_available("synergy")) {
    stop(
      "musyc_fit_fn requires the Python 'synergy' package through 'reticulate'; ",
      "install it with reticulate::py_install(\"synergy\", pip = TRUE)"
    )
  }
}

# Evaluate expr with the given seed and restore the global RNG state afterwards
#' @keywords internal
.with_local_seed <- function(seed, expr) {
  old_seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  })
  set.seed(seed)
  expr
}

# Versions of the Python runtime a MuSyC fit depends on, recorded in the fit
# provenance each time the fit runs (see gDRutils::set_fit_fn_info)
#' @keywords internal
.musyc_runtime <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE) ||
      !reticulate::py_available(initialize = FALSE)) {
    return(list(python = NA_character_))
  }
  py_metadata <- reticulate::import("importlib.metadata")
  version_of <- function(pkg) {
    tryCatch(py_metadata$version(pkg), error = function(e) NA_character_)
  }
  list(
    python = reticulate::py_config()$version,
    synergy = version_of("synergy"),
    numpy = version_of("numpy"),
    scipy = version_of("scipy")
  )
}
