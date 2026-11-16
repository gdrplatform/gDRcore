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
#' The surface (with \code{gamma12 = gamma21 = 1}) is fitted in R by bounded
#' least squares, \code{\link[stats]{nlminb}} with the analytic Jacobian; the
#' settings follow the Python \code{synergy} package, which an earlier version
#' of this function called.
#'
#' The whole grid is fitted, single-agent edges included, with three
#' additions to the \code{synergy} defaults: the untreated corner is set
#' to 1 (added when absent), \code{E} is bounded from -1 for GR and from 0 for
#' RV up to 1.5, and \code{C1}/\code{C2} stay within a factor of 1000 of the
#' tested doses.  \code{h}, \code{C} and \code{alpha} are fitted on a log
#' scale.  The fit starts from the \code{synergy} initial guess (Hill fits of
#' the two single-agent edges) and from \code{n_restarts} random initial
#' values; the fit with the best \code{r2}, computed on the measured points,
#' is kept, preferring fits that converged within the iteration limit.
#'
#' \code{alpha > 1} means synergistic potency, \code{beta > 0} synergistic
#' efficacy.  \code{alpha} is weakly constrained when a drug has little
#' effect, so extreme values of it are not a measure of synergy.
#'
#' @param avg_dt \code{data.table} for one (drug1 x drug2 x cell line x
#'   normalization_type) combination.  Must contain columns \code{Concentration},
#'   \code{Concentration_2}, \code{x}, and \code{normalization_type}.
#' @param n_restarts integer; number of restarts from random initial values.
#' @param seed integer; seed for the random initial values.  The global random
#'   number generator is left as it was.
#'
#' @return Named list with \code{normalization_type}, the parameters
#'   \code{E0}, \code{E1}, \code{E2}, \code{E3}, \code{h1}, \code{h2},
#'   \code{C1}, \code{C2}, \code{alpha12}, \code{alpha21}, the efficacy synergy
#'   \code{beta}, \code{r2} and \code{converged}.  Numeric fields are
#'   \code{NA_real_} when no fit succeeds.
#'
#' @examples
#' mae <- gDRutils::get_synthetic_data("finalMAE_combo_matrix_small")
#' combo_se <- mae[[gDRutils::get_supported_experiments("combo")]]
#' combo_se_out <- gDRutils::apply_fit(
#'   combo_se, musyc_fit_fn, "combination",
#'   output_assay = "musyc_params",
#'   summary_fn = musyc_summary_fn, summary_assay = "musyc_summary",
#'   fit_source = "musyc"
#' )
#'
#' @seealso \code{\link{musyc_summary_fn}}
#' @export
musyc_fit_fn <- function(avg_dt, n_restarts = 10L, seed = 1L) {
  checkmate::assert_data_table(avg_dt)
  checkmate::assert_count(n_restarts)
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

  e_lower <- if (norm_type == "RV") 0 else -1
  dose1 <- conc1[conc1 > 0]
  dose2 <- conc2[conc2 > 0]
  limits <- .musyc_fit_limits(.musyc_bounds(e_lower, dose1, dose2))
  corners <- .musyc_corners(obs, e_lower)
  # h1, h2, C1, C2, alpha12, alpha21
  log_low <- log10(c(0.1, 0.1, min(dose1), min(dose2), 0.01, 0.01))
  log_high <- log10(c(10, 10, max(dose1), max(dose2), 100, 100))
  random_p0 <- .with_local_seed(seed, lapply(seq_len(n_restarts), function(i) {
    c(corners, 10^stats::runif(length(log_low), log_low, log_high))
  }))

  starts <- c(list(.musyc_default_p0(obs)), random_p0)
  best <- .musyc_best_fit(lapply(starts, function(p0) .musyc_fit_once(obs, limits, p0)))

  if (is.null(best$params)) {
    return(na_result)
  }
  c(
    list(normalization_type = norm_type),
    best$params,
    list(beta = best$beta, r2 = best$r2, converged = best$converged)
  )
}

musyc_fit_fn <- gDRutils::set_fit_fn_info(
  musyc_fit_fn,
  model = "MuSyC (multi-dimensional synergy of combinations)",
  description = paste(
    "Fits the MuSyC dose-response surface to the whole grid by bounded least squares,",
    "keeps the best of several starts, and returns its parameters, the efficacy",
    "synergy beta and the R-squared of the fit."
  ),
  reference = .ref_musyc
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
#'   \code{log2(alpha12)} and \code{log2(alpha21)}, each clipped to
#'   \code{[-10, 10]}), \code{mean_r2},
#'   \code{any_converged} and \code{synergy_call} (\code{"synergistic"},
#'   \code{"antagonistic"}, \code{"additive"} or \code{"unknown"}, from
#'   \code{mean_beta} with a dead zone of 0.1).  \code{synergy_call} reflects
#'   efficacy only: a pair with strong potency synergy (\code{alpha > 1}) and
#'   \code{beta} near 0 is called \code{"additive"}, so read it together with
#'   \code{mean_log2_alpha}.
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
  # alpha is weakly constrained for a drug with little effect, so one extreme
  # value would dominate the mean
  log2_alpha <- log2(c(fit_dt[["alpha12"]][converged], fit_dt[["alpha21"]][converged]))
  mean_log2_alpha <- mean(pmin(pmax(log2_alpha, -10), 10), na.rm = TRUE)
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
    "Averages beta, clipped log2 alpha and R-squared of the MuSyC fits over the normalization types",
    "and calls the pair synergistic, antagonistic or additive from the mean beta."
  )
)



#' predict_musyc_surface
#'
#' Response of a fitted MuSyC surface (\code{gamma12 = gamma21 = 1}) at the
#' given concentrations, e.g. to compare a \code{\link{musyc_fit_fn}} fit
#' with the measured grid.
#'
#' @param conc,conc_2 numeric vectors of the same length; concentrations of
#'   drug 1 and drug 2.
#' @param params named list or one-row \code{data.table} with \code{E0},
#'   \code{E1}, \code{E2}, \code{E3}, \code{h1}, \code{h2}, \code{C1},
#'   \code{C2}, \code{alpha12} and \code{alpha21}, as returned by
#'   \code{musyc_fit_fn}.
#'
#' @return Numeric vector of responses, one per concentration pair.
#'
#' @examples
#' params <- list(E0 = 1, E1 = 0.4, E2 = 0.3, E3 = -0.2, h1 = 1.2, h2 = 0.9,
#'                C1 = 0.5, C2 = 0.2, alpha12 = 3, alpha21 = 0.5)
#' predict_musyc_surface(c(0, 0.1, 1), c(0, 1, 1), params)
#'
#' @seealso \code{\link{musyc_fit_fn}}
#' @export
predict_musyc_surface <- function(conc, conc_2, params) {
  checkmate::assert_numeric(conc, lower = 0)
  checkmate::assert_numeric(conc_2, lower = 0, len = length(conc))
  checkmate::assert_names(names(params), must.include = .musyc_param_names)
  values <- vapply(.musyc_param_names, function(name) as.numeric(params[[name]][1]), numeric(1))
  .musyc_surface(c(values[1:4], log(values[5:10])), conc, conc_2)
}


#' musyc_available
#'
#' Whether \code{\link{musyc_fit_fn}} can run.  Always \code{TRUE}: the fit
#' no longer needs Python.  Kept so that code written for the Python version
#' keeps working.
#'
#' @return Logical scalar, \code{TRUE}.
#'
#' @examples
#' musyc_available()
#'
#' @seealso \code{\link{musyc_fit_fn}}
#' @export
musyc_available <- function() {
  TRUE
}

.musyc_param_names <- c("E0", "E1", "E2", "E3", "h1", "h2", "C1", "C2", "alpha12", "alpha21")

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
.musyc_grid <- function(d1, d2, x) {
  obs <- data.table::data.table(d1 = d1, d2 = d2, E = x)
  obs <- obs[!is.na(obs$E), ]
  obs <- obs[!(obs$d1 == 0 & obs$d2 == 0), ]
  obs <- rbind(data.table::data.table(d1 = 0, d2 = 0, E = 1), obs)
  obs[order(obs$d1, obs$d2), ]
}

# E bounds from the normalization, C within a factor of 1000 of the tested
# doses; h and alpha keep the synergy defaults
.musyc_bounds <- function(e_lower, dose1, dose2) {
  list(
    E0_bounds = c(0.99, 1.01),
    E1_bounds = c(e_lower, 1.5),
    E2_bounds = c(e_lower, 1.5),
    E3_bounds = c(e_lower, 1.5),
    C1_bounds = c(min(dose1) / 1000, max(dose1) * 1000),
    C2_bounds = c(min(dose2) / 1000, max(dose2) * 1000)
  )
}

# The bounds on the fitting scale: E linear; h, C and alpha as logs, unbounded
# unless given
.musyc_fit_limits <- function(bounds) {
  one <- function(name, side) {
    b <- bounds[[paste0(name, "_bounds")]]
    if (is.null(b)) {
      return(if (side == 1L) -Inf else Inf)
    }
    if (startsWith(name, "E")) b[side] else log(b[side])
  }
  list(
    lower = vapply(.musyc_param_names, one, numeric(1), side = 1L),
    upper = vapply(.musyc_param_names, one, numeric(1), side = 2L)
  )
}

# Initial E0..E3 from the grid corners, clipped to [e_lower, 1]; a missing
# corner falls back to the median response
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

# (alpha * d / C)^h on the log scale, 0 at dose 0
.musyc_ratio_pow <- function(log_d, log_c, h, log_alpha = 0) {
  out <- exp(pmin(pmax(h * (log_d + log_alpha - log_c), -200), 200))
  out[!is.finite(log_d)] <- 0
  out
}

# MuSyC surface with gamma = 1 at doses d1, d2 for theta = (E0, E1, E2, E3,
# log h1, log h2, log C1, log C2, log alpha12, log alpha21), and optionally its
# Jacobian with respect to theta
.musyc_surface <- function(theta, d1, d2, jacobian = FALSE) {
  e0 <- theta[1]
  e1 <- theta[2]
  e2 <- theta[3]
  e3 <- theta[4]
  h1 <- exp(theta[5])
  h2 <- exp(theta[6])
  log_d1 <- suppressWarnings(log(d1))
  log_d2 <- suppressWarnings(log(d2))

  x1 <- .musyc_ratio_pow(log_d1, theta[7], h1)
  x2 <- .musyc_ratio_pow(log_d2, theta[8], h2)
  y1 <- .musyc_ratio_pow(log_d1, theta[7], h1, theta[10])
  y2 <- .musyc_ratio_pow(log_d2, theta[8], h2, theta[9])

  u_num <- 2 + y1 + y2
  a1_num <- x1 * (2 + y1) + x2 * y1
  a2_num <- x1 * y2 + x2 * (2 + y2)
  den <- u_num + x1 * (2 + y1) * (1 + y2) + x2 * (2 + y2) * (1 + y1)
  num <- u_num * (e0 - e3) + a1_num * (e1 - e3) + a2_num * (e2 - e3)
  if (!jacobian) {
    return(e3 + num / den)
  }

  d_e <- function(d_num, d_den) (d_num * den - num * d_den) / den^2
  d_x1 <- d_e((2 + y1) * (e1 - e3) + y2 * (e2 - e3), (2 + y1) * (1 + y2))
  d_x2 <- d_e(y1 * (e1 - e3) + (2 + y2) * (e2 - e3), (2 + y2) * (1 + y1))
  d_y1 <- d_e((e0 - e3) + (x1 + x2) * (e1 - e3), 1 + x1 * (1 + y2) + x2 * (2 + y2))
  d_y2 <- d_e((e0 - e3) + (x1 + x2) * (e2 - e3), 1 + x1 * (2 + y1) + x2 * (1 + y1))
  # d (ratio power) / d log h, 0 at dose 0
  d_log_h <- function(value, h, log_d, log_ratio) ifelse(is.finite(log_d), value * h * log_ratio, 0)

  cbind(
    u_num / den,
    a1_num / den,
    a2_num / den,
    1 - (u_num + a1_num + a2_num) / den,
    d_x1 * d_log_h(x1, h1, log_d1, log_d1 - theta[7]) +
      d_y1 * d_log_h(y1, h1, log_d1, log_d1 + theta[10] - theta[7]),
    d_x2 * d_log_h(x2, h2, log_d2, log_d2 - theta[8]) +
      d_y2 * d_log_h(y2, h2, log_d2, log_d2 + theta[9] - theta[8]),
    -h1 * (d_x1 * x1 + d_y1 * y1),
    -h2 * (d_x2 * x2 + d_y2 * y2),
    h2 * d_y2 * y2,
    h1 * d_y1 * y1
  )
}

# Bounded least squares on residual function f with Jacobian jac, by nlminb
# with the gradient and the Gauss-Newton Hessian. Returns the best point
# evaluated: on singular or false convergence nlminb's par can be worse.
.musyc_least_squares <- function(f, jac, start, lower, upper) {
  finite_or_zero <- function(m) {
    m[!is.finite(m)] <- 0
    m
  }
  best <- list(par = start, objective = Inf)
  objective <- function(p) {
    r <- f(p)
    value <- if (all(is.finite(r))) sum(r^2) else 1e30
    if (value < best$objective) {
      best <<- list(par = p, objective = value)
    }
    value
  }
  gradient <- function(p) finite_or_zero(2 * crossprod(jac(p), f(p))[, 1])
  hessian <- function(p) finite_or_zero(2 * crossprod(jac(p)))
  fit <- stats::nlminb(start, objective, gradient = gradient, hessian = hessian,
                       lower = lower, upper = upper,
                       control = list(iter.max = 1000L, eval.max = 1500L))
  # like synergy, the fit failed to converge only when it ran out of iterations
  c(best, list(converged = !grepl("limit reached", fit$message, fixed = TRUE)))
}

# Four-parameter Hill fit of a single-agent edge as synergy does it: E0, Emax,
# log h and log C on doses scaled by exp(median(log d)), starting from E at the
# lowest and highest dose, h = 1 and C = median(d)
.musyc_hill_fit <- function(d, E) {
  scale <- exp(stats::median(log(d[d > 0])))
  log_d <- suppressWarnings(log(d / scale))
  start <- c(stats::median(E[d == min(d)]), stats::median(E[d == max(d)]), 0,
             log(stats::median(d / scale)))
  if (!all(is.finite(start))) {
    return(NULL)
  }
  model <- function(p) {
    xh <- .musyc_ratio_pow(log_d, p[4], exp(p[3]))
    p[1] + (p[2] - p[1]) * xh / (1 + xh)
  }
  jac <- function(p) {
    h <- exp(p[3])
    xh <- .musyc_ratio_pow(log_d, p[4], h)
    frac <- xh / (1 + xh)
    d_frac <- (p[2] - p[1]) * frac / (1 + xh)
    cbind(1 - frac, frac, d_frac * ifelse(is.finite(log_d), h * (log_d - p[4]), 0), -d_frac * h)
  }
  fit <- .musyc_least_squares(function(p) model(p) - E, jac, start, rep(-Inf, 4), rep(Inf, 4))
  if (!all(is.finite(fit$par))) {
    return(NULL)
  }
  list(E0 = fit$par[1], Emax = fit$par[2], h = exp(fit$par[3]), C = exp(fit$par[4]) * scale)
}

# synergy's initial guess: E, h and C from Hill fits of the two single-agent
# edges, E3 at the highest doses, alpha = 1
.musyc_default_p0 <- function(obs) {
  edge1 <- obs$d2 == min(obs$d2)
  edge2 <- obs$d1 == min(obs$d1)
  hill1 <- .musyc_hill_fit(obs$d1[edge1], obs$E[edge1])
  hill2 <- .musyc_hill_fit(obs$d2[edge2], obs$E[edge2])
  if (is.null(hill1) || is.null(hill2)) {
    return(NULL)
  }
  e3 <- obs$E[obs$d1 == max(obs$d1) & obs$d2 == max(obs$d2)]
  e3 <- if (length(e3)) stats::median(e3) else min(obs$E)
  c((hill1$E0 + hill2$E0) / 2, hill1$Emax, hill2$Emax, e3,
    hill1$h, hill2$h, hill1$C, hill2$C, 1, 1)
}

.musyc_fit_once <- function(obs, limits, p0) {
  failed <- list(params = NULL, beta = NA_real_, r2 = NA_real_, converged = FALSE)
  if (is.null(p0)) {
    return(failed)
  }
  start <- pmin(pmax(c(p0[1:4], log(p0[5:10])), limits$lower), limits$upper)
  fit <- tryCatch(
    .musyc_least_squares(
      function(theta) .musyc_surface(theta, obs$d1, obs$d2) - obs$E,
      function(theta) .musyc_surface(theta, obs$d1, obs$d2, jacobian = TRUE),
      start, limits$lower, limits$upper
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || !all(is.finite(fit$par))) {
    return(failed)
  }

  theta <- fit$par
  params <- stats::setNames(as.list(c(theta[1:4], exp(theta[5:10]))), .musyc_param_names)
  # r2 on the measured points; the untreated corner is set, not measured
  measured <- !(obs$d1 == 0 & obs$d2 == 0)
  predicted <- .musyc_surface(theta, obs$d1[measured], obs$d2[measured])
  ss_res <- sum((obs$E[measured] - predicted)^2)
  ss_tot <- sum((obs$E[measured] - mean(obs$E[measured]))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  if (!is.finite(r2)) {
    return(failed)
  }
  strongest <- min(params$E1, params$E2)
  list(
    params = params,
    beta = (strongest - params$E3) / (params$E0 - strongest),
    r2 = r2,
    converged = fit$converged
  )
}

# The fit with the best r2 among the converged ones, as in synergy, where a
# start that runs out of iterations fails; among all fits when none converged
.musyc_best_fit <- function(fits) {
  r2 <- vapply(fits, function(fit) fit$r2, numeric(1))
  converged <- vapply(fits, function(fit) isTRUE(fit$converged), logical(1))
  pool <- which(!is.na(r2) & converged)
  if (length(pool) == 0L) {
    pool <- which(!is.na(r2))
  }
  if (length(pool) == 0L) fits[[1]] else fits[[pool[which.max(r2[pool])]]]
}

# Evaluate expr with the given seed and restore the global RNG state afterwards
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
