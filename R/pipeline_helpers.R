create_SE_step <- function(inl, 
                           exp,
                           data_type, 
                           readout, 
                           control_mean_fxn, 
                           nested_identifiers,
                           override_untrt_controls, 
                           add_raw_data) {
  df <- inl$df_list[[exp]]
  if (is.null(df)) {
 
    msg1 <- sprintf(
      "It's impossible to run pipeline from the first 
      step ('%s') for experiment: '%s'. ",
      get_pipeline_steps()[1], exp
    )
    msg2 <- sprintf(
      "Consider running the pipeline from the second ('%s') step", 
      get_pipeline_steps()[2]
    )
    stop(c(msg1, msg2))

  }

  se <- create_SE(
    df_ = df,
    data_type = data_type,
    readout = readout,
    nested_identifiers = nested_identifiers,
    nested_confounders = inl$nested_confounders,
    override_untrt_controls = override_untrt_controls
  )
  invisible(se)
}

run_pipeline_step <- function(run_vars, 
                              step, 
                              se,
                              step_fun, 
                              step_args, 
                              if_paste_warnings = FALSE,
                              if_read_intermediate_data = TRUE) {
  if (!run_vars$partial_run || !do_skip_step(step, run_vars$start_from)) {
    se <- purrr::quietly(do.call)(step_fun, step_args)

    if (!is.null(run_vars$data_dir)) {
      save_intermediate_data(run_vars$data_dir, step, run_vars$exp, se$result)
    }
    if (if_paste_warnings) {
      paste_warnings(se$warnings)
    }
  } else {
    if (if_read_intermediate_data && 
        is_preceding_step(step, run_vars$start_from)) {
      se$result <- read_intermediate_data(run_vars$data_dir, step, run_vars$exp)
    }
  }

  invisible(se)
}

#' get pipeline steps
#' 
#' @return vector with steps
#' 
#' @keywords internal
get_pipeline_steps <-
  function() {
    c("create_SE",
      "normalize_SE",
      "average_SE",
      "fit_SE")
  }

#' check if the given step can be skipped if partial run is chosen
#' 
#' @param current_step, string with the step to be evaluated
#' @param start_from string indicating the pipeline step from which partial 
#' run should be launched
#' @param steps charvect with all available steps
#' 
#' @keywords internal
#' @return logical
#' 
do_skip_step <-
  function(current_step, start_from, steps = get_pipeline_steps()) {
    
    checkmate::assert_choice(current_step, steps)
    checkmate::assert_choice(start_from, steps, null.ok = TRUE)
    checkmate::assert_character(steps)
    
    if (is.null(start_from)) {
      FALSE
    } else {
      c_idx <- which(steps %in% current_step)
      s_idx <- which(steps %in% start_from)
      c_idx < s_idx
    }
  }

#' check if the given step is preceding the step chosen in the partial run
#' 
#' @param current_step, string with the step to be evaluated
#' @param start_from string indicating the pipeline step from which 
#' partial run should be launched
#' @param steps charvect with all available steps
#' 
#' @keywords internal
#' @return logical
#' 
is_preceding_step <-
  function(current_step, start_from, steps = get_pipeline_steps()) {
   
    checkmate::assert_choice(current_step, steps)
    checkmate::assert_choice(start_from, steps)
    checkmate::assert_character(steps)
    
    c_idx <- which(steps %in% current_step)
    s_idx <- which(steps %in% start_from)
    s_idx - c_idx == 1
  }

#' save intermediate data for the given experiment and step to qs file
#' 
#' @param path string with the save directory for the qs file 
#' @param step, string with the step name
#' @param experiment string with the experiment name
#' @param se output se 
#' 
#' @keywords internal
#' 
#' @return \code{NULL}
#' 
save_intermediate_data <- function(path, step, experiment, se) {
  
  checkmate::assert_directory(path, "rw")
  checkmate::assert_string(step)
  checkmate::assert_string(experiment)
  
  fpath <- file.path(path, paste0(experiment, "__", step, ".qs"))
  qs::qsave(se, fpath)
}

#' read intermediate data for the given experiment and step to qs file
#' 
#' @param path string with the input directory of the qs file 
#' @param step, string with the step name
#' @param experiment string with the experiment name
#' 
#' @keywords internal
#' @return se
#' 
read_intermediate_data <- function(path, step, experiment) {
  
  checkmate::assert_directory(path, "r")
  checkmate::assert_string(step)
  checkmate::assert_string(experiment)
  
  fpath <- file.path(path, paste0(experiment, "__", step, ".qs"))
  qs::qread(fpath)
}

#' @keywords internal
paste_warnings <- function(list, sep = "\n") {
  pasted <- paste0(list, sep = sep)
  warning(pasted, call. = FALSE)
}

#' @keywords internal
.clear_rownames <- function(x) {
  lapply(x, function(y) {
    if (inherits(y, c("DFrame", "data.frame"))) {
      y <- y[do.call(order, y), ]
    }
    rownames(y) <- NULL
    y
  })
}
