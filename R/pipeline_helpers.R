create_SE_step <- function(inl, 
                           exp,
                           data_type, 
                           readout, 
                           control_mean_fxn, 
                           nested_identifiers,
                           override_untrt_controls) {
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

#' @keywords internal
prepareData <- function(cell_lines, drugs, conc = 10 ^ (seq(-3, 1, 0.5))) {
  df_layout <- merge(cell_lines, drugs, by = NULL)
  df_layout <- gDRtestData::add_data_replicates(df_layout)
  gDRtestData::add_concentration(df_layout, conc)
}

#' @keywords internal
prepareMergedData <- function(cell_lines, drugs, noise = 0.1) {
  df <- prepareData(cell_lines, drugs)
  gDRtestData::generate_response_data(df, noise)
}

#' @keywords internal
prepareComboMergedData <- function(cell_lines, 
                                   drugs, 
                                   drugsIdx1 = 2:4,
                                   drugsIdx2 = c(26, 26, 26), 
                                   concentration = c(0, .2, 1), 
                                   noise = 0.1, 
                                   modifyDf2 = FALSE) {
  df_layout <- prepareData(cell_lines, drugs[drugsIdx1, ])
  
  df_2 <- cbind(drugs[drugsIdx2, ], Concentration = concentration)
  colnames(df_2) <- paste0(colnames(df_2), "_2")
  df_layout_2 <- merge(df_layout, df_2, by = NULL)
  if (modifyDf2) {
    df_layout_2 <- df_layout_2[!(df_layout_2$Concentration == 0 & df_layout_2$Concentration_2 > 0), ]
  }
  
  gDRtestData::generate_response_data(df_layout_2, noise)
}

#' @keywords internal
prepareCodilutionData <- function(df, df_layout) {
  colnames(df) <- paste0(colnames(df), "_2")
  df_2 <- cbind(df_layout, df)
  df_2 <- df_2[df_2$DrugName != df_2$DrugName_2, ]
  rows <- df_2$Concentration_2 > 0
  cols <- c("Concentration", "Concentration_2")
  df_2[rows, cols] <- df_2[rows, cols] / 2
  
  df_2
}

#' @keywords internal
changeColNames <- function(df, drugs, suffix) {
  cols <- colnames(df) %in% c(colnames(drugs), "Concentration")
  colnames(df)[cols] <- paste0(colnames(df)[cols], suffix)
  
  df
}

#' @keywords internal
save_rds <- function(rdsObj, rdsName) {
  saveRDS(
    rdsObj,
    file.path(system.file("testdata", package = "gDRtestData"), rdsName), 
    compress = "gzip"
  )
}
