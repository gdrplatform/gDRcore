#' Deprecated function(s) in \pkg{gDR}
#' Overall_function
#' 
#' Process DR data using set of manifest, template and results files
#' 
#' @param manifest_file a file path for the manifest file
#' @param template_file a file path(s) for the template file(s)
#' @param results_file a file path(s) for the results file(s)
#' @param output_files 
#' @param selected_keys a vector of keys that should be included in the normalization
#' @param key_values a list of key values that should be included in the normalization
#' @param instrument character
#' @return a list with the raw, normalized, averaged, and metrics dataframes
#' 
#' @examples
#' \dontrun{
#' Overall_function(ManifestFile, TemplateFile, ResultsFile, OutputFile)
#' }
#' @import gneDB
#' @import reshape2
#' @import dplyr
Overall_function <-
  function(manifest_file,
           template_file,
           results_file,
           output_files, # TODO: investigate what that is
           selected_keys = NULL,
           key_values = NULL,
           instrument = "EnVision") {
    # output_files should contain file names for :
    #   QC_file, raw_result, process_results, metrics_results
    .Deprecated(msg = "The Overall_function is deprecated. Use separate functions to processing DR data")
    futile.logger::flog.info("Report from gDR pipeline")
    # Assertions:
    checkmate::assert_string(manifest_file)
    checkmate::assert_character(template_file)
    checkmate::assert_character(results_file)
    checkmate::assert_string(output_files)
    checkmate::assert_vector(selected_keys, null.ok = TRUE)
    checkmate::assert_list(key_values, null.ok = TRUE)
    checkmate::assert_string(instrument)
    
    
    lData <- load_data(manifest_file,
                       template_file,
                       results_file,
                       instrument)
    df_raw_data <-
      merge_data(lData$manifest, lData$treatments, lData$data)
    
    #returning two SEs with 'df_raw_data' assay and creating MAE from it
    untreatedSE <-
      gDR::createSE(df_raw_data, data_type = "untreated")
    treatedSE <-
      gDR::createSE(df_raw_data, data_type = "treated")
    rawMAE <- gDR::createMAE_raw(untreatedSE, treatedSE)
    
    # output_QC_byPlate(df_raw_data, output_files['QC_file']) # TODO: check column/row bias
    
    # TODO: this function could be overload to have an MAE as input; ok for now
    Keys <- identify_keys(df_raw_data) # may be manually changed
    if (!is.null(selected_keys)) {
      Keys[names(selected_keys)] <- selected_keys[names(selected_keys)]
    }
    
    df_normalized <-
      normalize_data(df_raw_data, Keys, key_values)
    df_averaged <- average_replicates(df_normalized, Keys$Trt)
    
    df_metrics <- calculate_DRmetrics(df_averaged, Keys$DoseResp)
    
    return(
      list(
        raw = df_raw_data,
        normalized = df_normalized,
        averaged = df_averaged,
        metrics = df_metrics
      )
    )
  }

#' normalize_data
#' 
#' Normalize raw DR data
#'
#' @param df_raw_data a dataframe with raw data
#' @param selected_keys a vector of keys that should be included in the normalization (NULL by default)
#' @param key_values a list of values for keys that should be included in the normalization (NULL by default)
#'
#' @return a dataframe with normalization
#' @export
#'
#' @examples
#' \dontrun{
#' normalize_data(df_raw_data)
#' }
normalize_data <-
  function(df_raw_data,
           selected_keys = NULL,
           key_values = NULL) {
    .Deprecated(msg = "The normalize_data is deprecated. Use normalize_SE function with SE object to normalize DR data")
    # Assertions
    stopifnot(inherits(df_raw_data, "data.frame"))
    checkmate::assert_vector(selected_keys, null.ok = TRUE)
    checkmate::assert_list(key_values, null.ok = TRUE)
    
    # average technical replicates and assign the right controls to each treated well
    
    # remove unused columns but keep barcodes to normalize by plate
    df_normalized <- df_raw_data[, setdiff(colnames(df_raw_data),
                                           c("Template", gDRutils::get_identifier("WellPosition")))]
    
    # Identify keys for assigning the controls
    Keys <- identify_keys(df_normalized)
    if (!is.null(selected_keys)) {
      Keys[names(selected_keys)] <- selected_keys[names(selected_keys)]
    }
    
    df_normalized$CorrectedReadout <-
      pmax(df_normalized$ReadoutValue -
             df_normalized$BackgroundValue,
           1)
    
    # enforced key values for end points (override selected_keys)
    Keys$untrt_Endpoint <- setdiff(Keys$untrt_Endpoint, names(key_values))
    endpoint_value_filter <- array(TRUE, dim(df_raw_data)[1])
    if (!is.null(key_values) & length(key_values) > 0) {
      for (i in 1:length(key_values)) {
        if (is.numeric(key_values[i])) {
          endpoint_value_filter <- endpoint_value_filter &
            (df_normalized[, names(key_values)[i] ] == key_values[i] &
               !is.na(df_normalized[, names(key_values)[i] ]))
        } else {
          endpoint_value_filter <- endpoint_value_filter &
            (df_normalized[, names(key_values)[i] ] %in% key_values[i])
        }}}
    # get the untreated controls at endpoint and perform interquartile mean
    df_end_untrt <- df_normalized[df_normalized[, gDRutils::get_identifier("duration")] > 0 & endpoint_value_filter &
                                    apply(df_normalized[, agrep("Concentration", colnames(df_normalized)), drop = FALSE] == 0, 1, all), ]
    df_end_mean <- aggregate(df_end_untrt[, "CorrectedReadout"],
                             by = as.list(df_end_untrt[, Keys$untrt_Endpoint]), function(x) mean(x, trim = .25))
    colnames(df_end_mean)[dim(df_end_mean)[2]] <- "UntrtReadout"
    
    # get the untreated controls at Day 0 and perform interquartile mean
    df_day0 <- df_normalized[df_normalized[,gDRutils::get_identifier("duration")] == 0 &
                               apply(df_normalized[, agrep("Concentration", colnames(df_normalized)), drop = FALSE] == 0, 1, all), ]
    df_day0_mean <- aggregate(df_day0[,"CorrectedReadout"],
                              by = as.list(df_day0[,Keys$Day0]), function(x) mean(x, trim = .25))
    colnames(df_day0_mean)[dim(df_day0_mean)[2]] <- "Day0Readout"
    
    df_controls <- merge(df_end_mean, df_day0_mean[, setdiff(colnames(df_day0_mean),
                                                             c(gDRutils::get_identifier("duration"), "Barcode"))], all.x = TRUE)
    if (length(setdiff(Keys$untrt_Endpoint, Keys$Day0)) > 0) {
      futile.logger::flog.warn("Not all control conditions found on the day 0 plate, dispatching values for field: %s",
                               paste(setdiff(Keys$untrt_Endpoint, Keys$Day0), collapse = " ; "))
    }
    # identify missing values in the Day0 that needs to be matched (usually for co-treatments)
    df_controls_NA <- which(is.na(df_controls$Day0Readout))
    
    if (length(df_controls_NA) > 0) {
      dispatched <- NULL
      for (i in df_controls_NA) {
        matches <- t(apply(df_day0_mean[, setdiff(Keys$Day0,
                                                  c(gDRutils::get_identifier("duration"), "Barcode"))], 1,
                           function(x)
                             df_controls[i, setdiff(Keys$Day0, c(gDRutils::get_identifier("duration"), "Barcode")),
                                         drop = FALSE] == x))
        colnames(matches) <-
          setdiff(Keys$Day0, c(gDRutils::get_identifier("duration"), "Barcode"))
        # try to find a good match for the day 0 (enforce same cell line)
        idx <-
          rowSums(matches) * matches[, gDRutils::get_identifier("cellline")]
        if (all(idx == 0)) {
          next
        }
        match_idx <- which.max(idx)
        mismatch <- df_day0_mean[match_idx, setdiff(Keys$Day0,
                                                    c(gDRutils::get_identifier("duration"), "Barcode"))] !=
          df_controls[i, setdiff(Keys$Day0, c(gDRutils::get_identifier("duration"), "Barcode"))]
        dispatched <-
          c(dispatched, colnames(mismatch)[mismatch])
        df_controls[i, "Day0Readout"] <-
          df_day0_mean[match_idx, "Day0Readout"]
      }
      futile.logger::flog.warn("Not all control conditions found on the day 0 plate")
      if(length(dispatched) > 0) {
        futile.logger::flog.warn(
          "dispatching values for mismatches in field: ",
          paste(unique(dispatched), collapse = " ; ")
        )
      } else {
        futile.logger::flog.warn("some Day0 are not being matched")
      }
    }
    
    df_to_norm <-
      df_normalized[df_normalized[, gDRutils::get_identifier("duration")] > 0 &
                      (apply(df_normalized[, agrep("Concentration", colnames(df_normalized)), drop = FALSE] != 0, 1, any) |
                         !endpoint_value_filter),]
    
    df_to_norm_conditions <-
      unique(df_to_norm[, intersect(colnames(df_to_norm),
                                    colnames(df_controls))])
    
    # if missing barcodes --> dispatch for similar conditions
    if (!all(df_to_norm_conditions$Barcode %in% df_controls$Barcode)) {
      futile.logger::flog.warn("Not all control conditions found at the end of treatment,
                               dispatching values for plates: %s",
                               paste(
                                 setdiff(df_to_norm_conditions$Barcode, df_controls$Barcode),
                                 collapse = " ; "
                               ))
      
      df_ctrl_mean <-
        aggregate(df_controls[, c("UntrtReadout", "Day0Readout")],
                  by <- as.list(subset(
                    df_controls,
                    select = -c(UntrtReadout, Day0Readout, Barcode)
                  )), mean)
      df_controls <- rbind(df_controls, merge(df_ctrl_mean,
                                              df_to_norm_conditions[!(df_to_norm_conditions$Barcode %in%
                                                                        df_controls$Barcode), ]))
    }
    
    df_normalized <- merge(df_to_norm, df_controls)
    
    df_normalized$RelativeViability <-
      round(df_normalized$CorrectedReadout /
              df_normalized$UntrtReadout,
            4)
    df_normalized$GRvalue <- round(2 ** (
      log2(df_normalized$CorrectedReadout / df_normalized$Day0Readout) /
        log2(df_normalized$UntrtReadout / df_normalized$Day0Readout)
    ), 4) - 1
    
    df_normalized$DivisionTime <-
      round(df_normalized[, gDRutils::get_identifier("duration")] /
              log2(df_normalized$UntrtReadout / df_normalized$Day0Readout), 4)
    
    
    if (any(is.na(df_normalized$Day0Readout))) {
      # need to use the reference doubling Time if day 0 missing
      InferedIdx <- is.na(df_normalized$Day0Readout)
      filtered <-
        df_normalized$ReferenceDivisionTime > (df_normalized[, gDRutils::get_identifier("duration")] * 2) |
        is.na(df_normalized$ReferenceDivisionTime)
      futile.logger::flog.warn(
        "Missing day 0 information --> calculate GR value based on reference doubling time")
      
      if (any(filtered & InferedIdx)) {
        futile.logger::flog.warn("Filtering %s conditions because of too short assay: %s",
                                 paste(unique(df_normalized$CellLineName[filtered &
                                                                           InferedIdx]), collpase = " ; ")
        )
      }
      
      InferedIdx <- !filtered & InferedIdx
      # calculate GR values using formula from https://www.nature.com/articles/nbt.3882
      df_normalized$GRvalue[InferedIdx] <-
        round(2 ^ (1 + (
          log2(pmin(1.25,
                    df_normalized[InferedIdx, "RelativeViability"])) /
            (df_normalized[, gDRutils::get_identifier("duration")][InferedIdx] /
               df_normalized$ReferenceDivisionTime[InferedIdx])
        )), 4) - 1
    }
    
    df_normalized <-
      cbind(df_normalized[, 1:(which(colnames(df_normalized) == "ReadoutValue") - 1)],
            df_normalized[, c("GRvalue", "RelativeViability", "DivisionTime")],
            df_normalized[, which(colnames(df_normalized) == "ReadoutValue"):(dim(df_normalized)[2] - 3)])
    df_normalized <- Order_result_df(df_normalized)
    futile.logger::flog.info("df normalized")
    return(df_normalized)
  }

#' average_replicates
#'
#' @param df_normalized a dataframe with normalized values of DR data
#' @param TrtKeys a vector of keys used for averaging (NULL by default)
#'
#' @return a dataframe with averaged replicated of DR data
#' @export
#'
#' @examples
#' \dontrun{
#' average_replicates(df_normalized)
#' }
average_replicates <- function(df_normalized, TrtKeys = NULL) {
  .Deprecated(msg = "The average_replicates is deprecated. Use average_SE function with SE object to average DR data")
  # Assertions
  stopifnot(inherits(df_normalized, "data.frame"))
  checkmate::assert_vector(TrtKeys, null.ok = TRUE)
  
  if (is.null(TrtKeys)) {
    TrtKeys <- identify_keys(df_normalized)$Trt
  }
  df_averaged <-
    aggregate(
      df_normalized[, c(
        "GRvalue",
        "RelativeViability",
        "CorrectedReadout",
        "UntrtReadout",
        "Day0Readout",
        "DivisionTime",
        "ReferenceDivisionTime"
      )],
      by = as.list(df_normalized[, TrtKeys]),
      FUN = function(x)
        mean(x, na.rm = TRUE)
    )
  df_std <-
    aggregate(
      df_normalized[, c("GRvalue", "RelativeViability")],
      by = as.list(df_normalized[, TrtKeys]),
      FUN = function(x)
        sd(x, na.rm = TRUE)
    )
  colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")] <-
    paste0("std_",
           colnames(df_std)[colnames(df_std) %in% c("GRvalue", "RelativeViability")])
  df_averaged <- merge(df_averaged, df_std, by = TrtKeys)
  
  #reorganize column order:
  df_averaged <- Order_result_df(df_averaged)
  
  futile.logger::flog.info("df averaged")
  return(df_averaged)
}

#' calculate_DRmetrics
#' 
#' Calculate metrics for DR data
#'
#' @param df_averaged a SummarizedExperiment with averaged and normalized assays
#' @param DoseRespKeys a vector of dose response keys used for calculation of metrics (NULL by default)
#' @param studyConcThresh a numeric with study concentration threshold (4 by default)
#'
#' @return a dataframe with with metrics
#'
#' @importFrom dplyr arrange_at group_by_at left_join summarise
#' @export 
#' @examples
#' \dontrun{
#' calculate_DRmetrics(df_averaged)
#' }
calculate_DRmetrics <-
  function(df_averaged,
           DoseRespKeys = NULL,
           studyConcThresh = 4) {
    .Deprecated(msg = "The calculate_DRmetrics is deprecated. Use metrics_SE function with SE object to calculate metrics of DR data")
    
    # Assertions
    stopifnot(inherits(df_averaged, "data.frame"))
    checkmate::assert_vector(DoseRespKeys, null.ok = TRUE)
    checkmate::assert_number(studyConcThresh)
    
    df_a <- df_averaged
    colnames(df_a)[colnames(df_a) == gDRutils::get_identifier("drugname")] <-
      "DrugName"
    
    if (is.null(DoseRespKeys)) {
      DoseRespKeys <- identify_keys(df_a)$DoseResp
    } else {
      DoseRespKeys [DoseRespKeys == gDRutils::get_identifier("drugname")] <-
        "DrugName"
    }
    DoseRespKeys <- setdiff(DoseRespKeys, "Concentration")
    DoseRespKeys <- c(DoseRespKeys, "DivisionTime")
    DoseRespKeys <- intersect(DoseRespKeys, colnames(df_a))
    
    df_a$log10Concentration <- log10(df_a$Concentration)
    
    metrics <-
      names(ICGRlogisticFit(c(-7, -6, -5, -4), c(1, .9, .8, .7), c(1, .9, .8, .7)))
    # dummy call to get variable names
    
    # define set of key for merging control and study data
    mergeKeys <-
      setdiff(DoseRespKeys, c(gDRutils::get_identifier("drug"), "DrugName"))
    
    # get avereage GRvalue ("GR_0") for control data
    controlSets <-
      df_a %>%
      dplyr::filter(DrugName %in% gDRutils::get_identifier("untreated_tag")) %>%
      dplyr::group_by_at(mergeKeys) %>%
      dplyr::summarise(GR_0 = mean(GRvalue),
                       e_0 = mean(RelativeViability))
    
    # get study data
    studySets <-
      df_a %>% dplyr::filter(!DrugName %in% gDRutils::get_identifier("untreated_tag"))
    
    # join study and control data
    # i.e. get  reference (average control) GRvalue ("GR_0") for study data
    fSets <-
      dplyr::left_join(studySets, controlSets, by = mergeKeys)
    # for study sets with no reference GRvalue, assing GRValue0 to 1
    fSets[is.na(fSets$GR_0), "GR_0"] <- 1
    fSets[is.na(fSets$e_0), "e_0"] <- 1
    
    #group study data by "DoseRespKeys"
    gSets <-
      fSets %>% dplyr::group_by_at(DoseRespKeys) %>% group_split()
    
    # filter to have at least 4 records with non-NA RelativeViability
    gSets <-
      gSets[lapply(gSets, function(x)
        sum(!is.na(x$RelativeViability))) >= studyConcThresh]
    
    futile.logger::flog.info(
      "Metadata variables for dose response curves: %s (%d groups)",
      paste(setdiff(
        DoseRespKeys, c(
          gDRutils::get_identifier("drug"),
          gDRutils::get_identifier("cellline"),
          paste(gDRutils::get_identifier("drug"), "_", 2:10)
        )
      ),
      collapse = " "),
      length(gSets)
    )
    
    
    #iterate over study groups
    resL <- lapply(1:length(gSets), function(x) {
      # the 'DoseRespKeys' columns in given grup are identical for each entry
      # let's get the first record then
      repCols <- as.vector(gSets[[x]][1, DoseRespKeys])
      #get selected columns ("metrics") from GRlogisticFit output
      # (if at least 4 records with non-NA RelativeViability)
      if (sum(!is.na(gSets[[x]]$RelativeViability)) >= studyConcThresh) {
        grLogCols <-
          ICGRlogisticFit(
            gSets[[x]]$log10Concentration,
            gSets[[x]]$RelativeViability,
            gSets[[x]]$GRvalue,
            e_0 = gSets[[x]]$e_0[1],
            GR_0 = gSets[[x]]$GR_0[1]
          )[metrics]
      } else {
        grLogCols <- rep(NA, length(metrics))
        names(grLogCols) <- metrics
        grLogCols$N_conc <-
          sum(!is.na(gSets[[x]]$RelativeViability))
      }
      cbind(repCols, t(grLogCols))
    })
    
    #return final data.frame
    resDf <- do.call(rbind, resL)
    resDf <- resDf [resDf$N_conc >= studyConcThresh,]
    resDf <- resDf %>% dplyr::arrange_at(DoseRespKeys)
    colnames(resDf)[colnames(resDf) == "DrugName"] <-
      gDRutils::get_identifier("drugname")
    resDf <- Order_result_df(resDf)
    return(resDf)
  }
