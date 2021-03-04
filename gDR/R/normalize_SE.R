#' normalize_SE
#'
#' Normalize raw DR data
#'
#' @param df_raw_data a dataframe with raw data
#' @param selected_keys a vector of keys that should be included in the normalization (NULL by default)
#' @param key_values a list of values for keys that should be included in the normalization (NULL by default)
#' @param discard_keys a vector of keys that should be discarded (NULL by default)
#' @param control_mean_fct the function used for averaging controls (trimmed arithmetic mean with trim = .25 by default)
#' @param ndigit_rounding an integer with number of digits of rounding (4 by default)
#'
#' @return a SummarizedExperiment object with normalized assays
#' @export
#'
normalize_SE <- function(df_raw_data,
                  selected_keys = NULL,
                  key_values = NULL,
                  discard_keys = NULL,
                  control_mean_fct = function(x) mean(x, trim = .25), # used for averaging controls
                  ndigit_rounding = 4 # rounding of normalized response values
                ) {
    .Deprecated(msg = "see normalize_SE2 for similar, but not identical functionality")

    # Assertions
    stopifnot(inherits(df_raw_data, "data.frame"))
    checkmate::assert_vector(selected_keys, null.ok = TRUE)
    checkmate::assert_list(key_values, null.ok = TRUE)
    checkmate::assert_vector(discard_keys, null.ok = TRUE)
    checkmate::assert_function(control_mean_fct, null.ok = TRUE)
    checkmate::assert_number(ndigit_rounding)
    checkmate::assertFALSE('RelativeViability' %in% colnames(df_raw_data))

    # average technical replicates and assign the right controls to each treated well
    Keys <- identify_keys(df_raw_data)
    if (!is.null(key_values)) {
      checkmate::assert_true(all(names(key_values) %in% unlist(Keys)))
    } 
    Keys$discard_keys <- discard_keys
    if (!is.null(selected_keys)) {
      Keys[names(selected_keys)] <- selected_keys[names(selected_keys)]
    }
    if (!is.null(discard_keys)) {
      Keys$DoseResp <- setdiff(Keys$DoseResp, discard_keys)
    }

    # adding 'masked = F' if missing from df_raw_data
    if (!(gDRutils::get_identifier('masked_tag') %in% colnames(df_raw_data))) {
      df_raw_data[, gDRutils::get_identifier('masked_tag')] <- FALSE
    }

    # remove background value to readout (at least 1e-10 to avoid artefactual normalized values)
    df_raw_data$CorrectedReadout <- pmax(df_raw_data$ReadoutValue -
                    df_raw_data$BackgroundValue, 1e-10)
    # creates the DataFrameMatrix and fill with the treated/untreated data
    normSE <- gDR::create_SE(df_raw_data, data_type = "treated", discard_keys = discard_keys)
    SummarizedExperiment::assayNames(normSE) <- "Normalized"
    ctrlSE <- gDR::create_SE(df_raw_data, data_type = "untreated", discard_keys = discard_keys)

    # enforced key values for end points (override selected_keys) --> for rows of the SE
    Keys$untrt_Endpoint <- setdiff(Keys$untrt_Endpoint, names(key_values))
    row_endpoint_value_filter <- rep(TRUE, nrow(ctrlSE))
    if (!is.null(key_values) & length(key_values) > 0) {
      for (i in which(names(key_values) %in% names(SummarizedExperiment::rowData(ctrlSE)))) {
	if (is.numeric(key_values[i])) {
	  row_endpoint_value_filter <- row_endpoint_value_filter &
	      (SummarizedExperiment::rowData(ctrlSE)[, names(key_values)[i] ] == key_values[i] &
		      !is.na(SummarizedExperiment::rowData(ctrlSE)[, names(key_values)[i] ]))
	} else {
	  row_endpoint_value_filter <- row_endpoint_value_filter &
	      (SummarizedExperiment::rowData(ctrlSE)[, names(key_values)[i] ] %in% key_values[i])
	}
      }
    }

    # perform the mapping for normalization
    # first the rows
    # matching the reference end point without any treatment
    row_maps_end <- map_SE(normSE, ctrlSE, row_endpoint_value_filter, Keys)
    names(row_maps_end) <- rownames(normSE)

    # matching the reference end point with the same co-treatment (all the same but conc=0/Gnumber="vehicle")
    row_maps_cotrt <- lapply(rownames(normSE), function(x) {
      ref_metadata_idx <- intersect(Keys$ref_Endpoint,
                            names(SummarizedExperiment::rowData(ctrlSE)))
      names(ref_metadata_idx) = ref_metadata_idx

      rownames(ctrlSE)[which(apply(as.matrix(IRanges::LogicalList(
            lapply(ref_metadata_idx, function(y) # matching the metadata
              SummarizedExperiment::rowData(ctrlSE)[,y] ==
                  (SummarizedExperiment::rowData(normSE)[x, y])
            ))),
            2, all))]
        })
    names(row_maps_cotrt) <- rownames(normSE)

    # reassess the cases without a match to find equivalent drug and concentration (only 2 drugs)
    # test if one can use one of the treatment as a reference
    if ('Gnumber_2' %in% colnames(SummarizedExperiment::rowData(normSE))) {
      for (rnames in names(row_maps_cotrt)[sapply(row_maps_cotrt, length) == 0L]) {
        # loop through the rows without co-treatment matched
        ref_metadata_idx <- setdiff(intersect(Keys$ref_Endpoint,
                                              names(SummarizedExperiment::rowData(ctrlSE))),
                                    c('Gnumber_2', "DrugName_2", 'Concentration_2'))
        names(ref_metadata_idx) <- ref_metadata_idx

        ref_match <- apply(as.matrix(c(IRanges::LogicalList(
          lapply(ref_metadata_idx, function(y) # matching the metadata
                unlist(SummarizedExperiment::rowData(normSE)[, y, drop = FALSE] ==
                    (SummarizedExperiment::rowData(normSE)[rnames, y, drop = FALSE]))
              )), # matching the drugs with mapping from Gnumber to Gnumber_2
              list( Gnumber = SummarizedExperiment::rowData(normSE)$Gnumber ==
                SummarizedExperiment::rowData(normSE)[rnames,'Gnumber_2']),
              list( Gnumber_2 = SummarizedExperiment::rowData(normSE)$Gnumber_2 %in% 
                gDRutils::get_identifier('untreated_tag') ))),
              2, all)
        if (any(ref_match)) {
          row_maps_cotrt[rnames] <- rownames(normSE)[ref_match]
        }
      }
    }

    # matching the reference at time 0 (if available)
    row_maps_T0 <- map_SE(normSE, ctrlSE, row_endpoint_value_filter, Keys, T0 = TRUE)
    names(row_maps_T0) <- rownames(normSE)

    # mapping for columns; 1 to 1 unless overridden by key_values
    col_maps <- array(colnames(ctrlSE), dimnames = list(colnames(normSE)))
    if (any(names(key_values) %in% names(SummarizedExperiment::colData(normSE)))) {
      col_maps[] <- colnames(ctrlSE)[
	which(key_values[names(key_values) %in% names(SummarizedExperiment::colData(normSE))] ==
	    SummarizedExperiment::colData(ctrlSE)[, names(SummarizedExperiment::colData(ctrlSE)) %in% names(key_values)])]
    }

    # creates the DataFrameMatrix for controls
    SummarizedExperiment::assay(normSE, "Controls", withDimnames=FALSE) <- 
      matrix(lapply(seq_len(prod(dim(normSE))), function(x) S4Vectors::DataFrame()), nrow = nrow(normSE), ncol = ncol(normSE))

    # temporary optimization (use 'normSE_n' and 'normSE_c' to avoid using 'assay<-` in for loops)
    # TODO: refactor this part of code once we switch to DataFrameMatrix class
    ctrl_original <- SummarizedExperiment::assay(gDR::aapply(ctrlSE, function(x) x[!x$masked,]))
    # need to keep original data for the case in which reference is such that Gnumber == Gnumber_2
    normSE_n <- normSE_original <- SummarizedExperiment::assay(normSE, "Normalized")
    normSE_c <- SummarizedExperiment::assay(normSE, "Controls")

    for (i in rownames(normSE_n)) {
        for (j in colnames(normSE_n)) {
            if (nrow(normSE_original[[i, j]]) == 0) next # skip if no data
            # get all the control endpoint data
            df_end <- do.call(rbind,
                    lapply(row_maps_end[[i]], function(x) ctrl_original[[x, col_maps[j]]]))
            if (nrow(df_end) == 0) {
              futile.logger::flog.warn(
                  "Missing control data.
                  Treatment Id: %s
                  Cell_line Id: %s",
                  i,
                  j
                )
              next
            }
            df_end <- df_end[, c("CorrectedReadout",
                    intersect(Keys$untrt_Endpoint, colnames(df_end))), drop = FALSE]
            colnames(df_end)[1] <- "UntrtReadout"
            if (ncol(df_end)>1) {
              df_end <- aggregate(df_end[, 1, drop = FALSE],
                by = as.list(df_end[, -1, drop = FALSE]),
                function(x) control_mean_fct(x))
            } else {
              df_end <- DataFrame(UntrtReadout = control_mean_fct(df_end$UntrtReadout))
            }
            # reference co-treatment is not always present
            if (i %in% names(row_maps_cotrt) && length(row_maps_cotrt[[i]])>0) {
              if (all(row_maps_cotrt[[i]] %in% rownames(ctrlSE))) {
                # get all the co-treatment reference endpoint data
                df_ref <- do.call(rbind,
                        lapply(row_maps_cotrt[[i]], function(x) ctrl_original[[x, col_maps[j]]]))
                df_ref <- df_ref[, c("CorrectedReadout",
                        intersect(Keys$ref_Endpoint, colnames(df_ref))), drop = F]
                colnames(df_ref)[1] <- "RefReadout"
                if (ncol(df_ref)>1) {
                  df_ref <- aggregate(df_ref[, 1, drop = FALSE],
                    by = as.list(df_ref[, -1, drop = FALSE]),
                    function(x) control_mean_fct(x))
                } else {
                  df_ref <- DataFrame(RefReadout = control_mean_fct(df_ref$RefReadout))
                }

                # check if all control have matching co-treated wells are on the same plate
                if (all(df_end$Barcode %in% df_ref$Barcode) && all(df_ref$Barcode %in% df_end$Barcode)) {
                  df_end <- merge(df_end, df_ref, # merge based on Barcode
                    by = intersect(colnames(df_end), c('Barcode', Keys$discard_keys)))
                } else {
                  futile.logger::flog.warn(
                      "Control data for the drug are propagated to other plates with co-drug controls.
                      Treatment Id: %s
                      Cell_line Id: %s",
                      i,
                      j
                    )
                  # merge (knowing that there will be NA)
                  df_end <- merge(df_end, df_ref, by = "Barcode", all = TRUE)
                  # propagate average values to the other plates
                  mean_UntrtReadout <- mean(df_end$UntrtReadout, na.rm = TRUE)
                  mean_RefReadout <- mean(df_end$RefReadout, na.rm = TRUE)
                  df_end$UntrtReadout[is.na(df_end$UntrtReadout)] <- mean_UntrtReadout
                  df_end$RefReadout[is.na(df_end$RefReadout)] <- mean_RefReadout
                }

                #gladkia: assert for control data
                if (nrow(df_end) == 0) {
                  stop(sprintf("Control dataframe failed.
                      Treatment Id: '%s'
                      Cell_line Id: %s",
                               i,
                               j))
                }
              } else if (all(row_maps_cotrt[[i]] %in% rownames(normSE))) {
                # case of the reference being with Gnumber == Gnumber_2
                ref_conc <- SummarizedExperiment::rowData(normSE)[i, 'Concentration_2']
                df_ref <- do.call(rbind,
                        lapply(row_maps_cotrt[[i]], function(x) {
                          if (any(normSE_original[[x, col_maps[j]]]$Concentration == ref_conc)) {
                            # the reference value with same concentration is found
                            normSE_original[[x, col_maps[j]]][
                              normSE_original[[x, col_maps[j]]]$Concentration == ref_conc,]
                          } else {
                            # the reference with proper concentration will be inferred --> need fits
                            ref_drc <- normSE_original[[x, col_maps[j]]]
                            if (length(unique(ref_drc$Concentration[!ref_drc$masked])) > 3) {
                              # tryCatch( 
                                # trycatch to write
                                drc_fit = drc::drm(
                                CorrectedReadout ~ Concentration,
                                data = ref_drc[!ref_drc$masked,],
                                fct = drc::LL.4(), # para = c(Hill, x_inf, x0, c50)
                                start = c(2, min(ref_drc$CorrectedReadout),
                                            max(ref_drc$CorrectedReadout),
                                            median(ref_drc$Concentration)),
                                lowerl = c(1e-5, min(ref_drc$CorrectedReadout)*.8, # wide range
                                            min(ref_drc$CorrectedReadout)*.9,
                                            min(ref_drc$Concentration)/1e3),
                                upperl =  c(12, max(ref_drc$CorrectedReadout)*1.1,
                                            max(ref_drc$CorrectedReadout)*1.2,
                                            max(ref_drc$Concentration)*1e3)
                              )
                              df_ref <- data.frame(Concentration = ref_conc,
                                  CorrectedReadout = predict(drc_fit,
                                          data.frame(Concentration = ref_conc)))
                        # )
                          } else {
                            df_ref <- data.frame(Concentration = ref_conc,
                                CorrectedReadout = NA)
                          }
                        }
                      }))
                
                df_ref <- df_ref[, c("CorrectedReadout", 
                                     intersect(Keys$ref_Endpoint, colnames(df_ref))), 
                                 drop = FALSE]
                colnames(df_ref)[1] <- "RefReadout"
                if (ncol(df_ref) > 1) {
                  df_ref <- aggregate(df_ref[, 1, drop = FALSE],
                    by = as.list(df_ref[, -1, drop = FALSE]),
                    function(x) control_mean_fct(x))
                } else {
                  df_ref <- DataFrame(RefReadout = control_mean_fct(df_ref$RefReadout))
                }

                # check if all control have matching co-treated wells are on the same plate
                if (all(df_end$Barcode %in% df_ref$Barcode) && all(df_ref$Barcode %in% df_end$Barcode)) {
                  df_end <- merge(df_end, df_ref,
                    by = intersect(colnames(df_end), c('Barcode', Keys$discard_keys)))
                } else {
                  futile.logger::flog.warn(
                      "Control data for the drug are propagated to other plates with co-drug controls.
                      Treatment Id: %s
                      Cell_line Id: %s",
                      i,
                      j
                    )
                  # propagate average values to the other plates
                  df_end <- merge(df_end, df_ref, by = "Barcode", all = TRUE)
                  mean_UntrtReadout <- mean(df_end$UntrtReadout, na.rm = TRUE)
                  mean_RefReadout <- mean(df_end$RefReadout, na.rm = TRUE)
                  df_end$UntrtReadout[is.na(df_end$UntrtReadout)] <- mean_UntrtReadout
                  df_end$RefReadout[is.na(df_end$RefReadout)] <- mean_RefReadout
                }
              } else {
                stop(sprintf("Reference failed.
                    Treatment Id: '%s'
                    Cell_line Id: %s",
                             i,
                             j))
              }
            } else if (i %in% names(row_maps_cotrt) && length(row_maps_cotrt[[i]])==0) {
              futile.logger::flog.warn("No reference condition found for
                  Treatment Id: '%s'
                  Cell_line Id: %s",
                           i,
                           j)
              df_end$RefReadout <- NA
            } else {
              df_end$RefReadout <- df_end$UntrtReadout
            }

            if (length(row_maps_T0[[i]]) > 0) {
              df_0 <- do.call(rbind,
                      lapply(row_maps_T0[[i]], function(x) ctrl_original[[x, col_maps[j]]]))
              df_0 <- df_0[, c("CorrectedReadout", intersect(Keys$Day0, colnames(df_0)))]
              colnames(df_0)[1] <- "Day0Readout"
              df_0 <- aggregate(df_0[, 1, drop = FALSE], by = as.list(df_0[, -1, drop = FALSE]),
                  function(x) control_mean_fct(x))

              if (!is.null(Keys$discard_keys) && all(Keys$discard_keys %in% colnames(df_0))) {
                df_ctrl <- merge(df_0[, setdiff(colnames(df_0), "Barcode")], df_end, all.y = TRUE, by = Keys$discard_keys)
              } else {
                df_ctrl <- merge(df_0[, setdiff(colnames(df_0), "Barcode")], df_end, all.y = TRUE)
                colnames(df_ctrl)[1] <- "Day0Readout"
              }
            } else {
              df_ctrl <- df_end
              df_ctrl$Day0Readout <- NA
            }

            # calculating the normalized response value for the control
            df_ctrl$RefRelativeViability <- round(df_ctrl$RefReadout/df_ctrl$UntrtReadout,
                ndigit_rounding)
            df_ctrl$RefGRvalue <- round(2 ^ (
                    log2(df_ctrl$RefReadout / df_ctrl$Day0Readout) /
                    log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout) ), ndigit_rounding) - 1
            df_ctrl$DivisionTime <- round(
                    SummarizedExperiment::rowData(normSE)[i,gDRutils::get_identifier("duration")] /
                        log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout), ndigit_rounding)


            #gladkia: assert for merged study/control data
            ctrl_bcodes <- sort(unique(df_ctrl$Barcode))
            trt_bcodes <-
              sort(unique(normSE_original[[i, j]]$Barcode))
            # check if all treated values have matching controls on the same plate
            if (!all(trt_bcodes %in% ctrl_bcodes)) {
              # if not, propagate to all plates
              futile.logger::flog.warn(
                  "Control data are averaged and propagated to treatment plates.
                      Treatment Id: %s (plates %s)
                      Control plates: %s",
                  i, paste(trt_bcodes, collapse = ", "),
                  paste(ctrl_bcodes, collapse = ", ")
                )
              data.table::setDF(data.table::rbindlist(list(df_ctrl, cbind(data.frame(Barcode = setdiff(trt_bcodes, ctrl_bcodes)),
                                                                          t(colMeans(df_ctrl[, setdiff(colnames(df_ctrl), "Barcode")])))), fill = TRUE))
            }

            # works with by = character(0) but changes the order of rows
            df_merged <- merge(
              data.frame(normSE_original[[i, j]]),
                    data.frame(df_ctrl),
                    by = intersect(colnames(df_ctrl), c('Barcode', Keys$discard_keys)),
                    all.x = T)

            # merge the data with the controls assuring that the order of the records is preseved
            # removing this line because failed when   by = character(0)
            # df_merged <- dplyr::left_join(
            #         data.frame(normSE_original[[i, j]]),
            #         data.frame(df_ctrl),
            #         by = intersect(colnames(df_ctrl), c('Barcode', Keys$discard_keys)))

            # calculate the normalized values
            df_merged$RelativeViability <- round(df_merged$CorrectedReadout / df_merged$UntrtReadout, ndigit_rounding)

            df_merged$GRvalue <- round(2 ^ (
              log2(df_merged$CorrectedReadout / df_merged$Day0Readout) /
                log2(df_merged$UntrtReadout / df_merged$Day0Readout)
            ), ndigit_rounding) - 1

            # use the reference doubling Time (ReferenceDivisionTime) for GRvalue if day 0 missing
            if (any(is.na(df_merged$Day0Readout)) ) {
                ref_div_time_col <- gDRutils::get_identifier("cellline_ref_div_time")
                cl_name_col <- gDRutils::get_identifier("cellline_name")
                if (!(ref_div_time_col %in% colnames(SummarizedExperiment::colData(normSE))) ||
                  is.na(SummarizedExperiment::colData(normSE)[j, ref_div_time_col]) ) {
                    # missing division time
                    futile.logger::flog.warn(paste(
                      "No day 0 information and no reference doubling time for cell line", SummarizedExperiment::colData(normSE)[j, cl_name_col],
                      "--> GR values are NA"))
                } else if (SummarizedExperiment::colData(normSE)[j, ref_div_time_col] >
                    1.5 * SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")]) {
                      # check if experiment is long enough relative to division time
                      futile.logger::flog.warn(paste( "Reference doubling time for cell line",
                      SummarizedExperiment::colData(normSE)[j,cl_name_col], "is",
                      SummarizedExperiment::colData(normSE)[j, ref_div_time_col],
                        "which is too long for GR calculation compared to assay duration (",
                      SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")],
                        "--> GR values are NA"))
                 } else {
                   # division time is correct for calculation
                  refDivisionTime <- as.numeric(SummarizedExperiment::colData(normSE)[j, ref_div_time_col])

                  futile.logger::flog.warn(paste(
                    "Missing day 0 information --> calculate GR value based on reference doubling time for", SummarizedExperiment::colData(normSE)[j, cl_name_col]))

                  df_merged$GRvalue <-
                  round(2 ^ (1 + (
                    log2(pmin(1.25, # capping to avoid artefacts
                              df_merged[, "RelativeViability"])) /
                      (SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")] / refDivisionTime)
                  )), ndigit_rounding) - 1

                  df_ctrl$RefGRvalue <-
                  round(2 ^ (1 + (
                    log2(pmin(1.25, # capping to avoid artefacts
                              df_ctrl[, "RefRelativeViability"])) /
                      (SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")] / refDivisionTime)
                  )), ndigit_rounding) - 1
                }
            }

            # more robust assignment in case the order of df_merged has changed
            normSE_n[[i, j]] <- merge(normSE_n[[i, j]],
                df_merged[, c(colnames(normSE_n[[i, j]]), 'GRvalue', 'RelativeViability')],
                by = colnames(normSE_n[[i, j]]))
            normSE_c[[i, j]] <- DataFrame(df_ctrl)
        }
    }
    metadata(normSE) <- c(metadata(normSE),
            list(df_raw_data = df_raw_data,
                Keys = Keys,
                row_maps = list(end = row_maps_end,
                                cotrt = row_maps_cotrt,
                                T0 = row_maps_T0)
                ))

    SummarizedExperiment::assay(normSE, 'Normalized') <- normSE_n
    SummarizedExperiment::assay(normSE, 'Controls') <- normSE_c

    return(normSE)
}


#' normalize_SE2
#'
#' Normalize drug response data from treated and untreated pairings.
#'
#' @param se \code{BumpyMatrix} object with assays \code{"RawTreated"} and \code{"Controls"}.
#' @param control_assay string containing the name of the assay representing the controls in the \code{se}.
#' @param raw_treated_assay string containing the name of the assay representing the raw treated data in the \code{se}.
#' @param ndigit_rounding integer specifying number of digits of rounding during calculations.
#' Defaults to \code{4}.
#'
#' @return \code{BumpyMatrix} object with a new assay named \code{"Normalized"} containing \code{DataFrame}s 
#' holding \code{RelativeViability} and \code{GRvalue}, as well as new assays named \code{RefRelativeViability}, \code{RefGRvalue}, and \code{DivisionTime} values.
#'
#' @export
#'
normalize_SE2 <- function(se, 
                          control_assay = "Controls", 
                          raw_treated_assay = "RawTreated", 
                          ndigit_rounding = 4) {

  # Assertions
  checkmate::assert_number(ndigit_rounding)

  refs <- SummarizedExperiment::assays(se)[[control_assay]]
  trt <- SummarizedExperiment::assays(se)[[raw_treated_assay]]

  discard_keys <- get_SE_keys(se, key_type = "discard_keys")
  trt_keys <- get_SE_keys(se, key_type = "Trt")

  norm_cols <- c("RelativeViability", "GRvalue", "DivisionTime")
  out <- vector("list", nrow(se) * ncol(se))

  ref_rel_viability <- ref_GR_value <- div_time <- matrix(NA, nrow = nrow(se), ncol = ncol(se), dimnames = dimnames(se))
  # Column major order, so go down first.
  cdata <- SummarizedExperiment::colData(se)
  rdata <- SummarizedExperiment::rowData(se)
  for (j in seq_along(colnames(se))) {
    cl_md <- cdata[j, ]
    cl_name <- cl_md[[gDRutils::get_identifier("cellline_name")]]
    ref_div_time <- cl_md[[gDRutils::get_identifier("cellline_ref_div_time")]]

    for (i in seq_along(rownames(se))) {
      duration <- rdata[i, gDRutils::get_identifier("duration")]

      ref_df <- refs[i, j][[1]]
      trt_df <- trt[i, j][[1]]

      if (length(trt_df) == 0L || nrow(trt_df) == 0L) {
	next # skip if no data
        # TODO: Does this still need to initialize an empty DFrame with appropriate colnames?
      }

      if (length(ref_df) == 0L || nrow(ref_df) == 0L) {
	futile.logger::flog.warn(
          sprintf("Missing control data. Treatment Id: '%s' Cell_line Id: '%s'", 
            rownames(se)[i], colnames(se)[j])
        )
	next
      }
      
      # Merge to ensure that the proper discard_key values are mapped.
      all_readouts_df <- merge(trt_df, 
        ref_df, 
	by = discard_keys,
	all.x = TRUE)

      normalized <- DataFrame(matrix(NA, nrow = nrow(trt_df), ncol = length(norm_cols)))
      colnames(normalized) <- c(norm_cols)

      # Normalized treated.
      normalized$RelativeViability <- round(all_readouts_df$CorrectedReadout/all_readouts_df$UntrtReadout, ndigit_rounding)
      normalized$GRvalue <- calculate_GR_value(rel_viability = normalized$RelativeViability, 
        corrected_readout = all_readouts_df$CorrectedReadout, 
        day0_readout = all_readouts_df$Day0Readout, 
        untrt_readout = all_readouts_df$UntrtReadout, 
        ndigit_rounding = ndigit_rounding, 
        duration = duration, 
        ref_div_time = ref_div_time, 
        cl_name = cl_name)

      # Carry over present treated keys.
      normalized <- cbind(all_readouts_df[intersect(c(trt_keys, gDRutils::get_identifier("masked_tag")), colnames(all_readouts_df))], normalized) 

      normalized$row_id <- rep(rownames(se)[i], nrow(trt_df))
      normalized$col_id <- rep(colnames(se)[j], nrow(trt_df))

      out[[nrow(se) * (j - 1) + i]] <- normalized

      ###########################
      ref_df <- all_readouts_df # TODO: keep this here for now in order to make comparisons between old and new method.
      # !!!!!!!!!!!!HOWEVER, THE ABOVE LINE SHOULD LATER BE DELETED!!!!!!!!
      ###########################

      ## Perform the calculations on all references.
      ## Then, take the mean to report the final reference normalized value.
      RV_vec <- ref_df$RefReadout/ref_df$UntrtReadout
      GR_vec <- calculate_GR_value(rel_viability = ref_rv_value, 
        corrected_readout = ref_df$RefReadout, 
        day0_readout = ref_df$Day0Readout, 
        untrt_readout = ref_df$UntrtReadout, 
        ndigit_rounding = ndigit_rounding, 
        duration = duration, 
        ref_div_time = ref_div_time, 
        cl_name = cl_name)
      ref_rel_viability[i, j] <- round(mean(RV_vec, na.rm=TRUE), ndigit_rounding)
      ref_GR_value[i, j] <- round(mean(GR_vec, na.rm = TRUE), ndigit_rounding)
      div_time[i, j] <- round(duration / log2(mean(ref_df$UntrtReadout/ref_df$Day0Readout, na.rm = TRUE)), ndigit_rounding)
    }
  }

  out <- DataFrame(do.call("rbind", out))
  norm <- BumpyMatrix::splitAsBumpyMatrix(out[!colnames(normalized) %in% c("row_id", "col_id")], 
    row = out$row_id, 
    col = out$col_id)

  SummarizedExperiment::assays(se)[["Normalized"]] <- norm

  SummarizedExperiment::assays(se)[["RefGRvalue"]] <- ref_GR_value
  SummarizedExperiment::assays(se)[["RefRelativeViability"]] <- ref_rel_viability
  SummarizedExperiment::assays(se)[["DivisionTime"]] <- div_time

  return(se)
}
