#' normalize_SE
#'
#' Normalize raw DR data
#'
#' @param df_raw_data a dataframe with raw data
#' @param selected_keys a vector of keys that should be included in the normalization (NULL by default)
#' @param key_values a list of values for keys that should be included in the normalization (NULL by default)
#' @param discard_keys a vector of keys that should be discarded (NULL by default)
#' @param control_mean_fct the function used for averaging controls (trimmed arithmetic mean with trim = .25 by default)
#' @param nDigits_rounding an integer with number of digits of rounding (4 by default)
#'
#' @return a SummarizedExperiment object with normalized assays
#' @export
#'
normalize_SE <- function(df_raw_data,
                  selected_keys = NULL,
                  key_values = NULL,
                  discard_keys = NULL,
                  control_mean_fct = function(x) mean(x, trim = .25), # used for averaging controls
                  nDigits_rounding = 4 # rounding of normalized response values
                ) {
    # Assertions
    stopifnot(inherits(df_raw_data, "data.frame"))
    checkmate::assert_vector(selected_keys, null.ok = TRUE)
    checkmate::assert_list(key_values, null.ok = TRUE)
    checkmate::assert_vector(discard_keys, null.ok = TRUE)
    checkmate::assert_function(control_mean_fct, null.ok = TRUE)
    checkmate::assert_number(nDigits_rounding)
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
    if ( !(gDRutils::get_identifier('masked_tag') %in% colnames(df_raw_data))) {
      df_raw_data[,gDRutils::get_identifier('masked_tag')] = FALSE
    }

    # remove background value to readout (at least 1e-10 to avoid artefactual normalized values)
    df_raw_data$CorrectedReadout = pmax(df_raw_data$ReadoutValue -
                    df_raw_data$BackgroundValue, 1e-10)
    # creates the DataFrameMatrix and fill with the treated/untreated data
    normSE <- gDR::createSE(df_raw_data, data_type = "treated", discard_keys = discard_keys)
    SummarizedExperiment::assayNames(normSE) = "Normalized"
    ctrlSE <- gDR::createSE(df_raw_data, data_type = "untreated", discard_keys = discard_keys)

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
            }}}

    # perform the mapping for normalization
    # first the rows
    # matching the reference end point without any treatment
    row_maps_end <- map_SE(normSE, ctrlSE, row_endpoint_value_filter, Keys)
    names(row_maps_end) <- rownames(normSE)

    # matching the reference end point with the same co-treatment (all the same but conc=0/Gnumber="vehicle")
    row_maps_cotrt <- lapply(rownames(normSE), function(x) {
      ref_metadata_idx = intersect(Keys$ref_Endpoint,
                            names(SummarizedExperiment::rowData(ctrlSE)))
      names(ref_metadata_idx) = ref_metadata_idx

      rownames(ctrlSE)[which(apply(as.matrix(  IRanges::LogicalList(
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
      for (rnames in names(row_maps_cotrt)[sapply(row_maps_cotrt, length)==0]) {
        # loop through the rows without co-treatment matched
        ref_metadata_idx = setdiff(intersect(Keys$ref_Endpoint,
                              names(SummarizedExperiment::rowData(ctrlSE))),
                            c('Gnumber_2', "DrugName_2", 'Concentration_2'))
        names(ref_metadata_idx) = ref_metadata_idx

        ref_match = apply(as.matrix(  c(IRanges::LogicalList(
          lapply(ref_metadata_idx, function(y) # matching the metadata
                unlist(SummarizedExperiment::rowData(normSE)[,y, drop=F] ==
                    (SummarizedExperiment::rowData(normSE)[rnames, y, drop=F]))
              )), # matching the drugs with mapping from Gnumber to Gnumber_2
              list( Gnumber = SummarizedExperiment::rowData(normSE)$Gnumber ==
                SummarizedExperiment::rowData(normSE)[rnames,'Gnumber_2']),
              list( Gnumber_2 = SummarizedExperiment::rowData(normSE)$Gnumber_2 %in% 
                gDRutils::get_identifier('untreated_tag') ))),
              2, all)
        if (any(ref_match)) row_maps_cotrt[rnames] = rownames(normSE)[which(ref_match)]
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
    SummarizedExperiment::assay(normSE, "Controls", withDimnames=FALSE) <- matrix(lapply(1:prod(dim(normSE)), function(x) S4Vectors::DataFrame()),
            nrow = nrow(normSE), ncol = ncol(normSE))

    # temporary optimization (use 'normSE_n' and 'normSE_c' to avoid using 'assay<-` in for loops)
    # TODO: refactor this part of code once we switch to DataFrameMatrix class
    ctrl_original = SummarizedExperiment::assay(gDR::aapply(ctrlSE, function(x) x[!x$masked,]))
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
              df_end = DataFrame(UntrtReadout = control_mean_fct(df_end$UntrtReadout))
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
                  df_ref = DataFrame(RefReadout = control_mean_fct(df_ref$RefReadout))
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
                ref_conc = SummarizedExperiment::rowData(normSE)[i,'Concentration_2']
                df_ref <- do.call(rbind,
                        lapply(row_maps_cotrt[[i]], function(x) {
                          if (any(normSE_original[[x, col_maps[j]]]$Concentration == ref_conc)) {
                            # the reference value with same concentration is found
                            normSE_original[[x, col_maps[j]]][
                              normSE_original[[x, col_maps[j]]]$Concentration == ref_conc,]
                          } else {
                            # the reference with proper concentration will be inferred --> need fits
                            ref_drc = normSE_original[[x, col_maps[j]]]
                            if (length(unique(ref_drc$Concentration[!ref_drc$masked]))>3) {
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
                              df_ref = data.frame(Concentration = ref_conc,
                                  CorrectedReadout = predict(drc_fit,
                                          data.frame(Concentration = ref_conc)))
                        # )
                          } else {
                            df_ref = data.frame(Concentration = ref_conc,
                                CorrectedReadout = NA)
                          }
                          }
                          }))
                
                df_ref <- df_ref[, c("CorrectedReadout",
                        intersect(Keys$ref_Endpoint, colnames(df_ref))), drop = F]
                colnames(df_ref)[1] <- "RefReadout"
                if (ncol(df_ref)>1) {
                  df_ref <- aggregate(df_ref[, 1, drop = FALSE],
                    by = as.list(df_ref[, -1, drop = FALSE]),
                    function(x) control_mean_fct(x))
                } else {
                  df_ref = DataFrame(RefReadout = control_mean_fct(df_ref$RefReadout))
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
              df_ctrl = df_end
              df_ctrl$Day0Readout = NA
            }

            # calculating the normalized response value for the control
            df_ctrl$RefRelativeViability <- round(df_ctrl$RefReadout/df_ctrl$UntrtReadout,
                nDigits_rounding)
            df_ctrl$RefGRvalue <- round(2 ** (
                    log2(df_ctrl$RefReadout / df_ctrl$Day0Readout) /
                    log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout) ), nDigits_rounding) - 1
            df_ctrl$DivisionTime <- round(
                    SummarizedExperiment::rowData(normSE)[i,gDRutils::get_identifier("duration")] /
                        log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout), nDigits_rounding)


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
              dplyr::bind_rows(df_ctrl, cbind(data.frame(Barcode = setdiff(trt_bcodes, ctrl_bcodes)),
                        t(colMeans(df_ctrl[, setdiff(colnames(df_ctrl), "Barcode")]))))
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
            df_merged$RelativeViability <- round(df_merged$CorrectedReadout / df_merged$UntrtReadout, nDigits_rounding)

            df_merged$GRvalue = round(2 ** (
              log2(df_merged$CorrectedReadout / df_merged$Day0Readout) /
                log2(df_merged$UntrtReadout / df_merged$Day0Readout)
            ), nDigits_rounding) - 1

            # use the reference doubling Time (ReferenceDivisionTime) for GRvalue if day 0 missing
            if ( any(is.na(df_merged$Day0Readout)) ) {

                if ( !(gDRutils::get_header('add_clid')[3] %in% colnames(SummarizedExperiment::colData(normSE))) ||
                  is.na(SummarizedExperiment::colData(normSE)[j, gDRutils::get_header('add_clid')[3]]) ) {
                    # missing division time
                    futile.logger::flog.warn(paste(
                      "No day 0 information and no reference doubling time for cell line", SummarizedExperiment::colData(normSE)[j,gDRutils::get_header('add_clid')[1]],
                      "--> GR values are NA"))
                } else if (SummarizedExperiment::colData(normSE)[j, gDRutils::get_header('add_clid')[3]] >
                    1.5 * SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")]) {
                      # check if experiment is long enough relative to division time
                      futile.logger::flog.warn(paste( "Reference doubling time for cell line",
                      SummarizedExperiment::colData(normSE)[j,gDRutils::get_header('add_clid')[1]], "is",
                      SummarizedExperiment::colData(normSE)[j, gDRutils::get_header('add_clid')[3]],
                        "which is too long for GR calculation compared to assay duration (",
                      SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")],
                        "--> GR values are NA"))
                 } else {
                   # division time is correct for calculation
                  refDivisionTime = as.numeric(SummarizedExperiment::colData(normSE)[j, gDRutils::get_header('add_clid')[3]])

                  futile.logger::flog.warn(paste(
                    "Missing day 0 information --> calculate GR value based on reference doubling time for", SummarizedExperiment::colData(normSE)[j,gDRutils::get_header('add_clid')[1]]))

                  df_merged$GRvalue <-
                  round(2 ^ (1 + (
                    log2(pmin(1.25, # capping to avoid artefacts
                              df_merged[, "RelativeViability"])) /
                      (SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")] / refDivisionTime)
                  )), nDigits_rounding) - 1

                  df_ctrl$RefGRvalue <-
                  round(2 ^ (1 + (
                    log2(pmin(1.25, # capping to avoid artefacts
                              df_ctrl[, "RefRelativeViability"])) /
                      (SummarizedExperiment::rowData(normSE)[i, gDRutils::get_identifier("duration")] / refDivisionTime)
                  )), nDigits_rounding) - 1
                }
            }

            # more robust assignment in case the order of df_merged has changed
            normSE_n[[i, j]] = merge(normSE_n[[i, j]],
                df_merged[, c(colnames(normSE_n[[i, j]]), 'GRvalue', 'RelativeViability')],
                by = colnames(normSE_n[[i, j]]))
            normSE_c[[i,j]] <- DataFrame(df_ctrl)
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
