normalize_SE_no0 <- function(df_raw_data, selected_keys = NULL,
                         key_values = NULL, discard_keys = NULL) {
  # average technical replicates and assign the right controls to each treated well
  Keys <- identify_keys(df_raw_data)
  Keys$discard_keys <- discard_keys
  if (!is.null(selected_keys)) {
    Keys[names(selected_keys)] <- selected_keys[names(selected_keys)]
  }
  if (!is.null(discard_keys)) {
    Keys$DoseResp <- setdiff(Keys$DoseResp, discard_keys)
  }
  # 
  # # the normalized SE only contains the treated conditions
  # if(length(Keys$Day0)==0){
  #   df_raw_data$CorrectedReadout <- pmax(df_raw_data$ReadoutValue - df_raw_data$BackgroundValue, 1)
  #   df_raw_data$UntrtReadout <- df_raw_data$CorrectedReadout
  #   
  #   df_ctrl <- aggregate(df_raw_data$CorrectedReadout[ df_raw_data$Concentration == 0],
  #                        by = list(clid = df_raw_data$clid[ df_raw_data$Concentration == 0]),
  #                        function(x) mean(x, trim = .25))
  #   df_raw_data$RelativeViability <- df_raw_data$CorrectedReadout / df_raw_data$UntrtReadout
  #   df_raw_data$GRvalue = round(2 ^ (1 + (
  #     log2(pmin(1.25,
  #               df_raw_data[, "RelativeViability"])) /
  #       (df_raw_data$Duration / df_raw_data$ReferenceDivisionTime)
  #   )), 4) - 1
    normSE <- gDR::createSE(df_raw_data, data_type = "treated", discard_keys = discard_keys) # in original function data_type = "treated"
    SummarizedExperiment::assayNames(normSE) = "Normalized"

    ctrlSE <- gDR::createSE(df_raw_data, data_type = "all", discard_keys = discard_keys) # we have no control (original data_type = "untreated" returns an error)
    
    # enforced key values for end points (override selected_keys) --> for rows of the SE
    Keys$untrt_Endpoint <- setdiff(Keys$untrt_Endpoint, names(key_values))
    row_endpoint_value_filter <- array(TRUE, nrow(ctrlSE))
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
    
    #TODO gladkia: move mapping code to the separate function
    # perform the mapping for normalization
    # first the rows
    # matching the reference end point without any treatment
    row_maps_end <- lapply(rownames(normSE), function(x) {
      # define matix with matching metadata
      match_mx <- c(
        (
          SummarizedExperiment::rowData(ctrlSE) == (SummarizedExperiment::rowData(normSE)[x,])
        )[intersect(Keys$untrt_Endpoint,
                    names(SummarizedExperiment::rowData(ctrlSE)))],
        IRanges::LogicalList(
          key_values = row_endpoint_value_filter,
          conc = apply(cbind(array(0, nrow(ctrlSE)), # padding to avoid empty df
                             SummarizedExperiment::rowData(ctrlSE)[, agrep("Concentration",
                                                                           colnames(SummarizedExperiment::rowData(ctrlSE))), drop = FALSE]), 1,
                       function(x)
                         all(x == 0))
        )
      )
      match_idx <- which(apply(as.matrix(match_mx), 2, all))
      if (length(match_idx) == 0) {
        # if not exact match, try to find best match
        futile.logger::flog.warn("Missing treated contols for: %s", x)
        idx <-
          apply(as.matrix(match_mx), 2, function(y)
            sum(y, na.rm = TRUE)) *
          match_mx[[get_identifier("duration")]]
        if (any(idx > 0)) {
          match_idx <- which.max(idx)
          futile.logger::flog.warn("Found partial match:",
                                   rownames(ctrlSE)[match_idx])
        } else {
          futile.logger::flog.warn("No partial match found")
        }
      }
      return(rownames(ctrlSE)[match_idx])
    })
    names(row_maps_end) <- rownames(normSE)
    
    # matching the reference end point with the same co-treatment (all the same but conc=0/Gnumber="vehicle")
    # row_maps_cotrt <- lapply(rownames(normSE), function(x)
    #   rownames(ctrlSE)[which(apply(as.matrix(
    #     (SummarizedExperiment::rowData(ctrlSE) == (SummarizedExperiment::rowData(normSE)[x, ]))[
    #       intersect(Keys$ref_Endpoint, names(SummarizedExperiment::rowData(ctrlSE)))]),
    #     2, all))])
    # names(row_maps_cotrt) <- rownames(normSE)
    row_maps_cotrt <- lapply(rownames(normSE), function(x)
      rownames(ctrlSE))
    names(row_maps_cotrt) <- rownames(normSE)
    
    # matching the reference at time 0 (if available)
    # row_maps_T0 <- lapply(rownames(normSE), function(x) {
    #   # define matix with matching metadata
    #   match_mx <- c(
    #     (SummarizedExperiment::rowData(ctrlSE) == (SummarizedExperiment::rowData(normSE)[x,]))[
    #       intersect(Keys$Day0, names(SummarizedExperiment::rowData(ctrlSE)))],
    #     IRanges::LogicalList(#key_values = row_endpoint_value_filter,
    #       T0 = SummarizedExperiment::rowData(ctrlSE)[, get_identifier("duration")] == 0,
    #       conc = apply(cbind(array(0, nrow(ctrlSE)),# padding to avoid empty df
    #                          SummarizedExperiment::rowData(ctrlSE)[, agrep("Concentration",
    #                                                                        colnames(SummarizedExperiment::rowData(ctrlSE))), drop = FALSE]), 1,
    #                    function(x) all(x == 0)) ))
    #   match_idx <- which(apply(as.matrix(match_mx), 2, all))
    #   if (length(match_idx) == 0) {
    #     # if not exact match, try to find best match
    #     futile.logger::flog.warn("Missing day 0 plate for: %s", x)
    #     idx <- apply(as.matrix(match_mx), 2, function(y) sum(y, na.rm = TRUE)) *
    #       match_mx[["T0"]]
    #     if (any(idx > 0)) {
    #       match_idx <- which.max(idx)
    #       futile.logger::flog.warn("Found partial match: %s",
    #                                rownames(ctrlSE)[match_idx])
    #     } else {
    #       futile.logger::flog.warn("No partial match found")
    #     }
    #   }
    #   return(rownames(ctrlSE)[match_idx])
    # })
    row_maps_T0 <- row_maps_cotrt
    names(row_maps_T0) <- rownames(normSE)
    
    # mapping for columns; 1 to 1 unless overridden by key_values
    col_maps <- array(colnames(ctrlSE), dimnames = list(colnames(normSE)))
    if (any(names(key_values) %in% names(SummarizedExperiment::colData(normSE)))) {
      col_maps[] <- colnames(ctrlSE)[
        which(key_values[names(key_values) %in% names(SummarizedExperiment::colData(normSE))] ==
                SummarizedExperiment::colData(ctrlSE)[, names(SummarizedExperiment::colData(ctrlSE)) %in% names(key_values)])]
    }
    
    # remove background value to readout
    normSE <- aapply(normSE, function(x) {
      x$CorrectedReadout <- pmax(x$ReadoutValue - x$BackgroundValue, 1)
      return(x)},
      "Normalized")
    ctrlSE <- aapply(ctrlSE, function(x) {
      x$CorrectedReadout <- pmax(x$ReadoutValue - x$BackgroundValue, 1)
      return(x)})
    
    SummarizedExperiment::assay(normSE, "Controls") <- matrix(lapply(1:prod(dim(normSE)), function(x) S4Vectors::DataFrame()),
                                                              nrow = nrow(normSE), ncol = ncol(normSE))
    
    # run through all conditions to assign controls and normalize the data
    # TODO: optimize (could that be replaced by a lapply?? or dyplr function??)
    
    # temporary optimization (use 'normSE_n' and 'normSE_c' to avoid using 'assay<-` in a foor loops)
    # TODO: refactor this part of code once we switch to DataFrameMatrix class
    normSE_n <- SummarizedExperiment::assay(normSE, "Normalized")
    normSE_c <- SummarizedExperiment::assay(normSE, "Controls")
    for (i in rownames(normSE_n)) {
      for (j in colnames(normSE_n)) {
        
        if (nrow(SummarizedExperiment::assay(normSE, "Normalized")[[i, j]]) == 0) next # skip if no data
        
        df_end <- do.call(rbind,
                          lapply(row_maps_end[[i]], function(x) SummarizedExperiment::assay(ctrlSE)[[x, col_maps[j]]]))
        df_end <- df_end[, c("CorrectedReadout",
                             intersect(Keys$untrt_Endpoint, colnames(df_end)))]
        colnames(df_end)[1] <- "UntrtReadout"
        df_end <- aggregate(df_end[, 1, drop = FALSE], by = as.list(df_end[, -1, drop = FALSE]),
                            function(x) mean(x, trim = .25))
        
        df_ref <- do.call(rbind,
                            lapply(row_maps_cotrt[[i]], function(x) SummarizedExperiment::assay(ctrlSE)[[x, col_maps[j]]]))
        df_ref <- df_ref[, c("CorrectedReadout",
                               intersect(Keys$ref_Endpoint, colnames(df_ref)))]
        colnames(df_ref)[1] <- "RefReadout"
          df_ref <- aggregate(df_ref[, 1, drop = FALSE], by = as.list(df_ref[, -1, drop = FALSE]),
                              function(x) mean(x, trim = .25))
          
                      df_end <- merge(df_end, df_ref, by = "Barcode", all = TRUE)
            mean_UntrtReadout <- mean(df_end$UntrtReadout, na.rm = TRUE)
            mean_RefReadout <- mean(df_end$RefReadout, na.rm = TRUE)
            df_end$UntrtReadout[is.na(df_end$UntrtReadout)] <- mean_UntrtReadout
            df_end$RefReadout[is.na(df_end$RefReadout)] <- mean_RefReadout

        df_0 <- do.call(rbind,
                        lapply(row_maps_T0[[i]], function(x) SummarizedExperiment::assay(ctrlSE)[[x, col_maps[j]]]))
        df_0 <- df_0[, c("CorrectedReadout", colnames(df_0))]
        colnames(df_0)[1] <- "Day0Readout"
        df_0 <- aggregate(df_0[, 1, drop = FALSE], by = as.list(df_0[, -1, drop = FALSE]),
                          function(x) mean(x, trim = .25))
        
          df_ctrl <- merge(df_0[, setdiff(colnames(df_0), "Barcode")], df_end, all.y = TRUE, by = Keys$discard_keys)
        
        df_ctrl$RefRelativeViability <- round(df_ctrl$RefReadout/df_ctrl$UntrtReadout, 4)
        
        df_ctrl$RefGRvalue <- round(2 ** (
          log2(df_ctrl$RefReadout / df_ctrl$Day0Readout) /
            log2(df_ctrl$UntrtReadout / df_ctrl$Day0Readout) ), 4) - 1
        # df_ctrl$RefGRvalue = round(2 ^ (1 + (
        #   log2(pmin(1.25,
        #             df_ctrl[, "RefRelativeViability"])) /
        #     (df_raw_data$Duration / df_raw_data$ReferenceDivisionTime)
        # )), 4) - 1
        # 
        
        
        df_ctrl$DivisionTime <- gneDB::annotateCLIDs(stringr::str_extract(j, "CL[0-9]*"))$doublingtime
        normSE_c[[i,j]] <- DataFrame(df_ctrl)
        
        #gladkia: assert for merged study/control data
        ctrl_bcodes <- sort(unique(df_ctrl$Barcode))
        trt_bcodes <-
          sort(unique(SummarizedExperiment::assay(normSE, "Normalized")[[i, j]]$Barcode))
        if (!all(trt_bcodes %in% ctrl_bcodes)) {
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
        
        # merge the data with the controls assuring that the order of the records is preseved
        df_merged <- data.frame(SummarizedExperiment::assay(normSE, "Normalized")[[i, j]])
        
        # calculate the normalized values
        normSE_n[[i, j]]$UntrtReadout <- normSE_n[[i, j]]$CorrectedReadout
        
        normSE_n[[i, j]]$RelativeViability <-
          round(normSE_n[[i, j]]$CorrectedReadout / normSE_n[[i, j]]$UntrtReadout, 4)
        

        normSE_n[[i, j]]$GRvalue = round(2 ^ (1 + (
          log2(pmin(1.25,
                    normSE_n[[i, j]]$RelativeViability)) /
            (unique(df_raw_data$Duration) / unique(df_ctrl$DivisionTime))
        )), 4) - 1
        
        }
      }
    metadata(normSE) <- c(metadata(normSE),
                          list(df_raw_data = df_raw_data,
                               Keys = Keys,
                               row_maps = list(end = row_maps_end, # empty list?2
                                               cotrt = row_maps_cotrt,
                                               T0 = row_maps_T0)
                          ))
    
    SummarizedExperiment::assay(normSE, 'Normalized') <- normSE_n
    SummarizedExperiment::assay(normSE, 'Controls') <- normSE_c
    
    return(normSE)
  }
  }
