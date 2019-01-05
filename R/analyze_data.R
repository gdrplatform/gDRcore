library(gCellGenomics) # best reference for cell line and drug names?

load_all_data = function(manifest_file, template_file, results_file, log_file) {

    manifest = load_manifest(manifest_file, log_file)
    stopifnot(all(unique(manifest$Template) %in% basename(template_file)))
    treatments = load_templates(template_file, log_file)
    data = load_results(results_file, log_file)

    # merge manifest and treatments first
    df_metadata = merge(manifest, treatments, by = 'Template')
    # sort out duplicate metadata columns
    duplicated_col = setdiff(intersect(colnames(manifest), colnames(treatments)), 'Template')
    for (m_col in duplicated_col) {
        df_metadata[,m_col] = df_metadata[,paste0(m_col,'.y')] # parse template values
        missing_idx = is.na(df_metadata[,m_col]) | df_metadata[,m_col] %in% c('', '-')
        # add manifest values when missing in template
        df_metadata[missing_idx,m_col] = df_metadata[missing_idx,paste0(m_col,'.x')]
        # check for conflicts
        double_idx = !(is.na(df_metadata[,paste0(m_col,'.x')]) |
                            df_metadata[,paste0(m_col,'.x')] %in% c('', '-')) &
                     !(is.na(df_metadata[,paste0(m_col,'.y')]) |
                            df_metadata[,paste0(m_col,'.y')] %in% c('', '-'))
        if (any(double_idx) &&
            any(df_metadata[,paste0(m_col,'.x')] != df_metadata[,paste0(m_col,'.y')], na.rm=T)) {
                WarnMsg = paste('Metadata field', m_col,
                    'found in both the manifest and some templates with inconsistent values;',
                    'values in template supersede the ones in the manifest')
                writeLines('Warning in load_all_data:', log_file)
                writeLines(WarnMsg, log_file)
                warning(WarnMsg)
            }
        df_metadata[,paste0(m_col,'.x')] = NULL
        df_metadata[,paste0(m_col,'.y')] = NULL
    }

    # clean up the metadata
    cleanedup_metadata = cleanup_metadata(df_metadata)
    stopifnot( dim(cleanedup_metadata)[1] == dim(df_metadata)[1] ) # should not happen

    df_merged = merge(cleanedup_metadata, data, by = c('Barcode', 'WellRow', 'WellColumn'))
    if (dim(df_merged)[1] != dim(data)[1]) {# need to identify issue and output relevant warning
        WarnMsg = 'Not all results have been matched with treatments; merged table is smaller than data table'
        writeLines('Warning in load_all_data:', log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
    }
    if (dim(df_merged)[1] != dim(df_metadata)[1]) {# need to identify issue and output relevant warning
        WarnMsg = 'Not all treatments have been matched with results; merged table is smaller than metadata table'
        writeLines('Warning in load_all_data:', log_file)
        writeLines(WarnMsg, log_file)
        warning(WarnMsg)
    }

    # remove wells not labeled
    df_raw_data = df_merged[ !is.na(df_merged$Gnumber), ]
    print(sprintf('%i well loaded, %i discarded for lack of annotation, %i data point selected',
                dim(data)[1], sum(is.na(df_merged$Gnumber)), dim(df_raw_data)[1]))

    # reorder the columns
    cols = c('CellLineName', 'Tissue', 'DrugName', 'Concentration',
        intersect(paste0(c('DrugName_','Concentration_'),sort(c(2:10,2:10))),colnames(df_raw_data)))
    cols = c(cols, setdiff(colnames(df_metadata), c('Barcode', 'Template', "WellRow", "WellColumn",
                        cols, 'Gnumber', paste0('Gnumber_', 2:10))))
    cols = c(cols, "ReadoutValue", "BackgroundValue")
    cols = c(cols, setdiff(colnames(df_raw_data), cols))
    df_raw_data = df_raw_data[, cols ]
    df_raw_data = df_raw_data[order(df_raw_data$CellLineName, df_raw_data$DrugName,
                            df_raw_data$Concentration, df_raw_data$Time), ]

    return(df_raw_data)
}






cleanup_metadata = function(df_metadata) {

    # clean up numberic fields
    df_metadata$Time = round(as.numeric(df_metadata$Time),6)
    # identify potential numeric fields and replace NA by 0 - convert strings in factors
    for (c in setdiff(1:dim(df_metadata)[2], c( agrep('Gnumber', colnames(df_metadata)),
            agrep('Concentration', colnames(df_metadata)),
            grep('Time|CLid|Barcode|WellRow|WellColumn|Template', colnames(df_metadata)) ))) {
        vals = unique(df_metadata[,c])
        if (is.numeric(vals)) {
            df_metadata[is.na(df_metadata[,c]),c] = 0
        } else if (is.character(vals)) {
            num_vals = as.numeric(vals)
            if (sum(is.na(num_vals))>2 || all(is.na(num_vals))) {
                df_metadata[,c] = factor(df_metadata[,c])
                WarnMsg = paste('Metadata field ', colnames(df_metadata)[c],
                                    ' converted to factors')
                writeLines('Warning in cleanup_metadata:', log_file)
                writeLines(WarnMsg, log_file)
                warning(WarnMsg)
            } else {
                is.na(df_metadata[,c]) = 0
                df_metadata[,c] = as.numeric(df_metadata[,c])
                WarnMsg = paste('Metadata field ', colnames(df_metadata)[c],
                                    ' converted to numeric values')
                writeLines('Warning in cleanup_metadata:', log_file)
                writeLines(WarnMsg, log_file)
                warning(WarnMsg)
            }
        }
    }

    # check that CLid are in the format 'CL####' and add common name
    gCLs = gCellGenomics::getSamples()[,c('clid', 'celllinename', 'tissue')]
    colnames(gCLs)[2:3] = c('CellLineName', 'Tissue')
    CLids = unique(df_metadata$CLid)
    bad_CL = !(CLids %in% gCLs$clid)
    if (any(bad_CL)) {
        ErrorMsg = paste('Cell line ID ', paste(CLids[bad_CL], collapse = ' ; '),
            ' not found in gCSI database')
        writeLines('Error in load_all_data:', log_file)
        writeLines(ErrorMsg, log_file)
        close(log_file)
        stop(ErrorMsg)
        }
    df_metadata = merge(df_metadata, gCLs, by.x='CLid', by.y='clid', all.x = T)

    # check that Gnumber_* are in the format 'G####' and add common name (or Vehicle or Untreated)
    untrt_flag = c('Vehicle', 'Untreated')
    for (i in agrep('Gnumber', colnames(df_metadata))) { # correct case issues
        for (w in untrt_flag) {
            df_metadata[grep(w, df_metadata[,i], ignore.case = T),i] = w
        }
    }
    gDrugs = gCellGenomics::getDrugs()[,c('drug', 'gcsi_drug_name')]
    gDrugs$drug = substr(gDrugs$drug, 1, 9)
    colnames(gDrugs)[2] = 'DrugName'
    gDrugs = rbind(data.frame(drug=untrt_flag, DrugName=untrt_flag),
                gDrugs)
    Gnbrs = unique(unlist(df_metadata[,agrep('Gnumber', colnames(df_metadata))]))
    bad_Gn = !(Gnbrs %in% gDrugs$drug) & !is.na(Gnbrs)
    if (any(bad_CL)) {
        ErrorMsg = paste('Drug ', paste(Gnbrs[bad_Gn], collapse = ' ; '),
            ' not found in gCSI database')
        writeLines('Error in load_all_data:', log_file)
        writeLines(ErrorMsg, log_file)
        close(log_file)
        stop(ErrorMsg)
        }
    df_metadata = merge(df_metadata, gDrugs, by.x='Gnumber', by.y='drug', all.x = T)
    # add info for columns Gnumber_*
    for (i in grep('Gnumber_\\d', colnames(df_metadata))) {
        df_metadata[ is.na(df_metadata[,i]), i] = 'Untreated' # set missing values to Untreated
        gDrugs_ = gDrugs
        colnames(gDrugs_)[2] = paste0(colnames(gDrugs_)[2], substr(colnames(df_metadata)[i], 8, 12))
        df_metadata = merge(df_metadata, gDrugs_, by.x=i, by.y='drug', all.x = T)
    }

    # clean up concentration fields
    for (i in agrep('Concentration', colnames(df_metadata))) {
        trt_n = ifelse(regexpr('_\\d', colnames(df_metadata)[i])>0,
                            substr(colnames(df_metadata)[i], 15,20), 1)
        G_col = ifelse(trt_n == 1, 'Gnumber', paste0('Gnumber_', trt_n))
        df_metadata[df_metadata[,G_col] %in% untrt_flag,i] = 0 # set all untreated to 0
        df_metadata[,i] = round(as.numeric(df_metadata[,i]),6) # avoid mismatch due to string truncation
    }


    return(df_metadata)
}
