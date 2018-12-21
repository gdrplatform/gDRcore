library(readxl) # openxlsx skip the first emprty rows


load_manifest = function (manifest_file, log_file) {
    # manifest_file is a string or a vector of strings
    # log_file is an open file to sink errors and warnings

    # read files
    manifest_data = lapply(manifest_file, function(x) {
        df = read_excel(x)
    } )

    # check default headers are in each df
    headersOK = sapply(manifest_data, function(x)
            all(c('Barcode', 'Time', 'Template') %in% colnames(x)))
    stopifnot(all(headersOK)) # --> need better error report

    cat_manifest_data = bind_rows(manifest_data)
    # check that barcodes are unique
    stopifnot(dim(cat_manifest_data)[1] == length(unique(cat_manifest_data$Barcode)))

    return(cat_manifest_data)
}


load_templates = function (template_file, log_file) {
    # manifest_file is a string or a vector of strings
    # log_file is an open file to sink errors and warnings

    # read sheets in files
    template_sheets = lapply(template_file, excel_sheets)

    for (iF in 1:length(template_file)) {
        # first check that the sheet names are ok
        # case of untreated plate
        if (length(template_sheets[[iF]])==1) {
            stopifnot(template_sheets[[iF]] == 'Gnumber')
            df = read_excel(template_file[[iF]], sheet = 'Gnumber', col_names = F)
            stopifnot(all(toupper(unlist(df)) %in% c('UNTREATED', 'VEHICLE')))
        } else {
        # normal case
            stopifnot('Gnumber' %in% template_sheets[[iF]])
            stopifnot('Concentration' %in% template_sheets[[iF]])
            # assess if multiple drugs and proper pairing
            n_drug = agrep('Gnumber', template_sheets[[iF]])
            n_conc = agrep('Concentration', template_sheets[[iF]])
            stopifnot(length(n_drug) == length(n_conc))
            if (length(n_drug)>1) {
                stopifnot(all(paste0('Gnumber_', 2:length(n_drug)) %in% template_sheets[[iF]]))
                stopifnot(all(paste0('Concentration_', 2:length(n_drug)) %in%
                    template_sheets[[iF]]))
            }
        }

        # read the different sheets and check for plate size
        # enforce range to avoid skipping empty rows at the beginning
        df = read_excel(template_file[[iF]], sheet = 'Gnumber', col_names = F, range = 'A1:AV32')
        # get the plate size
        n_row = 2**ceiling(log2(max(which(apply(!is.na(df), 1, any)))))
        n_col = 1.5*2**ceiling(log2(max(which(apply(!is.na(df), 2, any)))/1.5))
        n_row = max(n_row, n_col/1.5)
        n_col = max(1.5*n_row, n_col)
        plate_range = ifelse(n_col<26, paste0('A1:', LETTERS[n_col], n_row), 'A1:AV32')

        df_template = expand.grid(WellRow = LETTERS[1:n_row], WellColumn = 1:n_col)
        for (iS in template_sheets[[iF]]) {
            df = read_excel(template_file[[iF]], sheet = iS, col_names = F, range = plate_range)
            df$WellRow = LETTERS[1:n_row]
            df_melted = melt(df, id='WellRow', value.name = iS)
            colnames(df_melted)[colnames(df_melted) == 'variable'] = 'WellColumn'
            df_melted$WellColumn = gsub('X__', '', df_melted$WellColumn)
            df_template = merge(df_template, df_melted, by=c('WellRow', 'WellColumn'))
        }

    }


    df = read_excel(template_file[[iF]], sheet = 'Gnumber', col_names = F, range = 'A1:AV32')
                ##### this is an issue as the first rows, if empty, are discarded

    (x) {
        loadWorkbook(x)
        df = read.xlsx(x, skipEmptyCols=F)
    } )

    # check default headers are in each df
    headersOK = sapply(manifest_data, function(x)
            all(c('Barcode', 'Time', 'Template') %in% colnames(x)))
    stopifnot(all(headersOK)) # --> need better error report
