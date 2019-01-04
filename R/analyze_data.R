
load_all_data = function(manifest_file, template_file, results_file, log_file) {

    manifest = load_manifest(manifest_file, log_file)
    stopifnot(all(unique(manifest$Template) %in% basename(template_file)))
    treatments = load_templates(template_file, log_file)
    data = load_results(results_file, log_file)

    # need to merge manifest and treatments first and sort out duplicate metadata columns

    df_merge1 = merge(manifest, data, by = c('Barcode'))
    stopifnot( dim(df_merge1)[1] == dim(data)[1]) # need proper error message

    df_merge2 = merge(df_merge1, treatments, by = c('Template', 'WellRow', 'WellColumn'))
    stopifnot( dim(df_merge2)[1] == dim(data)[1]) # need proper error message



    # remove wells not labeled
    df_raw_data = df_merge2[ !is.na(df_merge2$Gnumber), ]
    print(sprintf('%i well loaded, %i discarded for lack of annotation, %i data point selected',
                dim(data)[1], sum(is.na(df_merge2$Gnumber)), dim(df_raw_data)[1]))

    return(df_raw_data)
}
