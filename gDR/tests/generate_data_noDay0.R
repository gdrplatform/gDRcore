
source('../R/analyze_data.R')
source('../R/GR_curve_fit.R')
# single agent ( matrix 2 drugs (9 doses) in triplicates on one plate
#   4 lines; last one without  ReferenceDivisionTime)
df_normalized = data.frame(
    CellLineName = as.vector(t(matrix(c('COV318', 'HCC2218', 'MCF-7', 'HCC1008'), 4, 60))),
            Tissue = as.vector(t(matrix(c('Ovary', 'Breast', 'Breast', 'Breast'), 4, 60))),
            Duration = 72,
            clid = as.vector(t(matrix(c('CL586056', 'CL131828', 'CL129757', 'CL586359'), 4, 60))),
            Barcode = 'Plate1',
            DrugName = as.vector(t(matrix(c('Palbociclib', 'Trametinib'), 12, 10))),
            Concentration = as.vector(matrix(c(0, 10**(seq(-3,1,.5))), 24, 10)),
            GRvalue = 1,
            RelativeViability = 1,
            ReadoutValue = 100,
            BackgroundValue = 1,
            Gnumber = as.vector(t(matrix(c('G02001876', 'G02442104'), 12, 10))),
            ReferenceDivisionTime = as.vector(t(matrix(c(57,85, 30, NA), 4, 60))), # from gneDB
            CorrectedReadout = 99,
            UntrtReadout = 99,
            Day0Readout = NA)

# apply the drug effect
ec50 = matrix(c(.08, 2, .5, .2, .4, 1, .04, .3),2,4)
colnames(ec50) = c('COV318', 'HCC2218', 'MCF-7', 'HCC1008')
rownames(ec50) = c('Palbociclib', 'Trametinib')

e_max = matrix(c(.62, .9, .2, .6, .1, .3, .5, .3),2,4)
colnames(e_max) = c('COV318', 'HCC2218', 'MCF-7', 'HCC1008')
rownames(e_max) = c('Palbociclib', 'Trametinib')

df_normalized$ReadoutValue = 100 * apply( df_normalized, 1, function(x)
    e_max[x['DrugName'],x['CellLineName']] + (1-e_max[x['DrugName'],x['CellLineName']])*
        (ec50[x['DrugName'],x['CellLineName']]**2 / (as.numeric(x['Concentration'])**2 +
            ec50[x['DrugName'],x['CellLineName']]**2)))

levels(df_normalized$Gnumber) = c(levels(df_normalized$Gnumber), 'Vehicle')
df_normalized$Gnumber[df_normalized$Concentration == 0] = 'Vehicle'
levels(df_normalized$DrugName) = c(levels(df_normalized$DrugName), 'Vehicle')
df_normalized$DrugName[df_normalized$Concentration == 0] = 'Vehicle'
# add some noise
set.seed(1)
df_normalized$ReadoutValue = df_normalized$ReadoutValue + runif(dim(df_normalized)[1], -2, 3)

# recalculate the normalized values

df_normalized$CorrectedReadout = df_normalized$ReadoutValue - df_normalized$BackgroundValue

df_ctrl = aggregate(df_normalized$CorrectedReadout[ df_normalized$Concentration == 0],
        by = list(clid = df_normalized$clid[ df_normalized$Concentration == 0]),
            function(x) mean(x, trim = .25))
for (cl in df_ctrl$clid) {
    df_normalized$UntrtReadout[df_normalized$clid == cl] = df_ctrl$x[df_ctrl$clid == cl]
}

df_normalized$GRvalue = round(2 ^ (1 + (
          log2(pmin(1.25,
                    df_normalized[, "RelativeViability"])) /
            (df_normalized$Duration / df_normalized$ReferenceDivisionTime)
        )), 4) - 1
df_normalized$RelativeViability = df_normalized$CorrectedReadout / df_normalized$UntrtReadout


# decompose de table into original files
df_raw_data = df_normalized[order(df_normalized$Barcode),]
df_raw_data$Template = 'Template_trt.tsv'
df_raw_data$WellRow = sort(rep(LETTERS[3:14],20))
df_raw_data$WellColumn = 3:22
levels(df_raw_data$Barcode) = paste0('Plate',1:8)

df_manifest = unique(df_raw_data[, c('Barcode', 'Template', 'Duration', 'clid')])
df_data = df_raw_data[,c('Barcode', 'WellRow', 'WellColumn', 'ReadoutValue', 'BackgroundValue')]
df_treatment = unique(df_raw_data[, c('Template', 'WellRow', 'WellColumn',
        'Gnumber', 'Concentration')])

dir.create('../inst/testdata/data11')
write.table(df_manifest, '../inst/testdata/data11/Manifest_data11.tsv',
        sep='\t', quote=F, row.names=F)
write.table(df_data, '../inst/testdata/data11/RawData_data11.tsv',
        sep='\t', quote=F, row.names=F)
for (trt_f in unique(df_treatment$Template)) {
        write.table(df_treatment[df_treatment$Template==trt_f,
                setdiff(colnames(df_treatment), 'Template')],
            paste0('../inst/testdata/data11/', trt_f),
                sep='\t', quote=F, row.names=F)
}
write.table(df_normalized, '../inst/testdata/data11/calculated_normalized_data.tsv',
            sep='\t', quote=F, row.names=F)
