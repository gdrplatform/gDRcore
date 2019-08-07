
source('../R/analyze_data.R')
source('../R/GR_curve_fit.R')
# drug combination ( matrix of 2 x 3 drugs (8 doses by 5 doses including 0) in triplicates, 2 lines)
df_normalized = data.frame( CellLineName = as.vector(t(matrix(c('COV318', 'HCC2218'), 2, 720))),
                    Tissue = as.vector(t(matrix(c('Ovary', 'Breast'), 2, 720))),
                    Duration = 72,
                    clid = as.vector(t(matrix(c('CL586056', 'CL131828'), 2, 720))),
                    Barcode = as.vector(t(matrix(paste0('Plate', 1:6), 6, 240))),
                    DrugName = as.vector(t(matrix(c('Palbociclib', 'Trametinib'), 12, 120))),
                    DrugName_2 = as.vector(t(matrix(c('GDC-0032', 'GDC-0941', 'GDC-0077'), 36,40))),
                    Concentration = as.vector(matrix(c(0, 10**(seq(-2.5,.5,.5))), 180, 8)),
                    Concentration_2=array(as.vector(t(matrix(c(0,10**(seq(-1,.5,.5))),5,8))),1440),
                    GRvalue = 1,
                    RelativeViability = 1,
                    DivisionTime = as.vector(t(matrix(c(58,90), 2, 720))),
                    ReadoutValue = 100,
                    BackgroundValue = 1,
                    Gnumber = as.vector(t(matrix(c('G02001876', 'G02442104'), 12, 120))),
                    Gnumber_2 = as.vector(t(matrix(c('G00505032','G00033829','G02967907'), 36,40))),
                    ReferenceDivisionTime = as.vector(t(matrix(c(60,85), 2, 720))),
                    CorrectedReadout = 99,
                    UntrtReadout = 99,
                    Day0Readout = as.vector(t(matrix(c(42,57), 2, 720))))

# apply the drug effect
GEC50 = matrix(c(.08, 2, .5, .2),2,2)
colnames(GEC50) = c('COV318', 'HCC2218')
rownames(GEC50) = c('Palbociclib', 'Trametinib')

GRmax = matrix(c(.62, .9, .75, .6),2,2)
colnames(GRmax) = c('COV318', 'HCC2218')
rownames(GRmax) = c('Palbociclib', 'Trametinib')

df_normalized$ReadoutValue = 100 * apply( df_normalized, 1, function(x)
    GRmax[x['DrugName'],x['CellLineName']] + (1-GRmax[x['DrugName'],x['CellLineName']])*
        (GEC50[x['DrugName'],x['CellLineName']]**2 / (as.numeric(x['Concentration'])**2 +
            GEC50[x['DrugName'],x['CellLineName']]**2)))

# apply the 2nd drug
GEC50_2 = matrix(c(.03, .2, 2,   .02, 3, .5),3,2)
colnames(GEC50_2) = c('COV318', 'HCC2218')
rownames(GEC50_2) = c('GDC-0032', 'GDC-0941', 'GDC-0077')

GRmax_2 = matrix(c(.75, .85, 1,  .7, .95, .8),3,2)
colnames(GRmax_2) = c('COV318', 'HCC2218')
rownames(GRmax_2) = c('GDC-0032', 'GDC-0941', 'GDC-0077')

df_normalized$ReadoutValue = df_normalized$ReadoutValue * apply( df_normalized, 1, function(x)
    GRmax_2[x['DrugName_2'],x['CellLineName']] + (1-GRmax_2[x['DrugName_2'],x['CellLineName']])*
        (GEC50_2[x['DrugName_2'],x['CellLineName']]**2 / (as.numeric(x['Concentration_2'])**2 +
            GEC50_2[x['DrugName_2'],x['CellLineName']]**2)))

# add some synergy/antagonism
idx = df_normalized$CellLineName == 'HCC2218' & df_normalized$DrugName == 'Trametinib' &
                    df_normalized$DrugName_2 == 'GDC-0941'
df_normalized$ReadoutValue[idx] = df_normalized$ReadoutValue[idx] +
                            .2 * (df_normalized$ReadoutValue[idx] - 100)


idx = df_normalized$CellLineName == 'HCC2218' & df_normalized$DrugName == 'Trametinib' &
                                df_normalized$DrugName_2 == 'GDC-0077'
df_normalized$ReadoutValue[idx] = df_normalized$ReadoutValue[idx] -
                            .1 * (df_normalized$ReadoutValue[idx] - 100)

idx = df_normalized$CellLineName == 'COV318' & df_normalized$DrugName == 'Palbociclib'
df_normalized$ReadoutValue[idx] = df_normalized$ReadoutValue[idx] +
                            .1 * (df_normalized$ReadoutValue[idx] - 100)

idx = df_normalized$CellLineName == 'COV318' & df_normalized$DrugName_2 != 'GDC-0032'
df_normalized$ReadoutValue[idx] = df_normalized$ReadoutValue[idx] -
                            .11 * (df_normalized$ReadoutValue[idx] - 100)

df_normalized$ReadoutValue = pmin(df_normalized$ReadoutValue, 100)

# add some noise
set.seed(1)
df_normalized$ReadoutValue = df_normalized$ReadoutValue + runif(dim(df_normalized)[1], -2, 3)

# recalculate the GR values
df_normalized$CorrectedReadout = df_normalized$ReadoutValue - df_normalized$BackgroundValue
df_normalized$GRvalue = 2 ** ( log2( df_normalized$CorrectedReadout / df_normalized$Day0Readout) /
            log2( df_normalized$UntrtReadout / df_normalized$Day0Readout)) - 1
df_normalized$RelativeViability = df_normalized$CorrectedReadout / df_normalized$UntrtReadout

levels(df_normalized$Gnumber) = c(levels(df_normalized$Gnumber), 'Vehicle')
df_normalized$Gnumber[df_normalized$Concentration == 0] = 'Vehicle'
levels(df_normalized$DrugName) = c(levels(df_normalized$DrugName), 'Vehicle')
df_normalized$DrugName[df_normalized$Concentration == 0] = 'Vehicle'
levels(df_normalized$Gnumber_2) = c(levels(df_normalized$Gnumber_2), 'Vehicle')
df_normalized$Gnumber_2[df_normalized$Concentration_2 == 0] = 'Vehicle'
levels(df_normalized$DrugName_2) = c(levels(df_normalized$DrugName_2), 'Vehicle')
df_normalized$DrugName_2[df_normalized$Concentration_2 == 0] = 'Vehicle'

Keys = identify_keys(df_normalized)
df_averaged = average_replicates(df_normalized[df_normalized$Concentration > 0 |
                                        df_normalized$Concentration_2 > 0, ] , Keys$Trt)
df_metrics = calculate_DRmetrics(df_averaged, Keys$DoseResp)

write.table(df_normalized, '../inst/testdata/test_combo1/df_normalized.tsv',
        sep='\t', quote=F, row.names=F)
write.table(df_averaged, '../inst/testdata/test_combo1/df_averaged.tsv',
        sep='\t', quote=F, row.names=F)
write.table(df_metrics, '../inst/testdata/test_combo1/df_metrics.tsv',
        sep='\t', quote=F, row.names=F)


# add a third co-treatment

df_normalized2 = df_normalized3 = df_normalized
df_normalized2$DrugName_3 = 'Vehicle'
df_normalized2$Gnumber_3 = 'Vehicle'
df_normalized2$Concentration_3 = 0
df_normalized3$DrugName_3 = 'Erlotinib'
df_normalized3$Gnumber_3 = 'G00022086'
df_normalized3$Concentration_3 = .5

df_normalized3$ReadoutValue = pmax(df_normalized3$ReadoutValue * .92 - 5, 4)
df_normalized3$ReadoutValue[df_normalized3$CellLineName == 'HCC2218'] =
        df_normalized3$ReadoutValue[df_normalized3$CellLineName == 'HCC2218'] * 1.02 + 2
set.seed(2)
df_normalized3$ReadoutValue = df_normalized3$ReadoutValue + runif(dim(df_normalized)[1], -1, 2)

df_normalized_3drugs = rbind(df_normalized2, df_normalized3)
Keys = identify_keys(df_normalized_3drugs)
df_averaged = average_replicates(df_normalized_3drugs[df_normalized_3drugs$Concentration > 0 |
            df_normalized_3drugs$Concentration_2 > 0 | df_normalized_3drugs$Concentration_3 > 0, ],
         Keys$Trt)
df_metrics = calculate_DRmetrics(df_averaged, Keys$DoseResp)

write.table(df_normalized, '../inst/testdata/test_combo_3d/df_normalized.tsv',
        sep='\t', quote=F, row.names=F)
write.table(df_averaged, '../inst/testdata/test_combo_3d/df_averaged.tsv',
        sep='\t', quote=F, row.names=F)
write.table(df_metrics, '../inst/testdata/test_combo_3d/df_metrics.tsv',
        sep='\t', quote=F, row.names=F)
