testthat::context("Test assay_to_df function")

testthat::test_that("Test assay_to_df works", {

seL_file <- system.file(package = "gDR", "testdata/seL.rds")
seL <- readRDS(seL_file)
se <- seL[[3]]
test1 = assay_to_df(se, "Avg_Controls")
test2 = assay_to_df(se, 4)

})