test_that("the synergy fits are announced to the gDRutils registry", {
  skip_if_not(
    exists("get_fit_fns", envir = asNamespace("gDRutils"), inherits = FALSE),
    "gDRutils predates the fit function registry"
  )
  combo <- gDRutils::get_fit_fns(data_type = "combination")
  expect_true(all(c("bliss_fit_fn", "hss_fit_fn", "musyc_fit_fn") %in%
                    combo$name))
  expect_true(all(combo[combo$name %in%
                          c("bliss_fit_fn", "hss_fit_fn", "musyc_fit_fn"),
                        ]$package == "gDRcore"))
})


test_that("a registered name resolves to the function itself", {
  skip_if_not(
    exists("get_fit_fn", envir = asNamespace("gDRutils"), inherits = FALSE),
    "gDRutils predates the fit function registry"
  )
  expect_identical(gDRutils::get_fit_fn("bliss_fit_fn"), bliss_fit_fn)
  expect_identical(gDRutils::get_fit_fn("hss_fit_fn"), hss_fit_fn)
})


test_that("MuSyC availability tracks musyc_available()", {
  skip_if_not(
    exists("get_fit_fns", envir = asNamespace("gDRutils"), inherits = FALSE),
    "gDRutils predates the fit function registry"
  )
  row <- gDRutils::get_fit_fns()[
    gDRutils::get_fit_fns()$name == "musyc_fit_fn", ]
  # Asserted against the predicate rather than against TRUE: whether Python
  # synergy is provisioned differs between a developer machine and CI, and the
  # point is that the registry reports it rather than what it reports.
  expect_identical(row$available, isTRUE(musyc_available()))
})


test_that("registration is skipped, not fatal, on an older gDRutils", {
  # The DESCRIPTION floor is below the release that added the registry, so
  # loading gDRcore must survive its absence.
  expect_type(gDRcore:::.fit_fn_registry_available(), "logical")
  expect_silent(gDRcore:::.register_fit_fns())
})
