test_that("paste_warnings works as expected", {
  fun1 <- function(x) {
    warning("warning 1")
    warning("warning 2")
    10
  }
  funOutput <- purrr::quietly(fun1)()
  expect_warning(paste_warnings(funOutput$warnings), regexp = "warning 1\nwarning 2")
})
