# tests/testthat/test-cpppResult.R
test_that("new_cppresults creates a cpppResult with expected fields", {
  x <- new_cppresults(CPPP = 0.4, repPPP = runif(5), obsPPP = 0.3)
  expect_s3_class(x, "cpppResult")
  expect_true(is.list(x))
  expect_equal(length(x$repPPP), 5)
  expect_true(x$CPPP >= 0 && x$CPPP <= 1)
})

test_that("new_cppresults enforces bounds when finite", {
  expect_error(new_cppresults(CPPP = 1.2), "in \\[0,1\\]")
  expect_error(new_cppresults(repPPP = c(-0.1, 0.2)), "in \\[0,1\\]")
  expect_error(new_cppresults(obsPPP = c(0.5, 1.1)), "in \\[0,1\\]")
})
