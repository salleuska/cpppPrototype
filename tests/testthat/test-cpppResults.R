# tests/testthat/test-cpppResult.R
test_that("new_cppresults creates a cpppResult with expected fields", {
  x <- new_cppresults(cppp = 0.4, ppp = runif(5), obs_ppp = 0.3)
  expect_s3_class(x, "cpppResult")
  expect_true(is.list(x))
  expect_equal(length(x$ppp), 5)
  expect_true(x$cppp >= 0 && x$cppp <= 1)
})

test_that("new_cppresults enforces bounds when finite", {
  expect_error(new_cppresults(cppp = 1.2), "in \\[0,1\\]")
  expect_error(new_cppresults(ppp = c(-0.1, 0.2)), "in \\[0,1\\]")
  expect_error(new_cppresults(obs_ppp = c(0.5, 1.1)), "in \\[0,1\\]")
})
