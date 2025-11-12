test_that("computeCppp works for basic cases", {
  pppCal <- c(0.1, 0.3, 0.8, 0.6)

  expect_equal(computeCppp(0.2, pppCal), 1/4)
  expect_equal(computeCppp(0.7, pppCal), 3/4)

  expect_true(computeCppp(0.2, pppCal) <= computeCppp(0.7, pppCal))

  expect_true(is.na(computeCppp(0.5, numeric())))
})
