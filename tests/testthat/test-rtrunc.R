test_that("rtrunc: samples within (linf, lsup]", {
  set.seed(42)
  r <- rtrunc("rnorm", n = 500, linf = 0, lsup = Inf)
  expect_true(all(r > 0))

  r2 <- rtrunc(rnorm, n = 500, linf = 3, lsup = 5, sd = 10)
  expect_true(all(r2 > 3))
  expect_true(all(r2 <= 5))
})

test_that("rtrunc: works with discrete distributions", {
  set.seed(42)
  r <- rtrunc(rpois, 1000, linf = 2, lsup = 4, lambda = 1)
  expect_true(all(r > 2))
  expect_true(all(r <= 4))
})

test_that("rtrunc: error when linf >= lsup", {
  expect_error(rtrunc("rnorm", n = 10, linf = 5, lsup = 3))
  expect_error(rtrunc("rnorm", n = 10, linf = 5, lsup = 5))
})

test_that("rtrunc: gamma distribution stays positive", {
  set.seed(42)
  r <- rtrunc("rgamma", n = 200, linf = 1, lsup = 5, shape = 2, rate = 1)
  expect_true(all(r > 1))
  expect_true(all(r <= 5))
})

test_that("rtrunc: works with function object", {
  set.seed(42)
  r <- rtrunc(rnorm, n = 100, linf = -2, lsup = 2)
  expect_true(all(r > -2))
  expect_true(all(r <= 2))
})

test_that("rtrunc: works with character function name", {
  set.seed(42)
  r <- rtrunc("rnorm", n = 100, linf = -2, lsup = 2)
  expect_true(all(r > -2))
  expect_true(all(r <= 2))
})
