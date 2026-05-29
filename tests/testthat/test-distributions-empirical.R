## Tests for empiricalD

test_that("dempiricalD: probabilities sum correctly", {
  vals <- 2:6
  prob <- c(10, 10, 70, 0, 10)
  # Normalized: 10/100, 10/100, 70/100, 0/100, 10/100
  d <- dempiricalD(2:6, values = vals, prob = prob)
  expect_equal(sum(d), 1, tolerance = 1e-10)
  expect_equal(d[3], 0.7, tolerance = 1e-10)  # value 4 has prob 70/100
})

test_that("dempiricalD: zero for values not in support", {
  d <- dempiricalD(1, values = 2:5, prob = c(1, 1, 1, 1))
  expect_equal(d, 0)
})

test_that("pempiricalD: CDF is monotone and valid", {
  vals <- 1:5
  p <- pempiricalD(1:6, values = vals)
  expect_true(all(diff(p[1:5]) >= 0))
  expect_equal(p[5], 1)
  expect_equal(p[6], 1)  # beyond max
})

test_that("qempiricalD: quantile function works", {
  vals <- c(10, 20, 30)
  prob <- c(1, 2, 1)  # probs 0.25, 0.50, 0.25
  q <- qempiricalD(0.5, values = vals, prob = prob)
  expect_equal(q, 20)
})

test_that("rempiricalD: samples within values", {
  set.seed(42)
  vals <- c(1, 2, 5, 10)
  r <- rempiricalD(1000, values = vals)
  expect_true(all(r %in% vals))
})

test_that("rempiricalD: equivalent to sample when prob given as vector", {
  set.seed(123)
  r1 <- rempiricalD(100, values = 1:3, prob = c(1, 2, 3))
  set.seed(123)
  r2 <- sample(1:3, 100, replace = TRUE, prob = c(1, 2, 3))
  expect_equal(r1, r2)
})

test_that("dempiricalD: empty input returns empty", {
  d <- dempiricalD(numeric(0), values = 1:3)
  expect_equal(d, numeric(0))
})

## Tests for empiricalC

test_that("dempiricalC: density is non-negative", {
  prob <- c(2, 3, 1, 6, 1)
  values <- 1:5
  x <- seq(0, 6, length.out = 200)
  d <- dempiricalC(x, min = 0, max = 6, values = values, prob = prob)
  expect_true(all(d >= 0, na.rm = TRUE))
})

test_that("dempiricalC: density integrates to approximately 1", {
  prob <- c(2, 3, 1, 6, 1)
  values <- 1:5
  x <- seq(0, 6, length.out = 2001)
  d <- dempiricalC(x, min = 0, max = 6, values = values, prob = prob)
  dx <- diff(x)[1]
  expect_equal(sum(d * dx, na.rm = TRUE), 1, tolerance = 0.02)
})

test_that("dempiricalC: zero outside [min, max]", {
  prob <- c(1, 1, 1)
  values <- 1:3
  expect_equal(dempiricalC(-1, min = 0, max = 4, values = values, prob = prob), 0)
  expect_equal(dempiricalC(5, min = 0, max = 4, values = values, prob = prob), 0)
})

test_that("pempiricalC and qempiricalC are inverses", {
  prob <- c(2, 3, 1, 6, 1)
  values <- 1:5
  p_vals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q_vals <- qempiricalC(p_vals, min = 0, max = 6, values = values, prob = prob)
  p_back <- pempiricalC(q_vals, min = 0, max = 6, values = values, prob = prob)
  expect_equal(p_back, p_vals, tolerance = 1e-5)
})

test_that("rempiricalC: samples within [min, max]", {
  set.seed(42)
  r <- rempiricalC(500, min = 0, max = 6, values = 1:5, prob = c(1,1,1,1,1))
  expect_true(all(r >= 0))
  expect_true(all(r <= 6))
})
