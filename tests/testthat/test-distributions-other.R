## Tests for rdirichlet/ddirichlet

test_that("rdirichlet: rows sum to 1", {
  set.seed(42)
  r <- rdirichlet(100, alpha = c(1, 2, 3, 4))
  row_sums <- rowSums(r)
  expect_equal(row_sums, rep(1, 100), tolerance = 1e-10)
})

test_that("rdirichlet: values are in [0, 1]", {
  set.seed(42)
  r <- rdirichlet(100, alpha = c(1, 2, 3))
  expect_true(all(r >= 0))
  expect_true(all(r <= 1))
})

test_that("ddirichlet: density is non-negative", {
  alpha <- c(2, 3, 5)
  set.seed(42)
  x <- rdirichlet(10, alpha)
  d <- ddirichlet(x, alpha)
  expect_true(all(d >= 0))
})

test_that("rdirichlet: vectorized alpha (matrix)", {
  set.seed(42)
  alpha <- matrix(c(1,2,3,4, 4,3,2,1), nrow=2, byrow=TRUE)
  r <- rdirichlet(2, alpha)
  expect_equal(nrow(r), 2)
  expect_equal(rowSums(r), c(1,1), tolerance = 1e-10)
})

## Tests for dmultinomial/rmultinomial

test_that("rmultinomial: rows sum to size", {
  set.seed(42)
  r <- rmultinomial(10, size = 100, prob = c(1, 2, 7))
  expect_equal(rowSums(r), rep(100, 10))
})

test_that("dmultinomial: probabilities are valid", {
  x <- matrix(c(10, 20, 70), nrow = 1)
  prob <- c(1, 2, 7)
  d <- dmultinomial(x, prob = prob)
  # Should match dmultinom
  d_ref <- dmultinom(c(10, 20, 70), prob = prob)
  expect_equal(d, d_ref, tolerance = 1e-10)
})

test_that("rmultinomial: vectorized size works", {
  set.seed(42)
  r <- rmultinomial(3, size = c(10, 100, 1000), prob = c(1, 2, 7))
  expect_equal(rowSums(r), c(10, 100, 1000))
})

## Tests for rmultinormal/dmultinormal

test_that("rmultinormal: correct dimensions", {
  set.seed(42)
  sigma <- as.vector(diag(2))
  r <- rmultinormal(10, mean = c(0, 0), sigma = sigma)
  expect_equal(dim(r), c(10, 2))
})

test_that("dmultinormal: matches dmvnorm", {
  mean <- c(0, 0)
  sigma <- as.vector(diag(2))
  x <- matrix(c(1, 0), nrow = 1)
  d1 <- dmultinormal(x, mean = mean, sigma = sigma)
  d2 <- mvtnorm::dmvnorm(x, mean = mean, sigma = diag(2))
  expect_equal(as.numeric(d1), d2, tolerance = 1e-10)
})

## Tests for betasubjective
# dbetasubj(x, min, mode, mean, max, log=FALSE)

test_that("dbetasubj: density is non-negative", {
  x <- seq(0.1, 0.9, by = 0.1)
  d <- dbetasubj(x, min = 0, mode = 0.4, mean = 0.45, max = 1)
  expect_true(all(d >= 0))
})

test_that("pbetasubj and qbetasubj are inverses", {
  p_vals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q_vals <- qbetasubj(p_vals, min = 0, mode = 0.4, mean = 0.45, max = 1)
  p_back <- pbetasubj(q_vals, min = 0, mode = 0.4, mean = 0.45, max = 1)
  expect_equal(p_back, p_vals, tolerance = 1e-5)
})

test_that("rbetasubj: samples within support", {
  set.seed(42)
  r <- rbetasubj(500, min = 0, mode = 0.4, mean = 0.45, max = 1)
  expect_true(all(r >= 0))
  expect_true(all(r <= 1))
})

## Tests for MQI (minimum quantile information)
# dmqi(x, mqi, mqi.quantile=c(0.05, 0.5, 0.95), ...)
# mqi is a vector of quantile values (e.g., c(q05, q50, q95))

test_that("dmqi: density is non-negative", {
  set.seed(42)
  mqi_vals <- c(1, 3, 8)  # 5th, 50th, 95th percentile values
  x <- seq(0.5, 10, by = 0.2)
  d <- dmqi(x, mqi = mqi_vals)
  expect_true(all(d >= 0, na.rm = TRUE))
})

test_that("pmqi and qmqi are inverses", {
  mqi_vals <- c(1, 3, 8)
  p_vals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q_vals <- qmqi(p_vals, mqi = mqi_vals)
  p_back <- pmqi(q_vals, mqi = mqi_vals)
  expect_equal(p_back, p_vals, tolerance = 1e-5)
})

test_that("rmqi: generates numeric values", {
  set.seed(42)
  mqi_vals <- c(1, 3, 8)
  r <- rmqi(200, mqi = mqi_vals)
  expect_true(is.numeric(r))
  expect_equal(length(r), 200)
})
