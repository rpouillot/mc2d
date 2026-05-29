test_that("dlnormb: default params give meanlog=0, sdlog=1", {
  # Default mean=exp(0.5), sd=sqrt(exp(2)-exp(1)) gives meanlog=0, sdlog=1
  expect_equal(dlnormb(1), dlnorm(1), tolerance = 1e-10)
  expect_equal(dlnormb(1), dnorm(0), tolerance = 1e-10)
})

test_that("dlnormb: density is non-negative", {
  x <- seq(0.1, 20, by = 0.1)
  d <- dlnormb(x, mean = 5, sd = 3)
  expect_true(all(d >= 0))
  expect_equal(dlnormb(0, mean = 5, sd = 3), 0)
})

test_that("dlnormb: log argument works", {
  d <- dlnormb(2, mean = 5, sd = 3)
  d_log <- dlnormb(2, mean = 5, sd = 3, log = TRUE)
  expect_equal(d_log, log(d), tolerance = 1e-10)
})

test_that("plnormb and qlnormb are inverses", {
  x_vals <- c(1, 2, 5, 10)
  p_vals <- plnormb(x_vals, mean = 5, sd = 3)
  x_back <- qlnormb(p_vals, mean = 5, sd = 3)
  expect_equal(x_back, x_vals, tolerance = 1e-6)
})

test_that("plnormb: lower.tail and log.p work", {
  p1 <- plnormb(5, mean = 5, sd = 3)
  p2 <- plnormb(5, mean = 5, sd = 3, lower.tail = FALSE)
  expect_equal(p1 + p2, 1)
  expect_equal(plnormb(5, mean = 5, sd = 3, log.p = TRUE), log(p1))
})

test_that("rlnormb: sample mean and sd close to specified values", {
  set.seed(42)
  r <- rlnormb(1e5, mean = 3, sd = 1)
  expect_equal(mean(r), 3, tolerance = 0.05)
  expect_equal(sd(r), 1, tolerance = 0.05)
})

test_that("dlnormb: density integrates to approximately 1", {
  x <- seq(0.001, 50, length.out = 5000)
  d <- dlnormb(x, mean = 5, sd = 3)
  dx <- diff(x)[1]
  expect_equal(sum(d * dx), 1, tolerance = 0.01)
})
