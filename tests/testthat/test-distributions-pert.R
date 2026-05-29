test_that("dpert: density is non-negative and zero outside support", {
  x <- seq(3, 10, length.out = 100)
  d <- dpert(x, min = 3, mode = 5, max = 10)
  expect_true(all(d >= 0))
  expect_equal(dpert(2.9, min = 3, mode = 5, max = 10), 0)
  expect_equal(dpert(10.1, min = 3, mode = 5, max = 10), 0)
  expect_equal(dpert(numeric(0), min = 3, mode = 5, max = 10), numeric(0))
})

test_that("dpert: density integrates to approximately 1", {
  x <- seq(3, 10, length.out = 2001)
  d <- dpert(x, min = 3, mode = 5, max = 10)
  dx <- diff(x)[1]
  expect_equal(sum(d * dx), 1, tolerance = 0.01)
})

test_that("ppert: CDF is monotone, starts at 0, ends at 1", {
  x <- seq(3, 10, length.out = 100)
  p <- ppert(x, min = 3, mode = 5, max = 10)
  expect_true(all(diff(p) >= 0))
  expect_equal(ppert(3, min = 3, mode = 5, max = 10), 0)
  expect_equal(ppert(10, min = 3, mode = 5, max = 10), 1)
})

test_that("ppert: lower.tail and log.p work correctly", {
  p1 <- ppert(5, min = 3, mode = 5, max = 10)
  p2 <- ppert(5, min = 3, mode = 5, max = 10, lower.tail = FALSE)
  expect_equal(p1 + p2, 1)
  expect_equal(ppert(5, min = 3, mode = 5, max = 10, log.p = TRUE), log(p1))
})

test_that("qpert and ppert are inverses", {
  p_vals <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  q_vals <- qpert(p_vals, min = 3, mode = 5, max = 10)
  p_back <- ppert(q_vals, min = 3, mode = 5, max = 10)
  expect_equal(p_back, p_vals, tolerance = 1e-6)
})

test_that("rpert: samples within support", {
  set.seed(42)
  r <- rpert(1000, min = 3, mode = 5, max = 10)
  expect_true(all(r >= 3))
  expect_true(all(r <= 10))
})

test_that("dpert: alternative mean parametrization", {
  # mean = (min + 4*mode + max) / 6 for shape=4 PERT
  mn <- (3 + 4 * 5 + 10) / 6
  d1 <- dpert(6, min = 3, mode = 5, max = 10)
  d2 <- dpert(6, min = 3, mean = mn, max = 10)
  expect_equal(d1, d2, tolerance = 1e-6)
})

test_that("dpert/rpert: error when both mode and mean given", {
  expect_error(dpert(5, min = 3, mode = 5, mean = 5, max = 10))
  expect_error(rpert(10, min = 3, mode = 5, mean = 5, max = 10))
})

test_that("dpert: shape parameter changes distribution", {
  # Shape=6 (more peaked) vs shape=4 (default)
  d4 <- dpert(5, min = 3, mode = 5, max = 10, shape = 4)
  d6 <- dpert(5, min = 3, mode = 5, max = 10, shape = 6)
  expect_false(isTRUE(all.equal(d4, d6)))
})

test_that("dpert: min == max returns NaN (degenerate)", {
  expect_true(is.nan(dpert(2, min = 2, mode = 2, max = 2)))
})

test_that("dpert: log argument works", {
  d <- dpert(5, min = 3, mode = 5, max = 10)
  d_log <- dpert(5, min = 3, mode = 5, max = 10, log = TRUE)
  expect_equal(d_log, log(d), tolerance = 1e-10)
})
