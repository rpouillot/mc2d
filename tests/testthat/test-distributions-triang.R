test_that("dtriang: density is non-negative and zero outside support", {
  x <- seq(3, 10, length.out = 100)
  d <- dtriang(x, min = 3, mode = 6, max = 10)
  expect_true(all(d >= 0))
  expect_equal(dtriang(2.9, min = 3, mode = 6, max = 10), 0)
  expect_equal(dtriang(10.1, min = 3, mode = 6, max = 10), 0)
  expect_equal(dtriang(numeric(0), min = 3, mode = 6, max = 10), numeric(0))
})

test_that("dtriang: density integrates to approximately 1", {
  x <- seq(3, 10, length.out = 2001)
  d <- dtriang(x, min = 3, mode = 6, max = 10)
  dx <- diff(x)[1]
  expect_equal(sum(d * dx), 1, tolerance = 0.01)
})

test_that("dtriang: NaN when min == mode == max", {
  expect_true(is.nan(suppressWarnings(dtriang(2, min = 2, mode = 2, max = 2))))
})

test_that("ptriang: CDF is monotone, starts at 0, ends at 1", {
  x <- seq(3, 10, length.out = 100)
  p <- ptriang(x, min = 3, mode = 6, max = 10)
  expect_true(all(diff(p) >= 0))
  expect_equal(ptriang(3, min = 3, mode = 6, max = 10), 0)
  expect_equal(ptriang(10, min = 3, mode = 6, max = 10), 1)
})

test_that("ptriang: lower.tail and log.p work", {
  p1 <- ptriang(6, min = 3, mode = 6, max = 10)
  p2 <- ptriang(6, min = 3, mode = 6, max = 10, lower.tail = FALSE)
  expect_equal(p1 + p2, 1)
  expect_equal(ptriang(6, min = 3, mode = 6, max = 10, log.p = TRUE), log(p1))
})

test_that("qtriang and ptriang are inverses", {
  p_vals <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  q_vals <- qtriang(p_vals, min = 3, mode = 6, max = 10)
  p_back <- ptriang(q_vals, min = 3, mode = 6, max = 10)
  expect_equal(p_back, p_vals, tolerance = 1e-6)
})

test_that("rtriang: samples within support", {
  set.seed(42)
  r <- rtriang(1000, min = 3, mode = 6, max = 10)
  expect_true(all(r >= 3))
  expect_true(all(r <= 10))
})

test_that("dtriang: mode parametrization via mean", {
  # For triangular, mean = (min + mode + max) / 3
  mn <- (3 + 6 + 10) / 3
  d1 <- dtriang(5, min = 3, mode = 6, max = 10)
  d2 <- dtriang(5, min = 3, mean = mn, max = 10)
  expect_equal(d1, d2, tolerance = 1e-6)
})

test_that("dtriang: error when both mode and mean given", {
  expect_error(dtriang(5, min = 3, mode = 6, mean = 6, max = 10))
})

test_that("dtriang: symmetric triangular has density peaking at mode", {
  d_mode <- dtriang(0, min = -1, mode = 0, max = 1)
  d_off <- dtriang(0.5, min = -1, mode = 0, max = 1)
  expect_true(d_mode > d_off)
})
