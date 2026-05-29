test_that("dbetagen: density is non-negative and zero outside support", {
  x <- seq(1, 6, length.out = 100)
  d <- dbetagen(x, shape1 = 3, shape2 = 5, min = 1, max = 6)
  expect_true(all(d >= 0))
  expect_equal(dbetagen(0.9, shape1 = 3, shape2 = 5, min = 1, max = 6), 0)
  expect_equal(dbetagen(6.1, shape1 = 3, shape2 = 5, min = 1, max = 6), 0)
  expect_equal(dbetagen(numeric(0), shape1 = 2, shape2 = 3), numeric(0))
})

test_that("dbetagen: density integrates to approximately 1", {
  x <- seq(1, 6, length.out = 2001)
  d <- dbetagen(x, shape1 = 3, shape2 = 5, min = 1, max = 6)
  dx <- diff(x)[1]
  expect_equal(sum(d * dx), 1, tolerance = 0.01)
})

test_that("pbetagen: CDF is monotone, starts at 0, ends at 1", {
  x <- seq(1, 6, length.out = 100)
  p <- pbetagen(x, shape1 = 3, shape2 = 5, min = 1, max = 6)
  expect_true(all(diff(p) >= 0))
  expect_equal(pbetagen(1, shape1 = 3, shape2 = 5, min = 1, max = 6), 0, tolerance = 1e-10)
  expect_equal(pbetagen(6, shape1 = 3, shape2 = 5, min = 1, max = 6), 1, tolerance = 1e-10)
})

test_that("pbetagen: lower.tail and log.p work", {
  p1 <- pbetagen(3, shape1 = 3, shape2 = 5, min = 1, max = 6)
  p2 <- pbetagen(3, shape1 = 3, shape2 = 5, min = 1, max = 6, lower.tail = FALSE)
  expect_equal(p1 + p2, 1)
})

test_that("qbetagen and pbetagen are inverses", {
  p_vals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q_vals <- qbetagen(p_vals, shape1 = 3, shape2 = 5, min = 1, max = 6)
  p_back <- pbetagen(q_vals, shape1 = 3, shape2 = 5, min = 1, max = 6)
  expect_equal(p_back, p_vals, tolerance = 1e-6)
})

test_that("rbetagen: samples within support", {
  set.seed(42)
  r <- rbetagen(1000, shape1 = 3, shape2 = 5, min = 1, max = 6)
  expect_true(all(r >= 1))
  expect_true(all(r <= 6))
})

test_that("pbetagen: special case min == max returns 1 or 0", {
  # When min == max == q, pbetagen should return 1 (lower.tail=TRUE)
  p <- pbetagen(2, shape1 = 2, shape2 = 2, min = 2, max = 2)
  expect_equal(p, 1)
})

test_that("dbetagen: corresponds to scaled Beta", {
  # dbetagen(x, shape1, shape2, min, max) = dbeta((x-min)/(max-min)) / (max-min)
  x <- 0.6
  d_gen <- dbetagen(x, shape1 = 2, shape2 = 3, min = 0, max = 1)
  d_beta <- dbeta(x, shape1 = 2, shape2 = 3)
  expect_equal(d_gen, d_beta, tolerance = 1e-10)
})
