test_that("dbern: correct probabilities at 0 and 1", {
  expect_equal(dbern(0, prob = 0.3), 0.7)
  expect_equal(dbern(1, prob = 0.3), 0.3)
  expect_equal(dbern(0, prob = 0.0), 1.0)
  expect_equal(dbern(1, prob = 1.0), 1.0)
})

test_that("dbern: zero for x outside {0,1}", {
  expect_equal(dbern(2, prob = 0.5), 0)
  expect_equal(dbern(-1, prob = 0.5), 0)
})

test_that("dbern: log argument works", {
  d <- dbern(1, prob = 0.3)
  d_log <- dbern(1, prob = 0.3, log = TRUE)
  expect_equal(d_log, log(d), tolerance = 1e-10)
})

test_that("pbern: CDF is correct", {
  expect_equal(pbern(0, prob = 0.3), 0.7)
  expect_equal(pbern(1, prob = 0.3), 1.0)
  expect_equal(pbern(-0.5, prob = 0.3), 0.0)
})

test_that("pbern: lower.tail works", {
  p1 <- pbern(0, prob = 0.3)
  p2 <- pbern(0, prob = 0.3, lower.tail = FALSE)
  expect_equal(p1 + p2, 1)
})

test_that("pbern: log.p works", {
  p <- pbern(1, prob = 0.3)
  expect_equal(pbern(1, prob = 0.3, log.p = TRUE), log(p))
})

test_that("qbern: quantile function works correctly", {
  # q(p) = 0 if p <= 1-prob, 1 if p > 1-prob
  expect_equal(qbern(0.5, prob = 0.3), 0)  # 0.5 <= 0.7
  expect_equal(qbern(0.8, prob = 0.3), 1)  # 0.8 > 0.7
  expect_equal(qbern(0.7, prob = 0.3), 0)  # exactly at boundary -> 0
})

test_that("rbern: returns only 0 and 1", {
  set.seed(42)
  r <- rbern(1000, prob = 0.4)
  expect_true(all(r %in% c(0, 1)))
  # Mean should be approximately prob
  expect_equal(mean(r), 0.4, tolerance = 0.08)
})

test_that("rbern: vectorized prob works", {
  set.seed(42)
  r <- rbern(n = 3, prob = c(0, 0.5, 1))
  expect_equal(r[1], 0)
  expect_equal(r[3], 1)
})
