## MC workflow integration tests

test_that("mc: creates mc object from mcnodes", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(10); ndunc(5)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  x <- mcstoc(runif)
  y <- mcdata(3, type = "0")
  z <- x * y
  m <- mc(x, y, z, name = c("x", "y", "z"))
  expect_s3_class(m, "mc")
  expect_equal(length(m), 3)
  expect_equal(names(m), c("x", "y", "z"))
})

test_that("evalmcmod: evaluates mcmodel correctly", {
  mod <- mcmodel({
    x <- mcstoc(runif, type = "V")
    y <- mcdata(2, type = "0")
    z <- x * y
    mc(x, y, z)
  })

  res <- evalmcmod(mod, nsv = 50, nsu = 1, seed = 42)
  expect_s3_class(res, "mc")
  expect_equal(length(res), 3)
  expect_equal(dim(res$x), c(50, 1, 1))
})

test_that("summary.mc: returns correct structure", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(100); ndunc(50)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  set.seed(42)
  xVU <- mcstoc(rnorm, type = "VU")
  m <- mc(xVU)
  s <- summary(m)
  expect_s3_class(s, "summary.mc")
  expect_equal(length(s), 1)
  # VU node should have both rows and columns
  expect_true(nrow(s[[1]]) > 1)
})

test_that("quantile.mc: returns list of quantiles", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(100); ndunc(50)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  set.seed(42)
  xV <- mcstoc(runif, type = "V")
  m <- mc(xV)
  q <- quantile(m)
  expect_true(is.list(q))
  expect_equal(length(q), 1)
})

test_that("tornado: correlation values in [-1, 1]", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(100); ndunc(1)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  set.seed(42)
  x1 <- mcstoc(rnorm, type = "V")
  x2 <- mcstoc(rnorm, type = "V")
  x3 <- x1 + x2
  m <- mc(x1, x2, x3)
  tor <- tornado(m, output = 3)
  expect_s3_class(tor, "tornado")
  vals <- unlist(tor$value)
  expect_true(all(vals >= -1 & vals <= 1))
})

test_that("tornadounc: correlation in uncertainty dimension", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(50); ndunc(100)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  set.seed(42)
  x1 <- mcstoc(rnorm, type = "U")
  x2 <- mcstoc(rnorm, type = "U")
  x3 <- x1 + x2
  m <- mc(x1, x2, x3)
  tor <- tornadounc(m, output = 3)
  expect_s3_class(tor, "tornadounc")
  vals <- unlist(tor$value)
  expect_true(all(vals >= -1 & vals <= 1, na.rm = TRUE))
})

test_that("ndvar/ndunc: get and set work", {
  old_v <- ndvar()
  old_u <- ndunc()
  on.exit({ ndvar(old_v); ndunc(old_u) })

  ndvar(42)
  expect_equal(ndvar(), 42)
  ndunc(17)
  expect_equal(ndunc(), 17)
})

test_that("mcapply: apply function over mc dimensions", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(10); ndunc(5)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  set.seed(42)
  xVU <- mcstoc(rnorm, type = "VU")
  # Sum over uncertainty dimension
  result <- mcapply(xVU, "unc", sum)
  expect_s3_class(result, "mcnode")
  expect_equal(typemcnode(result), "V")
  expect_equal(dim(result), c(10, 1, 1))
})

test_that("mcratio: returns matrix with expected columns", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(50); ndunc(30)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  set.seed(42)
  xVU <- mcstoc(rnorm, type = "VU", mean = 5, sd = 1)
  m <- mc(xVU)
  r <- mcratio(m, na.rm = TRUE)
  expect_true(is.matrix(r))
  expect_true("VariabilityR" %in% colnames(r))
  expect_true("UncertaintyR" %in% colnames(r))
})

test_that("lhs: latin hypercube has correct dimensions", {
  set.seed(42)
  p <- lhs(runif, nsv = 10, nsu = 5)
  expect_equal(length(p), 50)
  # Values should be in (0,1)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("cornode: induces rank correlation", {
  set.seed(42)
  x1 <- rnorm(1000)
  x2 <- rnorm(1000)
  mat <- cbind(x1, x2)
  target <- matrix(c(1, 0.8, 0.8, 1), ncol = 2)
  mat_c <- cornode(mat, target = target)
  cor_result <- cor(mat_c, method = "spearman")
  expect_equal(cor_result[1, 2], 0.8, tolerance = 0.05)
})
