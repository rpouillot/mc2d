## mcnode creation and basic properties

test_that("mcdata: creates node with correct dimensions", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(10); ndunc(5)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  x0 <- mcdata(1, type = "0")
  expect_equal(dim(x0), c(1, 1, 1))
  expect_equal(typemcnode(x0), "0")

  xV <- mcdata(1:10, type = "V")
  expect_equal(dim(xV), c(10, 1, 1))
  expect_equal(typemcnode(xV), "V")

  xU <- mcdata(1:5, type = "U")
  expect_equal(dim(xU), c(1, 5, 1))
  expect_equal(typemcnode(xU), "U")

  xVU <- mcdata(1:50, type = "VU")
  expect_equal(dim(xVU), c(10, 5, 1))
  expect_equal(typemcnode(xVU), "VU")
})

test_that("mcdata: multivariate node dimensions", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(5); ndunc(3)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  xVM <- mcdata(1:10, type = "V", nvariates = 2)
  expect_equal(dim(xVM), c(5, 1, 2))
  expect_equal(typemcnode(xVM), "V")
})

test_that("dimmcnode: returns named vector", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(7); ndunc(4)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  xV <- mcdata(1:7, type = "V")
  d <- dimmcnode(xV)
  expect_equal(d["nsv"], c(nsv = 7L))
  expect_equal(d["nsu"], c(nsu = 1L))
})

test_that("dimmc: returns max dimensions across mc", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(10); ndunc(5)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  x1 <- mcdata(1:10, type = "V")
  x2 <- mcdata(1:50, type = "VU")
  m <- mc(x1, x2)
  d <- dimmc(m)
  expect_equal(d["nsv"], c(nsv = 10L))
  expect_equal(d["nsu"], c(nsu = 5L))
})

test_that("is.mcnode and is.mc: type checking", {
  x <- mcdata(1, type = "0")
  expect_true(is.mcnode(x))
  expect_false(is.mc(x))

  m <- mc(x)
  expect_true(is.mc(m))
  expect_false(is.mcnode(m))

  expect_false(is.mcnode(42))
  expect_false(is.mc(42))
})

test_that("Ops.mcnode: arithmetic preserves class and type", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(5); ndunc(3)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  xV <- mcdata(1:5, type = "V")
  result <- xV * 2
  expect_s3_class(result, "mcnode")
  expect_equal(typemcnode(result), "V")

  x0 <- mcdata(3, type = "0")
  result2 <- xV + x0
  expect_s3_class(result2, "mcnode")
  expect_equal(typemcnode(result2), "V")
})

test_that("Ops.mcnode: V + U = VU", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(5); ndunc(3)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  xV <- mcdata(1:5, type = "V")
  xU <- mcdata(1:3, type = "U")
  result <- xV + xU
  expect_equal(typemcnode(result), "VU")
})

test_that("is.na.mcnode: returns logical mcnode", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(3); ndunc(1)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  x <- mcdata(c(1, NA, 3), type = "V")
  na_x <- is.na(x)
  expect_s3_class(na_x, "mcnode")
  expect_true(is.logical(unmc(na_x)))
  expect_equal(as.vector(unmc(na_x)), c(FALSE, TRUE, FALSE))
})

test_that("is.finite.mcnode: handles Inf and NaN", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(4); ndunc(1)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  x <- mcdata(c(1, Inf, NaN, NA), type = "V")
  fin_x <- is.finite(x)
  expect_equal(as.vector(unmc(fin_x)), c(TRUE, FALSE, FALSE, FALSE))

  inf_x <- is.infinite(x)
  expect_equal(as.vector(unmc(inf_x)), c(FALSE, TRUE, FALSE, FALSE))
})

test_that("mcstoc: creates stochastic node of correct type", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(10); ndunc(5)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  set.seed(42)
  xV <- mcstoc(rnorm, type = "V")
  expect_s3_class(xV, "mcnode")
  expect_equal(typemcnode(xV), "V")
  expect_equal(dim(xV), c(10, 1, 1))

  xVU <- mcstoc(runif, type = "VU")
  expect_equal(typemcnode(xVU), "VU")
  expect_equal(dim(xVU), c(10, 5, 1))
})

test_that("pmin.mcnode and pmax.mcnode work", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(10); ndunc(1)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  set.seed(42)
  x <- mcstoc(rnorm, type = "V", nsv = 10)
  y <- pmin(x, 0)
  expect_true(all(unmc(y) <= 0))

  z <- pmax(x, 0)
  expect_true(all(unmc(z) >= 0))
})

test_that("unmc: removes mc class and returns array/list", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(5); ndunc(1)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  x <- mcdata(1:5, type = "V")
  u <- unmc(x, drop = TRUE)
  expect_false(inherits(u, "mcnode"))
  expect_equal(as.numeric(u), 1:5)
})

test_that("outm: changes outm attribute", {
  old_v <- ndvar(); old_u <- ndunc()
  ndvar(3); ndunc(1)
  on.exit({ ndvar(old_v); ndunc(old_u) })

  x <- mcdata(1:6, type = "V", nvariates = 2)
  x2 <- outm(x, "mean")
  expect_equal(attr(x2, "outm"), "mean")
})
