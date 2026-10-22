## Regression tests for the bugs fixed in mc2d 0.2.3

test_that("Ops.mcnode accepts arrays for V, U and 0 nodes", {
  with_dims(7, 5, {
    set.seed(1)
    for (ty in c("0", "V", "U", "VU")) {
      x <- mcstoc(rnorm, type = ty, nvariates = 3)
      d <- dim(x)
      expect_identical(x + array(1, d), x + 1)
      expect_identical(array(2, d) * x, 2 * x)
      expect_identical(x + array(1, c(d[1:2], 1)), x + 1)
      expect_identical(x + matrix(1, d[1], d[2]), x + 1)
    }
    x0 <- mcdata(1, "0")
    expect_error(x0 + matrix(1, 1, 3), "array size")
  })
})

test_that("Ops.mcnode handles a vector outm", {
  with_dims(7, 5, {
    m <- mcstoc(rnorm, "V", outm = c("min", "max"))
    b <- mcdata(2, "V", outm = "mean")
    e <- mcdata(2, "V", outm = "each")
    expect_identical(attr(m + b, "outm"), c("min", "max"))
    expect_identical(attr(b + m, "outm"), "mean")
    expect_identical(attr(m + e, "outm"), "each")
    expect_identical(unclass(m + b), unclass(m) + 2)
  })
})

test_that("mcdata: vector outm and 1 x k matrix for a '0' node", {
  with_dims(7, 5, {
    x <- mcdata(1, type = "V", outm = c("min", "max"))
    expect_identical(attr(x, "outm"), c("min", "max"))
    expect_error(mcdata(matrix(1:5, nrow = 1), type = "0"), "array size")
  })
})

test_that("mcstoc: rtrunc and lhs with the default function or pkg::fun", {
  with_dims(50, 3, {
    set.seed(1)
    x <- mcstoc(type = "V", rtrunc = TRUE, linf = .5)
    expect_true(all(x > .5 & x <= 1))
    y <- mcstoc(stats::rnorm, type = "V", rtrunc = TRUE, linf = 0)
    expect_true(all(y > 0))
    z <- mcstoc(type = "VU", lhs = TRUE)
    # one value per stratum in each column
    expect_true(all(apply(unclass(z)[, , 1], 2, function(u) all(sort(floor(u * 50)) == 0:49))))
  })
})

test_that("pmin / pmax with a single mcnode and incompatible variates", {
  with_dims(7, 5, {
    x <- mcstoc(rnorm, "V")
    expect_identical(pmin(x), x)
    expect_identical(pmax(x), x)
    expect_error(pmin(mcstoc(rnorm, "VU", nvariates = 3), mcdata(1, "V", nvariates = 2)), "Incompatible")
  })
})

test_that("pmin / pmax fast path gives the same results as the R loop", {
  with_dims(11, 4, {
    set.seed(2)
    a <- mcstoc(rnorm, "VU", nvariates = 2); b <- mcstoc(rnorm, "V"); u <- mcstoc(rnorm, "U", nvariates = 2)
    a[c(1, 5, 30)] <- NA; a[7] <- NaN; b[c(1, 2)] <- NaN; b[3] <- NA; u[2] <- NA
    old <- function(..., na.rm = FALSE, min = TRUE) {
      argsd <- lapply(list(...), function(z) as.vector(unclass(mcdata(z, "VU", nvariates = 2))))
      mmm <- argsd[[1L]]; has.na <- FALSE
      for (each in argsd[-1L]) {
        nas <- cbind(is.na(mmm), is.na(each))
        if (has.na || (has.na <- any(nas))) { mmm[nas[, 1L]] <- each[nas[, 1L]]; each[nas[, 2L]] <- mmm[nas[, 2L]] }
        change <- if (min) mmm > each else mmm < each
        change <- change & !is.na(change)
        mmm[change] <- each[change]
        if (has.na && !na.rm) mmm[nas[, 1L] | nas[, 2L]] <- NA
      }
      mmm
    }
    for (narm in c(FALSE, TRUE)) {
      expect_identical(as.vector(unclass(pmin(a, b, u, na.rm = narm))), old(a, b, u, na.rm = narm))
      expect_identical(as.vector(unclass(pmax(b, a, u, na.rm = narm))), old(b, a, u, na.rm = narm, min = FALSE))
    }
  })
})

test_that("mcapply on a dimension of length 1", {
  with_dims(7, 5, {
    x <- mcstoc(rnorm, "V")
    expect_identical(unclass(mcapply(x, "unc", sum)), unclass(x))
    expect_identical(unclass(mcapply(x, "variates", exp)), exp(unclass(x)))
    u <- mcstoc(rnorm, "U")
    expect_identical(unclass(mcapply(u, "variates", sum)), unclass(u))
  })
})

test_that("mcprobtree: local variables, list given as a variable, single branch", {
  with_dims(7, 5, {
    set.seed(3)
    a <- mcstoc(rnorm, "VU", mean = 10); b <- mcdata(0, "VU")
    set.seed(4); r1 <- mcprobtree(c(.4, .6), list("1" = a, "2" = b), type = "VU")
    vals <- list("1" = a, "2" = b)
    set.seed(4); r2 <- mcprobtree(c(.4, .6), vals, type = "VU")
    expect_identical(r1, r2)
    expect_true(all(r1 == 0 | r1 == a))
    expect_identical(unclass(mcprobtree(1, list("5" = a), type = "VU")), unclass(a))
  })
})

test_that("cornode rejects a 'V' node not in first position", {
  with_dims(20, 3, {
    expect_error(cornode(mcstoc(rnorm, "VU"), mcstoc(rnorm, "V"), target = .5), "first position")
  })
})

test_that("is.mcnode and outm", {
  with_dims(7, 5, {
    x <- mcdata(1:5, "V", nsv = 5); attr(x, "type") <- "0"
    expect_false(is.mcnode(x))
    y <- mcdata(1, "V"); attr(y, "type") <- "X"
    expect_false(is.mcnode(y))
    m <- mc(a = mcdata(1, "V"), b = mcdata(2, "V"))
    expect_identical(attr(outm(m, "none", which.node = 2)$b, "outm"), "none")
    expect_identical(attr(outm(m, "mean", which.node = "a")$a, "outm"), "mean")
    expect_identical(attr(outm(mcdata(1, "V"), c("min", "max")), "outm"), c("min", "max"))
  })
})

test_that("rmultinomial: rows in the order of the parameters", {
  set.seed(1)
  r <- rmultinomial(6, size = c(10, 1000), prob = c(1, 1))
  expect_identical(rowSums(r), c(10, 1000, 10, 1000, 10, 1000))
  r3 <- rmultinomial(7, size = c(1, 10, 100), prob = c(.2, .8))
  expect_identical(rowSums(r3), c(1, 10, 100, 1, 10, 100, 1))
  with_dims(7, 5, {
    s <- mcstoc(rpois, "V", lambda = 50)
    m <- mcstoc(rmultinomial, "VU", nvariates = 3, size = s, prob = c(.2, .3, .5))
    expect_identical(as.vector(apply(unclass(m), c(1, 2), sum)), rep(as.vector(unclass(s)), 5))
    su <- mcstoc(rpois, "U", lambda = 50)
    mu <- mcstoc(rmultinomial, "VU", nvariates = 3, size = su, prob = c(.2, .3, .5))
    expect_identical(as.vector(apply(unclass(mu), c(1, 2), sum)), rep(as.vector(unclass(su)), each = 7))
  })
})

test_that("dmultinomial / ddirichlet / dbetasubj / dlnormb / qmqi edge cases", {
  expect_identical(dmultinomial(c(1, 1), prob = c(0, 1)), 0)
  expect_identical(dmultinomial(c(1, 1), prob = c(0, 1), log = TRUE), -Inf)
  expect_equal(ddirichlet(c(0, 1), c(1, 1)), 1)
  expect_identical(dbetasubj(c(-1, 6), 0, 1, 2, 5), c(0, 0))
  expect_true(is.nan(suppressWarnings(dlnormb(1, mean = -1, sd = 1))))
  expect_equal(qmqi(c(-0.1, .2, .6, .97, .3), c(40, 50, 60))[2:5], qmqi(c(.2, .6, .97, .3), c(40, 50, 60)))
  expect_true(is.nan(qmqi(-0.1, c(40, 50, 60))))
})

test_that("empirical distributions: support, single value, NA", {
  set.seed(2)
  r <- rempiricalC(2, 0, 10, values = c(1, 9, 9, 9, 9))
  expect_length(r, 2)
  expect_identical(rempiricalD(3, values = 10), c(10, 10, 10))
  expect_identical(suppressWarnings(qempiricalD(c(NA, .5), 1:3)), c(NA, 2))
})

test_that("pert / triang / betagen: own warnings are emitted; recycling of min", {
  expect_warning(rpert(3, 0, mode = 2, 1), "mode < min or mode > max")
  expect_warning(rbetagen(2, 1, 1, min = 1, max = 0), "NaN in rbetagen")
  expect_equal(dpert(0, -1, 0, 1, shape = 2000), dbeta(.5, 1001, 1001) / 2)
  expect_equal(dtriang(0.5, min = c(0, .25), mode = .5, max = 1), c(dtriang(.5, 0, .5, 1), dtriang(.5, .25, .5, 1)))
  expect_equal(ptriang(0.5, min = c(0, .25), mode = .5, max = 1), c(ptriang(.5, 0, .5, 1), ptriang(.5, .25, .5, 1)))
})

test_that("internal ggarrange (replaces ggpubr) gives a ggplot object", {
  utils::data("total", package = "mc2d", envir = environment())
  p <- on_null_device(ggplotmc(total))
  expect_true(inherits(p, "ggplot"))
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f))
  ggplot2::ggsave(f, p, width = 6, height = 5)
  expect_true(file.exists(f))
  expect_true(inherits(on_null_device(gghist(total)), "ggplot"))
})

test_that("plot.mccut / plot.tornadounc / plot with a list of xlim", {
  utils::data("total", package = "mc2d", envir = environment())
  on_null_device({
    expect_no_error(plot(tornadounc(total, 10, use = "complete.obs")))
    expect_no_error(plot(mc(xV, xVU), xlim = list(c(0, 1), c(0, 2))))
    expect_true(inherits(ggtornadounc(tornadounc(total, 10, use = "complete.obs")), "ggplot"))
  })
})

test_that("evalmccut with nsu = 1 and mixed U / VU nodes in the first block", {
  mod <- mcmodelcut({
    { a <- mcstoc(rnorm, "U"); b <- mcstoc(rnorm, "VU") }
    { d <- a + b; res <- mc(a, b, d) }
    { list(node = res$d, s = summary(res)) }
  })
  x1 <- suppressWarnings(capture.output(r1 <- evalmccut(mod, nsv = 10, nsu = 1, seed = 1)))
  expect_s3_class(r1, "mccut")
  capture.output(r2 <- evalmccut(mod, nsv = 10, nsu = 4, seed = 1))
  expect_identical(dim(r2$node), c(10L, 4L, 1L))
})

test_that("compatibility with 0.2.2: class 'ggarrange' and names from mc objects", {
  utils::data("total", package = "mc2d", envir = environment())
  on_null_device({
    expect_true(inherits(ggplotmc(total), "ggarrange"))
    p <- ggplotmc(total, griddim = c(1, 2))           # more graphs than cells: list of pages
    expect_identical(class(p), c("list", "ggarrange"))
    expect_false(inherits(p[[1]], "ggarrange"))
  })
  with_dims(7, 5, {
    m1 <- mc(x = mcdata(1, "V")); m2 <- mc(x = mcdata(2, "V"))
    expect_identical(names(mc(m1, m2)), c("x", "x"))   # names of mc elements are kept
    a <- mcdata(1, "V")
    expect_identical(names(mc(a, a)), c("a", "a1"))    # names of mcnode arguments are made unique
  })
})
