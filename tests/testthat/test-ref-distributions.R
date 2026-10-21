## Characterization tests: probability distributions (d, p, q, r).
## See helper-ref.R. Comments tagged "BUG:" document a behaviour of
## mc2d <= 0.2.2 that is known to be wrong and is fixed in 0.2.3.

pgrid <- c(0, 1e-10, .01, .25, .5, .75, .99, 1, NA, -0.1, 1.1)
lpgrid <- log(c(1e-10, .01, .25, .5, .75, .99, 1))

check_family <- function(prefix, d, p, q, r, xs, args, rn = 20) {
  expect_ref(paste0(prefix, "-d"), do.call(d, c(list(xs), args)))
  expect_ref(paste0(prefix, "-dlog"), do.call(d, c(list(xs), args, log = TRUE)))
  expect_ref(paste0(prefix, "-p"), do.call(p, c(list(xs), args)))
  expect_ref(paste0(prefix, "-pupper"), do.call(p, c(list(xs), args, lower.tail = FALSE)))
  expect_ref(paste0(prefix, "-plog"), do.call(p, c(list(xs), args, log.p = TRUE)))
  expect_ref(paste0(prefix, "-plogupper"), do.call(p, c(list(xs), args, lower.tail = FALSE, log.p = TRUE)))
  expect_ref(paste0(prefix, "-q"), do.call(q, c(list(pgrid), args)))
  expect_ref(paste0(prefix, "-qupper"), do.call(q, c(list(pgrid), args, lower.tail = FALSE)))
  expect_ref(paste0(prefix, "-qlog"), do.call(q, c(list(lpgrid), args, log.p = TRUE)))
  expect_ref(paste0(prefix, "-qlogupper"), do.call(q, c(list(lpgrid), args, log.p = TRUE, lower.tail = FALSE)))
  expect_ref(paste0(prefix, "-r"), do.call(r, c(list(rn), args)))
  expect_ref(paste0(prefix, "-r0"), do.call(r, c(list(0), args)))
}

test_that("bern", {
  check_family("bern", dbern, pbern, qbern, rbern, c(-1, 0, .5, 1, 2, NA), list(prob = .3))
  check_family("bern-vec", dbern, pbern, qbern, rbern, c(0, 1), list(prob = c(.1, .5, .9)))
  expect_ref("bern-bad", list(dbern(1, prob = 2), rbern(3, prob = -1), qbern(.5, prob = NA)))
})

test_that("betagen", {
  xs <- c(-1, 0, .5, 1, 2, 3, 5, NA)
  check_family("betagen", dbetagen, pbetagen, qbetagen, rbetagen, xs,
               list(shape1 = 2, shape2 = 3, min = 1, max = 4))
  check_family("betagen-ncp", dbetagen, pbetagen, qbetagen, rbetagen, xs,
               list(shape1 = 2, shape2 = 3, min = 1, max = 4, ncp = 1))
  check_family("betagen-vec", dbetagen, pbetagen, qbetagen, rbetagen, xs,
               list(shape1 = c(1, 2), shape2 = 3, min = c(0, 1, 2), max = 4))
  # BUG: warnings are suppressed by options(warn = -1)
  expect_ref("betagen-bad", list(rbetagen(2, 1, 1, min = 1, max = 0), dbetagen(.5, 1, 1, min = 1, max = 0),
                                 qbetagen(.5, 1, 1, min = 1, max = 0), pbetagen(.5, 1, 1, min = 1, max = 0)))
  expect_ref("betagen-equal", list(rbetagen(3, 1, 1, min = 1, max = 1), qbetagen(.5, 1, 1, min = 1, max = 1)))
})

test_that("betasubj", {
  xs <- c(-1, 0, 1, 2, 3, 4, 5, 6, NA)
  args <- list(min = 0, mode = 2, mean = 2.2, max = 5)
  check_family("betasubj", dbetasubj, pbetasubj, qbetasubj, rbetasubj, xs, args)
  check_family("betasubj-vec", dbetasubj, pbetasubj, qbetasubj, rbetasubj, xs,
               list(min = 0, mode = c(1, 2), mean = c(1.5, 2.2), max = 5))
  # BUG: density outside the support is NaN instead of 0
  expect_ref("betasubj-outside", dbetasubj(c(-1, 6), 0, 1, 2, 5))
  expect_ref("betasubj-bad", list(dbetasubj(1, 0, 3, 2, 1), rbetasubj(2, 0, 2, 2, 5), qbetasubj(.5, 0, 6, 2, 5)))
})

test_that("empiricalC", {
  xs <- c(-1, 0, 1, 2.5, 4, 5, 6, NA)
  check_family("empC", dempiricalC, pempiricalC, qempiricalC, rempiricalC, xs,
               list(min = 0, max = 5, values = c(1, 2, 3, 4), prob = c(1, 2, 3, 4)))
  check_family("empC-noprob", dempiricalC, pempiricalC, qempiricalC, rempiricalC, xs,
               list(min = 0, max = 5, values = c(1, 2, 2, 4)))
  check_family("empC-vecmin", dempiricalC, pempiricalC, qempiricalC, rempiricalC, xs,
               list(min = c(-1, 0, .5), max = c(5, 6), values = c(1, 2, 3, 4)))
  vm <- matrix(c(1, 2, 3, 1.5, 2.5, 3.5), nrow = 2, byrow = TRUE)
  pm <- matrix(c(1, 1, 1, 1, 2, 3), nrow = 2, byrow = TRUE)
  check_family("empC-matrix", dempiricalC, pempiricalC, qempiricalC, rempiricalC, xs,
               list(min = 0, max = 5, values = vm, prob = pm))
  # BUG: support truncated to n values
  expect_ref("empC-truncate", rempiricalC(2, 0, 10, values = c(1, 9, 9, 9, 9)))
  # BUG: grouping of duplicated values uses as.character
  expect_ref("empC-ties", dempiricalC(0.25, 0, 1, values = c(0.1 + 0.2, 0.3, 0.5)))
  expect_ref("empC-bad", list(dempiricalC(1, 2, 1, values = 1.5), rempiricalC(2, 0, 5, values = 6),
                              qempiricalC(.5, 0, 5, values = c(1, 2), prob = c(-1, 1))))
  expect_ref("empC-large", rempiricalC(1000, 0, 5, values = c(1, 2, 3, 4)))
  expect_ref("empC-large-vecmin", rempiricalC(200, min = seq(-1, 0, length = 200), max = 5, values = c(1, 2, 3, 4)))
})

test_that("empiricalD", {
  xs <- c(0, 1, 1.5, 2, 3, 5, 9, NA)
  check_family("empD", dempiricalD, pempiricalD, qempiricalD, rempiricalD, xs,
               list(values = c(1, 2, 3, 5), prob = c(1, 2, 3, 4)))
  check_family("empD-noprob", dempiricalD, pempiricalD, qempiricalD, rempiricalD, xs,
               list(values = c(1, 2, 2, 5)))
  vm <- matrix(c(1, 2, 3, 2, 3, 5), nrow = 2, byrow = TRUE)
  pm <- matrix(c(1, 1, 1, 1, 2, 3), nrow = 2, byrow = TRUE)
  expect_ref("empD-matrix", list(dempiricalD(xs[-8], vm, pm), pempiricalD(xs[-8], vm, pm),
                                 qempiricalD(pgrid[1:8], vm, pm), rempiricalD(20, vm, pm)))
  # BUG: a single support value fails (sample(x) with length(x) == 1)
  expect_ref("empD-single", rempiricalD(3, values = 10))
  expect_ref("empD-single-1", rempiricalD(3, values = 1))
  # BUG: NA probability
  expect_ref("empD-qNA", qempiricalD(NA_real_, 1:3))
  expect_ref("empD-q0", qempiricalD(0, 1:3, prob = c(0, 1, 1)))
  expect_ref("empD-large", rempiricalD(1000, values = c(1, 2, 3, 5), prob = c(1, 2, 3, 4)))
  expect_ref("empD-bad", list(dempiricalD(1, 1:3, prob = c(-1, 1, 1)), rempiricalD(2, 1:3, prob = 1:2)))
})

test_that("empiricalC: varying min / max", {
  set.seed(6)
  mi <- runif(60, -2, 0.9); ma <- runif(60, 4.1, 7)
  v <- c(1, 2, 2, 3, 4); pr <- c(1, 2, 3, 4, 5)
  xs <- c(-3, -1, 0.5, 1, 1.5, 2, 3.5, 4, 5, 8, mi[1:5], ma[1:5])
  xx <- sample(xs, 60, TRUE)
  pp <- c(0, 1e-12, runif(56), 1 - 1e-12, 1)
  for (k in list(list(mi, ma), list(mi, 5), list(-1, ma), list(mi[1:20], ma[1:30]), list(mi[1:7], ma[1:3]))) {
    nm <- paste(lengths(k), collapse = "x")
    expect_ref(paste0("empC2-d-", nm), list(dempiricalC(xx, k[[1]], k[[2]], v, pr), dempiricalC(xx, k[[1]], k[[2]], v),
                                            dempiricalC(xx[1:13], k[[1]], k[[2]], v, pr, log = TRUE)))
    expect_ref(paste0("empC2-p-", nm), list(pempiricalC(xx, k[[1]], k[[2]], v, pr), pempiricalC(xx, k[[1]], k[[2]], v),
                                            pempiricalC(xx[1:13], k[[1]], k[[2]], v, pr, lower.tail = FALSE)))
    expect_ref(paste0("empC2-q-", nm), list(qempiricalC(pp, k[[1]], k[[2]], v, pr), qempiricalC(pp, k[[1]], k[[2]], v),
                                            qempiricalC(pp[1:13], k[[1]], k[[2]], v, pr, lower.tail = FALSE)))
    expect_ref(paste0("empC2-r-", nm), rempiricalC(60, k[[1]], k[[2]], v, pr))
  }
  expect_ref("empC2-one-value", list(qempiricalC(pp, mi, ma, 2), dempiricalC(xx, mi, ma, 2), pempiricalC(xx, mi, ma, 2)))
  expect_ref("empC2-NA", list(qempiricalC(c(NA, .5, NaN), mi[1:3], ma[1:3], v), dempiricalC(c(NA, 2, NaN), mi[1:3], ma[1:3], v),
                              pempiricalC(c(NA, 2, NaN), mi[1:3], ma[1:3], v), qempiricalC(c(NA, NA), mi[1:2], ma[1:2], v)))
  expect_ref("empC2-big", rempiricalC(2000, min = runif(2000, -1, 0), max = 10, values = 1:5))
})

test_that("empiricalD: quantiles and multi-row parameters", {
  p0 <- c(0, 1e-10, .01, .1, .25, .3, .5, .6, .75, .99, 1, -0.1, 1.1)
  v <- c(1, 2, 3, 5); pr <- c(1, 2, 3, 4)
  expect_ref("empD2-q", list(qempiricalD(p0, v, pr), qempiricalD(p0, v), qempiricalD(p0, c(5, 1, 3, 1)),
                             qempiricalD(p0, v, pr, lower.tail = FALSE), qempiricalD(log(p0[2:11]), v, pr, log.p = TRUE)))
  vm <- matrix(c(1, 2, 3, 2, 3, 5, 0, 1, 1), nrow = 3, byrow = TRUE)
  pm <- matrix(c(1, 1, 1, 1, 2, 3, 0, 1, 2), nrow = 3, byrow = TRUE)
  set.seed(5)
  xx <- sample(c(0, 1, 1.5, 2, 3, 5, 9), 30, replace = TRUE)
  pp <- runif(30)
  expect_ref("empD2-matrix-d", list(dempiricalD(xx, vm, pm), dempiricalD(xx, vm), dempiricalD(xx, vm[1, ], pm),
                                    dempiricalD(xx[1:31 %% 31 > 0][1:29], vm, pm)))
  expect_ref("empD2-matrix-p", list(pempiricalD(xx, vm, pm), pempiricalD(xx, vm), pempiricalD(xx, vm[1, ], pm),
                                    pempiricalD(xx[1:29], vm, pm, lower.tail = FALSE)))
  expect_ref("empD2-matrix-q", list(qempiricalD(pp, vm, pm), qempiricalD(pp, vm), qempiricalD(pp, vm[1, ], pm),
                                    qempiricalD(pp[1:29], vm, pm)))
  bigv <- matrix(round(runif(200 * 4, 0, 10)), ncol = 4)
  bigp <- matrix(runif(200 * 4), ncol = 4)
  expect_ref("empD2-many-rows", list(dempiricalD(bigv[, 2], bigv, bigp), pempiricalD(bigv[, 3] + .5, bigv, bigp),
                                     qempiricalD(runif(200), bigv, bigp), rempiricalD(200, bigv, bigp)))
  expect_ref("empD2-long", list(dempiricalD(sample(0:6, 5000, TRUE), 1:5, 5:1), pempiricalD(runif(5000, 0, 6), 1:5, 5:1),
                                qempiricalD(runif(5000), 1:5, 5:1)))
})

test_that("lnormb", {
  xs <- c(-1, 0, .5, 1, 2, 10, NA)
  check_family("lnormb", dlnormb, plnormb, qlnormb, rlnormb, xs, list())
  check_family("lnormb-par", dlnormb, plnormb, qlnormb, rlnormb, xs, list(mean = 2, sd = c(1, 3)))
  # BUG: negative mean accepted; precision with tiny sd/mean
  expect_ref("lnormb-bad", list(rlnormb(5, mean = -3, sd = 1), qlnormb(.5, mean = 1, sd = 1e-9),
                                dlnormb(1, mean = 1, sd = -1)))
})

test_that("mqi", {
  xs <- c(30, 40, 45, 50, 55, 60, 70, NA)
  args <- list(mqi = c(40, 50, 60))
  check_family("mqi", dmqi, pmqi, qmqi, rmqi, xs, args)
  check_family("mqi-opts", dmqi, pmqi, qmqi, rmqi, xs,
               list(mqi = c(40, 50, 60), mqi.quantile = c(.1, .5, .9), k = .2, intrinsic = c(0, 100)))
  check_family("mqi-real", dmqi, pmqi, qmqi, rmqi, xs, list(mqi = c(40, 50, 60), realization = 65))
  # BUG: an invalid p corrupts the other results
  expect_ref("mqi-q-bad", qmqi(c(-0.1, .2, .6, .97, .3), c(40, 50, 60)))
  expect_ref("mqi-bad", list(dmqi(1, c(60, 50, 40)), rmqi(2, c(40, 50))))
})

test_that("pert", {
  xs <- c(-2, -1, -.5, 0, .5, 1, 2, NA)
  check_family("pert", dpert, ppert, qpert, rpert, xs, list())
  check_family("pert-par", dpert, ppert, qpert, rpert, xs, list(min = -1, mode = .5, max = 1, shape = 6))
  check_family("pert-mean", dpert, ppert, qpert, rpert, xs, list(min = -1, mean = .2, max = 1))
  check_family("pert-vec", dpert, ppert, qpert, rpert, xs, list(min = c(-1, -2), mode = c(0, .5, -.5), max = 1))
  check_family("pert-equal", dpert, ppert, qpert, rpert, xs, list(min = 0, mode = 0, max = 0))
  check_family("pert-modemin", dpert, ppert, qpert, rpert, xs, list(min = 0, mode = 0, max = 1))
  # BUG: warnings are suppressed; density underflow for large shape
  expect_ref("pert-bad", list(rpert(3, 0, mode = 2, 1), dpert(0, -1, 0, 1, shape = 2000),
                              qpert(.5, 0, 2, 1), ppert(.5, 0, 2, 1), dpert(.5, 0, 2, 1)))
  expect_ref("pert-shape0", list(rpert(5, 0, .5, 1, shape = 0), dpert(.5, 0, .5, 1, shape = 0)))
  expect_ref("pert-large", rpert(1000, min = 0, mode = 1, max = 4))
  expect_ref("pert-large-vec", rpert(1000, min = runif(1000), mode = 2, max = 4))
})

test_that("triang", {
  xs <- c(-2, -1, -.5, 0, .5, 1, 2, NA)
  check_family("triang", dtriang, ptriang, qtriang, rtriang, xs, list())
  check_family("triang-par", dtriang, ptriang, qtriang, rtriang, xs, list(min = -1, mode = .5, max = 1))
  check_family("triang-mean", dtriang, ptriang, qtriang, rtriang, xs, list(min = -1, mean = .2, max = 1))
  check_family("triang-vec", dtriang, ptriang, qtriang, rtriang, xs, list(min = c(-1, -2), mode = c(0, .5, -.5), max = 1))
  check_family("triang-equal", dtriang, ptriang, qtriang, rtriang, xs, list(min = 0, mode = 0, max = 0))
  check_family("triang-modemax", dtriang, ptriang, qtriang, rtriang, xs, list(min = 0, mode = 1, max = 1))
  # BUG: min is not part of the recycling length
  expect_ref("triang-recycle", list(dtriang(0.5, min = c(0, .25), mode = .5, max = 1),
                                    ptriang(0.5, min = c(0, .25), mode = .5, max = 1)))
  expect_ref("triang-bad", list(rtriang(3, 0, 2, 1), dtriang(.5, 0, 2, 1), qtriang(.5, 0, 2, 1)))
  expect_ref("triang-large", rtriang(1000, min = 0, mode = 1, max = 4))
  expect_ref("triang-large-vec", rtriang(1000, min = runif(1000), mode = 2, max = 4))
})

test_that("dirichlet", {
  x <- rbind(c(.2, .3, .5), c(.1, .1, .8), c(0, .5, .5), c(.5, .6, -.1))
  expect_ref("dirichlet-d", list(ddirichlet(x, c(1, 2, 3)), ddirichlet(x[1, ], c(1, 2, 3)),
                                 ddirichlet(x, rbind(c(1, 2, 3), c(2, 2, 2)))))
  # BUG: 0 * log(0) gives NaN
  expect_ref("dirichlet-d-zero", ddirichlet(c(0, 1), c(1, 1)))
  expect_ref("dirichlet-r", list(rdirichlet(10, c(1, 2, 3)), rdirichlet(4, rbind(c(1, 2, 3), c(5, 5, 5)))))
  expect_ref("dirichlet-r0", rdirichlet(0, c(1, 2)))
  expect_ref("dirichlet-small", rdirichlet(20, c(.001, .001)))
  expect_ref("dirichlet-bad", list(rdirichlet(2, c(-1, 1)), ddirichlet(c(.5, .5), c(-1, 1))))
})

test_that("multinomial", {
  x <- rbind(c(1, 2, 7), c(3, 3, 4), c(0, 0, 10))
  expect_ref("multinomial-d", list(dmultinomial(x, prob = c(.2, .3, .5)),
                                   dmultinomial(x, prob = c(.2, .3, .5), log = TRUE),
                                   dmultinomial(x[1, ], prob = c(.2, .3, .5)),
                                   dmultinomial(x, prob = rbind(c(.2, .3, .5), c(0, .5, .5)))))
  expect_ref("multinomial-r", list(rmultinomial(5, 10, c(.2, .3, .5)),
                                   rmultinomial(4, c(10, 100), c(.2, .3, .5)),
                                   rmultinomial(4, 10, rbind(c(.2, .3, .5), c(.5, .5, 0)))))
  # BUG: rows in the wrong order when parameters vary
  expect_ref("multinomial-order", rmultinomial(6, size = c(10, 1000), prob = c(1, 1)))
  expect_ref("multinomial-order3", rmultinomial(6, size = c(10, 100, 1000), prob = c(1, 1)))
  expect_ref("multinomial-bad", list(rmultinomial(2, -1, c(.5, .5)), dmultinomial(c(-1, 2), prob = c(.5, .5))))
})

test_that("multinormal", {
  s <- matrix(c(1, .5, .5, 2), 2)
  x <- rbind(c(0, 0), c(1, 2), c(-1, .5))
  expect_ref("multinormal-d", list(dmultinormal(x, c(0, 1), as.vector(s)),
                                   dmultinormal(x, c(0, 1), as.vector(s), log = TRUE),
                                   dmultinormal(c(0, 0), c(0, 1), as.vector(s)),
                                   dmultinormal(x, rbind(c(0, 1), c(1, 1)), rbind(as.vector(s), c(2, 0, 0, 2)))))
  expect_ref("multinormal-r", list(rmultinormal(5, c(0, 1), as.vector(s)),
                                   rmultinormal(4, rbind(c(0, 1), c(10, 10)), as.vector(s)),
                                   rmultinormal(4, c(0, 1), rbind(as.vector(s), c(2, 0, 0, 2))),
                                   rmultinormal(3, c(0, 1), as.vector(s), method = "chol"),
                                   rmultinormal(3, c(0, 1), as.vector(s), method = "svd")))
  # BUG: dimension of the result with a single variate
  expect_ref("multinormal-nv1", rmultinormal(3, 0, matrix(1:3, ncol = 1)))
  expect_ref("multinormal-bad", rmultinormal(3, c(0, 1), c(1, 2, 3)))
})

test_that("distribution invariants (run everywhere)", {
  ps <- c(.01, .1, .5, .9, .99)
  expect_equal(ppert(qpert(ps, 0, 1, 4), 0, 1, 4), ps, tolerance = 1e-10)
  expect_equal(ptriang(qtriang(ps, 0, 1, 4), 0, 1, 4), ps, tolerance = 1e-10)
  expect_equal(pbetagen(qbetagen(ps, 2, 3, 1, 4), 2, 3, 1, 4), ps, tolerance = 1e-10)
  expect_equal(plnormb(qlnormb(ps, 2, 1), 2, 1), ps, tolerance = 1e-10)
  expect_equal(pbetasubj(qbetasubj(ps, 0, 2, 2.2, 5), 0, 2, 2.2, 5), ps, tolerance = 1e-10)
  expect_equal(pempiricalC(qempiricalC(ps, 0, 5, 1:4), 0, 5, 1:4), ps, tolerance = 1e-10)
  expect_equal(pmqi(qmqi(ps, c(40, 50, 60)), c(40, 50, 60)), ps, tolerance = 1e-10)
  expect_equal(integrate(dpert, -1, 1)$value, 1, tolerance = 1e-6)
  expect_equal(integrate(dtriang, -1, 1)$value, 1, tolerance = 1e-6)
  expect_equal(integrate(dempiricalC, 0, 5, min = 0, max = 5, values = 1:4)$value, 1, tolerance = 1e-4)
  expect_equal(integrate(dbetasubj, 0, 5, min = 0, mode = 2, mean = 2.2, max = 5)$value, 1, tolerance = 1e-6)
  expect_equal(integrate(dmqi, 30, 70, mqi = c(40, 50, 60))$value, 1, tolerance = 1e-4)
  expect_equal(sum(dempiricalD(1:4, 1:4, prob = 1:4)), 1)
  set.seed(1)
  expect_true(all(rpert(1000, 0, 1, 4) >= 0 & rpert(1000, 0, 1, 4) <= 4))
  expect_true(all(rtriang(1000, 0, 1, 4) >= 0))
  expect_true(all(rempiricalD(100, c(2, 4, 8)) %in% c(2, 4, 8)))
  r <- rempiricalC(1000, 0, 5, values = 1:4)
  expect_true(all(r >= 0 & r <= 5))
  expect_true(all(abs(rowSums(rdirichlet(10, c(1, 2, 3))) - 1) < 1e-12))
  expect_true(all(rowSums(rmultinomial(10, 7, c(.2, .3, .5))) == 7))
  expect_equal(dim(rmultinormal(4, c(0, 1), c(1, 0, 0, 1))), c(4L, 2L))
})

test_that("reference: r* with parameter lengths that are not multiples, d/p/q on matrices and mcnodes", {
  expect_ref("rgen-coprime", list(
    rtriang(6, min = c(0, 10), mode = c(1, 11, 21), max = 30),
    rpert(6, min = c(0, 10), mode = c(5, 15, 25), max = 30),
    rpert(12, min = c(0, 10), mode = c(5, 15, 25), max = 30, shape = c(2, 4, 6, 8)),
    rbetagen(6, 2, 3, min = c(0, 10), max = c(20, 30, 40))))
  m <- matrix(c(4, 5, 7, 8), 2, dimnames = list(c("a", "b"), NULL))
  expect_ref("dpq-matrix", list(
    dtriang(m, 3, 6, 10), ptriang(m, 3, 6, 10), qtriang(m / 10, 3, 6, 10),
    dpert(m, 3, 6, 10), ppert(m, 3, 6, 10), qpert(m / 10, 3, 6, 10),
    dbetagen(m, 2, 3, 3, 10), pbetagen(m, 2, 3, 3, 10), qbetagen(m / 10, 2, 3, 3, 10)))
  expect_ref("dpq-mcnode", with_dims(5, 3, {
    x <- mcstoc(runif, "VU", min = 3, max = 10)
    list(dtriang(x, 3, 6, 10), ptriang(x, 3, 6, 10), qtriang(x / 10, 3, 6, 10),
         dpert(x, 3, 6, 10), ppert(x, 3, 6, 10))
  }))
})
