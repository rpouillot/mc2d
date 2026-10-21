## Characterization tests: summaries, quantiles, convergence, tornado, models.
## See helper-ref.R. Comments tagged "BUG:" document a behaviour of
## mc2d <= 0.2.2 that is known to be wrong and is fixed in 0.2.3.

load_total <- function(env = parent.frame()) {
  utils::data("total", package = "mc2d", envir = env)
}

small_model <- function() {
  nodes <- make_nodes(nsv = 31, nsu = 11)
  with_dims(31, 11, {
    set.seed(2)
    a <- nodes$V; b <- nodes$U; c2 <- nodes$VU
    d <- mcstoc(rpois, "VU", lambda = abs(a * b) + 1)
    e <- a + b * c2 + d
    e[3, 2, 1] <- NA
    mc(a = a, b = b, c2 = c2, d = d, e = e)
  })
}

test_that("summary and print of mc / mcnode", {
  load_total()
  expect_ref("summary-total", summary(total))
  expect_ref("summary-total-probs", summary(total, probs = c(0, .1, .5, 1), lim = c(.05, .95)))
  expect_ref("summary-total-digits", summary(total, digits = 5))
  for (nm in c("x0", "x0M", "xU", "xUM", "xV", "xVM", "xVU", "xVUM", "xVUM2", "xVUM3")) {
    node <- get(nm)
    expect_ref(paste0("summary-", nm), summary(node))
    expect_ref(paste0("summary-print-", nm), capture.output(print(summary(node))))
    expect_ref(paste0("print-", nm), capture.output(print(node)))
  }
  expect_ref("print-total", capture.output(print(total)))
  expect_ref("print-total-digits", capture.output(print(total, digits = 2)))
  expect_ref("print-summary-total", capture.output(print(summary(total))))
  m <- small_model()
  # BUG: the NA warning is issued once per node
  expect_ref("summary-na", summary(m))
  expect_ref("summary-outm-fun", summary(mc(mcdata(array(1:60, c(5, 4, 3)), "VU", nsv = 5, nsu = 4,
                                                       nvariates = 3, outm = c("mean", "sum", "max")))))
  expect_ref("summary-logical", summary(mcdata(c(TRUE, FALSE, TRUE), "V", nsv = 3)))
})

test_that("quantile", {
  load_total()
  expect_ref("quantile-total", quantile(total))
  expect_ref("quantile-total-probs", quantile(total, probs = c(.1, .9), lim = c(.25, .75)))
  for (nm in c("x0", "xU", "xV", "xVU", "xVUM", "xVUM3"))
    expect_ref(paste0("quantile-", nm), quantile(get(nm)))
  # BUG: a single probability fails for VU nodes
  expect_ref("quantile-oneprob", quantile(xVU, probs = .5))
  expect_ref("quantile-na", quantile(small_model(), na.rm = TRUE))
})

test_that("converg", {
  load_total()
  on_null_device({
    expect_ref("converg-var", converg(xVU, margin = "var"))
    expect_ref("converg-unc", converg(xVU, margin = "unc"))
    expect_ref("converg-V", converg(xV))
    expect_ref("converg-U", converg(xU))
    expect_ref("converg-iter", converg(xVU, margin = "var", iter = 5))
    expect_ref("converg-probs", converg(xVU, probs = c(.1, .9), lim = c(.1, .9)))
    expect_ref("converg-log", converg(abs(xVU) + 1, log = TRUE))
    expect_ref("converg-mc", converg(total, node = 2))
    # BUG: nvariates > 1 is out of bounds
    expect_ref("converg-nvariates", converg(xVUM, margin = "var", nvariates = 2))
    expect_ref("converg-probs01", converg(xV, probs = c(0, 1)))
  })
})

test_that("mcratio", {
  load_total()
  expect_ref("mcratio-total", mcratio(total, na.rm = TRUE))
  expect_ref("mcratio-total-p", mcratio(total, pcentral = .4, pvar = .9, punc = .9, na.rm = TRUE))
  expect_ref("mcratio-node", mcratio(xVU))
  expect_ref("mcratio-outm", mcratio(xVUM3))
  expect_ref("mcratio-na", mcratio(small_model()))
  expect_ref("mcratio-none", mcratio(mc(xVUM2)))
})

test_that("tornado and tornadounc", {
  load_total()
  m <- small_model()
  for (me in c("spearman", "pearson", "kendall")) {
    expect_ref(paste0("tornado-total-", me), tornado(total, method = me, use = "complete.obs"))
    expect_ref(paste0("tornado-small-", me), tornado(m, method = me, use = "complete.obs"))
    expect_ref(paste0("tornado-small-pairwise-", me), tornado(m, method = me, use = "pairwise.complete.obs"))
  }
  expect_ref("tornado-small-allobs", tornado(m))
  expect_ref("tornado-output", tornado(total, output = 3, use = "complete.obs"))
  expect_ref("tornado-lim", tornado(total, lim = c(.1, .9), use = "complete.obs"))
  expect_ref("tornado-print", capture.output(print(tornado(total, use = "complete.obs"))))
  expect_ref("tornado-outnone", tornado(mc(xVU, xVUM2), output = 2))
  expect_ref("tornado-baduse", tornado(total, use = "foo"))
  expect_ref("tornado-total-allobs", tornado(total))
  expect_ref("tornadounc-total", tornadounc(total, 3))
  expect_ref("tornadounc-total-10", tornadounc(total, 10, use = "complete.obs"))
  expect_ref("tornadounc-quant", tornadounc(total, 10, quant = c(.25, .75), method = "pearson"))
  expect_ref("tornadounc-print", capture.output(print(tornadounc(total, 10, use = "complete.obs"))))
  expect_ref("tornadounc-small", tornadounc(m, output = 5, use = "complete.obs"))
  # BUG: pairwise with NA mis-shapes the data
  expect_ref("tornadounc-small-pairwise", tornadounc(m, output = 5, use = "pairwise.complete.obs"))
  expect_ref("tornadounc-default", tornadounc(xVU, xV))
})

test_that("evalmcmod, mcmodel", {
  utils::data("ec", package = "mc2d", envir = environment())
  expect_ref("evalmcmod-ec1", evalmcmod(ec$modEC1, nsv = 50, nsu = 10, seed = 666), seed = NULL)
  expect_ref("evalmcmod-ec2", evalmcmod(ec$modEC2, nsv = 50, nsu = 10, seed = 666), seed = NULL)
  mod <- mcmodel({
    x <- mcstoc(runif, type = "V")
    y <- mcdata(2, type = "0")
    z <- x * y
    mc(x, y, z)
  })
  expect_ref("mcmodel", deparse(mod), seed = NULL)
  expect_ref("evalmcmod-mod", evalmcmod(mod, nsv = 20, nsu = 1, seed = 1), seed = NULL)
  expect_ref("mcmodel-expr", deparse(mcmodel(expression({ x <- 1; mc(x) }), is.expr = TRUE)), seed = NULL)
  expect_ref("mcmodel-print", capture.output(print(mod)), seed = NULL)
  expect_ref("ndvar-after-eval", c(ndvar(), ndunc()), seed = NULL)
})

modEC3 <- mcmodelcut({
  {
    cook <- mcstoc(rempiricalD, type = "V", values = c(0, 1/5, 1/50), prob = c(0.027, 0.373, 0.6))
    serving <- mcstoc(rgamma, type = "V", shape = 3.93, rate = 0.0806)
    conc <- mcstoc(rnorm, type = "U", mean = 10, sd = 2)
    r <- mcstoc(runif, type = "U", min = 5e-04, max = 0.0015)
  }
  {
    expo <- conc * cook * serving
    dose <- mcstoc(rpois, type = "VU", lambda = expo)
    risk <- 1 - (1 - r)^dose
    res <- mc(conc, cook, serving, expo, dose, r, risk)
  }
  {
    list(
      sum = summary(res),
      plot = plot(res, draw = FALSE),
      minmax = lapply(res, range),
      tor = tornado(res, use = "complete.obs"),
      node = res$risk
    )
  }
})

test_that("evalmccut, mcmodelcut and methods", {
  expect_ref("mcmodelcut", deparse(modEC3), seed = NULL)
  x <- evalmccut(modEC3, nsv = 41, nsu = 9, seed = 666)
  expect_ref("evalmccut", x, seed = NULL)
  expect_ref("evalmccut-summary", summary(x), seed = NULL)
  expect_ref("evalmccut-print", capture.output(print(x)), seed = NULL)
  expect_ref("evalmccut-print-summary", capture.output(print(summary(x))), seed = NULL)
  expect_ref("evalmccut-tornadounc", tornadounc(x), seed = NULL)
  expect_ref("evalmccut-tornadounc-co", tornadounc(x, use = "complete.obs"), seed = NULL)
  on_null_device({
    expect_ref("evalmccut-plot", plot(x), seed = NULL)
  })
  expect_ref("evalmccut-ind", evalmccut(modEC3, nsv = 11, nsu = 3, seed = 1, ind = "j"), seed = NULL)
  # BUG: fails when nsu == 1
  expect_ref("evalmccut-nsu1", evalmccut(modEC3, nsv = 11, nsu = 1, seed = 1), seed = NULL)
  expect_ref("ndvar-after-cut", c(ndvar(), ndunc()), seed = NULL)
})

test_that("summary invariants (run everywhere)", {
  load_total()
  s <- summary(xVU)
  expect_s3_class(s, "summary.mc")
  q <- quantile(xV, probs = c(.25, .75))
  expect_true(is.list(q) || is.numeric(q))
  t1 <- tornado(total, use = "complete.obs")
  expect_s3_class(t1, "tornado")
  vals <- unlist(t1$value)
  expect_true(all(is.na(vals) | abs(vals) <= 1))
})
