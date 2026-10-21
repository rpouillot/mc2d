## Characterization tests: mcnode construction, arithmetic and utilities.
## See helper-ref.R. Comments tagged "BUG:" document a behaviour of
## mc2d <= 0.2.2 that is known to be wrong and is fixed in 0.2.3.

types <- c("0", "V", "U", "VU")

test_that("mcdata: every type from scalars, vectors, matrices and arrays", {
  with_dims(7, 5, {
    for (ty in types) {
      expect_ref(paste0("mcdata-scalar-", ty), mcdata(3, type = ty))
      expect_ref(paste0("mcdata-scalarM-", ty), mcdata(3, type = ty, nvariates = 3))
      n <- switch(ty, "0" = 1, V = 7, U = 5, VU = 35)
      expect_ref(paste0("mcdata-vector-", ty), mcdata(seq_len(n), type = ty))
      expect_ref(paste0("mcdata-vectorM-", ty), mcdata(seq_len(3 * n), type = ty, nvariates = 3))
      expect_ref(paste0("mcdata-logical-", ty), mcdata(rep(c(TRUE, FALSE), length.out = n), type = ty))
      expect_ref(paste0("mcdata-NA-", ty), mcdata(NA, type = ty))
      for (om in list("none", "each", "mean", "min"))
        expect_ref(paste0("mcdata-outm-", om, "-", ty), mcdata(1, type = ty, outm = om))
    }
    # BUG: a vector of outm functions errors with R >= 4.3
    expect_ref("mcdata-outm-minmax", mcdata(1, type = "V", outm = c("min", "max")))
    expect_ref("mcdata-matrix-VU", mcdata(matrix(1:35, 7), type = "VU"))
    expect_ref("mcdata-matrix-V", mcdata(matrix(1:7, 7), type = "V"))
    expect_ref("mcdata-array-VUM", mcdata(array(1:105, c(7, 5, 3)), type = "VU", nvariates = 3))
    expect_ref("mcdata-array-VM", mcdata(array(1:21, c(7, 1, 3)), type = "V", nvariates = 3))
    expect_ref("mcdata-array-UM", mcdata(array(1:15, c(1, 5, 3)), type = "U", nvariates = 3))
    # BUG: a 1 x k matrix is accepted for a "0" node
    expect_ref("mcdata-matrix-0-bad", mcdata(matrix(1:5, nrow = 1), type = "0"))
    expect_ref("mcdata-badlength", mcdata(1:3, type = "V"))
    expect_ref("mcdata-character", mcdata("a", type = "V"))
    expect_ref("mcdata-nsv-nsu", mcdata(1:6, type = "VU", nsv = 3, nsu = 2))
  })
})

test_that("mcdata: from mcnodes of every type", {
  nodes <- make_nodes()
  with_dims(7, 5, {
    for (nm in names(nodes)) for (ty in types) {
      nv <- dim(nodes[[nm]])[3]
      expect_ref(paste0("mcdata-node-", nm, "-to-", ty),
                 mcdata(nodes[[nm]], type = ty, nvariates = nv))
      expect_ref(paste0("mcdata-node-", nm, "-to-", ty, "-nv3"),
                 mcdata(nodes[[nm]], type = ty, nvariates = 3))
    }
  })
})

test_that("mcdatanocontrol", {
  with_dims(7, 5, {
    for (ty in types)
      expect_ref(paste0("mcdatanocontrol-", ty), mcdatanocontrol(1:35, type = ty))
    expect_ref("mcdatanocontrol-recycle", mcdatanocontrol(1:3, type = "VU"))
    expect_ref("mcdatanocontrol-nv", mcdatanocontrol(1:105, type = "VU", nvariates = 3))
    expect_ref("mcdatanocontrol-dims", mcdatanocontrol(1:6, type = "V", nsv = 6, nsu = 1))
  })
})

test_that("mcstoc: every type and distribution", {
  with_dims(7, 5, {
    for (ty in types) {
      expect_ref(paste0("mcstoc-runif-", ty), mcstoc(runif, type = ty))
      expect_ref(paste0("mcstoc-default-", ty), mcstoc(type = ty))
      expect_ref(paste0("mcstoc-rnorm-", ty), mcstoc(rnorm, type = ty, mean = 2, sd = 3))
      expect_ref(paste0("mcstoc-rnorm-vecpar-", ty), mcstoc(rnorm, type = ty, mean = 1:3))
      expect_ref(paste0("mcstoc-rpois-", ty), mcstoc(rpois, type = ty, lambda = 4))
      expect_ref(paste0("mcstoc-rbinom-", ty), mcstoc(rbinom, type = ty, size = 10, prob = .3))
      expect_ref(paste0("mcstoc-rempiricalD-", ty),
                 mcstoc(rempiricalD, type = ty, values = c(1, 5, 9), prob = c(1, 2, 3)))
      expect_ref(paste0("mcstoc-nv3-", ty), mcstoc(rnorm, type = ty, nvariates = 3))
      expect_ref(paste0("mcstoc-rmultinomial-", ty),
                 mcstoc(rmultinomial, type = ty, nvariates = 3, size = 10, prob = c(.2, .3, .5)))
      expect_ref(paste0("mcstoc-rdirichlet-", ty),
                 mcstoc(rdirichlet, type = ty, nvariates = 3, alpha = c(1, 2, 3)))
      expect_ref(paste0("mcstoc-seed-", ty), mcstoc(rnorm, type = ty, seed = 99), seed = NULL)
      expect_ref(paste0("mcstoc-outm-", ty), mcstoc(rnorm, type = ty, outm = "mean"))
      expect_ref(paste0("mcstoc-outm-none-", ty), mcstoc(rnorm, type = ty, outm = "none"))
      expect_ref(paste0("mcstoc-outm-minmax-", ty), mcstoc(rnorm, type = ty, outm = c("min", "max")))
      expect_ref(paste0("mcstoc-rtrunc-", ty),
                 mcstoc(rnorm, type = ty, rtrunc = TRUE, linf = -1, lsup = 1))
      expect_ref(paste0("mcstoc-rtrunc-lsup-", ty),
                 mcstoc(rlnorm, type = ty, rtrunc = TRUE, lsup = 2))
      expect_ref(paste0("mcstoc-lhs-", ty), mcstoc(rnorm, type = ty, lhs = TRUE))
      expect_ref(paste0("mcstoc-lhs-rtrunc-", ty),
                 mcstoc(rnorm, type = ty, lhs = TRUE, rtrunc = TRUE, linf = 0))
      # BUG: default func with rtrunc / lhs fails (the function name is not found)
      expect_ref(paste0("mcstoc-default-rtrunc-", ty), mcstoc(type = ty, rtrunc = TRUE, linf = .5))
      expect_ref(paste0("mcstoc-default-lhs-", ty), mcstoc(type = ty, lhs = TRUE))
    }
    expect_ref("mcstoc-nsv-nsu", mcstoc(rnorm, type = "VU", nsv = 3, nsu = 4))
    expect_ref("mcstoc-nsample", mcstoc(sample, type = "V", x = 1:3, replace = TRUE, nsample = "size"))
    expect_ref("mcstoc-rtrunc-error", mcstoc(rnorm, type = "V", rtrunc = TRUE, linf = 2, lsup = 1))
    expect_ref("mcstoc-rtrunc-multi", mcstoc(rmultinomial, type = "V", nvariates = 3, size = 10,
                                               prob = c(.2, .3, .5), rtrunc = TRUE, linf = 1))
    expect_ref("mcstoc-array-param", mcstoc(rnorm, type = "V", mean = array(1, c(7, 1, 1))))
    expect_ref("mcstoc-bad-outm", mcstoc(rnorm, type = "V", outm = "notAFunction"))
    expect_ref("mcstoc-bad-length", mcstoc(rmultinomial, type = "V", nvariates = 2, size = 10,
                                             prob = c(.2, .3, .5)))
  })
})

test_that("mcstoc: mcnodes as parameters", {
  nodes <- make_nodes()
  with_dims(7, 5, {
    for (ty in types) for (nm in names(nodes)) {
      p <- abs(nodes[[nm]])
      nv <- dim(p)[3]
      expect_ref(paste0("mcstoc-param-", nm, "-in-", ty),
                 mcstoc(rnorm, type = ty, mean = p, sd = 0.1, nvariates = nv))
      expect_ref(paste0("mcstoc-param-", nm, "-in-", ty, "-two"),
                 mcstoc(runif, type = ty, min = p, max = p + 1, nvariates = nv))
      expect_ref(paste0("mcstoc-param-", nm, "-in-", ty, "-rtrunc"),
                 mcstoc(rnorm, type = ty, mean = p, rtrunc = TRUE, linf = p - 0.5, nvariates = nv))
      if (nv == 1) {
        expect_ref(paste0("mcstoc-param-", nm, "-in-", ty, "-lhs"),
                   mcstoc(rnorm, type = ty, mean = p, lhs = TRUE))
        expect_ref(paste0("mcstoc-param-", nm, "-in-", ty, "-nv3"),
                   mcstoc(rnorm, type = ty, mean = p, nvariates = 3))
      }
    }
    expect_ref("mcstoc-param-multinomial",
               mcstoc(rmultinomial, type = "VU", nvariates = 3, size = mcstoc(rpois, "V", lambda = 10),
                      prob = c(.2, .3, .5)))
    expect_ref("mcstoc-param-multinomial-U",
               mcstoc(rmultinomial, type = "VU", nvariates = 3, size = mcstoc(rpois, "U", lambda = 10),
                      prob = c(.2, .3, .5)))
  })
})

test_that("Ops.mcnode: all pairs of types, univariate and multivariate", {
  nodes <- make_nodes()
  for (a in names(nodes)) for (b in names(nodes)) {
    e1 <- nodes[[a]]; e2 <- nodes[[b]]
    expect_ref(paste0("ops-", a, "-plus-", b), e1 + e2)
    expect_ref(paste0("ops-", a, "-times-", b), e1 * e2)
    expect_ref(paste0("ops-", a, "-gt-", b), e1 > e2)
    expect_ref(paste0("ops-", a, "-pow-", b), abs(e1)^e2)
  }
})

test_that("Ops.mcnode: mcnode with numbers, vectors, arrays; unary ops", {
  nodes <- make_nodes()
  with_dims(7, 5, {
    for (nm in names(nodes)) {
      x <- nodes[[nm]]
      d <- dim(x)
      expect_ref(paste0("ops-num-", nm), list(x + 1, 2 * x, x / 3, 3 - x, x %% 1, x %/% 1, x == x, -x, !(x > 0)))
      expect_ref(paste0("ops-vec-", nm), list(x + seq_len(d[1] * d[2]), seq_len(prod(d)) * x))
      expect_ref(paste0("ops-badvec-", nm), x + 1:2)
      expect_ref(paste0("ops-array-", nm), list(x + array(1, d), array(2, d) * x))
      expect_ref(paste0("ops-matrix-", nm), x + matrix(1, d[1], d[2]))
      expect_ref(paste0("ops-array1-", nm), x + array(1, c(d[1], d[2], 1)))
      expect_ref(paste0("ops-badarray-", nm), x + array(1, c(2, 2, 2, 2)))
      expect_ref(paste0("ops-char-", nm), x + "a")
      expect_ref(paste0("ops-logical-", nm), (x > 0) & (x < 1))
    }
    # outm combination rules
    a <- mcdata(1, "V", outm = "none"); b <- mcdata(2, "V", outm = "mean")
    c2 <- mcdata(3, "V", outm = "each")
    expect_ref("ops-outm", list(attr(a + b, "outm"), attr(b + a, "outm"), attr(b + c2, "outm"),
                                attr(a + a, "outm")))
    # BUG: a vector outm errors with R >= 4.3
    m <- mcstoc(rnorm, "V", outm = c("min", "max"))
    expect_ref("ops-outm-vector", m + b)
    expect_ref("ops-incompatible", mcdata(1, "VU", nsv = 7, nsu = 5) + mcdata(1, "VU", nsv = 6, nsu = 5))
    expect_ref("ops-incompatible-nv", nodes$VUM + mcdata(1, "V", nvariates = 2))
  })
})

test_that("pmin / pmax", {
  nodes <- make_nodes()
  for (a in names(nodes)) for (b in names(nodes)) {
    expect_ref(paste0("pmin-", a, "-", b), pmin(nodes[[a]], nodes[[b]]))
    expect_ref(paste0("pmax-", a, "-", b), pmax(nodes[[a]], nodes[[b]]))
  }
  x <- nodes$VU; x[c(1, 5, 9)] <- NA
  y <- nodes$V; y[2] <- NA
  expect_ref("pmin-na", list(pmin(x, y), pmin(x, y, na.rm = TRUE), pmax(x, y), pmax(x, y, na.rm = TRUE)))
  expect_ref("pmin-three", list(pmin(nodes$V, nodes$U, nodes$VU), pmax(nodes$V, nodes$U, nodes$VU)))
  expect_ref("pmin-number", list(pmin(nodes$VU, 3), pmax(nodes$VU, 3), pmin(nodes$V, 0.5)))
  # BUG: a single argument returns a list
  expect_ref("pmin-single", list(pmin(nodes$V), pmax(nodes$V)))
  expect_ref("pmin-default", list(pmin(1:5, 3), pmax(1:5, 3), pmin(c(1, NA, 3), 2, na.rm = TRUE)))
  expect_ref("pmin-incompatible", pmin(nodes$VUM, mcdata(1, "V", nsv = 7, nvariates = 2)))
})

test_that("mcapply", {
  nodes <- make_nodes()
  for (nm in names(nodes)) {
    x <- nodes[[nm]]
    for (mg in c("all", "var", "unc", "variates")) {
      expect_ref(paste0("mcapply-", nm, "-", mg, "-sum"), mcapply(x, mg, sum))
      expect_ref(paste0("mcapply-", nm, "-", mg, "-exp"), mcapply(x, mg, exp))
      expect_ref(paste0("mcapply-", nm, "-", mg, "-sort"), mcapply(x, mg, sort))
      expect_ref(paste0("mcapply-", nm, "-", mg, "-cumsum"), mcapply(x, mg, cumsum))
    }
  }
  expect_ref("mcapply-dots", mcapply(nodes$VU, "var", quantile, probs = .5))
})

test_that("mc: building and naming", {
  nodes <- make_nodes()
  a <- nodes$V; b <- nodes$VU; d <- nodes$U
  expect_ref("mc-basic", mc(a, b, d))
  expect_ref("mc-name", mc(a, b, name = c("first", "second")))
  expect_ref("mc-named-args", mc(x = a, y = b))
  expect_ref("mc-nested", mc(mc(a, b), d))
  expect_ref("mc-devname", mc(mc(a, b), d, devname = TRUE))
  expect_ref("mc-duplicated", mc(a, a))
  expect_ref("mc-attr", list(attr(mc(a), "type"), attr(mc(a, d), "type"), attr(mc(nodes[["0"]]), "type")))
  expect_ref("mc-incompatible", mc(a, mcdata(1, "V", nsv = 3)))
  expect_ref("mc-bad", mc(a, 1:3))
  expect_ref("mc-empty", mc())
})

test_that("mcprobtree", {
  with_dims(7, 5, {
    c1 <- mcstoc(rnorm, type = "VU", mean = 10, sd = 2)
    c2 <- mcstoc(runif, type = "VU", min = -6, max = -5)
    c3 <- mcdata(0, type = "VU")
    w <- mcstoc(rempiricalD, type = "VU", values = 1:3, prob = c(.75, .20, .05))
    expect_ref("mcprobtree-node", mcprobtree(w, list("1" = c1, "2" = c2, "3" = c3), type = "VU"))
    expect_ref("mcprobtree-weights",
               mcprobtree(c(.75, .20, .05), list("1" = c1, "2" = c2, "3" = c3), type = "VU"))
    expect_ref("mcprobtree-seed",
               mcprobtree(c(.5, .5), list("1" = c1, "2" = c2), type = "VU", seed = 3), seed = NULL)
    v1 <- mcstoc(rnorm, "V"); u1 <- mcstoc(rnorm, "U"); z1 <- mcdata(5, "0")
    expect_ref("mcprobtree-mixed", mcprobtree(c(.3, .3, .4), list("1" = v1, "2" = u1, "3" = z1), type = "VU"))
    expect_ref("mcprobtree-V", mcprobtree(c(.3, .7), list("1" = v1, "2" = z1), type = "V"))
    expect_ref("mcprobtree-U", mcprobtree(c(.3, .7), list("1" = u1, "2" = z1), type = "U"))
    expect_ref("mcprobtree-numbers", mcprobtree(c(.3, .7), list("1" = 1, "2" = 2), type = "V"))
    vals <- list("1" = c1, "2" = c2)
    # BUG: a list given as a variable is not handled
    expect_ref("mcprobtree-listvar", mcprobtree(c(.5, .5), vals, type = "VU"))
    # BUG: a single branch fails
    expect_ref("mcprobtree-single", mcprobtree(1, list("5" = c1), type = "VU"))
    expect_ref("mcprobtree-badprob", mcprobtree(c(.5, -.5), list("1" = c1, "2" = c2), type = "VU"))
  })
})

test_that("mcprobtree with global variables", {
  ge <- globalenv()
  nm <- c("mc2dt_c1", "mc2dt_c2", "mc2dt_c3", "mc2dt_w", "mc2dt_v1", "mc2dt_u1", "mc2dt_z1")
  on.exit(rm(list = nm, envir = ge))
  with_dims(7, 5, {
    set.seed(10)
    assign("mc2dt_c1", mcstoc(rnorm, type = "VU", mean = 10, sd = 2), envir = ge)
    assign("mc2dt_c2", mcstoc(runif, type = "VU", min = -6, max = -5), envir = ge)
    assign("mc2dt_c3", mcdata(0, type = "VU"), envir = ge)
    assign("mc2dt_w", mcstoc(rempiricalD, type = "VU", values = 1:3, prob = c(.75, .20, .05)), envir = ge)
    assign("mc2dt_v1", mcstoc(rnorm, "V"), envir = ge)
    assign("mc2dt_u1", mcstoc(rnorm, "U"), envir = ge)
    assign("mc2dt_z1", mcdata(5L, "0"), envir = ge)
    expect_ref("mcprobtreeG-node", mcprobtree(mc2dt_w, list("1" = mc2dt_c1, "2" = mc2dt_c2, "3" = mc2dt_c3), type = "VU"))
    expect_ref("mcprobtreeG-weights",
               mcprobtree(c(.75, .20, .05), list("1" = mc2dt_c1, "2" = mc2dt_c2, "3" = mc2dt_c3), type = "VU"))
    expect_ref("mcprobtreeG-seed",
               mcprobtree(c(.5, .5), list("1" = mc2dt_c1, "2" = mc2dt_c2), type = "VU", seed = 3), seed = NULL)
    expect_ref("mcprobtreeG-mixed", mcprobtree(c(.3, .3, .4), list("1" = mc2dt_v1, "2" = mc2dt_u1, "3" = mc2dt_z1), type = "VU"))
    expect_ref("mcprobtreeG-V", mcprobtree(c(.3, .7), list("1" = mc2dt_v1, "2" = mc2dt_z1), type = "V"))
    expect_ref("mcprobtreeG-U", mcprobtree(c(.3, .7), list("1" = mc2dt_u1, "2" = mc2dt_z1), type = "U"))
    expect_ref("mcprobtreeG-0", mcprobtree(c(.3, .7), list("1" = mc2dt_z1, "2" = mcdata(2L, "0")), type = "0"))
    expect_ref("mcprobtreeG-integer", mcprobtree(c(.3, .7), list("1" = mc2dt_z1, "2" = mcdata(2L, "V")), type = "V"))
    expect_ref("mcprobtreeG-calls", mcprobtree(c(.3, .7), list("1" = mcstoc(rnorm, "V"), "2" = mcstoc(runif, "V")), type = "V"))
    expect_ref("mcprobtreeG-nv", mcprobtree(c(.3, .7), list("1" = mcstoc(rnorm, "V", nvariates = 2), "2" = mc2dt_z1),
                                            type = "V", nvariates = 2))
    expect_ref("mcprobtreeG-outm", mcprobtree(c(.3, .7), list("1" = mc2dt_v1, "2" = mc2dt_z1), type = "V", outm = "mean"))
    expect_ref("mcprobtreeG-badname", mcprobtree(c(.3, .7), list("a" = mc2dt_v1, "2" = mc2dt_z1), type = "V"))
    expect_ref("mcprobtreeG-badtype", mcprobtree(c(.3, .7), list("1" = mc2dt_u1, "2" = mc2dt_z1), type = "V"))
  })
})

test_that("cornode", {
  set.seed(1)
  mat <- cbind(rnorm(50), rnorm(50), rnorm(50))
  corr <- matrix(c(1, .5, .2, .5, 1, .2, .2, .2, 1), ncol = 3)
  expect_ref("cornode-matrix", cornode(mat, target = corr))
  expect_ref("cornode-matrix-result", cornode(mat, target = corr, result = TRUE))
  expect_ref("cornode-matrix-outrank", cornode(mat, target = corr, outrank = TRUE))
  expect_ref("cornode-matrix-seed", cornode(mat, target = corr, seed = 4), seed = NULL)
  expect_ref("cornode-bad-target", cornode(mat, target = matrix(2, 3, 3)))
  with_dims(30, 4, {
    a <- mcstoc(rnorm, "V"); b <- mcstoc(rnorm, "V"); c2 <- mcstoc(rnorm, "VU"); d <- mcstoc(rnorm, "VU")
    expect_ref("cornode-V", cornode(a, b, target = 0.7))
    expect_ref("cornode-V-result", cornode(a = a, b = b, target = -0.5, result = TRUE))
    expect_ref("cornode-VU", cornode(c2, d, target = 0.6))
    expect_ref("cornode-V-VU", cornode(a, c2, d, target = matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3)))
    # BUG: a non-first "V" is not rejected
    expect_ref("cornode-VU-V", cornode(c2, a, target = 0.5))
    expect_ref("cornode-U", cornode(mcstoc(rnorm, "U"), mcstoc(rnorm, "U"), target = .5))
  })
})

test_that("extractvar, addvar, unmc", {
  nodes <- make_nodes()
  for (nm in c("0M", "VM", "UM", "VUM")) {
    expect_ref(paste0("extractvar-", nm), list(extractvar(nodes[[nm]], 1), extractvar(nodes[[nm]], 2:3)))
  }
  expect_ref("extractvar-outm", extractvar(mcstoc(rnorm, "V", nvariates = 2, outm = "mean", nsv = 7), 1))
  expect_ref("extractvar-bad", extractvar(nodes$VUM, 4))
  expect_ref("addvar", list(addvar(nodes$V, nodes$V), addvar(nodes$VU, nodes$VUM)))
  expect_ref("addvar-bad", addvar(nodes$V, nodes$U))
  for (nm in names(nodes))
    expect_ref(paste0("unmc-", nm), list(unmc(nodes[[nm]]), unmc(nodes[[nm]], drop = FALSE)))
  expect_ref("unmc-mc", list(unmc(mc(nodes$V, nodes$VU)), unmc(mc(nodes$V, nodes$VU), drop = FALSE)))
})

test_that("predicates, dimensions, types and outm", {
  nodes <- make_nodes()
  m <- mc(nodes$V, nodes$U, nodes$VU)
  expect_ref("is-mcnode", lapply(nodes, is.mcnode))
  bad <- nodes$V; attr(bad, "type") <- "U"
  bad2 <- nodes$V; class(bad2) <- NULL
  bad3 <- mcdata(1, "0"); dim(bad3) <- c(1, 1, 1)
  expect_ref("is-mcnode-bad", list(is.mcnode(bad), is.mcnode(bad2), is.mcnode(1:3), is.mcnode(bad3)))
  expect_ref("is-mc", list(is.mc(m), is.mc(nodes$V), is.mc(list(1)), is.mc(mc(nodes$V))))
  x <- nodes$VU; x[2] <- NA; x[3] <- NaN; x[4] <- Inf
  expect_ref("is-na", list(is.na(x), is.nan(x), is.finite(x), is.infinite(x)))
  expect_ref("dimmc", list(dimmc(m), lapply(nodes, dimmcnode)))
  expect_ref("typemcnode", list(lapply(nodes, typemcnode), lapply(nodes, typemcnode, index = TRUE)))
  expect_ref("outm-get", list(outm(nodes$V), outm(m, which.node = 2)))
  expect_ref("outm-set", list(outm(nodes$V, "none"), outm(nodes$V, "mean")))
  # BUG: outm on an mc object does not change the node
  expect_ref("outm-set-mc", outm(m, "none", which.node = 2))
  # BUG: a vector of functions errors with R >= 4.3
  expect_ref("outm-set-vector", outm(nodes$V, c("min", "max")))
  expect_ref("outm-set-bad", outm(nodes$V, "notAFunction"))
})

test_that("ndvar / ndunc", {
  ov <- ndvar(); ou <- ndunc()
  on.exit({ ndvar(ov); ndunc(ou) })
  expect_ref("ndvar-ndunc", list(ndvar(123), ndvar(), ndunc(45), ndunc()), seed = NULL)
  expect_ref("ndvar-bad", ndvar(-1), seed = NULL)
})

test_that("lhs and rtrunc", {
  expect_ref("lhs-runif", lhs(nsv = 10, nsu = 3))
  expect_ref("lhs-rnorm", lhs(distr = "rnorm", nsv = 10, nsu = 3, mean = 5))
  expect_ref("lhs-nvariates", lhs(distr = "runif", nsv = 6, nsu = 2, nvariates = 2))
  expect_ref("lhs-function", lhs(distr = rnorm, nsv = 5, nsu = 2))
  expect_ref("rtrunc-norm", rtrunc(rnorm, 20, linf = -1, lsup = 1))
  expect_ref("rtrunc-string", rtrunc("rnorm", 20, linf = 0, mean = 2))
  expect_ref("rtrunc-lnorm", rtrunc(rlnorm, 20, lsup = 1))
  expect_ref("rtrunc-pois", rtrunc(rpois, 20, linf = 2, lambda = 3))
  expect_ref("rtrunc-vector", rtrunc(rnorm, 10, linf = 1:10, lsup = 2:11))
  expect_ref("rtrunc-bad", rtrunc(rnorm, 10, linf = 2, lsup = 1))
  expect_ref("rtrunc-unif-out", rtrunc(runif, 10, linf = 2, lsup = 4, max = 1))
  # BUG: default distr fails
  expect_ref("rtrunc-default", rtrunc(n = 10, linf = .5))
})

test_that("characterization of invariants (run everywhere)", {
  nodes <- make_nodes()
  for (nm in names(nodes)) {
    x <- nodes[[nm]]
    expect_s3_class(x, "mcnode")
    expect_true(is.mcnode(x))
    y <- x + x
    expect_identical(attributes(y)[c("dim", "class", "type")], attributes(x)[c("dim", "class", "type")])
    expect_equal(unclass(y), unclass(x) * 2, ignore_attr = TRUE)
  }
  expect_identical(typemcnode(nodes$V + nodes$U), "VU")
  expect_identical(dim(nodes$V * nodes$UM), c(7L, 5L, 3L))
  expect_identical(dim(nodes$V * nodes$U), c(7L, 5L, 1L))
})
