## Characterization tests: base and ggplot2 graphics.
## See helper-ref.R. Comments tagged "BUG:" document a behaviour of
## mc2d <= 0.2.2 that is known to be wrong and is fixed in 0.2.3.

test_that("base graphics: plot, hist, spaghetti", {
  utils::data("total", package = "mc2d", envir = environment())
  on_null_device({
    expect_ref("plot-total", plot(total))
    expect_ref("plot-xVUM3", plot(xVUM3, lim = c(0.025, 0.975)))
    expect_ref("plot-xV", plot(xV))
    expect_ref("plot-xU", plot(xU))
    expect_ref("plot-x0", plot(x0))
    expect_ref("plot-nodraw", plot(total, draw = FALSE))
    expect_ref("plot-plotmc", plot(quantile(total)))
    expect_ref("plot-stat", plot(xVU, stat = "mean", lim = c(.1, .9)))
    expect_ref("plot-ylim", plot(xVU, ylim = c(0, 1), xlim = c(-1, 1)))
    # BUG: a list of xlim is rejected
    expect_ref("plot-xlim-list", plot(mc(xV, xVU), xlim = list(c(0, 1), c(0, 2))))
    expect_ref("hist-total", hist(total))
    expect_ref("hist-node", hist(xVUM))
    expect_ref("hist-xVUM3", hist(xVUM3))
    expect_ref("hist-options", hist(xVU, griddim = c(1, 1), xlab = "a", main = "b", breaks = 5))
    expect_ref("spaghetti-mc", spaghetti(mc(xVUM)))
    expect_ref("spaghetti-node", spaghetti(xVU))
    expect_ref("spaghetti-maxlines", spaghetti(xVU, maxlines = 5))
    expect_ref("spaghetti-xVUM3", spaghetti(xVUM3))
    expect_ref("plot-tornado", plot(tornado(total, use = "complete.obs")))
    expect_ref("plot-tornado-mean", plot(tornado(total, use = "complete.obs"), stat = "mean"))
    expect_ref("plot-tornadounc", plot(tornadounc(total, 10, use = "complete.obs"), 1, 1))
    # BUG: the default stat ("median") does not match the row names
    expect_ref("plot-tornadounc-default", plot(tornadounc(total, 10, use = "complete.obs")))
  })
})

test_that("ggplot2 graphics", {
  utils::data("total", package = "mc2d", envir = environment())
  gg <- function(p) if (inherits(p, "ggplot")) gg_layers(p) else lapply(p, function(z)
    if (inherits(z, "ggplot")) gg_layers(z) else class(z))
  on_null_device({
    expect_ref("ggplotmc-total", gg(ggplotmc(total)))
    expect_ref("ggplotmc-node", gg(ggplotmc(xVU)))
    expect_ref("ggplotmc-V", gg(ggplotmc(xV)))
    expect_ref("ggplotmc-xVUM3", gg(ggplotmc(xVUM3)))
    expect_ref("ggplotmc-lim", gg(ggplotmc(xVU, lim = c(.1, .9))))
    expect_ref("gghist-V", gg(gghist(xV)))
    expect_ref("gghist-VUM", gg(gghist(xVUM)))
    expect_ref("gghist-which", gg(gghist(xVUM, which = 1)))
    expect_ref("gghist-total", gg(gghist(total)))
    expect_ref("gghist-xVUM3", gg(gghist(xVUM3)))
    expect_ref("ggspaghetti-node", gg(ggspaghetti(xVU)))
    expect_ref("ggspaghetti-mc", gg(ggspaghetti(mc(xVUM))))
    expect_ref("ggspaghetti-maxlines", gg(ggspaghetti(xVU, maxlines = 5)))
    expect_ref("ggtornado", gg(ggtornado(tornado(total, use = "complete.obs"))))
    expect_ref("ggtornado-unc", gg(ggtornado(tornadounc(total, 10, use = "complete.obs"), which = 1)))
    expect_ref("ggtornadounc", gg(ggtornadounc(tornadounc(total, 10, use = "complete.obs"), stat = 1)))
    expect_ref("ggtornadounc-default", gg(ggtornadounc(tornadounc(total, 10, use = "complete.obs"))))
  })
})

test_that("graphics run without error (run everywhere)", {
  utils::data("total", package = "mc2d", envir = environment())
  on_null_device({
    expect_no_error(plot(total))
    expect_no_error(hist(xVU))
    expect_no_error(spaghetti(xVU))
    expect_no_error(plot(tornado(total, use = "complete.obs")))
  })
  on_null_device(expect_s3_class(ggplotmc(xVU), "ggplot"))
  on_null_device(expect_s3_class(gghist(xV), "ggplot"))
})
