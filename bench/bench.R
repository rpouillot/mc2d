## Timing of mc2d hot paths on large nodes (not part of the package build).
## Usage (from the package root):
##   Rscript bench/bench.R [label]
## Prints a table of elapsed times (seconds) and appends it to bench/results.csv

suppressMessages(pkgload::load_all(".", quiet = TRUE))
label <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(label)) label <- "run"

nsv <- 1e5; nsu <- 200
ndvar(nsv); ndunc(nsu)
set.seed(1)
xV   <- mcstoc(rnorm, "V")
xU   <- mcstoc(rnorm, "U")
xVU  <- mcstoc(rnorm, "VU")
xVUM <- mcstoc(rnorm, "VU", nvariates = 3)
xVUo <- mcstoc(rnorm, "VU", nvariates = 3, outm = "mean")
m    <- mc(xV = xV, xU = xU, xVU = xVU, out = xV + xU * xVU)

tm <- function(expr) {
  gc(FALSE)
  unname(system.time(expr)["elapsed"])
}

res <- c(
  "mcstoc VU (rnorm)"           = tm(mcstoc(rnorm, "VU")),
  "mcstoc VU, U param"          = tm(mcstoc(rnorm, "VU", mean = xU)),
  "mcstoc VU nv=3, U param"     = tm(mcstoc(rnorm, "VU", mean = mcstoc(rnorm, "U", nvariates = 3), nvariates = 3)),
  "Ops V + U"                   = tm(xV + xU),
  "Ops VU * U"                  = tm(xVU * xU),
  "Ops VUM + U"                 = tm(xVUM + xU),
  "Ops VUM + V nv=3"            = tm(xVUM + mcstoc(rnorm, "V", nvariates = 3)),
  "pmin VU, U"                  = tm(pmin(xVU, xU)),
  "pmax VUM, V"                 = tm(pmax(xVUM, xV)),
  "mcdata U -> VU"              = tm(mcdata(xU, "VU")),
  "mcdata V -> VU nv=3"         = tm(mcdata(mcstoc(rnorm, "V", nvariates = 3), "VU", nvariates = 3)),
  "mcapply var sum"             = tm(mcapply(xVU, "var", sum)),
  "summary mc"                  = tm(summary(m)),
  "summary outm mean nv=3"      = tm(summary(xVUo)),
  "quantile mc"                 = tm(quantile(m)),
  "plot(draw=FALSE)"            = tm(plot(m, draw = FALSE)),
  "tornado spearman"            = tm(tornado(m)),
  "tornadounc"                  = tm(tornadounc(m)),
  "print mc"                    = tm(capture.output(print(m))),
  "is.mc"                       = tm(is.mc(m)),
  "unmc"                        = tm(unmc(m)),
  "rpert 1e7"                   = tm(rpert(1e7, 0, 1, 4)),
  "rtriang 1e7"                 = tm(rtriang(1e7, 0, 1, 4)),
  "rbetagen 1e7"                = tm(rbetagen(1e7, 2, 3, 1, 4)),
  "qpert 1e6 vec param"         = tm(qpert(runif(1e6), runif(1e6), 1, 4)),
  "dempiricalD 1e5"             = tm(dempiricalD(sample(1:5, 1e5, TRUE), 1:5)),
  "qempiricalD 1e5"             = tm(qempiricalD(runif(1e5), 1:5)),
  "rempiricalC 2e4 vec min"     = tm(rempiricalC(2e4, min = runif(2e4, -1, 0), max = 10, values = 1:5)),
  "ddirichlet 2e5"              = tm(ddirichlet(rdirichlet(2e5, c(1, 2, 3)), c(1, 2, 3))),
  "dmultinomial 2e5"            = tm(dmultinomial(rmultinomial(2e5, 10, c(.2, .3, .5)), prob = c(.2, .3, .5))),
  "rmultinormal 2e4 vec sigma"  = tm(rmultinormal(2e4, c(0, 0), cbind(runif(2e4, 1, 2), 0, 0, 1)))
)

## converg is quadratic in the current version: smaller node
ndvar(5000); ndunc(50)
xVc <- mcstoc(rnorm, "VU")
pdf(NULL)
res <- c(res,
  "converg var (nsv=5000)"     = tm(converg(xVc, margin = "var")),
  "spaghetti (nsv=1e5)"        = tm(spaghetti(xVU)),
  "hist (nsv=1e5)"             = tm(hist(xVUo))
)
dev.off()

## evalmccut
modcut <- mcmodelcut({
  { a <- mcstoc(rnorm, "V"); b <- mcstoc(rnorm, "U") }
  { c1 <- a + b; res <- mc(a, b, c1) }
  { list(sum = summary(res), plot = plot(res, draw = FALSE), node = res$c1) }
})
res <- c(res, "evalmccut nsv=1e4 nsu=200" = tm(capture.output(evalmccut(modcut, nsv = 1e4, nsu = 200, seed = 1))))

out <- data.frame(label = label, R = paste(R.version$major, R.version$minor, sep = "."),
                  op = names(res), seconds = round(res, 3), row.names = NULL)
print(out[, c("op", "seconds")], row.names = FALSE)
f <- file.path("bench", "results.csv")
write.table(out, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
