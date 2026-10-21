## Helpers for the characterization ("golden reference") tests.
##
## expect_ref() evaluates an expression (optionally after set.seed()), captures
## its value, its warnings and its error (if any) and compares the result with
## a reference stored in tests/testthat/_ref/<name>.rds using expect_identical().
## When the reference file does not exist yet, it is created.
## These strict comparisons are skipped on CRAN (platform dependent rounding).

ref_dir <- function() testthat::test_path("_ref")

capture_ref <- function(expr) {
  warns <- character()
  val <- withCallingHandlers(
    tryCatch(expr,
             error = function(e) structure(list(message = conditionMessage(e)),
                                           class = "ref_error")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  # Base R recycling warnings differ between R versions (R-devel adds
  # "object length is not a multiple of subscript length"): not compared.
  warns <- unique(warns[!grepl("is not a multiple of", warns)])
  list(value = val, warnings = warns)
}

expect_ref <- function(name, expr, seed = 1234) {
  testthat::skip_on_cran()
  # The references are not shipped in the package tarball (.Rbuildignore)
  if (!dir.exists(ref_dir())) testthat::skip("no reference directory")
  if (!is.null(seed)) set.seed(seed)
  obs <- capture_ref(expr)
  f <- file.path(ref_dir(), paste0(name, ".rds"))
  if (!file.exists(f)) {
    dir.create(ref_dir(), showWarnings = FALSE)
    saveRDS(obs, f, version = 2)
    testthat::succeed(paste("reference created:", name))
    return(invisible(obs))
  }
  ref <- readRDS(f)
  # Developer helper: names of the references that differ
  log <- Sys.getenv("MC2D_REF_LOG")
  if (nzchar(log) && !identical(obs, ref)) {
    cat(name, "\n", file = log, append = TRUE)
    nd <- file.path(dirname(log), "new")
    dir.create(nd, showWarnings = FALSE)
    saveRDS(obs, file.path(nd, paste0(name, ".rds")), version = 2)
  }
  testthat::expect_identical(obs, ref, label = name)
  invisible(obs)
}

## Evaluate code with temporary ndvar() / ndunc()
with_dims <- function(nsv, nsu, code) {
  ov <- ndvar(); ou <- ndunc()
  on.exit({ ndvar(ov); ndunc(ou) })
  ndvar(nsv); ndunc(nsu)
  code
}

## A set of nodes of every type, univariate and multivariate (nvariates = 3)
make_nodes <- function(nsv = 7, nsu = 5, seed = 1) {
  with_dims(nsv, nsu, {
    set.seed(seed)
    list(
      "0"   = mcstoc(runif, type = "0", min = 1, max = 2),
      V     = mcstoc(rnorm, type = "V"),
      U     = mcstoc(rnorm, type = "U", mean = 2),
      VU    = mcstoc(rnorm, type = "VU", mean = 3),
      "0M"  = mcstoc(runif, type = "0", nvariates = 3, min = 1, max = 2),
      VM    = mcstoc(rnorm, type = "V", nvariates = 3),
      UM    = mcstoc(rnorm, type = "U", nvariates = 3, mean = 2),
      VUM   = mcstoc(rnorm, type = "VU", nvariates = 3, mean = 3)
    )
  })
}

## Summarise a ggplot object by the coordinates of its built layers
gg_layers <- function(p) {
  b <- ggplot2::ggplot_build(p)
  keep <- c("x", "y", "xmin", "xmax", "ymin", "ymax", "xend", "yend", "group", "PANEL")
  lapply(b$data, function(d) {
    d <- as.data.frame(d)[intersect(names(d), keep)]
    if ("PANEL" %in% names(d)) d$PANEL <- as.integer(d$PANEL)
    d
  })
}

## Evaluate a base-graphics call on a null device
on_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  code
}
