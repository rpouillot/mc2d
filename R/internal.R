# Internal helpers (not exported).
# They avoid full copies of large arrays and R-level loops on large mcnodes.

# Set the attributes of an mcnode in a single step (no copy of the data).
# The data is recycled only if its length differs from prod(dim).
.as_mcnode <- function(x, dim, type, outm) {
  n <- prod(dim)
  if (length(x) != n) x <- rep_len(as.vector(x), n)
  attributes(x) <- list(dim = as.integer(dim), class = "mcnode", type = type, outm = outm)
  x
}

# Same as array(x, dim) but without a copy when x is a plain vector of the right length
.arr <- function(x, dim) {
  if (is.null(attributes(x)) && length(x) == prod(dim)) {
    dim(x) <- dim
    x
  } else array(x, dim = dim)
}

# Expand the first two dimensions of a "U" (by row), "V" (by column) or "0" node
# to nsv x nsu. Returns a plain (nsv*nsu) x nvariates matrix.
# Same result as apply(x, 3, matrix, nrow = nsv, ncol = nsu, byrow = (type == "U"))
.recycle_node <- function(x, from, nsv, nsu) {
  d <- dim(x)
  nvar <- if (length(d) >= 3) d[3] else 1L
  x <- as.vector(x)
  x <- switch(from,
    "U"  = rep(x, each = nsv),
    "V"  = if (nsu == 1) x else {           # columns of the nsv x nvar matrix, each one repeated nsu times
      dim(x) <- c(nsv, nvar)
      x[, rep(seq_len(nvar), each = nsu), drop = FALSE]
    },
    "0"  = rep(x, each = nsv * nsu),
    "VU" = x)
  dim(x) <- c(nsv * nsu, nvar)
  x
}

# Apply the "outm" function 'fname' over the variates of a node (third dimension).
# Same result as: node <- apply(x, c(1, 2), get(fname, mode = "function")); dim(node) <- c(dim(node), 1)
# with vectorized versions of "sum", "min", "max" (identical results) and "mean"
# (rowMeans: may differ from mean() in the last bit) for double nodes.
.outm_reduce <- function(x, fname) {
  d <- dim(x)
  nv <- d[3]
  res <- NULL
  if(is.double(x) && fname %in% c("sum", "mean", "min", "max") &&
     identical(get(fname, mode = "function"), get(fname, envir = baseenv(), mode = "function"))) {
    res <- switch(fname,
      sum  = rowSums(x, dims = 2),
      mean = rowMeans(x, dims = 2),
      min  = , max = {
        sl <- lapply(seq_len(nv), function(k) x[, , k])
        r <- do.call(if(fname == "min") base::pmin else base::pmax, c(sl, list(na.rm = TRUE)))
        nas <- Reduce("|", lapply(sl, function(s) is.na(s) & !is.nan(s)))
        nans <- Reduce("|", lapply(sl, is.nan))
        r[nans] <- NaN                       # NaN if any NaN ...
        r[nas] <- NA                         # ... but NA if any NA, as min() / max()
        r
      })
  }
  if(is.null(res)) res <- apply(x, c(1, 2), get(fname, mode = "function"))
  dim(res) <- c(d[1:2], 1L)
  res
}

# Quantiles of each column of a 2D/3D array (first variate if 3D).
# Same results as apply(x[, , 1, drop = FALSE], 2, quantile, probs = probs, ...)
# without the full copy made by apply() / aperm().
.col_quantile <- function(x, probs, ...) {
  d <- dim(x)
  nr <- d[1]
  vapply(seq_len(d[2]), function(j)
    stats::quantile(x[(j - 1L) * nr + seq_len(nr)], probs = probs, ...),
    numeric(length(probs)))
}

# Name of a function given as argument (a name, a string, or pkg::name).
# 'expr' is the result of substitute(arg). Returns 'default' if not found.
.fun_name <- function(expr, default) {
  if (is.null(expr) || (is.symbol(expr) && identical(as.character(expr), ""))) return(default)
  if (is.character(expr)) return(expr[1])
  if (is.symbol(expr)) return(as.character(expr))
  if (is.call(expr) && as.character(expr[[1]])[1] %in% c("::", ":::")) return(as.character(expr[[3]]))
  default
}

# Check the outm argument: "none", "each" or a vector of names of functions.
.check_outm <- function(outm) {
  is.character(outm) && length(outm) > 0 &&
    (identical(outm, "none") || identical(outm, "each") ||
       all(vapply(outm, exists, NA, mode = "function")))
}

# Number of graphs for the nodes of an mc object: nvariates for outm "each" (or NULL),
# 0 for "none", one per function otherwise. 'loutm': list of outm, 'nvar': nvariates of each node.
.n_outm_panels <- function(loutm, nvar) {
  sum(mapply(function(o, nv) if(is.null(o) || o[1] == "each") nv else if(o[1] == "none") 0 else length(o),
             loutm, nvar))
}

# Arrange a list of ggplot objects on a grid (nrow x ncol, filled by row).
# Replaces ggpubr::ggarrange (to avoid the dependency). Returns a ggplot object,
# or a list of ggplot objects (one per page) if there are more plots than cells.
# As ggpubr::ggarrange, "ggarrange" is added to the class (compatibility with mc2d <= 0.2.2).
.ggarrange <- function(plotlist, ncol, nrow, tag = TRUE) {
  plotlist <- plotlist[!vapply(plotlist, is.null, NA)]
  ncell <- ncol * nrow
  if (length(plotlist) > ncell) {
    pages <- split(plotlist, ceiling(seq_along(plotlist) / ncell))
    res <- lapply(pages, .ggarrange, ncol = ncol, nrow = nrow, tag = FALSE)   # pages not tagged, as in ggpubr
    class(res) <- c(class(res), "ggarrange")
    return(res)
  }
  grobs <- lapply(plotlist, ggplot2::ggplotGrob)
  grobs <- c(grobs, rep(list(grid::nullGrob()), ncell - length(grobs)))
  tab <- gtable::gtable_matrix("mc2d_arrange",
                               matrix(grobs, nrow = nrow, ncol = ncol, byrow = TRUE),
                               widths = grid::unit(rep(1, ncol), "null"),
                               heights = grid::unit(rep(1, nrow), "null"))
  res <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::annotation_custom(tab, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    ggplot2::theme_void()
  if (tag) class(res) <- c(class(res), "ggarrange")
  res
}
