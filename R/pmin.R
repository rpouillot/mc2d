#' Parallel Maxima and Minima for mcnodes
#'
#' Returns the parallel maxima and minima of the input values, extended to \samp{mcnode} objects.
#' The same rules as in \code{\link{Ops.mcnode}} are applied for type determination.
#'
#' \samp{pmax} and \samp{pmin} take one or more \samp{mcnode} and/or vectors as arguments and
#' return a \samp{mcnode} of adequate type and size giving the "parallel" maxima or minima.
#' Note that the first element of \samp{...} should be an \samp{mcnode}.
#'
#' @param ... one or more \samp{mcnode}s or \samp{mcnode}s and compatible vectors.
#'   The first element must be an \samp{mcnode}.
#' @param na.rm a logical indicating whether missing values should be removed.
#' @return an \samp{mcnode} of adequate type and dimension.
#' @seealso \code{\link{min}}, \code{\link{Ops.mcnode}}
#' @keywords utilities
#' @name pmin.mcnode
#' @examples
#' ndvar(10)
#' ndunc(21)
#' x <- mcstoc(rnorm, "V")
#' pmin(x, 0)
#' y <- mcdata(rep(c(-1, 1), length=ndunc()), "U")
#' unclass(pmin(x, y))
#' @export
pmin <- function (..., na.rm = FALSE)
 UseMethod("pmin")

#' @rdname pmin.mcnode
#' @export
pmax <- function (..., na.rm = FALSE)
 UseMethod("pmax")

#' @rdname pmin.mcnode
#' @export
pmin.default <- function(..., na.rm = FALSE) base::pmin(..., na.rm = na.rm)

#' @rdname pmin.mcnode
#' @export
pmax.default <- function(..., na.rm = FALSE) base::pmax(..., na.rm = na.rm)

#' @rdname pmin.mcnode
#' @export
pmin.mcnode <- function(..., na.rm = FALSE)
  .pminmax_mcnode(list(...), na.rm = na.rm, min = TRUE)

#' @rdname pmin.mcnode
#' @export
pmax.mcnode <- function(..., na.rm = FALSE)
  .pminmax_mcnode(list(...), na.rm = na.rm, min = FALSE)

# Internal: parallel minima (min = TRUE) or maxima (min = FALSE) of mcnodes
.pminmax_mcnode <- function(argsd, na.rm = FALSE, min = TRUE)
{
    largsd <- length(argsd)
    if(largsd == 1) return(argsd[[1]])

    outm <- attr(argsd[[1]],"outm")

    typemc <- lapply(argsd, attr, which = "type")

    type <- unlist(typemc)
    type <- if(all(type == "0")) "0" else
        if(all(type %in% c("0","V"))) "V" else
            if(all(type %in% c("0","U"))) "U" else "VU"

# evaluate the minimal common to build the minimal recycling level...
    maxdim1 <- max(unlist(lapply(argsd,function(x) dim(x)[1])))
    maxdim2 <- max(unlist(lapply(argsd,function(x) dim(x)[2])))
    dim3 <- unlist(lapply(argsd,function(x) dim(x)[3]))
    maxdim3 <- max(dim3)
    if(any(dim3 != 1 & dim3 != maxdim3)) stop("Incompatible mcnode dimensions")

#### A function to deal mcnodes as arguments, mostly copied from Ops.mcnode

    LAFUNC <- function(argsd,typethismc){

        if(!is.null(typethismc)){                                    #mcnode as arguments

            dimm <- dim(argsd)
            if ((typethismc == "V" && dimm[1] != maxdim1) ||
                (typethismc == "U" && dimm[2] != maxdim2) ||
                (typethismc == "VU" && (dimm[1] != maxdim1 || dimm[2] != maxdim2)))
                stop("Nodes of incompatible dimensions") # incompatible dimension

            if(maxdim3 > 1)  #at least one multivariate node as parameter, need recycling on the third dimension
                argsd <- .recycle_node(argsd, typethismc, maxdim1, maxdim2)   # recycling as matrix (maxdim1*maxdim2) x nvariates
            else {
                dim(argsd) <- NULL    # as vector
                if(typethismc == "U" && maxdim1 != 1) argsd <- rep(argsd, each = maxdim1)     #recycling U as vector nsv*nsu
            }
        }
        else if(is.array(argsd)) stop("Array prohibited in pmin as parameter. Use an mcnode or a vector instead")
        else if(is.vector(argsd)){
            l <- length(argsd)
            if(l != 1 && l != maxdim1*maxdim2 && l!= maxdim1*maxdim2*maxdim3) stop("The vector size is not compatible
                    with the node dimension. Length should be 1 or n=",maxdim1*maxdim2," or n=",maxdim1*maxdim2*maxdim3)
        }
        return(as.vector(argsd))
    }

####

    argsd <- mapply(LAFUNC, argsd=argsd, typethismc=typemc, SIMPLIFY=FALSE)

    if(all(vapply(argsd, is.double, NA))) {
      # Fast path (C code of base::pmin / base::pmax), with the NA rules of the R code below
      mmm <- do.call(if(min) base::pmin else base::pmax, c(argsd, list(na.rm = TRUE)))
      n <- length(mmm)
      nas <- lapply(argsd, function(x) rep_len(is.na(x), n))
      anyna <- Reduce("|", nas)
      if(any(anyna)) {
        if(!na.rm) mmm[anyna] <- NA
        else {
          allna <- Reduce("&", nas)
          if(any(allna)) mmm[allna] <- rep_len(argsd[[largsd]], n)[allna]   # all NA: the last NA value
        }
      }
    }
    else {
      # Copied from function pmin
      mmm <- argsd[[1L]]
      has.na <- FALSE
      for (each in argsd[-1L]) {
        if (length(mmm) < (m <- length(each)))
          mmm <- rep(mmm, length.out = m)
        else if (length(each) < (m <- length(mmm)))
          each <- rep(each, length.out = m)
        nas <- cbind(is.na(mmm), is.na(each))
        if (has.na || (has.na <- any(nas))) {
          mmm[nas[, 1L]] <- each[nas[, 1L]]
          each[nas[, 2L]] <- mmm[nas[, 2L]]
        }
        change <- if(min) mmm > each else mmm < each
        change <- change & !is.na(change)
        mmm[change] <- each[change]
        if (has.na && !na.rm)
          mmm[nas[, 1L] | nas[, 2L]] <- NA
      }
    }

    .as_mcnode(mmm, c(maxdim1, maxdim2, maxdim3), type, outm)
}
