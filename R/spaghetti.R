#' @title Spaghetti Plot of mc/mcnode Object
#' @description Use plot to draw spaghetti plots for the mc/mcnode objects.
#' @param x mc/mcnode object  
#' @param griddim a vector of two integers, indicating the size of the grid of the graph. If NULL, the grid is calculated to produce a "nice" graph.
#' @param xlab vector of labels for the x-axis. If NULL, use the name of the node.
#' @param ylab vector of labels for the y-axis.
#' @param main vector of main titles of the graph.
#' @param maxlines the maximum number of ecdf to draw.
#' @param ... further arguments to be passed to plot.stepfun()
#' @aliases spaghetti.mc
#' @aliases spaghetti.mcnode
#' @examples 
#' data(total)
#' spaghetti(mc(xVUM))
#' spaghetti(xVUM)
#' 
#' @rdname spaghetti
#' @export
spaghetti <- function(x,...){
  UseMethod("spaghetti")
}

#' @rdname spaghetti
#' @method spaghetti mc
#' @export
spaghetti.mc <- function(x, 
                         griddim = NULL, 
                         xlab = names(x),
                         ylab = "F(n)", 
                         main = "",  
                         maxlines = 100,
                         ...)
{
  beau <- function(n){
    nc <- round(sqrt(n))
    nr <- ceiling(n/nc)
    c(nc,nr)
  }
  
  l <- length(x)
  main <- rep(main,l)
  xlab <- rep(xlab,l)
  ylab <- rep(ylab,l)
  
  
  loutm <- lapply(x,attr,which="outm")
  dimm <- sapply(x,dim)
  n <- .n_outm_panels(loutm, dimm[3,])
  
  
  if(is.null(griddim)) griddim <- beau(n)
  if(prod(griddim) < n) {
    op <- par(mfrow=griddim,ask=TRUE,mar=c(5,4,.2,.2))
  } else {
    op <- par(mfrow=griddim,mar=c(5,4,.2,.2))
  }
  
  
  try({  #to restore par in case of error
    
    for(i in 1:l){
      
      if(is.null(loutm[[i]])) loutm[[i]] <- "each"
      if(loutm[[i]][1] == "none") next                                             # Pass outm == none
      for(j in loutm[[i]]){                                                     # loop on the nb of stat, j is the name of the function
        
        xi <- x[[i]]
        if(j == "each"){
          nvar <- dim(x[[i]])[3]
          if(nvar==1) xlab2 <- xlab[i] else xlab2 <- paste(xlab[i],1:nvar,sep="")
        } else {
          xi <- .outm_reduce(x[[i]], j)                                         # apply the function on the original node
          nvar <- 1                                                             #1 dimension now for this stat
          xlab2 <- paste(j,xlab[i])                                             #change the name with the name of the stat
        }
        
        if(is.logical(xi)) xi <- unclass(xi) * 1                                  # unclass to avoid Ops.mcnode
        
        for(k in 1:nvar){
          
          value <- as.matrix(unmc(xi[,,k])) # nvariates data
          ncol.value <- min(maxlines, dim(value)[2])
          rg <- suppressWarnings(range(value, na.rm = TRUE))
          xlim <- if(is.finite(rg[1]) && rg[1] < 0) rg else c(0, rg[2])        # from 0 for positive values
          plot.stepfun(value[!is.na(value[, 1]), 1],
                       xlim=xlim,
                       main = main[i],
                       xlab = xlab2[k],
                       ylab = ylab[i],
                       do.points = FALSE)
          if(ncol.value>1){
            for (p in 2: ncol.value){
              plot.stepfun(value[!is.na(value[, p]), p], do.points = FALSE, add=T)}}}
      }
    }
  })
  par(op)
  return(invisible())
}

#' @rdname spaghetti
#' @method spaghetti mcnode
#' @export
spaghetti.mcnode <- function(x, ...)
{ 
  nom <- deparse(substitute(x), width.cutoff = 500L, nlines=1)
  x <- list(x)
  names(x) <- nom
  spaghetti.mc(x, ...)
}


