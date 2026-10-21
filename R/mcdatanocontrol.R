#' Build mcnode Objects Without Dimension Control
#'
#' \samp{mcdatanocontrol} is a fast but unchecked version of \code{\link{mcdata}} which
#' forces the dimension of data to be \samp{(nsv x nsu x nvariates)} and sets the attributes
#' and class without any validation. Useful when your model is tested and speed is important.
#'
#' @rdname mcnode
#' @export
mcdatanocontrol <- function(data, type=c("V","U","VU","0"), nsv=ndvar(), nsu=ndunc(), nvariates=1, outm="each")
{
  type <- match.arg(type)
  n <- nsv*nsu*nvariates
  dim(data) <- NULL
  if(length(data) != n) data[1:n] <- data                 # recycle only when needed (avoids a copy)
  dim(data) <- c(nsv,nsu,nvariates)
  class(data) <- "mcnode"
  attr(data,which="type") <- type
  attr(data,which="outm") <- outm
  return(data)
  }
