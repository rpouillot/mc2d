#' The Discrete Empirical Distribution
#'
#' Density, distribution function, quantile function and random generation for a discrete
#' empirical distribution. This function is vectorized to accept different sets of
#' \samp{values} or \samp{prob}.
#'
#' If \samp{prob} is missing, the discrete distribution is obtained directly from the
#' vector of \samp{values}, otherwise \samp{prob} is used to weight the values.
#' \samp{prob} is normalized before use. Thus, \samp{prob} may be the count of each
#' \samp{values}. \samp{prob} values should be non negative and their sum should not be 0.
#'
#' \samp{values} and/or \samp{prob} may vary: in that case, \samp{values} and/or \samp{prob}
#' should be sent as matrices, the first row being used for the first element of \samp{x},
#' \samp{q}, \samp{p} or the first random value, the second row for the second element, ...
#' Recycling is permitted if the number of rows of \samp{prob} and \samp{values} are equal
#' or if the number of rows of \samp{prob} and/or \samp{values} is one.
#'
#' \samp{rempiricalD(n, values, prob)} with \samp{values} and \samp{prob} as vectors is
#' equivalent to \samp{sample(x=values, size=n, replace=TRUE, prob=prob)}.
#'
#' @note In the future, the functions should be written for non numerical values.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values. If \samp{length(n) > 1}, the length is taken to be
#'   the number required.
#' @param values vector or matrix of numerical values. See details.
#' @param prob optional vector or matrix of counts or probabilities. See details.
#' @param log,log.p logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.
#' @param lower.tail logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]},
#'   otherwise, \samp{P[X > x]}.
#' @return \samp{dempiricalD} gives the density, \samp{pempiricalD} gives the distribution
#' function, \samp{qempiricalD} gives the quantile function and \samp{rempiricalD} generates
#' random deviates.
#' @seealso \code{\link{sample}}, \code{\link{empiricalC}}.
#' @keywords distribution
#' @name empiricalD
#' @examples
#' dempiricalD(1:6, 2:6, prob=c(10, 10, 70, 0, 10))
#' pempiricalD(1:6, 2:6, prob=c(10, 10, 70, 0, 10))
#' qempiricalD(seq(0, 1, 0.1), 2:6, prob=c(10, 10, 70, 0, 10))
#' table(rempiricalD(10000, 2:6, prob=c(10, 10, 70, 0, 10)))
#'
#' ## Varying values
#' (values <- matrix(1:10, ncol=5))
#' dempiricalD(c(1, 1), values)
#' @export
dempiricalD <- function(x,values,prob=NULL,log=FALSE)
{
  if(length(x) == 0) return(numeric(0))
  if(is.vector(values)) values <- matrix(values,nrow=1)
  if(is.null(prob)) prob <- array(1,dim=dim(values))
  if(is.vector(prob)) prob <- matrix(prob,nrow=1)

  if(!is.numeric(values) || !is.numeric(prob) || !is.numeric(x)) stop("Non numeric arguments in a mathematical function")

  if(any(is.na(prob)) || any(prob < 0) || any(rowSums(prob) == 0))
    stop("Nas in prob, or sum(prob)=0 or negative values of prob")

  ncv <- ncol(values)
  ncp <- ncol(prob)
  if(ncv != ncp) stop("prob and values should be of same length or have the same number of columns.")

  l <- length(x)
  nrv <- nrow(values)
  nrp <- nrow(prob)
  mnr <- min(max(nrv,nrp),l)
  if(nrv != 1 && nrp != 1 && nrv != nrp)
    stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

  prob   <- lapply(1:min(nrp,mnr), function(x) prob[x,])
  values <- lapply(1:min(nrv,mnr), function(x) values[x,])

  prob <- mapply(function(pr,val) tapply(pr,val,sum), prob, values, SIMPLIFY=FALSE)
  prob <- lapply(prob, function(y) y /sum(y))
  values <- lapply(values,function(y) sort(unique(y)))

  res <- .empD_find(x, values, prob, "d", mnr)   # vectorized search
  if(is.null(res)) {
    find <- Vectorize(function(x,val,pr) {
      if(is.na(x)) return(NA)
      quel <- x==val
      if(any(quel)) return(pr[quel])
      return(0)}
      ,vectorize.args = "x")
  
    # to gain time and memory
    n1 <- l %/% mnr
    n2 <- l %% mnr
    l2 <- l-n2
    x1 <- x[1:l2]
  
    res <- vector(mode = "numeric", length = l)
    res[1:l2] <- mapply(find,x1,values,prob)
    if(n2!=0){
      x2 <- x[(l2+1):l]
      res[(l2+1):l] <- mapply(find,x2,values[1:n2],prob[1:n2])
      }
  
  }
  if(any(is.na(res))) warning("NaN in dempiricalD")
  if(log) res <- log(res)
  return(res)}

#' @rdname empiricalD
#' @export
pempiricalD <- function(q,values,prob=NULL,lower.tail = TRUE, log.p = FALSE)
{
  if(length(q) == 0) return(numeric(0))
  if(is.vector(values)) values <- matrix(values,nrow=1)
  if(is.null(prob)) prob <- array(1,dim=dim(values))
  if(is.vector(prob)) prob <- matrix(prob,nrow=1)

  if(!is.numeric(values) || !is.numeric(prob) || !is.numeric(q)) stop("Non numeric arguments in a mathematical function")

  if(anyNA(prob) || any(prob < 0) || any(rowSums(prob)==0))
    stop("Nas in prob, or sum(prob)=0 or negative values of prob")

  ncv <- ncol(values)
  ncp <- ncol(prob)
  if(ncv != ncp) stop("prob and values should be of same length or have the same number of columns.")

  l <- length(q)
  nrv <- nrow(values)
  nrp <- nrow(prob)
  mnr <- min(max(nrv,nrp),l)
  if(nrv != 1 && nrp != 1 && nrv != nrp)
    stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

  values <- lapply(1:min(nrv,mnr),function(y) values[y,])
  prob <- lapply(1:min(nrp,mnr),function(y) prob[y,])

  prob <- mapply(function(pr,val) tapply(pr,val,sum),prob,values,SIMPLIFY=FALSE)
  prob <- lapply(prob, function(x) cumsum(x)/sum(x))
  values <- lapply(values,function(y) sort(unique(y)))

  # the "-1" combined with the "<" allows duplicated values
  res <- .empD_find(q, values, prob, "p", mnr)   # vectorized search
  if(is.null(res)) {
    find <- Vectorize(function(x,val,pr){
          if(is.na(x)) return(NA)
          if(x < val[1]) return(0)
          return(max(pr[x >= val]))},vectorize.args = "x")
  
    # to gain time and memory
    n1 <- l %/% mnr
    n2 <- l %% mnr
    l2 <- l-n2
    q1 <- q[1:l2]
  
    res <- vector("numeric",length=l)
    res[1:l2] <- mapply(find,q1,values,prob)
    if(n2!=0){
      q2 <- q[(l2+1):l]
      res[(l2+1):l] <- mapply(find,q2,values[1:n2],prob[1:n2])
      }
  
  }
  if(any(is.na(res))) warning("NaN in pempiricalD")
  if(!lower.tail) res <- 1-res
  if(log.p) res <- log(res)
  return(res)}

#' @rdname empiricalD
#' @export
qempiricalD <- function(p,values,prob=NULL,lower.tail = TRUE, log.p = FALSE)
{
  if(length(p) == 0) return(numeric(0))
  if(is.vector(values)) values <- matrix(values,nrow=1)
  if(is.null(prob)) prob <- array(1,dim=dim(values))
  if(is.vector(prob)) prob <- matrix(prob,nrow=1)

  if(!is.numeric(values) || !is.numeric(prob) || !is.numeric(p)) stop("Non numeric arguments in a mathematical function")

  if(anyNA(prob) || any(prob < 0) || any(rowSums(prob)==0))
    stop("Nas in prob, or sum(prob)=0 or negative values of prob")

  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p

  ncv <- ncol(values)
  ncp <- ncol(prob)
  if(ncv != ncp) stop("prob and values should be of same length or have the same number of columns.")

  l <- length(p)
  nrv <- nrow(values)
  nrp <- nrow(prob)
  mnr <- min(max(nrv,nrp),l)
  if(nrv != 1 && nrp != 1 && nrv != nrp)
    stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

  values <- lapply(1:min(nrv,mnr),function(y) values[y,])
  prob <- lapply(1:min(nrp,mnr),function(y) prob[y,])

  prob <- mapply(function(pr,val) tapply(pr,val,sum),prob,values,SIMPLIFY=FALSE)
  prob <- lapply(prob, function(x) cumsum(x)/sum(x))
  values <- lapply(values,function(y) sort(unique(y)))

  res <- .empD_find(p, values, prob, "q", mnr)   # vectorized search
  if(is.null(res)) {
    find <- Vectorize(function(x,val,pr){
      if(x > 1 || x < 0) return(NaN)
      return(min(val[x <= pr]))},vectorize.args = "x")
  
    # to gain time and memory
    n1 <- l %/% mnr
    n2 <- l %% mnr
    l2 <- l-n2
    p1 <- p[1:l2]
  
    res <- vector("numeric",length=l)
    res[1:l2] <- mapply(find,p1,values,prob)
    if(n2!=0){
      p2 <- p[(l2+1):l]
      res[(l2+1):l] <- mapply(find,p2,values[1:n2],prob[1:n2])
      }
  
  }
  if(any(is.na(res))) warning("NaN in qempiricalD")
  return(res)
}

#' @rdname empiricalD
#' @export
rempiricalD <- function(n,values,prob=NULL)
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) cannot be negative in rempiricalD")

  if(is.vector(values)) values <- matrix(values,nrow=1)
  if(is.null(prob))     prob <- rep(1,length(values[1,]))
  if(is.vector(prob))   prob <- matrix(prob,nrow=1)

  if(!is.numeric(values) || !is.numeric(prob) || !is.numeric(n)) stop("Non numeric arguments in a mathematical function")

  ncv <- ncol(values)
  ncp <- ncol(prob)
  if(ncv!=ncp) stop("Prob and values should be of same length or have the same number of columns.")

  nrv <- nrow(values)
  nrp <- nrow(prob)
  mnr <- min(max(nrv,nrp),n)
  if(nrv!= 1 && nrp!=1 && nrv!=nrp) stop("values/prob should be vector(s), matrix(es) of 1 row or matrix(es) of the same number of rows")

  values <- lapply(1:min(nrv,mnr), function(y) values[y,])
  prob <- lapply(1:min(nrp,mnr), function(y) prob[y,])

  n1 <- n %/% mnr
  n2 <- n%% mnr
  # To gain time, send sample on n1
  # same random draws as sample(x, ...) when length(x) > 1, and correct for a single value
  samplex <- function(x, size, replace, prob) x[sample.int(length(x), size = size, replace = replace, prob = prob)]
  res <- mapply(samplex, x=values, prob=prob, MoreArgs=list(size=n1,replace=TRUE))
  dim(res) <- c(n1,mnr)
  res <- t(res)
  dim(res) <- NULL
  # then the remaining n2, if any
  if(n2!=0)
    res <- c(res,mapply(samplex,x=values[1:n2],prob=prob[1:n2],MoreArgs=list(size=1,replace=TRUE)))
  return(as.vector(res))
}

# Internal: vectorized version of the element by element search used in
# dempiricalD / pempiricalD / qempiricalD. 'values' is the list of the sorted unique values
# per row, 'prob' the list of (cumulative for "p" and "q") probabilities per row,
# as computed in these functions. Element i of x uses values[[(i-1) %% length(values) + 1]]
# and prob[[(i-1) %% length(prob) + 1]] (the recycling of mapply).
# Returns NULL when the fast path does not apply (the original code is then used).
.empD_find <- function(x, values, prob, what = c("d", "p", "q"), mnr)
{
  what <- match.arg(what)
  l <- length(x)
  nv <- length(values)
  np <- length(prob)
  if(l %% mnr != 0 && (nv != mnr || np != mnr)) return(NULL)
  if(any(lengths(prob) != lengths(values)[(seq_len(np) - 1L) %% nv + 1L])) return(NULL)
  if(nv == 1 && np == 1) {
    val <- values[[1]]
    pr <- as.vector(prob[[1]])
    res <- switch(what,
      d = { r <- pr[match(x, val)]; r[is.na(r)] <- 0; r },
      p = c(0, pr)[findInterval(x, val) + 1L],
      q = { i <- findInterval(x, pr, left.open = TRUE) + 1L
            r <- val[i]
            r[!is.na(i) & i > length(val)] <- Inf
            r })
  }
  else {
    tomat <- function(lst, fill) {
      K <- max(lengths(lst))
      m <- matrix(fill, nrow = length(lst), ncol = K)
      for(i in seq_along(lst)) m[i, seq_along(lst[[i]])] <- lst[[i]]
      m
    }
    V <- tomat(values, NA_real_)
    P <- tomat(lapply(prob, as.vector), NA_real_)
    K <- ncol(V)
    rv <- (seq_len(l) - 1L) %% nv + 1L
    rp <- (seq_len(l) - 1L) %% np + 1L
    res <- numeric(l)
    if(what == "d") {
      for(k in seq_len(K)) {
        hit <- which(x == V[rv, k])
        res[hit] <- P[cbind(rp[hit], k)]
      }
    }
    else if(what == "p") {
      idx <- integer(l)
      for(k in seq_len(K)) {
        ge <- x >= V[rv, k]
        ge[is.na(ge)] <- FALSE
        idx <- idx + ge
      }
      pos <- which(idx > 0)
      res[pos] <- P[cbind(rp[pos], idx[pos])]
    }
    else {
      idx <- rep(1L, l)
      for(k in seq_len(K)) {
        gt <- x > P[rp, k]
        gt[is.na(gt)] <- FALSE
        idx <- idx + gt
      }
      lr <- lengths(values)[rv]
      ok <- idx <= lr
      res[ok] <- V[cbind(rv[ok], idx[ok])]
      res[!ok] <- Inf
    }
  }
  if(what == "q") {
    if(any(is.infinite(res) & !is.na(x) & x >= 0 & x <= 1)) warning("no non-missing arguments to min; returning Inf")
    res[!is.na(x) & (x > 1 | x < 0)] <- NaN
  }
  res[is.na(x)] <- NA
  as.vector(res)
}
