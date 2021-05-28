#' The Log Normal Distribution parameterized through its mean and standard deviation.
#' 
#' Density, distribution function, quantile function and random generation for a log normal distribution whose 
#' arithmetic mean equals to `mean` and standard deviation equals to `sd`.
#' 
#' This function calls the corresponding density, distribution function, quantile function and random generation 
#' from the log normal (see \code{\link[stats]{Lognormal}}) after evaluation of \eqn{meanlog = log(mean^2 / sqrt(sd^2+mean^2))} and
#' \eqn{sqrt{(log(1+sd^2/mean^2))}}  
#' 
#' @name Lognormalb
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param mean the mean of the distribution.
#' @param sd the standard deviation of the distribution.
#' @param log,log.p logical. if `TRUE` probabilities `p` are given as `log(p)`.
#' @param lower.tail logical. if `TRUE`, probabilities are \eqn{P[X \le x]}, otherwise,  \eqn{P[X > x]}.
#' @keywords distribution
#' @seealso code{\link[stats]{Lognormal}}
#' 
#' @return `dlnormb`  gives the density, `plnormb` gives the distribution function, 
#' `qlnormb` gives the quantile function, and `rlnormb` generates random deviates.
#' The length of the result is determined by `n` for `rlnorm`, and is the maximum of the lengths 
#' of the numerical arguments for the other functions. 
#' The numerical arguments other than `n` are recycled to the length of the result. 
#' Only the first elements of the logical arguments are used.
#' The default `mean` and `sd` are chosen to provide a distribution close to a lognormal with 
#' `meanlog = 0` and `sdlog = 1`.
#' 
#' @exampleS 
#' x <- rlnormb(1E5,3,6)
#' mean(x) 
#' sd(x)
#' dlnormb(1) == dnorm(0)
#' dlnormb(1) == dlnorm(1)
#' 
dlnormb <- function (x, mean = exp(0.5), sd = sqrt(exp(2)-exp(1)), log = FALSE)
{ 
  gmean <- log(mean^2 / sqrt(sd^2+mean^2))
  gstd <- sqrt(log(1+sd^2/mean^2))
  
  return(dlnorm(x, meanlog = gmean, sdlog = gstd, log = log))
}

#' @rdname Lognormalb
plnormb <- function (q, mean = exp(0.5), sd = sqrt(exp(2)-exp(1)), lower.tail = TRUE, log.p = FALSE)
{ 
  gmean <- log(mean^2 / sqrt(sd^2+mean^2))
  gstd <- sqrt(log(1+sd^2/mean^2))
  
  return(plnorm(q, meanlog = gmean, sdlog = gstd, lower.tail = lower.tail, log.p = log.p))
}

#' @rdname Lognormalb
qlnormb <- function (p, mean = exp(0.5), sd = sqrt(exp(2)-exp(1)), lower.tail = TRUE, log.p = FALSE)
{ 
  gmean <- log(mean^2 / sqrt(sd^2+mean^2))
  gstd <- sqrt(log(1+sd^2/mean^2))
  
  return(qlnorm(p, meanlog = gmean, sdlog = gstd, lower.tail = lower.tail, log.p = log.p))
}

#' @rdname Lognormalb
rlnormb <- function (n, mean = exp(0.5), sd = sqrt(exp(2)-exp(1)))
{ 
  gmean <- log(mean^2 / sqrt(sd^2+mean^2))
  gstd <- sqrt(log(1+sd^2/mean^2))
  
  return(rlnorm(n, meanlog = gmean, sdlog = gstd))
}

