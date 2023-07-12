\name{multinormal}
\alias{multinormal}
\alias{rmultinormal}
\alias{dmultinormal}
\title{
The Vectorized Multivariate Random Deviates
}
\description{
This function is the vectorized version of the \samp{rmvnorm} from
the \samp{mvtnorm} library. It provides a random number generator for
the multivariate normal distribution with varying vectors of means
and varying covariance matrixes.
}
\usage{
rmultinormal(n, mean, sigma, method=c("eigen", "svd", "chol"))
dmultinormal(x, mean, sigma, log=FALSE)
}
\examples{
## including equivalence with dmvnorm
## mean and sigma as vectors
(mean <- c(10, 0))
(sigma <- matrix(c(1, 2, 2, 10), ncol=2))
sigma <- as.vector(sigma)
(x <- matrix(c(9, 8, 1, -1), ncol=2))
round(rmultinormal(10, mean, sigma))
dmultinormal(x, mean, sigma) 
## Eq
dmvnorm(x, mean, matrix(sigma, ncol=2)) 

## mean as matrix
(mean <- matrix(c(10, 0, 0, 10), ncol=2))
round(rmultinormal(10, mean, sigma))
dmultinormal(x, mean, sigma)
## Eq
dmvnorm(x[1, ], mean[1, ], matrix(sigma, ncol=2)) 
dmvnorm(x[2, ], mean[2, ], matrix(sigma, ncol=2)) 

## sigma as matrix
(mean <- c(10, 0))
(sigma <- matrix(c(1, 2, 2, 10, 10, 2, 2, 1), nrow=2, byrow=TRUE))
round(rmultinormal(10, mean, sigma))
dmultinormal(x, mean, sigma) 
## Eq
dmvnorm(x[1, ], mean, matrix(sigma[1, ], ncol=2)) 
dmvnorm(x[2, ], mean, matrix(sigma[2, ], ncol=2)) 

## mean and sigma as matrix
(mean <- matrix(c(10, 0, 0, 10), ncol=2))
(sigma <- matrix(c(1, 2, 2, 10, 10, 2, 2, 1), nrow=2, byrow=TRUE))
round(rmultinormal(10, mean, sigma))
dmultinormal(x, mean, sigma) 
## Eq
dmvnorm(x[1, ], mean[1, ], matrix(sigma[1, ], ncol=2)) 
dmvnorm(x[2, ], mean[2, ], matrix(sigma[2, ], ncol=2)) 

(mean <- c(10, 0))
(sigma <- matrix(c(1, 2, 2, 10, 10, 2, 2, 1), nrow=2, byrow=TRUE))
x <- rmultinormal(1000, mean, sigma)
plot(x)

}
\arguments{
\item{x}{Vector or matrix of quantiles. If x is a matrix, each row is
taken to be a quantile.}
\item{n}{Number of observations. If \samp{length(n) > 1}, the length
is taken to be the number required.}
\item{mean}{Vector or matrix of means. If a matrix, each row is taken
to be a quantile. Default is a vector of 0 of convenient length.}
\item{sigma}{Covariance vector corresponding to the coercion of the
covariance matrix into a vector (if unique for all \samp{n} or
\samp{x}) or array of covariance vectors (if varying according to
\samp{n} or \samp{x}). default is a diagonal matrix of convenient
size.}
\item{method}{Matrix decomposition used to determine the matrix root
of sigma, possible methods are eigenvalue decomposition ("eigen",
default), singular value decomposition ("svd"), and Cholesky
decomposition ("chol").}
\item{log}{Logical; if \samp{TRUE}, densities d are given as log(d).}
}
\details{
\samp{rmvnorm(n, m, s)} is equivalent to \samp{rmultinormal(n, m,
as.vector(s))}. \samp{dmvnorm(x, m, s)} is equivalent to
\samp{dmultinormal(x, m, as.vector(s))}.

If \samp{mean} and/or \samp{sigma} is a matrix, the first random
deviate will use the first row of \samp{mean} and/or \samp{sigma},
the second random deviate will use the second row of \samp{mean}
and/or \samp{sigma}, ... recycling being permitted by raw. If
\samp{mean} is a vector of length \samp{l} or is a matrix with
\samp{l} columns, \samp{sigma} should be a vector of length \samp{l x
l} or a matrix of number of \samp{l x 2} columns.
}
\note{
The use of a varying sigma may be very time consuming.
}
\keyword{distribution}
