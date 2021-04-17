\name{empiricalD}
\alias{empiricalD}
\alias{dempiricalD}
\alias{pempiricalD}
\alias{qempiricalD}
\alias{rempiricalD}
\title{
The Discrete Empirical Distribution
}
\description{
Density, distribution function and random generation for a discrete
empirical distribution. This function is vectorized to accept
different sets of \samp{values} or \samp{prob}.
}
\usage{
dempiricalD(x, values, prob=NULL, log=FALSE)
pempiricalD(q, values, prob=NULL, lower.tail=TRUE, log.p=FALSE)
qempiricalD(p, values, prob=NULL, lower.tail=TRUE, log.p=FALSE)
rempiricalD(n, values, prob=NULL)
}
\examples{
dempiricalD(1:6, 2:6, prob=c(10, 10, 70, 0, 10))
pempiricalD(1:6, 2:6, prob=c(10, 10, 70, 0, 10))
qempiricalD(seq(0, 1, 0.1), 2:6, prob=c(10, 10, 70, 0, 10))
table(rempiricalD(10000, 2:6, prob=c(10, 10, 70, 0, 10)))

## Varying values
(values <- matrix(1:10, ncol=5))
## the first x apply to the first row : p = 0.2
## the second x to the second one: p = 0
dempiricalD(c(1, 1), values)


##Use with mc2d
##Non Parameteric Bootstrap
val <- c(100, 150, 170, 200)
pr <- c(6, 12, 6, 6)
out <- c("min", "mean", "max")
##First Bootstrap in the uncertainty dimension
(x <- mcstoc(rempiricalD, type = "U", outm = out, nvariates = 30, values = val, prob = pr))
##Second one in the variability dimension
mcstoc(rempiricalD, type = "VU", values = x)



}
\arguments{
\item{x, q}{Vector of quantiles.}
\item{p}{Vector of probabilities.}
\item{n}{Number of random values. If length(n) \eqn{>} 1, the length
is taken to be the number required.}
\item{values}{Vector or matrix of numerical values. See details.}
\item{prob}{Optionnal vector or matrix of count or probabilities. See
details.}
\item{log, log.p}{logical; if \samp{TRUE}, probabilities \samp{p} are
given as \samp{log(p)}.}
\item{lower.tail}{logical; if \samp{TRUE} (default), probabilities
are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.}
}
\details{
If \samp{prob} is missing, the discrete distribution is obtained
directly from the vector of \samp{values}, otherwise \samp{prob} is
used to weight the values. \samp{prob} is normalized before use.
Thus, \samp{prob} may be the count of each \samp{values}. \samp{prob}
values should be non negative and their sum should not be 0.

 \samp{values} and/or \samp{prob} may vary: in that case,
\samp{values} and/or \samp{prob} should be sent as matrixes, the
first row being used for the first element of \samp{x}, \samp{q},
\samp{p} or the first random value, the second row for the second
element of \samp{x}, \samp{q}, \samp{p} or random value, ...
Recycling is permitted if the number of rows of \samp{prob} and
\samp{values} are equal or if the number of rows of \samp{prob}
and/or \samp{values} are one.

 \samp{rempiricalD(n, values, prob)} with \samp{values} and
\samp{prob} as vectors is equivalent to \samp{sample(x=values,
size=n, replace=TRUE, prob=prob)}.
}
\value{
\samp{dempiricalD} gives the density, \samp{pempiricalD} gives the
distribution function, \samp{qempiricalD} gives the quantile function
and \samp{rempiricalD} generates random deviates.
}
\seealso{
\code{\link{sample}}. \code{\link{empiricalC}}.
}
\note{
In the future, the fonctions should be written for non numerical
values.
}
\keyword{distribution}
