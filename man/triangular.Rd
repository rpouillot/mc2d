\name{triangular}
\alias{triangular}
\alias{dtriang}
\alias{ptriang}
\alias{qtriang}
\alias{rtriang}
\title{
The Triangular Distribution
}
\description{
Density, distribution function, quantile function and random
generation for the triangular distribution with minimum equal to
\samp{min}, mode equal \samp{mode} and maximum equal to \samp{max}.
}
\usage{
dtriang(x, min=-1, mode=0, max=1, log=FALSE)
ptriang(q, min=-1, mode=0, max=1, lower.tail=TRUE, log.p=FALSE)
qtriang(p, min=-1, mode=0, max=1, lower.tail=TRUE, log.p=FALSE)
rtriang(n, min=-1, mode=0, max=1)
}
\examples{
curve(dtriang(x, min=3, mode=5, max=10), from = 2, to = 11)
##no density when min == mode == max
dtriang(c(1, 2, 3), min=2, mode=2, max=2)



}
\arguments{
\item{x,q}{vector of quantiles.}
\item{p}{vector of probabilities.}
\item{n}{number of observations. If length(n) \eqn{>} 1, the length
is taken to be the number required.}
\item{min}{vector of minima.}
\item{mode}{vector of modes.}
\item{max}{vector of maxima.}
\item{log, log.p}{logical; if \samp{TRUE}, probabilities \samp{p} are
given as \samp{log(p)}.}
\item{lower.tail}{logical; if \samp{TRUE} (default), probabilities
are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.}
}
\details{
For the case of u := min == mode == max, there is no density in that
case and dtriang will return NaN (the error condition) (Similarity
with dunif).
}
\value{
\samp{dtriang} gives the density, \samp{ptriang} gives the
distribution function, \samp{qtriang} gives the quantile function,
and \samp{rtriang} generates random deviates.
}
\keyword{distribution}
