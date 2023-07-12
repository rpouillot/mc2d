\name{betagen}
\alias{betagen}
\alias{dbetagen}
\alias{pbetagen}
\alias{qbetagen}
\alias{rbetagen}
\title{
The Generalised Beta Distribution
}
\description{
Density, distribution function, quantile function and random
generation for the Beta distribution defined on the \samp{[min, max]}
domain with parameters \samp{shape1} and \samp{shape2} ( and optional
non-centrality parameter \samp{ncp}).
}
\usage{
dbetagen(x, shape1, shape2, min=0, max=1, ncp=0, log=FALSE)
pbetagen(q, shape1, shape2, min=0, max=1, ncp=0, lower.tail=TRUE,
	  log.p=FALSE)
qbetagen(p, shape1, shape2, min=0, max=1, ncp=0, lower.tail=TRUE,
	  log.p=FALSE)
rbetagen(n, shape1, shape2, min=0, max=1, ncp=0)
}
\examples{
curve(dbetagen(x, shape1=3, shape2=5, min=1, max=6), from = 0, to = 7)
curve(dbetagen(x, shape1=1, shape2=1, min=2, max=5), from = 0, to = 7, lty=2, add=TRUE)
curve(dbetagen(x, shape1=.5, shape2=.5, min=0, max=7), from = 0, to = 7, lty=3, add=TRUE)



}
\arguments{
\item{x,q}{Vector of quantiles.}
\item{p}{Vector of probabilities.}
\item{n}{Number of observations. If \samp{length(n) > 1}, the length
is taken to be the number required.}
\item{shape1, shape2}{Positive parameters of the Beta distribution.}
\item{min}{Vector of minima.}
\item{max}{Vector of maxima.}
\item{ncp}{Non-centrality parameter of the Beta distribution.}
\item{log, log.p}{Logical; if \samp{TRUE}, probabilities \samp{p} are
given as \samp{log(p)}.}
\item{lower.tail}{Logical; if \samp{TRUE} (default), probabilities
are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.}
}
\details{
\deqn{x \sim betagen(shape1, shape2, min, max, ncp)}{x ~
betagen(shape1, shape2, min, max, ncp)} if
\deqn{\frac{x-min}{max-min}\sim
beta(shape1,shape2,ncp)}{(x-min)/(max-min)~beta(shape1,shape2,ncp)}
These functions use the \code{\link{Beta}} distribution functions
after correct parameterization.
}
\value{
\samp{dbetagen} gives the density, \samp{pbetagen} gives the
distribution function, \samp{qbetagen} gives the quantile function,
and \samp{rbetagen} generates random deviates.
}
\seealso{
\code{\link{Beta}}
}
\keyword{distribution}
