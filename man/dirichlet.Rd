\name{dirichlet}
\alias{dirichlet}
\alias{ddirichlet}
\alias{rdirichlet}
\title{
The Dirichlet Distribution
}
\description{
Density function and random generation from the Dirichlet
distribution.
}
\usage{
ddirichlet(x, alpha)
rdirichlet(n, alpha)
}
\examples{
dat <- c(1, 10, 100, 1000, 1000, 100, 10, 1)
(alpha <- matrix(dat, nrow=4, byrow=TRUE))
round(x <- rdirichlet(4, alpha), 2)
ddirichlet(x, alpha)

## rdirichlet used with mcstoc
mcalpha <- mcdata(dat, type="V", nsv=4, nvariates=2)
(x <- mcstoc(rdirichlet, type="V", alpha=mcalpha, nsv=4, nvariates=2))
unclass(x)
x <- mcstoc(rdirichlet, type="VU", alpha=mcalpha, nsv=4, nsu=10, nvariates=2)
unclass(x)

}
\arguments{
\item{x}{A vector containing a single deviate or a matrix containing
one random deviate per row.}
\item{alpha}{A vector of shape parameters, or a matrix of shape
parameters by rows. Recycling (by row) is permitted.}
\item{n}{Number of random vectors to generate. If length(n) \eqn{>}
1, the length is taken to be the number required.}
}
\details{
The Dirichlet distribution is the multidimensional generalization of
the beta distribution. The original code was adapted to provide a
kind of "vectorization" used in multivariates \samp{mcnode}.
}
\value{
\samp{ddirichlet} gives the density. \samp{rdirichlet} returns a
matrix with \samp{n} rows, each containing a single Dirichlet random
deviate.
}
\author{
Code is adapted from \samp{MCMCpack}. It originates from Greg's
Miscellaneous Functions (gregmisc).
}
\seealso{
\code{\link{Beta}}
}
\keyword{distribution}
