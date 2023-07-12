\name{lhs}
\alias{lhs}
\title{
Random Latin Hypercube Sampling
}
\description{
Creates a Latin Hypercube Sample (LHS) of the specified distribution.
}
\usage{
lhs(distr="runif", nsv=ndvar(), nsu=ndunc(), nvariates=1, ...)
}
\examples{
ceiling(lhs(runif, nsu=10, nsv=10)*10)
}
\arguments{
\item{distr}{The function for generating random sample or its name.
If \samp{distr} is "rdist", the function "qdist" must be the quantile
function of this distribution with argument \samp{p} as a vector of
probabilities, as all univariates distributions of the \samp{stat}
library.}
\item{nsv}{The number of rows of the final matrix.}
\item{nsu}{The number of columns of the final matrix}
\item{nvariates}{The number of variates}
\item{\dots}{All arguments to be passed to \samp{distr} except the
size of the sample.}
}
\value{
A \samp{nsv x nsu} matrix of random variates.
}
\author{
adapted from a code of Rob Carnell (library \samp{lhs})
}
\seealso{
\code{\link{mcstoc}}
}
\note{
The resulting lhs is in fact a latin hypersquare sampling: the lhs is
provided only in the first 2 dimensions.

 It is not possible to send truncated distribution with
\code{\link{rtrunc}}. Use \code{\link{mcstoc}} for this purpose, with
\samp{lhs=TRUE} and \samp{rtrunc=TRUE}.

 The \dots arguments will be recycled.
}
\keyword{design}
