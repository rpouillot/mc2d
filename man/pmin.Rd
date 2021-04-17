\name{pmin}
\alias{pmin}
\alias{pmax}
\alias{pmin.default}
\alias{pmax.default}
\alias{pmin.mcnode}
\alias{pmax.mcnode}
\title{
Maxima and Minima for mcnodes
}
\description{
Returns the parallel maxima and minima of the input values.
}
\usage{
\method{pmin}{mcnode}(..., na.rm=FALSE)
\method{pmax}{mcnode}(..., na.rm=FALSE)
}
\examples{
ndvar(10);ndunc(21)
x <- mcstoc(rnorm, "V")
pmin(x, 0)
y <- mcdata(rep(c(-1, 1), length=ndunc()), "U")
unclass(pmin(x, y))





}
\arguments{
\item{\dots}{One or more \samp{mcnodes}s or one or more
\samp{mcnode}s and vector(s) of compatible size. Note that one
\samp{mcnode} must be at the first place.}
\item{na.rm}{a logical indicating whether missing values should be
removed.}
}
\details{
\samp{pmax} and \samp{pmin} take one or more \samp{mcnode} and/or
vectors as arguments and return a \samp{mcnode} of adequate type and
size giving the "parallel" maxima (or minima) of the \samp{mcnode}
and/or vectors. Note that the first element of ... should be an
\samp{mcnode}. The resulting type of \samp{mcnode} is variable
according to the elements that are passed. The same rules as in
\code{\link{Ops.mcnode}} are applied.
}
\value{
an \samp{mcnode} of adequate type and dimension.
}
\seealso{
\code{\link{min}}, \code{\link{Ops.mcnode}}
}
\keyword{utilities}
