\name{NA.mcnode}
\alias{NA.mcnode}
\alias{is.na.mcnode}
\alias{is.nan.mcnode}
\alias{is.finite.mcnode}
\alias{is.infinite.mcnode}
\title{
Finite, Infinite, NA and NaN Numbers in mcnode.
}
\description{
\samp{is.na}, \samp{is.nan}, \samp{is.finite} and \samp{is.infinite}
return a logical \samp{mcnode} of the same dimension as \samp{x}.
}
\usage{
\method{is.na}{mcnode}(x)
\method{is.nan}{mcnode}(x)
\method{is.finite}{mcnode}(x)
\method{is.infinite}{mcnode}(x)
}
\examples{
x <- log(mcstoc(rnorm, nsv=1001))
x
is.na(x)



}
\arguments{
\item{x}{A \samp{mcnode} object.}
}
\value{
A logical \samp{mcnode} object.
}
\seealso{
\code{\link{is.finite}}, \code{\link{NA}}
}
\keyword{NA}
