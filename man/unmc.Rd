\name{unmc}
\alias{unmc}
\title{
Unclasses the mc or the mcnode Object
}
\description{
Unclasses the \samp{mc} object in a list of arrays or the
\samp{mcnode} object in an array.
}
\usage{
unmc(x, drop=TRUE)
}
\examples{
data(total)
## A vector
unmc(total$xV, drop=TRUE)
## An array
unmc(total$xV, drop=FALSE)
}
\arguments{
\item{x}{A \samp{mc} or a \samp{mcnode} object.}
\item{drop}{Should the dimensions of size 1 be dropped (see
\code{\link{drop}}).}
}
\value{
if x is an \samp{mc} object: a list of arrays. If \samp{drop=TRUE}, a
list of vectors, matrixes and arrays. if x is an \samp{mcnode}
object: an array. If \samp{drop=TRUE}, a vector, matrix or array.
}
\keyword{manip}
