\name{dimmcnode}
\alias{dimmcnode}
\alias{dimmc}
\title{
Dimension of mcnode and mc Objects
}
\description{
Provides the dimension (i.e. the number of simulations in the
variability dimension, the number of simulations in the uncertainty
dimension and the maximum number of variates of a \samp{mcnode} or a
\samp{mc} object.
}
\usage{
dimmcnode(x)
dimmc(x)
}
\examples{
data(total)
dimmcnode(xVUM2)
dimmc(total)

}
\arguments{
\item{x}{a \samp{mcnode} or a \samp{mc} object.}
}
\value{
A vector of three scalars: the dimension of variability (1 for
\samp{"0"} and \samp{"U" mcnode}), the dimension of uncertainty (1
for \samp{"0"} and \samp{"V" mcnode}) and the number of variates (the
maximal number of variates for an \samp{mc} object.
}
\note{
This function does not test if the object is correctly built. See
\code{\link{is.mcnode}} and \code{\link{is.mc}} .
}
\keyword{utilities}
