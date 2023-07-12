\name{summary.mc}
\alias{summary.mc}
\alias{summary.mcnode}
\alias{print.summary.mc}
\alias{summary.mccut}
\title{
Summary of mcnode and mc Object
}
\description{
Provides a summary of a \samp{mcnode}, a \samp{mc} or a \samp{mccut}
object.
}
\usage{
\method{summary}{mc}(object, probs=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), lim=c(0.025,
	  0.975), ...)
\method{summary}{mcnode}(object, probs=c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1), lim=c(0.025,
	  0.975), digits=3, ...)
\method{print.summary}{mc}(x, digits=3, ...)
\method{summary}{mccut}(object, lim=c(0.025, 0.975), ...)
}
\examples{
data(total)
summary(xVUM3)
summary(total)



}
\arguments{
\item{object}{a \samp{mcnode} or a \samp{mc} object or a \samp{mccut}
object.}
\item{x}{A \samp{summary.mc} object as provided by the
\samp{summary.mc} function.}
\item{probs}{A vector of values used for the quantile function
(variability dimension).}
\item{digits}{Number of digits in the print.}
\item{lim}{A vector of values used for the quantile function
(uncertainty dimension).}
\item{\dots}{For generic functions consistency.}
}
\details{
The mean, the standard deviation and the \samp{probs} quantiles will
be evaluated in the variability dimension. The median, the mean and
the \samp{lim} quantiles will then be evaluated on these statistics
in the uncertainty dimension.

Multivariate nodes:

If the \samp{"outm"} attributes of the mcnode is "none", the node is
not evaluated, if it is "each" the variates are evaluated one by one,
if it is a function (e.g. "mean"), the function is applied on the
\samp{nvariates} dimension before providing a classical output.
}
\value{
a list.
}
\seealso{
\code{\link{mcnode}} for mcnode objects, \code{\link{mc}} for mc
objects, \code{\link{mccut}} for mccut objects,
\code{\link{quantile}}
}
\keyword{univar}
