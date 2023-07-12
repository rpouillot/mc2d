\name{quantile.mc}
\alias{quantile.mc}
\alias{quantile.mcnode}
\title{
Quantiles of a mc Object
}
\description{
Evaluates quantiles of a \samp{mc} object. This function is used by
\samp{plot.mc}
}
\usage{
\method{quantile}{mc}(x, probs=seq(0, 1, 0.01), lim=c(0.025, 0.975), na.rm=TRUE, ...)
\method{quantile}{mcnode}(x, ...)
}
\examples{
data(total)
quantile(total$xVUM3)
quantile(total)

}
\arguments{
\item{x}{a \samp{mc} objects}
\item{probs}{the quantiles to be calculated}
\item{na.rm}{TRUE or FALSE}
\item{lim}{a vector of numbers (between 0 and 1) indicating the
envelope. Maybe \samp{NULL} or empty.}
\item{\dots}{For generic method consistency.}
}
\details{
The quantiles are evaluated in the variability dimension. Then, the
median, the mean and the \samp{lim} quantiles are evaluated for each
of these quantiles.
}
\value{
A list of quantiles.
}
\seealso{
\code{\link{plot.mc}}, \code{\link{quantile}}.
}
\keyword{univar}
