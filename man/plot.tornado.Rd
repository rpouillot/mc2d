\name{plot.tornado}
\alias{plot.tornado}
\alias{plot.tornadounc}
\title{
Draws a Tornado chart.
}
\description{
Draws a Tornado chart as provided by \samp{tornado}.
}
\usage{
\method{plot}{tornado}(x, which=1, name=NULL, stat=c("median", "mean"), xlab="method",
	  ylab="", ...)
\method{plot}{tornadounc}(x, which=1, stat="median", name=NULL, xlab="method", ylab="", ...)
}
\examples{
data(ec)
x <- evalmcmod(ec$modEC2, nsv=100, nsu=100, seed=666)
tor <- tornado(x, 7)
plot(tor)

}
\arguments{
\item{x}{A \code{\link{tornado}} object or a \code{\link{tornadounc}}
object.}
\item{which}{Which output to print -for multivariates output-.}
\item{name}{Vector of name of input variables. If NULL, the name will
be given from the name of the elements.}
\item{stat}{The name of the statistics of the output to be
considered. For a \samp{tornado} object: "median" or "mean". For a
\samp{tornadounc} object: the value should match one row name of the
\samp{tornadounc} object. Alternatively, for a \samp{tornadounc}
object, the number of the row may be used.}
\item{xlab}{Label of the x axis. if "method", use the correlation
method used in the \samp{tornado} object.}
\item{ylab}{Label of the y axis.}
\item{\dots}{Further arguments to be passed to the \samp{plot}
function.}
}
\details{
A point is drawn at the estimate and the segment reflects the
uncertainty around this estimate.
}
\value{
NULL
}
\seealso{
\code{\link{tornado}}
}
\keyword{hplot}
