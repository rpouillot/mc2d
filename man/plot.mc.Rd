\name{plot.mc}
\alias{plot.mc}
\alias{plot.mcnode}
\alias{plot.plotmc}
\alias{plot.mccut}
\title{
Plots Results of a Monte Carlo Simulation
}
\description{
Plots the empirical cumulative distribution function of a
\samp{mcnode} or a \samp{mc} object ("0" and "V" nodes) or the
empirical cumulative distribution function of the estimate of a
\samp{mcnode} or \samp{mc} object ("U" and "VU" nodes).
}
\usage{
\method{plot}{mc}(x, prec=0.001, stat=c("median", "mean"), lim=c(0.025, 0.25, 0.75,
	  0.975), na.rm=TRUE, griddim=NULL, xlab=NULL, ylab="Fn(x)",
	  main="", draw=TRUE, paint=TRUE, xlim=NULL, ylim=NULL, ...)
\method{plot}{mcnode}(x, ...)
\method{plot}{plotmc}(x, ...)
\method{plot}{mccut}(x, stat=c("median", "mean"), lim=c(0.025, 0.25, 0.75, 0.975),
	  griddim=NULL, xlab=names(x), ylab="Fn(x)", main="",
	  draw=TRUE, ...)
}
\examples{
data(total)

plot(xVUM3)
## only one enveloppe corresponding to quantiles 0.025 and 0.975
plot(xVUM3, lim=c(0.025, 0.975)) 
## only one enveloppe not painted
plot(xVUM3, lim=c(0.025, 0.975), paint=FALSE) 

def.par <- par(no.readonly = TRUE)
par(mar=c(4, 4, 1, 1))
plot(total)
par(def.par)



}
\arguments{
\item{x}{a \samp{mcnode} or a \samp{mc} objects}
\item{prec}{the precision of the plot. 0.001 will provide an ecdf
from the 0.000, 0.001, .002, ..., 1.000 quantiles.}
\item{stat}{the function used for estimates (2D \samp{mc} or
\samp{mcnode}). By default the median.}
\item{lim}{a vector of numbers (between 0 and 1) indicating the
enveloppe (2D \samp{mc} or \samp{mcnode}) . Maybe \samp{NULL} or
empty.}
\item{na.rm}{Should NA values be discarded}
\item{griddim}{a vector of two integers, indicating the size of the
grid of the graph. If \samp{NULL}, the grid is calculated to produce
a "nice" graph.}
\item{xlab}{vector of labels for the x-axis. If \samp{NULL}, use the
name of the node.}
\item{ylab}{vector of labels for the y-axis.}
\item{main}{vector of main titles of the graph.}
\item{draw}{Should the plot be drawn?}
\item{paint}{Should the enveloppes be filled?}
\item{xlim}{x coordinate range. \samp{xlim} is either a vector of
length 2, used for each graph, or a list of vectors of length 2,
whose ith element is used for the ith graph. By default, the data
range is used as \samp{xlim}.}
\item{ylim}{y coordinate range. \samp{ylim} is either a vector of
length 2, used for each graph, or a list of vectors of length 2,
whose ith element is used for the ith graph. By default, the data
range is 0-1.}
\item{\dots}{further arguments to be passed to \samp{plot.stepfun}.}
}
\details{
\samp{plot.mcnode} is a user-friendly function that send the
\samp{mcnode} to \samp{plot.mc}.

 For \samp{"VU"} and \samp{"U"} \samp{mcnode}s, quantiles are
calculated using \code{\link{quantile.mc}} within each of the
\samp{nsu} simulations (i.e. by columns of each \samp{mcnode}). The
medians (but may be the means using \samp{stat="mean"}) calculated
from the \samp{nsu} values are plotted. The 0.025 and 0.975
quantiles, and the 0.25 and 0.75 quantiles (default values of
\samp{lim}) of these quantiles are used as the enveloppe.
}
\value{
A \samp{plot.mc} object, list of the quantiles used to plot the draw.
}
\references{
Cullen AC and Frey HC (1999) Probabilistic techniques in exposure
assessment. Plenum Press, USA, pp. 81-155.
}
\seealso{
\code{\link{ecdf}}, \code{\link{plot}}, \code{\link{quantile.mc}}
}
\keyword{hplot}
