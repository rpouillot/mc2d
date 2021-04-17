\name{hist.mc}
\alias{hist.mc}
\alias{hist.mcnode}
\title{
Histogram of a Monte Carlo Simulation
}
\description{
Shows histogram of a \samp{mcnode} or a \samp{mc} object.
}
\usage{
\method{hist}{mc}(x, griddim=NULL, xlab=names(x), ylab="Frequency", main="", ...)
\method{hist}{mcnode}(x, ...)
}
\examples{
data(total)
hist(xVUM3)
hist(total)

}
\arguments{
\item{x}{An \samp{mcnode} or an \samp{mc} object.}
\item{griddim}{A vector of two integers, indicating the size of the
grid of plots. If \samp{NULL}, the grid is calculated to produce a
"nice" graph.}
\item{xlab}{A vector of labels for the x-axis for drawn histograms
(those whose \samp{outm(x)!="none"}). May be recycled.}
\item{ylab}{A vector of labels for the y-axis for drawn histograms.
May be recycled.}
\item{main}{A vector of main title of histograms for drawn
histograms. May be recycled.}
\item{\dots}{Other arguments to be passed to all calls of
\samp{hist}.}
}
\note{
For Two-dimensional \samp{mc}, the histogram is based on all data
(variability and uncertainty) pooled together.
}
\keyword{hplot}
