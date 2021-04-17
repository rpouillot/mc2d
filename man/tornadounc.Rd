\name{tornadounc}
\alias{tornadounc}
\alias{tornadounc.mc}
\alias{tornadounc}
\alias{tornadounc.default}
\alias{print.tornadounc}
\alias{tornadounc.mccut}
\title{
Computes Correlation between Inputs and Output in a mc Object
(tornado) in the Uncertainty Dimension
}
\description{
Provides statistics for a tornado chart. Evaluates correlations
between output and inputs of a \samp{mc} object in the uncertainty
dimension.
}
\usage{
\method{tornadounc}{mc}(mc, output=length(mc), quant=c(0.5, 0.75, 0.975), use="all.obs",
	  method=c("spearman", "kendall", "pearson"), ...)
\method{tornadounc}{default}(mc, ...)
\method{print}{tornadounc}(x, ...)
\method{tornadounc}{mccut}(mc, output=length(mc), quant=c(0.5, 0.75, 0.975), use="all.obs",
	  method=c("spearman", "kendall", "pearson"), ...)
}
\examples{
data(total)
tornadounc(total, 3)
tornadounc(total, 4, use="complete")
tornadounc(total, 7, use="complete.obs")
tornadounc(total, 8, use="complete.obs")
(y <- tornadounc(total, 10, use="complete.obs"))
plot(y, 1, 1)




}
\arguments{
\item{mc}{a \samp{mc} object.}
\item{x}{a \samp{tornadounc} object.}
\item{output}{The rank or the name of the output to be considered.
Should be a \samp{"VU"} or a \samp{"U" type mcnode}. By default: the
last element of \samp{mc}.}
\item{quant}{The vector of quantiles used in the variability
dimension.}
\item{use}{An optional character string giving a method for computing
covariances in the presence of missing values. This must be (an
abbreviation of) one of the strings "all.obs", "complete.obs" or
"pairwise.complete.obs" (see \code{\link{cor}}).}
\item{method}{A character string indicating which correlation
coefficient (or covariance) is to be computed. One of "spearman"
(default), "kendall" or "pearson", can be abbreviated (see
\code{\link{cor}}). Warning : "pearson" is the default for
\code{\link{cor}}).}
\item{\dots}{Further arguments to be passed to the final print
function.}
}
\details{
The \samp{tornadounc.mc} function computes the spearman's rho
statistic between 

\itemize{
\item values (\samp{"U" type mcnode}) or statistics calculated in the
variability dimension (\samp{"VU" type mcnode}) of inputs and
\item values (\samp{"U" type mcnode}) or statistics calculated in the
variability dimension (\samp{"VU" type mcnode}) of one output.}
 The statistics are the mean, the median and the quantiles specified
by \samp{quant}.

It is useful to estimate a rank-based measure of association between
one set of random variable of a \samp{mc} object (the output) and the
others in the uncertainty dimension.

 \samp{tornadounc.mccut} may be applied on a \code{\link{mccut}}
object if a \samp{summary.mc} function was used in the third block of
the \code{\link{evalmccut}} call.

If output refers to a \samp{"0"} or \samp{"V" mcnode}, it is an
error.

If use is "all.obs", then the presence of missing observations will
produce an error. If use is "complete.obs" then missing values are
handled by casewise deletion. Finally, if use has the value
"pairwise.complete.obs" then the correlation between each pair of
variables is computed using all complete pairs of observations on
those variables.
}
\value{
An invisible object of class \samp{tornadounc}. A \samp{tornadounc}
object is a list of objects containing the following objects: 

\item{value}{a matrix of values of correlation coefficients. Each row
are the value or the statistics of inputs, each columns the value or
the statistics of outputs.}
\item{output}{the name of the output}
\item{method}{the method used}
\item{use}{the \samp{use} parameter}
}
\seealso{
\code{\link{cor}}.

 \code{\link{tornado}} for tornado in the variability dimension.

 \code{\link{plot.tornadounc}} to draw the results.
}
\keyword{univar}
