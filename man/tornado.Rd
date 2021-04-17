\name{tornado}
\alias{tornado}
\alias{print.tornado}
\title{
Computes Correlation between Inputs and Output in a mc Object
(tornado) in the Variability Dimension;
}
\description{
Provides statistics for a tornado chart. Evaluates correlations
between output and inputs of a \samp{mc} object.
}
\usage{
tornado(mc, output=length(mc), use="all.obs", method=c("spearman",
	  "kendall", "pearson"), lim=c(0.025, 0.975))
\method{print}{tornado}(x, ...)
}
\examples{
data(total)
tornado(total, 2, "complete.obs", "spearman", c(0.025, 0.975))
tornado(total, 4, "pairwise.complete.obs", "spearman", c(0.025, 0.975))
tornado(total, 6, "complete.obs", "kendall", c(0.025, 0.975))
tornado(total, 8, "complete.obs", "spearman", c(0.025, 0.975))
(y <- tornado(total, 10, "complete.obs", "spearman", c(0.025, 0.975)))
plot(y)

}
\arguments{
\item{mc}{a \code{\link{mc}} object or a \code{\link{mccut}} object.}
\item{x}{A \samp{tornado} object as provided by the \samp{tornado}
function.}
\item{output}{(for \samp{mc} objects only). The rank or the name of
the output to be considered. By default: the last element of the
\samp{mc}.}
\item{use}{(for \samp{mc} objects only). An optional character string
giving a method for computing covariances in the presence of missing
values. This must be (an abbreviation of) one of the strings
"all.obs", "complete.obs" or "pairwise.complete.obs" (see
\code{\link{cor}}).}
\item{method}{(for \samp{mc} objects only). A character string
indicating which correlation coefficient (or covariance) is to be
computed. One of "spearman" (default), "kendall" or "pearson", can be
abbreviated (see \code{\link{cor}}). Warning : the default is not the
same in \code{\link{cor}}.}
\item{lim}{A vector of quantiles used to compute the credible
interval in two-dimensional models.}
\item{\dots}{Further arguments to be passed to the final print
function.}
}
\details{
The tornado function computes the spearman's rho statistic. It is
used to estimate a rank-based measure of association between one set
of random variable of a \samp{mc} object (the output) and the others
(the inputs).

 \samp{tornado} may be applied on a \samp{mccut} object if a
\samp{tornado} function was used in the third block of the
\code{\link{evalmccut}} call.

If "output" refers to a \samp{"0" mcnode}, it is an error. If
"output" refers to a \samp{"V" mcnode}, correlations are only
provided for other \samp{"V" mcnode}s. If "output" refers to a
\samp{"U" mcnode}, correlations are only provided for other \samp{"U"
mcnode}s. If "output" refers to a \samp{"VU" mcnode}, correlations
are only provided for other \samp{"VU" mcnode}s and \samp{"V"
mcnode}s.

If use is "all.obs", then the presence of missing observations will
produce an error. If use is "complete.obs" then missing values are
handled by casewise deletion. Finally, if use has the value
"pairwise.complete.obs" then the correlation between each pair of
variables is computed using all complete pairs of observations on
those variables.
}
\value{
An invisible object of class tornado. A tornado object is a list of
objects containing the following objects: 

\item{value}{the value of correlation coefficients}
\item{output}{the name of the output}
\item{method}{the method used}
\item{use}{the use parameter}
}
\seealso{
\code{\link{cor}}.

 \code{\link{plot.tornado}} to draw the results.
}
\keyword{univar}
