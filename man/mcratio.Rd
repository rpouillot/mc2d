\name{mcratio}
\alias{mcratio}
\title{
Ratio of uncertainty and the variability
}
\description{
Provides measures of variability, uncertainty, and both combined for
an \samp{mc} or an \samp{mcnode} object.
}
\usage{
mcratio(x, pcentral=.5, pvar=.975, punc=.975, na.rm=FALSE)
}
\examples{
data(total)
mcratio(total, na.rm=TRUE)
}
\arguments{
\item{x}{an \samp{mc} or an \samp{mcnode} object}
\item{pcentral}{the quantile for the central tendency.}
 . 

\item{pvar}{the quantile for the measure of variability.}
\item{punc}{the quantile for the measure of uncertainty.}
\item{na.rm}{a logical value indicating whether NA values should be
stripped before the computation proceeds.}
}
\details{
The function evaluates three ratios for each \samp{mcnode}. Given: 

\describe{
\item{A}{the \samp{(100 * pcentral)}th percentile of uncertainty (by
default the median) for the \samp{(100 * pcentral)}th percentile of
variability}
\item{B}{the \samp{(100 * pcentral)}th percentile of uncertainty for
the \samp{(100 * pvar)}th percentile of variability}
\item{C}{the \samp{(100 * punc)}th percentile of uncertainty for the
\samp{(100 * pcentral)}th percentile of variability}
\item{D}{the \samp{(100 * punc)}th percentile of uncertainty for the
\samp{(100 * pvar)}th percentile of variability}}
The following ratio are estimated 

\itemize{
\item Variability Ratio: B / A
\item Uncertainty Ratio: C / A
\item Overall Uncertainty Ratio: D / A}
For multivariate nodes, the statistics are evaluate on each dimension
or on statistics according to the corresponding \samp{outm} value.
}
\value{
A matrix.
}
\references{
Ozkaynak, H., Frey, H.C., Burke, J., Pinder, R.W. (2009) "Analysis of
coupled model uncertainties in source-to-dose modeling of human
exposures to ambient air pollution: A PM2.5 case study", Atmospheric
environment, Volume 43, Issue 9, March 2009, Pages 1641-1649.
}
\keyword{distribution}
