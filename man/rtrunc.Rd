\name{rtrunc}
\alias{rtrunc}
\title{
Random Truncated Distributions
}
\description{
Provides samples from classical \R distributions and \samp{mc2d}
specific distributions truncated between \samp{linf} (excluded) and
\samp{lsup} (included).
}
\usage{
rtrunc(distr=runif, n, linf=-Inf, lsup=Inf, ...)
}
\examples{
rtrunc("rnorm", n=10, linf=0)
range(rtrunc(rnorm, n=1000, linf=3, lsup=5, sd=10))
## Discrete distributions
range(rtrunc(rpois, 1000, linf=2, lsup=4, lambda=1))
##Examples of rounding problems. 
##The first one will provide a warning while the results are unexpected, 
##The second will provide an error.
\dontrun{
table(rtrunc(rbinom, n=1000, size=10, prob=1-1E-20, lsup=9))
table(rtrunc(rbinom, n=1000, size=10, prob=1E-14, linf=0))
}
}
\arguments{
\item{distr}{A function providing random data or its name as
character. The function 'rdistr' should have a 'qdistr' form (with
argument 'p') and a 'pdistr' form (with argument 'q'). Example :
'rnorm' (has a 'qnorm' and a 'pnorm' form), 'rbeta', 'rbinom',
'rgamma', ...}
\item{n}{The size of the sample.}
 .

\item{linf}{A vector of lower bounds.}
\item{lsup}{A vector of upper bounds, with \samp{lsup < linf}
(strictly).}
\item{\dots}{All arguments to be passed to \samp{pdistr} and
\samp{qdistr}.}
}
\details{
The function 1) evaluates the \samp{p} values corresponding to
\samp{linf} and \samp{lsup} using \samp{pdistr}; 2) samples \samp{n}
values using \samp{runif(n, min=pinf, max=psup)}, and 3) takes the
\samp{n} corresponding quantiles from the specified distribution
using \samp{qdistr}.

All distributions (but sample) implemented in the stats library could
be used. The arguments in \dots should be named. Do not use 'log' or
'log.p' or 'lower.tail'. For discrete dictribution, rtrunc sample
within \samp{(linf, lsup]}. See example.
}
\value{
A vector of \samp{n} values.
}
\note{
The inversion of the quantile function leads to time consuming
functions for some distributions. WARNING: The method is flexible,
but can lead to problems linked to rounding errors in some extreme
situations. The function checks that the values are in the expected
range and returns an error if not. It also warns some extreme
situation that could lead to unexpected results. See Examples.
}
\keyword{distribution}
