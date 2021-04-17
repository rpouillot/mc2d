\name{empiricalC}
\alias{empiricalC}
\alias{dempiricalC}
\alias{pempiricalC}
\alias{qempiricalC}
\alias{rempiricalC}
\title{
The Continuous Empirical Distribution
}
\description{
Density, distribution function and random generation for a continuous
empirical distribution.
}
\usage{
dempiricalC(x, min, max, values, prob=NULL, log=FALSE)
pempiricalC(q, min, max, values, prob=NULL, lower.tail=TRUE, log.p=FALSE)
qempiricalC(p, min, max, values, prob=NULL, lower.tail=TRUE, log.p=FALSE)
rempiricalC(n, min, max, values, prob=NULL)
}
\examples{
prob <- c(2, 3, 1, 6, 1)
values <- 1:5
par(mfrow=c(1, 2))
curve(dempiricalC(x, min=0, max=6, values, prob), from=-1, to=7, n=1001)
curve(pempiricalC(x, min=0, max=6, values, prob), from=-1, to=7, n=1001)

## Varying values
(values <- matrix(1:10, ncol=5))
## the first x apply to the first row 
## the second x to the second one
dempiricalC(c(1, 1), values, min=0, max=11)


##Use with mc2d 
val <- c(100, 150, 170, 200)
pr <- c(6, 12, 6, 6)
out <- c("min", "mean", "max")
##First Bootstrap in the uncertainty dimension
##with rempirical D
(x <- mcstoc(rempiricalD, type = "U", outm = out, nvariates = 30, values = val, prob = pr))
##Continuous Empirical distribution in the variability dimension
mcstoc(rempiricalC, type = "VU", values = x, min=90, max=210)



}
\arguments{
\item{x, q}{Vector of quantiles.}
\item{p}{Vector of probabilities.}
\item{n}{Number of random values. If \samp{length(n) > 1}, the length
is taken to be the number required.}
\item{min}{A finite minimal value.}
\item{max}{A finite maximal value.}
\item{values}{Vector of numerical values.}
\item{prob}{Optionnal vector of count or probabilities.}
\item{log, log.p}{logical; if \samp{TRUE}, probabilities \samp{p} are
given as \samp{log(p)}.}
\item{lower.tail}{logical; if \samp{TRUE} (default), probabilities
are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.}
}
\details{
Given \eqn{p_{i}}{p_i}, the distribution value for \eqn{x_{i}}{x_i}
with \samp{i} the rank \eqn{i = 0, 1, 2, \ldots, N+1},
\eqn{x_{0}=min}{x_0 = min} and \eqn{x_{N+1}=max}{x_(N+1) = max} the
density is:
\deqn{f(x)=p_{i}+(\frac{x-x_{i}}{x_{i+1}-x_{i}})(p_{i+1}-p_{i})}{f(x)
= p_i + (p_(i+1) - p_i)/(x_(i+1) - x_i) for x_i<=x<x_(i+1)} The
\samp{p} values being normalized to give the distribution a unit
area.

\samp{min} and/or \samp{max} and/or \samp{values} and/or \samp{prob}
may vary: in that case, \samp{min} and/or \samp{max} should be
vector(s). \samp{values} and/or \samp{prob} should be matrixes, the
first row being used for the first element of \samp{x}, \samp{q},
\samp{p} or the first random value, the second row for the second
element of \samp{x}, \samp{q}, \samp{p} or random value, ...
Recycling is permitted if the number of elements of \samp{min} or
\samp{max} or the number of rows of \samp{prob} and \samp{values} are
equal or equals one.
}
\value{
\samp{dempiricalC} gives the density, \samp{pempiricalC} gives the
distribution function, \samp{qempiricalC} gives the quantile function
and \samp{rempiricalC} generates random deviates.
}
\seealso{
\code{\link{empiricalD}}
}
\keyword{distribution}
