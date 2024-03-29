\name{dmultinomial}
\alias{dmultinomial}
\alias{rmultinomial}
\title{
The Vectorized Multinomial Distribution
}
\description{
Generate multinomially distributed random number vectors and compute
multinomial probabilities.
}
\usage{
dmultinomial(x, size=NULL, prob, log=FALSE)
rmultinomial(n, size, prob)
}
\examples{
x <- c(100, 200, 700)
x1 <- matrix(c(100, 200, 700, 200, 100, 700, 700, 200, 100), byrow=TRUE, ncol=3)
p <- c(1, 2, 7)
p1 <- matrix(c(1, 2, 7, 2, 1, 7, 7, 2, 1), byrow=TRUE, ncol=3)
dmultinomial(x1, prob=p) 
## is equivalent to 
c(	dmultinom(x1[1, ], prob=p), 
	dmultinom(x1[2, ], prob=p), 
	dmultinom(x1[3, ], prob=p))

dmultinomial(x1, prob=p1, log=TRUE) 
## is equivalent to 
c(	dmultinom(x1[1, ], prob=p1[1, ], log=TRUE), 
	dmultinom(x1[2, ], prob=p1[2, ], log=TRUE), 
	dmultinom(x1[3, ], prob=p1[3, ], log=TRUE))

dmultinomial(x, prob=p1, log=TRUE)
## is equivalent to 
c(	dmultinom(x, prob=p1[1, ], log=TRUE), 
	dmultinom(x, prob=p1[2, ], log=TRUE), 
	dmultinom(x, prob=p1[3, ], log=TRUE))

prob <- c(1, 2, 7)
rmultinomial(4, 1000, prob)
rmultinomial(4, c(10, 100, 1000, 10000), prob)

## rmultinomial used with mcstoc
## (uncertain size and prob)
s <- mcstoc(rpois, "U", lambda=50)
p <- mcstoc(rdirichlet, "U", nvariates=3, alpha=c(4, 10, 20))
mcstoc(rmultinomial, "VU", nvariates=3, size=s, prob=p)

}
\arguments{
\item{x}{vector or matrix of length (or ncol) K of integers in
\samp{0:size}.}
\item{n}{number of random vectors to draw.}
\item{size}{a vector of integers, say N, specifying the total number
of objects that are put into K boxes in the typical multinomial
experiment. For \samp{dmultinom}, it defaults to \samp{sum(x)}. The
first element correspond to the vector \samp{prob} or the first row
of \samp{prob}, ...}
\item{prob}{Numeric non-negative vector of length K, or matrix of
size \samp{(x x K)} specifying the probability for the K classes; is
internally normalized to sum 1.}
\item{log}{Logical; if TRUE, log probabilities are computed.}
}
\details{
These functions are the vectorized versions of
\code{\link{rmultinom}} and \code{\link{dmultinom}}. Recycling is
permitted.
}
\keyword{distribution}
