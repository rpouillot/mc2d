\name{extractvar}
\alias{extractvar}
\alias{addvar}
\title{
Utilities for multivariate nodes
}
\description{
\samp{extractvar} extracts one variate from a multivariate node. 

 \samp{addvar} adds consistent \samp{mcnode}s to build a multivariate
\samp{mcnode} .
}
\usage{
extractvar(x, which=1)
addvar(...)
}
\examples{
x <- mcdata(0:3, "0", nvariates = 4)
y <- extractvar(x, c(1, 3)) 
y
addvar(x, y)

}
\arguments{
\item{x}{a multivariates \samp{mcnode}.}
\item{which}{a vector. which variate(s) should be extracted?}
\item{\dots}{ \samp{mcnode}s to be gathered in a multivariate
\samp{mcnode}. These \samp{mcnode}s should be of same type and
dimension.}
}
\details{
The \samp{outm} attribute of the output of \samp{addvar} will be the
one of the first element.
}
\value{
The new built \samp{mcnode}.
}
\seealso{
\code{\link{mcnode}} for \samp{mcnode} objects.
}
\keyword{methods}
