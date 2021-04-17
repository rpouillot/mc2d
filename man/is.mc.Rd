\name{is.mc}
\alias{is.mc}
\alias{is.mcnode}
\title{
Tests mc and mcnode Objects
}
\description{
\samp{is.mc} tests \samp{mc} objects and \samp{is.mcnode} tests
\samp{mcnode} objects.
}
\usage{
is.mc(x)
is.mcnode(x)
}
\examples{
data(total)
is.mcnode(xVU)
is.mcnode(total)
is.mc(total)

}
\arguments{
\item{x}{An \samp{mc} or a \samp{mcnode} object.}
}
\details{
\samp{is.mc} tests if \samp{x} is a list of \samp{mcnode}, each
elements being of compatible dimension. It tests if the class
\samp{"mc"} is affected to the object.

 \samp{is.mcnode} tests if \samp{x} is an array of numeric or
logical, if it has a "type" attribute and compatible dimensions, and
if the class \samp{"mcnode"} is affected to the object.
}
\value{
\samp{TRUE} or \samp{FALSE}
}
\keyword{utilities}
