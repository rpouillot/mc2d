\name{mc}
\alias{mc}
\title{
Monte Carlo Object
}
\description{
Creates \samp{mc} objects from \code{\link{mcnode}} or \samp{mc}
objects.
}
\usage{
mc(..., name=NULL, devname=FALSE)
}
\examples{
 
x <- mcstoc(runif)
y <- mcdata(3, type="0")
z <- x * y
(m <- mc(x, y, z, name=c('n1', 'n2', 'n3')))
mc(m, x, devname=TRUE)
}
\arguments{
\item{...}{\samp{mcnode} and/or \samp{mc} object(s) to be gathered in
a \samp{mc} object separated by a coma.}
\item{name}{Vector of character of the same length of the final
\samp{mc} object. If NULL, the name will be given from the name of
the elements. }
\item{devname}{Develop the name from the name of the \samp{mc}
objects, if any.}
}
\details{
A \samp{mc} object is a list of \code{\link{mcnode}} objects.
\samp{mcnode} objects must be of coherent dimensions.

If one of the arguments is a \samp{mc} object, the name of the
elements of this \samp{mc} object are used. \samp{devname = TRUE}
will develop the name, using as a prefix the name of the \samp{mc}
object.

 Finally, names are transformed to be unique.
}
\value{
An object of class \samp{mc}.
}
\seealso{
\code{\link{mcnode}}, the basic element of a \samp{mc} object.

 To evaluate \samp{mc} objects: \code{\link{mcmodel}},
\code{\link{evalmcmod}}, \code{\link{evalmccut}}

 Informations about an \samp{mc} object: \code{\link{is.mc}},
\code{\link{dimmc}}

 To study \samp{mc} objects: \code{\link{print.mc}},
\code{\link{summary.mc}}, \code{\link{plot.mc}},
\code{\link{converg}}, \code{\link{hist.mc}}, \code{\link{tornado}},
\code{\link{tornadounc.mc}}
}
\keyword{methods}
