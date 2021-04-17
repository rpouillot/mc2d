\name{typemcnode}
\alias{typemcnode}
\title{
Provides the Type of a mcnode Object
}
\description{
Provide the type of a \samp{mcnode} object.
}
\usage{
typemcnode(x, index=FALSE)
}
\examples{
data(total)
typemcnode(total$xVUM2)
}
\arguments{
\item{x}{a \samp{mcnode} object}
\item{index}{if \samp{TRUE} give the index of the type rather than
the type.}
}
\value{
\samp{"0", "V","U" or "VU"} or the corresponding index if
\samp{index=TRUE}.

 \samp{NULL} if none of this element is found.
}
\note{
This function does not test if the object is correct. See
\code{\link{is.mcnode}}.
}
\keyword{utilities}
