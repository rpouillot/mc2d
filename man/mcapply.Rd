\name{mcapply}
\alias{mcapply}
\title{
Apply Functions Over mc or mcnode Objects
}
\description{
Apply a function on all values or over a given dimension of an
\samp{mcnode} object. May be used for all \samp{mcnode} of an
\samp{mc} object.
}
\usage{
mcapply(x, margin=c("all", "var", "unc", "variates"), fun, ...)
}
\examples{
data(total)
xVUM
mcapply(xVUM, "unc", sum)
mcapply(xVUM, "var", sum)
mcapply(xVUM, "all", sum)
mcapply(xVUM, "variates", sum)
mcapply(total, "all", exp)
}
\arguments{
\item{x}{A \samp{mc} or a \samp{mcnode} object.}
\item{margin}{The dimension on which applying the function. Maybe
\samp{"all"} (default) to apply the function on all values,
\samp{"var"} to apply the function on the variability dimension,
\samp{"unc"} to apply the function on the uncertainty dimension, or
\samp{"variates"} to apply the function on the variates. Watch out:
do not use 'var' for 'variates'}
\item{fun}{The function to be applied. When applied to a vector of
length \samp{n}, \samp{fun} should return a vector of length \samp{n}
or \samp{1}.}
\item{\dots}{Optional arguments to \samp{fun}.}
}
\value{
If \samp{fun} returns a function of length \samp{n} or if
\samp{margin="all"}, the returned \samp{mcnode}s are of type and
dimension of \samp{x}. In other cases, the type of \samp{mcnode} is
changed.
}
\seealso{
\code{\link{apply}}, \code{\link{mc}}, \code{\link{mcnode}}.
}
\keyword{misc}
