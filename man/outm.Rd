\name{outm}
\alias{outm}
\title{
Output of Nodes
}
\description{
Changes the output of Nodes
}
\usage{
outm(x, value="each", which.node=1)
}
\examples{
data(total)
total$xVUM2
## since outm = NULL
summary(total$xVUM2) 
x <- outm(total$xVUM2, c("min"))
summary(x)
}
\arguments{
\item{x}{A \samp{mcnode} or a \samp{mc} object.}
\item{value}{The output of the \samp{mcnode} for multivariates nodes.
May be "each" (default) if output should be provided for each
variates considered independently, "none" for no output or a vector
of name of function(s) (as a character string) that will be applied
on the variates dimension before any output (ex: \samp{"mean"},
\samp{"median"}, \samp{c("min","max")}). The function should have no
other arguments and send one value per vector of values (ex. do not
use \samp{"range"}).}
\item{which.node}{which node should be changed in a \samp{mc} object}
}
\value{
\samp{x} with a modified \samp{outm} attribute.
}
\keyword{misc}
