\name{evalmcmod}
\alias{evalmcmod}
\title{
Evaluates a Monte-Carlo model
}
\description{
Evaluates a \code{\link{mcmodel}} object (or a valid expression)
using a specified number of simulations and with (or without) a
specified seed.
}
\usage{
evalmcmod(expr, nsv=ndvar(), nsu=ndunc(), seed=NULL)
}
\examples{
data(ec)
ec$modEC1
evalmcmod(ec$modEC1, nsv=100, nsu=100, seed=666)
}
\arguments{
\item{expr}{A model of class \code{\link{mcmodel}} or a valid
expression.}
\item{nsv}{The number of simulations in the dimension of variability
used in the evaluation.}
\item{nsu}{The number of simulations in the dimension of uncertainty
used in the evaluation.}
\item{seed}{The random seed used for the evaluation. If \samp{NULL}
the \samp{seed} is unchanged.}
}
\details{
The model is evaluated. The intermediate variables used to build the
\samp{mc} object are not stored.
}
\value{
The results of the evaluation. It should be a \samp{mc} object.
}
\seealso{
\code{\link{mcmodel}}

 \code{\link{evalmccut}} to evaluate high dimension Monte Carlo Model
in a loop.
}
\note{
The seed is set at the beginning of the evaluation. Thus, the
complete similarity of two evaluations with similar seed is not
certain, depending on the structure of your model.
}
\keyword{methods}
