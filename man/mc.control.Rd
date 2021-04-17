\name{mc.control}
\alias{mc.control}
\alias{ndvar}
\alias{ndunc}
\title{
Sets or Gets the Default Number of Simulations.
}
\description{
Sets or retrieves the default number of simulations.
}
\usage{
ndvar(n)
ndunc(n)
}
\examples{
(oldvar <- ndvar())
(oldunc <- ndunc())
mcstoc(runif, type="VU")
ndvar(12)
ndunc(21)
mcstoc(runif, type="VU")
ndvar(oldvar)
ndunc(oldunc)

}
\arguments{
\item{n}{Number of simulations.}
}
\details{
\samp{ndvar()} gets and \samp{ndvar(n)} sets the default number of
simulation in the 1D simulations or the number of simulation in the
variability dimension in the 2D simulations.

 \samp{ndunc()} gets and \samp{ndunc(n)} sets the number of
simulations in the uncertainty dimension in the 2D simulations.

 \samp{n} is rounded to its ceiling value.

 The default values when loaded are 1001 for \samp{ndvar} and 101 for
\samp{ndunc}.
}
\value{
The current value, AFTER modification if \samp{n} is present (!=
\samp{options}).
}
\keyword{misc}
