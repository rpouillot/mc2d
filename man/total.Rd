\name{total}
\docType{data}
\alias{total}
\alias{x0}
\alias{x0M}
\alias{xU}
\alias{xUM}
\alias{xV}
\alias{xVM}
\alias{xVU}
\alias{xVUM}
\alias{xVUM2}
\alias{xVUM3}
\title{An Example of all Kind of mcnode}
\description{
An example for each kind of \samp{mcnode}s. They are used in some \samp{mc2d} examples. They
have been built using the following code:

\code{ndvar(101)}
\code{ndunc(51)}

\code{x0 <- mcstoc(type="0")}

\code{xV <- mcstoc(type="V")}

\code{xU <- mcstoc(type="U")}

\code{xVU <- mcstoc(type="VU")}

\code{x0M <- mcstoc(type="0",nvariates=2)}

\code{xVM <- mcstoc(type="V",nvariates=2)}

\code{xUM <- mcstoc(type="U",nvariates=2)}

\code{xVUM <- mcstoc(type="VU",nvariates=2)}

\code{xVUM[c(1,12,35)] <- NA}
                              
\code{xVUM2 <- mcstoc(type="VU",nvariates=2,outm="none")}

\code{xVUM3 <- mcstoc(type="VU",nvariates=2,outm=c("mean","min"))}

\code{total <- mc(x0,xV,xU,xVU,x0M,xVM,xUM,xVUM,xVUM2,xVUM3)}
}
\usage{data(total)}
\format{
  Some \samp{mcnode} objects and one \samp{mc} object.
}
\source{None}
\references{
  None
}
\keyword{datasets}