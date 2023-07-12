\name{ec}
\docType{data}
\alias{ec}
\alias{modEC1}
\alias{modEC2}
\title{An example on Escherichia coli in ground beef}
\description{
The fictive example is as following:

A batch of ground beef is contaminated with \emph{E. coli},
with a mean concentration \samp{conc}.

Consumers may eat the beef "rare", "medium rare" or
"well cooked". If "rare", no bacteria is killed. If "medium rare", 1/5 of bacteria survive. If "well cooked", 1/50 of bacteria survive.

The serving size is variable.

The risk of infection follows an exponential model.

For the one-dimensional model, it is assumed that:

\code{conc <- 10}

\code{cook <- sample(n, x=c(1,1/5,1/50),replace=TRUE,prob=c(0.027,0.373,0.600))}

\code{serving <- rgamma(n, shape=3.93,rate=0.0806)}

\code{expo <- conc * cook * serving}

\code{dose <- rpois(n, lambda=expo)}

\code{risk <- 1-(1-0.001)^dose}

For the two-dimensional model, it is assumed moreover that the concentration and the \samp{r} parameter
of the dose response are uncertain.

\code{conc <- rnorm(n,mean=10,sd=2)}
  
\code{r <- runif(n ,min=0.0005,max=0.0015)}
  
}
\usage{data(ec)}
\format{
  A list of two expression to be passed in \code{mcmodel}
}
\source{Fictive example}
\references{
  None
}
\keyword{datasets}