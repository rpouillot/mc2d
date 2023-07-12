\name{Ops.mcnode}
\alias{Ops.mcnode}
\title{
Operations on mcnode Objects
}
\description{
This function alters the way operations are performed on
\samp{mcnode} objects for a better consistency of the theory.
}
\usage{
\method{Ops}{mcnode}(e1, e2)
}
\examples{
oldvar <- ndvar()
oldunc <- ndunc()
ndvar(30)
ndunc(20)

## Given
x0 <- mcdata(3, type="0")
xV <- mcdata(1:ndvar(), type="V")
xU <- mcdata(1:ndunc(), type="U")
xVU <- mcdata(1:(ndunc()*ndvar()), type="VU")
x0M <- mcdata(c(5, 10), type="0", nvariates=2)
xVM <- mcdata(1:(2*ndvar()), type="V", nvariates=2)
xUM <- mcdata(1:(2*ndunc()), type="U", nvariates=2)
xVUM <- mcdata(1:(2*(ndunc()*ndvar())), type="VU", nvariates=2)

## All possible combinations
## "0"
-x0
x0 + 3

## "V"
-xV
3 + xV
xV * (1:ndvar())
xV * x0
xV - xV

## "U"
-xU
xU + 3
(1:ndunc()) * xU
xU * x0
xU - xU

## Watch out the resulting type
xV + xU
xU + xV

## "VU"
-xVU
3 + xVU
(1:(ndunc()*ndvar())) * xVU
xVU + xV
x0 + xVU
xU + xVU
xVU - xVU

## Some Multivariates
x0M+3
xVM * (1:ndvar())
xVM - xV
xUM - xU
xVUM - xU
}
\arguments{
\item{e1}{An \samp{mcnode} object, a vector or an array.}
\item{e2}{An optional \samp{mcnode} object, a vector or a matrix
with at least one of both objects as an \samp{mcnode}.}
}
\details{
This method will be used for any of the Group \code{\link{Ops}}
functions.

The rules are as following (illustrated with a \samp{+} function and
ignoring the \samp{nvariates} dimension): 

\itemize{
\item \samp{0 + 0 = 0};
\item \samp{0 + V = V}: classical recycling of the scalar;
\item \samp{0 + U = U}: classical recycling of the scalar;
\item \samp{0 + VU = VU}: classical recycling of the scalar;
\item \samp{V + V = V}: if both of the same \samp{(nsv)} dimension;
\item \samp{V + U = VU}: the \samp{U} object will be recycled "by
row". The \samp{V} object will be recycled classically "by column";
\item \samp{V + VU = VU}: if the dimension of the \samp{V} is
\samp{(nsv)} and the dimension of the \samp{VU} is \samp{(nsv x
nsu)}. The \samp{V} object will be recycled classically "by column";
\item \samp{U + U = U}: if both of the same \samp{(nsu)} dimension;
\item \samp{U + VU = VU}: if the dimension of the \samp{U} is
\samp{(nsu)} and the dimension of the \samp{VU} is \samp{(nsv x
nsu)}. The \samp{U} object will be recycled "by row";
\item \samp{VU + VU = VU}: if the dimension of the \samp{VU} nodes is
\samp{(nsu x nsv)};}
A vector or an array may be combined with an \samp{mcnode} of size
\samp{(nsv x nsu)} if an \samp{mcnode} of this dimension may be built
from this vector/array using the \samp{mcdata} function. See
\code{\link{mcdata}} for the rules.

The \samp{outm} attribute is transferred as following: \samp{each +
each = each}; \samp{none + other = other}; \samp{other1 + other2 =
other1}. The \samp{outm} attribute of the resulting node may be
changed using the \code{\link{outm}} function.

For multivariate nodes, a recycling on the \samp{nvariates} dimension
is done if a \samp{(nsu x nsv x nvariates)} node is combined with a
\samp{(nsu x nsv x 1)} node.
}
\value{
The results as a \samp{mcnode} object.
}
\seealso{
\code{\link{mcdata}}, \code{\link{mcstoc}}
}
\keyword{utilities}
