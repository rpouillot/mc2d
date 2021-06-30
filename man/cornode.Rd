\name{cornode}
\alias{cornode}
\title{
Builds a Rank Correlation using the Iman and Conover Method.
}
\description{
This function builds a rank correlation structure between columns of
a matrix or between \samp{mcnode} objects using the Iman and Conover
method (1982).
}
\usage{
cornode(..., target, outrank=FALSE, result=FALSE, seed=NULL)
}
\examples{
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rnorm(1000)
mat <- cbind(x1, x2, x3)
## Target
(corr <- matrix(c(1, 0.5, 0.2, 0.5, 1, 0.2, 0.2, 0.2, 1), ncol=3))
## Before
cor(mat, method="spearman")
matc <- cornode(mat, target=corr, result=TRUE)
## The first row is unchanged
all(matc[, 1] == mat[, 1])

##Using mcnode and outrank
cook <- mcstoc(rempiricalD, values=c(0, 1/5, 1/50), prob=c(0.027, 0.373, 0.600), nsv=1000)
serving <- mcstoc(rgamma, shape=3.93, rate=0.0806, nsv=1000)
roundserv <- mcdata(round(serving), nsv=1000)
## Strong relation between roundserv and serving (of course)
cor(cbind(cook, roundserv, serving), method="spearman")

##The classical way to build the correlation structure 
matcorr <- matrix(c(1, 0.5, 0.5, 1), ncol=2)
matc <- cornode(cook=cook, roundserv=roundserv, target=matcorr)
## The structure between cook and roundserv is OK but ...
## the structure between roundserv and serving is lost
cor(cbind(cook=matc$cook, serv=matc$roundserv, serving), method="spearman")

##An alternative way to build the correlation structure
matc <- cornode(cook=cook, roundserv=roundserv, target=matcorr, outrank=TRUE)
## Rebuilding the structure
roundserv[] <- roundserv[matc$roundserv, , ]
serving[] <- serving[matc$roundserv, , ]
## The structure between cook and roundserv is OK and ...
## the structure between roundserv and serving is preserved
cor(cbind(cook, roundserv, serving), method="spearman")
}
\arguments{
\item{...}{A matrix (each of its \samp{n} columns but the first
one will be reordered) or \samp{n mcnode} objects (each elements but
the first one will be reordered).}
\item{target}{A scalar (only if \samp{n=2}) or a \samp{(n x n)}
matrix of correlation.}
\item{outrank}{Should the order be returned?}
\item{result}{Should the correlation eventually obtained be printed?}
\item{seed}{The random seed used for building the correlation. If
\samp{NULL} the \samp{seed} is unchanged.}
}
\details{
The arguments should be named.

The function accepts for \samp{data} a matrix or: 

\itemize{
\item some \samp{"V" mcnode} objects separated by a comma;
\item some \samp{"U" mcnode} objects separated by a comma;
\item some \samp{"VU" mcnode} objects separated by a comma. In that
case, the structure is built columns by colums (the first column of
each \samp{"VU" mcnode} will have a correlation structure, the second
ones will have a correlation structure, ....).
\item one \samp{"V" mcnode} as a first element and some \samp{"VU"
mcnode} objects, separated by a comma. In that case, the structure is
built between the \samp{"V" mcnode} and each column of the \samp{"VU"
mcnode} objects. The correlation result (\samp{result = TRUE}) is not
provided in that case.}
The number of variates of the elements should be equal.

\samp{target} should be a scalar (two columns only) or a real
symmetric positive-definite square matrix. Only the upper triangular
part of \samp{target} is used (see \code{\link{chol}}).

 The final correlation structure should be checked because it is not
always possible to build the target correlation structure.

 In a Monte-Carlo simulation, note that the order of the values
within each \samp{mcnode} will be changed by this function (excepted
for the first one of the list). As a consequence, previous links
between variables will be broken. The \samp{outrank} option may help
to rebuild these links (see the Examples).
}
\value{
If \samp{rank = FALSE}: the matrix or a list of rearranged
\samp{mcnode}s. 

 If \samp{rank = TRUE}: the order to be used to rearranged the matrix
or the \samp{mcnodes} to build the desired correlation structure.
}
\references{
Iman, R. L., & Conover, W. J. (1982). A distribution-free approach to inducing rank correlation among input variables. \emph{Communication in Statistics - Simulation and Computation}, 11(3), 311-334. 
}
\keyword{multivariate}