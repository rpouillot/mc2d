#<<BEGIN>>
dpert <- function(x,min=-1,mode=0,max=1,shape=4,log=FALSE)
#TITLE The (Modified) PERT Distribution
#NAME pert
#KEYWORDS distribution
#DESCRIPTION
#Density, distribution function, quantile function and random generation
#for the PERT (\emph{aka} Beta PERT) distribution with minimum equals to \samp{min}, mode equals to \samp{mode}
#and maximum equals to \samp{max}. 
#INPUTS
#{x,q}<<Vector of quantiles.>>
#{p}<<Vector of probabilities.>>
#{n}<<Number of observations. If length(n) > 1, the length is taken to be the number required.>>
#[INPUTS]
#{min}<<Vector of minima.>>
#{mode}<<Vector of modes.>>
#{max}<<Vector of maxima.>>
#{shape}<<Vector of scaling parameters. Default value: 4.>>
#{log, log.p}<<Logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.>>
#{lower.tail}<<Logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.>>
#DETAILS
#The PERT distribution is a beta distribution extended to the domain \samp{[min, max]} with mean 
#\deqn{\mu=\frac{min+shape\times mode+max}{shape+2}}{mu = (min + shape * mode + max)/(shape + 2)}
#
#The underlying beta distribution is specified by \eqn{\alpha_{1}}{shape1} and \eqn{\alpha_{2}}{shape2} defined as
#
#\deqn{\alpha_{1}=\frac{(\mu-min)(2\times mode-min-max)}{(mode-\mu)(max-min)}}{shape1=(mu - min)*(2 mode-min-max)/((mode-mu)*(max - min))}
#
#\deqn{\alpha_{2}=\frac{\alpha_{1}\times (max-\mu)}{mu-min}}{shape2=shape1*(max - mu)/(mu - min)}
#
#If \eqn{\mu=mode}{mu=mode}, \eqn{\alpha_{1}}{shape1} is set to \eqn{1+\nu/2}{1+shape/2}.
#
#David Vose (See reference) proposed a modified PERT distribution with a shape parameter different from 4. 
#
#The PERT distribution is frequently used (with the \link{triangular} distribution) to translate expert estimates
#of the min, max and mode of a random variable in a smooth parametric distribution. 
#REFERENCE
#Vose D. Risk Analysis - A Quantitative Guide (2nd and 3rd editions, John Wiley and Sons, 2000, 2008).  
#AUTHOR Regis Pouillot and Matthew Wiener
#EXAMPLE
#curve(dpert(x,min=3,mode=5,max=10,shape=6), from = 2, to = 11, lty=3)
#curve(dpert(x,min=3,mode=5,max=10), from = 2, to = 11, add=TRUE)
#curve(dpert(x,min=3,mode=5,max=10,shape=2), from = 2, to = 11, add=TRUE,lty=2)
#legend(x = 8, y = 2, c("Default","shape:2","shape:6"), lty=1:3)
#SEE ALSO
#\code{\link{Beta}}
#VALUE
#\samp{dpert} gives the density, \samp{ppert} gives the distribution function,
#\samp{qpert} gives the quantile function, and \samp{rpert} generates random deviates.

#CREATED 08-02-20
#MODIFIED 13-07-10
#--------------------------------------------
{
	if(length(x) == 0) return(numeric(0))
  
  min <- as.vector(min)
	mode <- as.vector(mode)
	max <- as.vector(max)
	shape <- as.vector(shape)
  
	a1 <- 1+shape*(mode-min)/(max-min)	
	a2 <- 1+shape*(max-mode)/(max-min)
	
	oldw <- options(warn = -1)
	d <- (x-min)^(a1-1) * (max-x)^(a2-1) / beta(a=a1,b=a2) / (max-min)^(a1+a2-1)	
	options(warn = oldw$warn)
	
  d[x < min | x > max] <- 0
	d[mode < min | max < mode] <- NaN
	d[shape <= 0] <- NaN
	if(log) d <- log(d)
	if(any(is.na(d))) warning("NaN in dpert")
  return(d)}

#<<BEGIN>>
ppert <- function(q,min=-1,mode=0,max=1,shape=4,lower.tail = TRUE, log.p = FALSE)
#ISALIAS dpert
#--------------------------------------------
{
	if(length(q) == 0) return(numeric(0))

  min <- as.vector(min)
	mode <- as.vector(mode)
	max <- as.vector(max)
	shape <- as.vector(shape)
	
	a1 <- 1+shape*(mode-min)/(max-min)	
	a2 <- 1+shape*(max-mode)/(max-min)
  
	oldw <- options(warn = -1)
	p <- pbeta(q=(q-min)/(max-min),shape1=a1,shape2=a2)
	options(warn = oldw$warn)
	
  p[q < min] <- 0
	p[q >= max] <- 1
	p[mode < min | max < mode] <- NaN
	p[shape <= 0] <- NaN
	
  if(!lower.tail) p <- 1-p
	if(log.p) p <- log(p)
	if(any(is.na(p))) warning("NaN in ppert")
  return(p)}

#<<BEGIN>>
qpert <- function(p,min=-1,mode=0,max=1,shape=4,lower.tail=TRUE,log.p=FALSE)
#ISALIAS dpert
#--------------------------------------------
{
    if (length(p) == 0) 
      return(numeric(0))
    min <- as.vector(min)
    mode <- as.vector(mode)
    max <- as.vector(max)
    shape <- as.vector(shape)
    
    lout <- max(length(p),length(min),length(mode),length(max),length(shape))
    min <- rep(min, length.out=lout)
    mode <- rep(mode, length.out=lout)
    max <- rep(max, length.out=lout)
    shape <- rep(shape, length.out=lout)
    
    if (log.p) 
      p <- exp(p)
    if (!lower.tail) 
      p <- 1 - p
    a1 <- 1 + shape * (mode - min)/(max - min)
    a2 <- 1 + shape * (max - mode)/(max - min)
    oldw <- options(warn = -1)
    q <- qbeta(p, shape1 = a1, shape2 = a2)
    options(warn = oldw$warn)
    q <- q * (max - min) + min
    minmodemax <- (abs(min - max) < (.Machine$double.eps^0.5))
    q <- ifelse(minmodemax, min, q)
    q[p < 0 | p > 1] <- NaN
    q[mode < min | max < mode] <- NaN
    q[shape <= 0] <- NaN
    if (any(is.na(q))) 
      warning("NaN in qpert")
    return(q)
  }



#<<BEGIN>>
rpert <- function(n,min=-1,mode=0,max=1,shape=4)
#ISALIAS dpert
#--------------------------------------------
{
  if (length(n) > 1) 
    n <- length(n)
  if (length(n) == 0 || as.integer(n) == 0) 
    return(numeric(0))
  n <- as.integer(n)
  
  min <- rep(as.vector(min),length.out=n)
  mode <- rep(as.vector(mode),length.out=n)
  max <- rep(as.vector(max),length.out=n)
  shape <- rep(as.vector(shape),length.out=n)
  
  a1 <- 1 + shape * (mode - min)/(max - min)
  a2 <- 1 + shape * (max - mode)/(max - min)
  oldw <- options(warn = -1)
  r <- rbeta(n, shape1 = a1, shape2 = a2) * (max - min) + min
  options(warn = oldw$warn)
  minmodemax <- (abs(min - max) < (.Machine$double.eps^0.5))
  r <- ifelse(minmodemax, min, r)
  if (any(is.na(r))) 
    warning("NaN in rpert")
  return(r)
}
