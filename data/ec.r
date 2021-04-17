ec <- list(
  modEC1 = expression({
    conc <- mcdata(10,"0")
    cook <- mcstoc(rempiricalD,values=c(0,1/5,1/50),prob=c(0.027,0.373,0.600))
    serving <- mcstoc(rgamma,shape=3.93,rate=0.0806)
    expo <- conc * cook * serving
    dose <- mcstoc(rpois, lambda=expo)
    risk <- 1-(1-0.001)^dose
   mc(conc,cook,serving,expo,dose,risk)
   }) ,

  modEC2 = expression({
    conc <- mcstoc(rnorm,type="U",mean=10,sd=2)
    cook <- mcstoc(rempiricalD,type="V",values=c(0,1/5,1/50),prob=c(0.027,0.373,0.600))
    serving <- mcstoc(rgamma,type="V",shape=3.93,rate=0.0806)
    expo <- conc * cook * serving
    dose <- mcstoc(rpois,type="VU",lambda=expo)
    r <- mcstoc(runif,type="U",min=0.0005,max=0.0015)
    risk <- 1-(1-r)^dose
    mc(conc,cook,serving,expo,dose,r,risk)
  })
)