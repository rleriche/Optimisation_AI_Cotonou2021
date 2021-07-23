# set of utilities for simple optimization program
# 
# Rodolphe Le Riche, CNRS LIMOS, july 2021
# 

# calculate f and gradient of f by forward finite difference
f.gradf <- function(x,f,h=1.e-8){
  d<-length(x)
  res <- list()
  res$gradf <- rep(NA,d)
  res$fofx <- f(x)
  for (i in 1:d){
    xp <- x
    xp[i] <- x[i]+h
    res$gradf[i] <- (f(xp)-res$fofx)/h
  }
  return(res)
}

# record points online
updateRec <- function(rec,x,f,t){
  rec$X <- rbind(rec$X,x)
  rec$F <- c(rec$F,f)
  rec$Time <- c(rec$Time,t)
  return(rec)
}

# L2 norm
l2norm <- function(x){
  return(sqrt(sum(x^2)))
}
