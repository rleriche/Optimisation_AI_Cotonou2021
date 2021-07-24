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
updateRec <- function(arec,x,f,t){
  if (is.null(arec$X)) {arec$X<-x} else {arec$X <- rbind(arec$X,x)}
  if (is.null(arec$F)) {arec$F<-f} else {arec$F <- c(arec$F,f)}
  if (is.null(arec$Time)) {arec$Time<-t} else {arec$Time <- c(arec$Time,t)}
  return(arec)
}

# L2 norm
l2norm <- function(x){
  return(sqrt(sum(x^2)))
}
