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

# plot contour of function when d==2
plot_contour <- function(LB,UB,f){
  no.grid <- 100
  x1 <- seq(LB[1], UB[1], length.out=no.grid)
  x2 <- seq(LB[2], UB[2], length.out=no.grid)
  x.grid <- expand.grid(x1, x2)
  z <- apply(x.grid, 1, f)
  z.grid <- matrix(z, no.grid)
  contour(x1, x2, z.grid, nlevels=20, xlab="x1", ylab="x2")
}

# increasing sequence: useful to plot first points and then fewer and fewer
inc.geom.seq <- function(from=1,to=1000,coef=1.4)
{
  s <- c(round(from))
  x <- from
  i <- 1
  ieff <- 1
  done <- FALSE
  while (!done){
    x <- x*coef
    sp <- round(x)
    if (sp != s[ieff]){
      s <- c(s,sp)
      ieff <- ieff+1
      if (sp>to) done<-TRUE
    }
    i <- i+1
  }
  s<-s[-ieff]
  return(s)
}