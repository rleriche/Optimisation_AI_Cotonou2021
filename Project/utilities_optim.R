# set of utilities for simple optimization program

# forward finite difference
ffinite.diff <- function(x,f,h=1.e-8){
  d<-length(x)
  df <- rep(NA,d)
  for (i in 1:d){
    xp <- x
    xp[i] <- x[i]+h
    df[i] <- (f(xp)-f(x))/h
  }
  return(df)
}

# project on bounds