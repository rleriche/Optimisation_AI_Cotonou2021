### line search functions
#
# Implementation notes: 
# * Some variables come from the calling environment:
#   - LB, UB : lower and upper bounds
#   - data used for recording calls to f: nbFun (nb of calls to f before routine is called)
#     printlevel, rec (the list used for recording all calls)


### Backtracking with Armijo (sufficient decrease) condition
#   s must satisfy : f(x+s*d) <= f(x) + stepSize* (suffDecFact * d^T*gradf)
#   suffDecFact in [0,1[  used in the sufficient decrease test 
#   decFact in ]0,1[   used to reduce stepSize
BacktrackLineSearch <- function(x,fofx,gradf,direction,f,
                                suffDecFact=0.1,decFact=0.5,initStepFact=1){
  res <- list()
  normGrad <- sqrt(sum(gradf^2))
  # calculate initial stepSize
  # either as initStepFact*norm of gradient (but this may fail in flat regions)
  # or as a fraction of domain diagonal
  stepSize <- max(initStepFact*normGrad,(l2norm(UB-LB)/100))
  decConst <- suffDecFact*(direction%*%gradf)
  maxloop <- 100 # max line search budget
  #
  xpp <- x+stepSize*direction
  # project on bounds
  xp <- ifelse(xpp < LB, LB, ifelse(xpp > UB, UB, xpp))
  fp <- f(xp)
  nloop <- 1
  if (printlevel>=2) {
    lrec <- rec
    lrec<-updateRec(rec=lrec,x=xp,f=fp,t=nbFun+nloop)
  }
  while ( (fp > fofx+stepSize*decConst) & (nloop<maxloop)) {
    stepSize <- stepSize*decFact
    xpp <- x+stepSize*direction
    # project on bounds
    xp <- ifelse(xpp < LB, LB, ifelse(xpp > UB, UB, xpp))
    if (l2norm(xpp-xp)<1.e-10){ # only evaluate if point is in bounds,
      # otherwise decrease stepSize
      fp <- f(xp)
      nloop <- nloop+1
      if (printlevel>=2) {lrec<-updateRec(rec=lrec,x=xp,f=fp,t=nbFun+nloop)}
    }
  }
  if (nloop >= maxloop){ 
    msg <- paste("nloop=",nloop," larger than maxloop=",maxloop)
    warning(msg)
  }
  res$xnew<-xp
  res$nFcalls <- nloop
  if (printlevel>=2) {res$rec <- lrec}
  return(res)
} ### end BacktrackLineSearch function
