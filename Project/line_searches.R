### line search functions
# 
# Rodolphe Le Riche, CNRS LIMOS, july 2021
# 

### Backtracking with Armijo (sufficient decrease) condition
#   find stepSize that satisfies : 
#     f(x+stepSize*direction) <= f(x) + stepSize* (suffDecFact * direction^T*gradf)
#   If direction is not a descend direction ( direction^T*gradf>0 ), 
#     turn around (multiply direction by -1).
# 
# INPUTS
# 
#   x, fofx, gradf : current point, its objective function, the function gradient at x
#   direction : direction in which to search. Not necessarily a descent direction (hence works with 
#     momentum and NAG).
#   suffDecFact : scalar in [0,1[  used in the sufficient decrease test, 
#     i.e. how much better than Taylor is demanded. 
#   decFact : scalar in ]0,1[ , used to reduce stepSize
#   initStepFact : multiplicative factor of gradient norm to determine initial stepSize
#   LB,UB : d-dimensional vectors of lower and upper bounds for the variables
#   nbFun : total number of calls to the objective function before the line search started
#
# OUTPUTS
#
#   res$xnew : point solution of the line search (or last iterate if failure)
#   res$nFcalls : number of calls to the objective function during line search
#   res$rec : complete record (list with $X, $F and $Time fields) of the points 
#     calculated during line search.
# 

BacktrackLineSearch <- function(x,fofx,gradf,direction,f,
                                suffDecFact=0.1,decFact=0.5,initStepFact=1,LB,UB,printlevel,nbFun)
{

  normGrad <- l2norm(gradf)
  # calculate initial stepSize
  # either as initStepFact*norm of gradient (but this may fail in flat regions)
  # or as a fraction of domain diagonal. Take the max of both initial step sizes.
  stepSize <- max(initStepFact*normGrad,(l2norm(UB-LB)/100))
  decConst <- suffDecFact*(direction%*%gradf)
  # if direction is not a descent direction, -direction is, use turnaround for this
  turnaround <- 1 
  if (decConst>0){
    turnaround <- -1
    decConst <- -decConst
  }
  maxloop <- 100 # max line search budget
  nloop <- 1 # initialize loop counter
  #
  res <- list()
  if (printlevel>=3) {lrec <- list()}
  fp <- .Machine$double.xmax # a very large number to get into the while loop
  
  while ( (fp > fofx+stepSize*decConst) & (nloop<maxloop)) {
    xpp <- x+stepSize*turnaround*direction
    # project on bounds
    xp <- ifelse(xpp < LB, LB, ifelse(xpp > UB, UB, xpp))
    if (l2norm(xpp-xp)<1.e-10){ # only evaluate if point is in bounds,
      # otherwise just decrease stepSize
      fp <- f(xp)
      nloop <- nloop+1
      if (printlevel>=3) {lrec<-updateRec(rec=lrec,x=xp,f=fp,t=nbFun+nloop)}
    }
    stepSize <- stepSize*decFact
  } # end while loop
  
  if (nloop >= maxloop){ 
    msg <- paste("nloop=",nloop," larger than maxloop=",maxloop, ", output last iterate")
    warning(msg)
  }
  res$xnew<-xp
  res$nFcalls <- nloop
  if (printlevel>=3) {res$rec <- lrec}
  return(res)
} ### end BacktrackLineSearch function

