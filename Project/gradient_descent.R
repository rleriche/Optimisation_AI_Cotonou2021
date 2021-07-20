#############################################
# Two gradient descent algorithms
# 1) Gradient descent where step size depends on gradient norm
# 2) Steepest descent where step size depends on a line search
#
# Rodolphe Le Riche
#
# Implementation notes: 
# * no attempt at making this code efficient, it is for teaching purpose.
# * the "descent" algorithm should be preferred over the "gradient" because
#   the line search robustifies a lot the search
# * delete the global variables if you change dimension and optimise quadratic function.
#############################################
rm(list=ls()) # clear environment
source('test_functions.R')
source('utilities_optim.R')

### problem definition
# search space
d<-2 # dimension
LB<-rep(-5,d) #lower bounds
UB<-rep(5,d) #upper bounds
fun<-rosen #function to minimize

### algorithm settings
xinit <- c(2.2,3.5) #rep(-4.9,d) # initial point
algo_type <- "descent" # choices are : "gradient" or "descent"
#   Both algorithms have minus the normalized gradient as search direction
#       x_{t+1} <- x_t + stepSize*direction
#   for "gradient" :  stepSize = stepFactor*normGrad
#   for "descent" : stepSize comes from a linesearch where 
#   sufficientDecreaseFactor is a control parameter
stepFactor <- 0.2 # step factor for "gradient" version
sufficientDecreaseFactor <- 0.1 # controls stepSize for "descent" version
#
printlevel <- 2 # controls how much is stored and printed
#                 =1 store best and minimal output
#                 =2 store all points and more outputs
stopBudget <- 100 # maximum number of function evaluations:
#   the maximum number of calls to the function is larger because of the 
#   finite differences scheme: nbFun = (d+1)*iter+nbFunLS
stopGradNorm <- 1.e-6 # stop when gradient norm / sqrt(d) is smaller than 

### initializations
x <- xinit
eval <- f.gradf(x=x,f=fun,h=1.e-8) #eval$fofx is the function 
                          # eval$gradf the associated gradient
normGrad <- sqrt(sum(eval$gradf^2))
iter <- 1
nbFun <- 1
nbFunLS <- 0 # calls to fun done during the line search
# recordings of best so far
Fbest <- eval$fofx
recBest <- list()
recBest$X <- matrix(x,nrow=1)
recBest$F <- Fbest
recBest$Time <- c(nbFun)
# below are recordings of the whole history. Memory consuming. 
# If memory issues, set printlevel<2
if (printlevel >= 2){
  rec <- list()
  rec$X <- matrix(x,nrow=1)
  rec$F <- eval$fofx
  rec$Time <- nbFun
}

### line search function
# backtracking with Armijo (sufficient decrease) condition
# s must satisfy : f(x+s*d) <= f(x) + stepSize* (suffDecFact * d^T*gradf)
#   suffDecFact in [0,1[  used in the sufficient decrease test 
#   decFact in ]0,1[   used to reduce stepSize
#   data for recording is copied from global environnement: rec, printlevel, nbFun
#   UB,LB also come from global environment (a bit ugly)
BacktrackLineSearch <- function(x,fofx,gradf,direction,f,suffDecFact=0.1,decFact=0.5,initStepFact=1){
  res <- list()
  normGrad <- sqrt(sum(gradf^2))
  stepSize <- initStepFact*normGrad
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
    fp <- f(xp)
    nloop <- nloop+1
    if (printlevel>=2) {lrec<-updateRec(rec=lrec,x=xp,f=fp,t=nbFun+nloop)}
  }
  if (nloop >= maxloop){ 
    msg <- paste("nloop=",nloop," larger than maxloop=",maxloop)
    warning(msg)
  }
  res$stepSize<-stepSize
  res$nFcalls <- nloop
  if (printlevel>=2) {res$rec <- lrec}
  return(res)
} ### end BacktrackLineSearch function


### run the algo
if (printlevel >=1) {cat("Start gradient search\n")}
while ((nbFun <= stopBudget) & ((normGrad/sqrt(d)) > stopGradNorm) ){
  
  direction <- -eval$gradf/normGrad # search direction
  # no line search, step size proportional to gradient norm
  if (algo_type=="gradient"){
  stepSize <- stepFactor*normGrad
  } 
  else if (algo_type=="descent"){
    lsres <- BacktrackLineSearch(x,eval$fofx,eval$gradf,direction,f=fun,suffDecFact=0.5)
    stepSize <- lsres$stepSize
    nbFun <- nbFun+lsres$nFcalls-1 # -1 for copy of last point, see below
    nbFunLS <- nbFunLS + lsres$nFcalls-1 
    if (printlevel>=2) {
      rec<-lsres$rec
      # delete last point because it will be recalculated 
      # (cheating to save code but the rec is correct then)
      rec<- shave1Rec(rec=rec)
    }
  }
  else {
    stop("unknown algo_type")
  }
  xcandidate <- x + stepSize*direction
  # project on bounds if necessary
  xnew <- ifelse(xcandidate < LB, LB, ifelse(xcandidate > UB, UB, xcandidate))
  # evaluate new point
  eval <- f.gradf(x=xnew,f=fun,h=1.e-8)
  normGrad <- sqrt(sum(eval$gradf^2))
  iter <- iter+1
  nbFun <- nbFun+1
  # make the step
  x <- xnew
  # bookkeeping
  if (eval$fofx < Fbest) {
    Fbest <- eval$fofx
    recBest <- updateRec(rec=recBest,x=xnew,f=Fbest,t=nbFun)
  }
  if (printlevel >= 1){
    if (printlevel >=2 ) {rec <- updateRec(rec=rec,x=xnew,f=eval$fofx,t=nbFun)}
    cat(" iteration : ",iter,"\r")
  } 

} # end while of main loop 
if (printlevel >= 1){
  cat("gradient search exited after ",iter," iterations, ",nbFun," fct evaluations\n")
  lrecBest <- dim(recBest$X)[1]
  cat('best x:',recBest$X[lrecBest,],'\n')
  cat('best f:',recBest$F[lrecBest],'\n')
}

### Vizualization
plot(x = recBest$Time,y=recBest$F,type = "l",xlab = "nb. iterations",ylab="f",col="red")
if (printlevel >= 2){
  lines(x = rec$Time,y = rec$F,col="blue")
}
if (d==2) { 
  # the code below is mainly a duplicate of what is in 3Dplots ... 
  no.grid <- 100
  x1 <- seq(LB[1], UB[1], length.out=no.grid)
  x2 <- seq(LB[2], UB[2], length.out=no.grid)
  x.grid <- expand.grid(x1, x2)
  z <- apply(x.grid, 1, fun)
  z.grid <- matrix(z, no.grid)
  # png(filename="./contour.png") # save the contour in the current directory
  contour(x1, x2, z.grid, nlevels=20, xlab="x1", ylab="x2")
  # with search points on top
  if (printlevel >=2){
    points(rec$X[,1], rec$X[,2], pch=20, col="blue")
    text(rec$X, labels=1:iter, pos=3, cex=1.0) # (un)comment for labeling (or not) nb of calls to f when points created
    points(recBest$X[,1], recBest$X[,2], pch=19, col="red") # bests so far in red (should resemble iterates)
  }
  # dev.off()
}
