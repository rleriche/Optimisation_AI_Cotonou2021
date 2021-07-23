#############################################
# A family of gradient-based descent algorithms
# Rodolphe Le Riche
#
# They differ in how the search direction is calculated
#   direction_type = (gradient, momentum, NAG) , where NAG= Nesterov Accelerated Gradient
# and they differ on whether a single step is taken or a line search performed
#   linesearch_type = (none, armijo)
#
# Implementation notes: 
# * no attempt at making this code efficient, it is for teaching purpose.
# * delete the global variables if you change dimension and optimize quadratic function.
#############################################
rm(list=ls()) # clear environment
source('test_functions.R')
source('utilities_optim.R')
source('line_searches.R')

### problem definition
# search space
d<-2 # dimension
LB<-rep(-5,d) #lower bounds
UB<-rep(5,d) #upper bounds
fun<-rosen #function to minimize

### algorithm settings
xinit <- rep(-4.9,d)#c(4.5,3.5) rep(-4.9,d) # initial point
direction_type <- "NAG" # choices are : "gradient", "momentum", "NAG"
linesearch_type <- "armijo" # choices are: "none", "armijo"
#   all algorithms have a search direction and a step size
#       x_{t+1} <- x_t + stepSize*direction
#   where direction is a vector of norm = 1 and stepSize a scalar
#   The direction and the stepSize are calculated in various ways:
#   step : vector specific to direction_type
#     e.g., direction_type == gradient , step = -gradf
#     and other formula for other direction_type. Then,
#     direction = step/l2norm(step) , 
#   Without linesearch, stepSize = stepFactor*l2norm(step)
#   Expl, for gradient&none :  stepSize = stepFactor*normGrad
#   when linesearch_type is not none : stepSize comes from a linesearch, cf line_searches.R file
stepFactor <- 0.1 # step factor when there is no line search, use depends on direction
#   stepfactor ~ 1/Lipschitz_constant : the steeper the function, the smaller stepfactor
beta <- 0.9 # momentum term
#
printlevel <- 2 # controls how much is stored and printed
#                 =1 store best and minimal output
#                 =2 store all points and more outputs
stopBudget <- 10000 # maximum number of function evaluations:
#   the maximum number of calls to the function is larger because of the 
#   finite differences scheme: nbFun = (d+1)*iter+nbFunLS
stopGradNorm <- 1.e-6 # stop when gradient norm / sqrt(d) is smaller than 

### initializations
x <- xinit
eval <- f.gradf(x=x,f=fun,h=1.e-8) #eval$fofx is the function 
                          # eval$gradf the associated gradient
normGrad <- l2norm(eval$gradf)
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


### run the algo
if (printlevel >=1) {
  startmsg <- paste0("Start ",direction_type," (direction) + ",linesearch_type," (line search) descent\n")
  cat(startmsg)
  }
while ((nbFun <= stopBudget) & ((normGrad/sqrt(d)) > stopGradNorm) ){
  
  if (direction_type=="gradient"){
    step <- -eval$gradf
  }
  else if (direction_type=="momentum"){
      if (iter <= 1){
        step <- -eval$gradf
      }
      else {
        step <- -eval$gradf + beta*previous_step
      }
  }
  else if (direction_type=="NAG"){
    if (iter <= 1){
      step <- -eval$gradf
    }
    else {
      # evaluate gradient at anticipated next point
      xnag <- x + beta*(x-xprevious)
      evalnag <- f.gradf(x=xnag,f=fun,h=1.e-8)
      nbFun <- nbFun+1 # this is a really basic accounting for the numerical cost
      # of the gradient evaluation. In fact, if done with finite differences, 
      # it costs d+1 evaluations. The cost of 1 is more in the spirit of having 
      # an analytical evaluation of the gradient such as with retropropagation in NN.
      # For the same reason we don't account for xnag in the best point so far 
      # (it might only be a gradient evaluation).
      step <- -evalnag$gradf + beta*previous_step
    }
  }
  else{stop("unknown direction_type")}
  # if the current point is near a boundary, the direction should be projected on that boundary
  tol <- 1.e-15
  violation <- ifelse(x>(LB+tol),0,-1) + ifelse(x<(UB-tol),0,1)
  step[which(violation*step>0)]<-0
  #
  rawStepSize <- l2norm(step)
  direction <- step/rawStepSize
  
  # no line search, step size proportional to gradient norm
  if (linesearch_type == "none"){
    stepSize <- stepFactor*rawStepSize
    xcandidate <- x + stepSize*direction
    # project on bounds if necessary
    xnew <- ifelse(xcandidate < LB, LB, ifelse(xcandidate > UB, UB, xcandidate))
    # evaluate new point
    eval <- f.gradf(x=xnew,f=fun,h=1.e-8)
    normGrad <- l2norm(eval$gradf)
    nbFun <- nbFun+1
    if (printlevel >=2 ) {rec <- updateRec(rec=rec,x=xnew,f=eval$fofx,t=nbFun)}
  } 
  else if (linesearch_type=="armijo"){
    lsres <- BacktrackLineSearch(x,eval$fofx,eval$gradf,direction,f=fun)
    xnew <- lsres$xnew
    nbFun <- nbFun+lsres$nFcalls
    nbFunLS <- nbFunLS + lsres$nFcalls
    # re-evaluate for the gradient and do not increment nbFun as the point has already been calculated in linesearch
    eval <- f.gradf(x=xnew,f=fun,h=1.e-8)
    normGrad <- l2norm(eval$gradf)
    if (printlevel>=2) {
      rec<-lsres$rec
    }
  }
  else {
    stop("unknown algo_type")
  }
  # make the step
  iter <- iter+1
  previous_step <- step
  # previous_step <- xnew-x
  xprevious <- x
  x <- xnew

  # bookkeeping
  if (eval$fofx < Fbest) {
    Fbest <- eval$fofx
    recBest <- updateRec(rec=recBest,x=xnew,f=Fbest,t=nbFun)
  }
  if (printlevel >= 1){
    cat(" iteration : ",iter,"\r")
  } 
} # end while of main loop

if (printlevel >= 1){
  cat("Search exited after ",iter," iterations, ",nbFun," fct evaluations\n")
  lrecBest <- dim(recBest$X)[1]
  cat('best x:',recBest$X[lrecBest,],'\n')
  cat('best f:',recBest$F[lrecBest],'\n')
}

### Vizualization
plot(x = recBest$Time,y=log(1+recBest$F),type = "l",xlab = "nb. fct. eval",ylab="log(1+f)",col="red",xlim=c(1,nbFun))
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
  if (printlevel >=2){
    # png(filename="./contour.png") # save the contour in the current directory
    contour(x1, x2, z.grid, nlevels=20, xlab="x1", ylab="x2")
    # with search points on top
    points(rec$X[,1], rec$X[,2], pch=20, col="blue")
    text(rec$X, labels=1:iter, pos=3, cex=1.0) # (un)comment for labeling (or not) nb of calls to f when points created
    points(recBest$X[,1], recBest$X[,2], pch=19, col="red") # bests so far in red (should resemble iterates)
  }
  # dev.off()
}
