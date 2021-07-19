#############################################
# Two gradient descent algorithms
# 1) Gradient descent where step size depends on gradient norm
# 2) Steepest descent where step size depends on a line search
#
# Rodolphe Le Riche
#
# Implementation notes: 
# * no attempt at making this code efficient, it is for teaching purpose. 
#############################################
source('test_functions.R')
source('utilities_optim.R')

### problem definition
# search space
d<-2 # dimension
LB<-c(-5,-5) #lower bounds
UB<-c(5,5) #upper bounds
fun<-quadratic #function to minimize

### algorithm settings
xinit <- c(-4,4.9) # initial point
algo_type <- "gradient" # choices are : "gradient" or "descent"
#   Both algorithms have minus the normalized gradient as search direction
#       x_{t+1} <- x_t + stepSize*direction
#   for "gradient" :  stepSize = stepFactor*normGrad
#   for "descent" : stepSize comes from a linesearch where 
#   sufficientDecreaseFactor is a control parameter
stepFactor <- 0.2 # step factor for "gradient" version
sufficientDecreaseFactor <- 0.1 # controls stepSize for "descent" version
#
printlevel <- 2 # controls how much is stored and printed
#                 =1 store best
#                 =2 store all points
stopBudget <- 10 # maximum number of iterations:
#   the maximum number of calls to the function is larger because of the 
#   finite differences scheme: nb_calls_to_fun = (d+1)*iter
stopGradNorm <- 1.e-6 # stop when gradient norm / sqrt(d) is smaller than 

### initializations
x <- xinit
eval <- f.gradf(x=x,f=fun,h=1.e-8) #eval$fofx is the function 
                          # eval$gradf the associated gradient
normGrad <- sqrt(sum(eval$gradf^2))
iter <- 1
# ...best are recordings of best so far
Fbest <- eval$fofx
recXbest <- matrix(x,nrow=1)
recFbest <- Fbest
recTime <- c(iter)
# below are recordings of the whole history. Memory consuming. 
# If memory issues, set printlevel<2
if (printlevel >= 2){
  recX <- matrix(x,nrow=1)
  recF <- eval$fofx
}

### line search function
# backtracking with Armijo (sufficient decrease) condition
# s must satisfy : f(x+s*d) <= f(x) + sufficientDecrease * s * d^T*gradf
BacktrackLineSearch <- function(x,fofx,gradf,direction,f,suffDecFact=0.1,decFact=0.5,initStepSize=1){
  stepSize <- initStepSize
  decConst <- suffDecFact*(direction%*%gradf)
  maxloop <- 100 # safety
  #
  xp <- x+stepSize*direction
  fp <- f(xp)
  nloop <- 1
  while ( (fp > fofx+stepSize*decConst) & (nloop<maxloop)) {
    stepSize <- stepSize*decFact
    xp <- x+stepSize*direction
    fp <- f(xp)
    nloop <- nloop+1
  }
  if (nloop >= maxloop){ 
    msg <- paste("nloop=",nloop," larger than maxloop=",maxloop)
    warning(msg)
  }
  return(stepSize)
}

### run the algo
if (printlevel >=2) {cat("Start gradient search\n")}
while ((iter <= stopBudget) & ((normGrad/sqrt(d)) > stopGradNorm) ){
  
  direction <- -eval$gradf/normGrad # search direction
  # no line search, step size proportional to gradient norm
  if (algo_type=="gradient"){
  stepSize <- stepFactor*normGrad
  } 
  else if (algo_type=="descent"){
    stop("not implemented yet")
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
  # make the step
  x <- xnew
  # bookkeeping
  if (eval$fofx < Fbest) {
    Fbest <- eval$fofx
    recXbest <- rbind(recXbest,xnew) 
    recFbest <- c(recFbest,Fbest)
    recTime <- c(recTime,iter)
  }
  if (printlevel >= 2){
    recX <- rbind(recX,xnew)
    recF <- c(recF,eval$fofx)
    cat(" iteration : ",iter,"\r")
  } 

} # end while of main loop 
if (printlevel >= 2){cat("gradient search exited after ",iter," iterations, start plotting\n")}

### Vizualization
plot(x = recTime,y=recFbest,type = "l",xlab = "nb. evaluations",ylab="f",col="red")
if (printlevel >= 2){
  lines(x = seq(1,iter),y = recF,col="blue")
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
    points(recX[,1], recX[,2], pch=20, col="blue")
    text(recX, labels=1:iter, pos=3, cex=1.0) # (un)comment for labeling (or not) nb of calls to f when points created
    points(recX[1,1], recX[1,2], pch=19, col="red") # initial point drawn in red
  }
  # dev.off()
}
