#############################################
# a constant step-size gradient algorithm
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
xinit <- c(3,2) # initial point
stepFactor <- 0.2 # constant step factor
# the maximum number of calls to the function is larger because of the 
# finite differences scheme
printlevel <- 2 # controls how much is stored and printed
#                 =1 store best
#                 =2 store all points
stopBudget <- 10 # maximum number of iterations
stopGradNorm <- 1.e-8 # stop when mean gradient norm is smaller than 

### initializations
x <- xinit
eval <- f.gradf(x=x,f=fun,h=1.e-8) #eval$fofx is the function 
                          # eval$gradf the associated gradient
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


### run the algo
while ((iter <= stopBudget) & ((norm(as.matrix(eval$gradf),type = "F")/sqrt(d)) > stopGradNorm) ){
  
  xcandidate <- x - stepFactor*eval$gradf
  # project on bounds if necessary
  xnew <- ifelse(xcandidate < LB, LB, ifelse(xcandidate > UB, UB, xcandidate))
  # evaluate new point
  eval <- f.gradf(x=xnew,f=fun,h=1.e-8)
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
  } 

} # end while of main loop 


