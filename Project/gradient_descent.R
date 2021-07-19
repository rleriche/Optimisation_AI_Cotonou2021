#############################################
# a constant step-size gradient algorithm
# Rodolphe Le Riche
#
# Implementation notes: 
# * no attempt at making this code efficient, it is for teaching purpose
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
stepFactor <- 0.8 # constant step factor
# the maximum number of calls to the function is larger because of the 
# finite differences scheme
printlevel <- 2 # controls how much is stored and printed
#                 =1 store best
#                 =2 store all points
stopBudget <- 100 # maximum number of iterations
stopGradNorm <- 1.e-8 # stop when mean gradient norm is smaller than 

### initializations
x <- xinit
Fbest <- fun(x)
iter <- 1
recXbest <- matrix(x,nrow=1)
recFbest <- Fbest
recTime <- c(iter)
if (printlevel >= 2){
  recX <- matrix(x,nrow=1)
  recF <- Fbest
}
gradf <- ffinite.diff(x=x,f=fun,h=1.e-8)

### run the algo
while ((iter <= stopBudget) | ((norm(as.matrix(df),type = "F")/sqrt(d)) > stopGradNorm) ){
  
  xstep <- 
  
}


