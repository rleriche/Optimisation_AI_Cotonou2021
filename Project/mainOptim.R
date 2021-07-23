# main program for optimization project
# 
# Rodolphe Le Riche, CNRS LIMOS, july 2021
# 
# For explanations, comments : cf. source files of the optimizers, 
# in particular gradient_descent.R
#
rm(list=ls()) # clear environment
source('test_functions.R')
source('gradient_descent.R')

### problem definition
# search space
pbFormulation <- list()
pbFormulation$fun<-rosen #function to minimize
d<-2
pbFormulation$d<-d # dimension
pbFormulation$LB<-rep(-5,d) #lower bounds
pbFormulation$UB<-rep(5,d) #upper bounds


### algorithm settings
optAlgoParam <- list()
optAlgoParam$xinit <- rep(-4.9,d)#c(4.5,3.5) rep(-4.9,d) # initial point
#
optAlgoParam$budget <- 3000
optAlgoParam$minGradNorm <- 1.e-6 
optAlgoParam$minStepSize <- 1.e-11 
#
optAlgoParam$direction_type <- "NAG" # choices are : "gradient", "momentum", "NAG"
optAlgoParam$linesearch_type <- "armijo" # choices are: "none", "armijo"
optAlgoParam$stepFactor <- 0.1 # step factor when there is no line search, 
optAlgoParam$beta <- 0.9 # momentum term for direction_type == "momentum" or "NAG" 
# 
printlevel <- 2 # controls how much is stored and printed, choices: 0,1,2

res<-gradient_descent(pbFormulation=pbFormulation,algoParam=optAlgoParam,printlevel=2)
  
  