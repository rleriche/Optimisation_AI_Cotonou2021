# main program for optimization project
# 
# Rodolphe Le Riche, CNRS LIMOS, july 2021
# 
# For explanations, comments : cf. source files of the optimizers, 
# in particular gradient_descent.R
#
# rm(list=ls()) # clear environment
source('test_functions.R')
source('utilities_optim.R')
source('line_searches.R')
source('gradient_descent.R')
source('restarted_descent.R')

### problem definition
# search space
pbFormulation <- list()
pbFormulation$fun<-rosen #function to minimize, examples in test_functions.R
d<-2
pbFormulation$d<-d # dimension
pbFormulation$LB<-rep(-5,d) #lower bounds
pbFormulation$UB<-rep(5,d) #upper bounds


### algorithm settings
optAlgoParam <- list()
optAlgoParam$xinit <- runif(n = d,min = pbFormulation$LB,max = pbFormulation$UB) #c(4.5,3.5) rep(-4.9,d) # initial point
#
optAlgoParam$budget <- 4000
optAlgoParam$minGradNorm <- 1.e-6 
optAlgoParam$minStepSize <- 1.e-11 
#
optAlgoParam$direction_type <- "momentum" # choices are : "gradient", "momentum", "NAG"
optAlgoParam$linesearch_type <- "armijo" # choices are: "none", "armijo"
optAlgoParam$stepFactor <- 0.1 # step factor when there is no line search, 
optAlgoParam$beta <- 0.9 # momentum term for direction_type == "momentum" or "NAG" 
# 
printlevel <- 4 # controls how much is stored and printed, choices: 0 to 4

# a single descent
res<-gradient_descent(pbFormulation=pbFormulation,algoParam=optAlgoParam,printlevel=printlevel)

# a restarted descent
# optAlgoParam$nb_restarts <- 4
# cres <- restarted_descent(pbFormulation=pbFormulation,algoParam = optAlgoParam,printlevel=printlevel)
  