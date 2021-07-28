###############################################
### main program of experiments with neural net
#   Rodolphe Le Riche, July 2021
###############################################
source("neural_net.R")

### get or generate data
source('test_functions.R')
fun<-quadratic
d<-2
LB <- rep(-5,d)
UB <- rep(5,d)
ntrain <- 50
ntest <- 50
ndata <- ntrain + ntest
set.seed(1) # unset seed to always have same data
rawX <-t(replicate(n = ndata,expr = runif(n = d,min = LB,max = UB)))
set.seed(Sys.time()) # unset seed, back to "random"
rawYtrue <- apply(X = rawX,MARGIN = 1,FUN = fun)
# normalize the data between -1 and 1
X <- normByRow(rawX)
# you can recover unnormalized data with, for expl : X <- unnormByRow(normIn)
Ytrue <- normByRow(as.matrix(rawYtrue))
itrain <- 1:ntrain
itest <- (ntrain+1):ndata
Xtrain<-X$dat[itrain,]
Xtest<-X$dat[itest,]
Ytrain<-matrix(data=Ytrue$dat[itrain,],ncol=1)
Ytest<-matrix(data=Ytrue$dat[itest,],ncol=1)
### end data collection

###############################################
### demonstrate NN utilities
###############################################

###  describe the network
n_x <- d
# hidden layer activation function names
hnames <- c("sigmoid","sigmoid","lrelu","lrelu")
# output layer activation function names
onames <- c("identity")

# random network generation
NN<-make_randomNN(n_x,hnames,onames,seednb = 2)

# evaluate NN predictions
# Ytest<-NNpred(NN,Xtest)

# evaluate NN MSE w.r.t. certain data
MSEtrain<-NNmse(NN = NN,X = Xtrain,Ytarget = Ytrain)
MSEtest<-NNmse(NN = NN,X = Xtest,Ytarget = Ytest)

###############################################
### optimize the NN
###############################################

# objective function wrapper, a bit ugly
# !!! NN , Xtrain, Ytrain passed as global variables
# !!! x now means part of the weights of the NN
fmse <- function(x){
  NNx <- xtoNN(x,NN)
  mse<-NNmse(NN = NNx,X = Xtrain,Ytarget = Ytrain)
  return(mse)
}

source('utilities_optim.R')
source('line_searches.R')
source('gradient_descent.R')
source('restarted_descent.R')

### problem definition
pbFormulation <- list()
pbFormulation$fun<-fmse #function to minimize
d<-17
pbFormulation$d<-d # dimension
pbFormulation$LB<-rep(-5,d) #lower bounds
pbFormulation$UB<-rep(5,d) #upper bounds


### algorithm settings
optAlgoParam <- list()
# optAlgoParam$xinit <- runif(n = d,min = -1,max = 1) # initial point
optAlgoParam$xinit <-c(as.vector(t(randomW(3,4))),as.vector(t(randomW(5,1))))
#
optAlgoParam$budget <- 1000
optAlgoParam$minGradNorm <- 1.e-6 
optAlgoParam$minStepSize <- 1.e-11 
#
optAlgoParam$direction_type <- "momentum" # choices are : "gradient", "momentum", "NAG"
optAlgoParam$linesearch_type <- "armijo" # choices are: "none", "armijo"
optAlgoParam$beta <- 0.9 # momentum term for direction_type == "momentum" or "NAG" 
# 
printlevel <- 4 # controls how much is stored and printed, choices: 0 to 4
# a single descent
res<-gradient_descent(pbFormulation=pbFormulation,algoParam=optAlgoParam,printlevel=printlevel)

# get the NN solution
NNopt<-xtoNN(x = res$xbest,NN = NN)
# see how it does in terms of test error
# 
Ynntest<-NNpred(NNopt,Xtest)
Ynntrain<-NNpred(NNopt,Xtrain)
msetrain<-NNmse(NN = NNopt,X = Xtrain,Ytarget = Ytrain)
msetest<-NNmse(NN = NNopt,X = Xtest,Ytarget = Ytest)
plot(Ytrain,Ynntrain)
lines(c(-1,1),c(-1,1))
title(paste("train RMSE=",toString(msetest)))
plot(Ytest,Ynntest)
lines(c(-1,1),c(-1,1))
title(paste("test RMSE=",toString(msetest)))
