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

