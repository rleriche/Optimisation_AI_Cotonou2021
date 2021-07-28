#######################################
#### utilities for a single hidden layer neural net
#  Rodolphe Le Riche, July 2021
#  
#  really simple, only nice thing so far is the 
#  ability of the hidden layer to have different 
#  types of activation functions
#######################################

### a list of activation functions

# identity() is part of base R
# logistic a.k.a. sigmoid 
sigmoid <- function(x){1 / (1 + exp(-x))}
# tanh() is part of base R
# relu, rectified linear unit function
relu <- function(x){ifelse(x < 0 , 0, x )}
# lrelu, leaky relu
lrelu <- function(x){ifelse(x < 0 , 0.01 *x , x )}
# gauss
gauss <- function(x){exp(-(x^2))}

# random initialization of weights according to LeCun and Bottou
# E(weight)=0 , V(weight)=1/number_of_inputs
randomW <- function(nin,nout){
  return(matrix(data = rnorm(n=nin*nout,mean=0,sd=sqrt(1/nin)),nrow = nout))
}

# associate activation function names to the functions
make_Fun <- function(hnames){
  hfun<-list()
  nh <- length(hnames)
  for (i in 1:nh){
    if (hnames[[i]]=="sigmoid") {hfun[[i]]<-sigmoid}
    else if (hnames[[i]]=="tanh") {hfun[[i]]<-tanh}
    else if (hnames[[i]]=="identity") {hfun[[i]]<-identity}
    else if (hnames[[i]]=="relu") {hfun[[i]]<-relu}
    else if (hnames[[i]]=="lrelu") {hfun[[i]]<-lrelu}
    else if (hnames[[i]]=="gauss") {hfun[[i]]<-gauss}
    else stop(paste("unknow activation function name",hnames[[i]]))
  }
  return(hfun)
}

# normalization routines
normByRow <- function(X){
  nr <- dim(X)[1]
  nc <- dim(X)[2]
  Xnorm <- matrix(nrow = nr,ncol = nc)
  minAndMax <- matrix(nrow = nc,ncol=2)
  for (i in 1:nc){
    zmin<-min(X[,i])
    zmax<-max(X[,i])
    minAndMax[i,]<-c(zmin,zmax)
    Xnorm[,i]<-2*(X[,i]-zmin)/(zmax-zmin)-1
  }
  res<-list()
  res$dat <- Xnorm
  res$minAndMax <- minAndMax
  return(res)
}

unnormByRow <- function(normDat){
  nr <- dim(normDat$dat)[1]
  nc <- dim(normDat$dat)[2]
  X <- matrix(nrow = nr,ncol = nc)
  for (i in 1:nc){
    zmin<-normDat$minAndMax[i,1]
    zmax<-normDat$minAndMax[i,2]
    X[,i]<-(normDat$dat[,i]+1)/2*(zmax-zmin)+zmin
  }
  return(X)
}


# make a network
make_randomNN<- function(n_x,hnames,onames,seednb=1){
  NN<-list()
  NN$n_x<-n_x
  NN$n_h <- length(hnames)
  NN$n_y <- length(onames)
  NN$hfun <- make_Fun(hnames) # list of hidden layer activation functions
  NN$ofun <- make_Fun(onames) # list of output layer activation functions
  # the +1's below are for the biases
  set.seed(seednb) # reproductible weights
  NN$W1 <- randomW(nin = (NN$n_x+1),nout=NN$n_h)
  NN$W2 <- randomW(nin = (NN$n_h+1),nout=NN$n_y)
  set.seed(Sys.time()) # unset seed, back to "random"  
  return(NN)
}

# forward propagation
NNpred<-function(NN,X){
  # It is customary in computer experiments to have X which is ndata*n_x , 
  # Y ndata*n_y . However, to work with a neural net, it is easier to 
  # transpose because it allows calculations through matrix products.
  # So, transpose X, and transpose the result.
  ndat <- dim(X)[1]
  if (dim(X)[2]!=NN$n_x) stop("incompatible dimensions between NN and data")
  XP <- cbind(X,rep(1,ndat)) # augmentation to account for biases
  Z1 <- NN$W1%*%t(XP)
  A1 <- matrix(nrow = NN$n_h,ncol = ndat)
  for (i in 1:NN$n_h){
    A1[i,]<-sapply(X = as.matrix(Z1[i,]),FUN = NN$hfun[[i]])
  }
  A1P <- rbind(A1,rep(1,ndat))
  Z2 <- NN$W2%*%A1P
  # do the transposition at the same time as the passing through output layer
  Y <- matrix(nrow = ndat,ncol = NN$n_y)
  for (i in 1:NN$n_y){
    Y[,i]<-sapply(X = as.matrix(Z2[i,]),FUN = NN$ofun[[i]])
  }
  return(Y)
}

# calculate root mean square error
# X is ndata*n_x  (n_x dimension of the input to the NN)
# Y is ndata*n_y
NNrmse<-function(NN,X,Ytarget){
  ndat<-dim(X)[1]
  Y<-NNpred(NN,X)
  rmse <- sqrt((sum((Ytarget-Y)^2))/ndat)
  return(rmse)
}

# transform x, the optimization variables, i.e., weights of the NN,
# into a NN. 
# Input :
#   NN : a model NN that will be overwritten
xtoNN <- function(x,NN){
  NNx <- NN
  # assign x to the NN
  # W line i is the connexion to neuron i of the next layer,
  #  the last column is the bias
  # pick 2 variables (useful for plots)
  # NNx$W1[1,1]<-x[1]
  # NNx$W2[1,5]<-x[2]
  # optimize the entire network
  NNx$W1<-matrix(data = x[1:12],nrow=4,byrow = T)
  NNx$W2<-matrix(data = x[13:17],nrow=1,byrow = T)
  return(NNx)
}

