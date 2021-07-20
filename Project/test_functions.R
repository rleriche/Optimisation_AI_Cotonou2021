## set of test functions

# Ackley function. Global minimum at glob_xstar
ackley <- function(xx, a=20, b=0.2, c=2*pi){
  glob_xstar <- rep(0,length(xx))
  xx <- xx - glob_xstar
  aa <- 6.4*xx
  d <- length(aa)
  sum <- -a*exp(-b*sqrt((1/d)*sum(aa*aa))) - exp((1/d)*sum(cos(c*aa)))+ a + exp(1)
  y <- sum 
  return(y)
}

## Rastrigin function. Global minimum at glob_xstar
rastrigin <- function(xx){
  glob_xstar <- rep(0,length(xx))
  scaling_fact <- 1.024
  xx <- xx - glob_xstar
  aa <- scaling_fact*xx
  d <- length(aa)
  sum <- sum(aa*aa - 10*cos(2*pi*aa))
  y <- sum + 10*d
  return(y)
}

## Schwefel function.
schwefel <- function(xx){
  glob_xstar <- rep(0,length(xx))
  xx <- xx - glob_xstar + 1
  aa <- 100*xx
  d <- length(aa)
  sum <- 418.9829*d - sum(aa*sin(sqrt(abs(aa))))
  y <- sum
  return(y)
}

## Sphere function. Global Minimum at glob_xstar
sphere <- function(xx){
  glob_xstar <- rep(0,length(xx))
  xx <- xx - glob_xstar
  #aa <- 1.024*xx
  aa <- xx
  d <- length(aa)
  sum <- sum(aa*aa)
  y <- sum
  return(y)
}

##### michalewicz function #########
#glob min in 2D: -1.83, 5D: -4.71 , 10D: -9.64
michalewicz <- function(xx, m=10){
  aa <- 0.1*pi*(xx + 5)
#   aa <- xx
  d <- length(aa)
  i <- 1:d
  sum <- sum(sin(aa)*sin((i*aa*aa)/pi)^(2*m))
  y <- -sum
  return(y + 2)
}

##### quadratic function #########"
quadratic <- function(xx){
  glob_xstar <- rep(0,length(xx))
  cond.no <- 5
  aa <- 1.024*xx
  d <- length(aa)
  xstar <- glob_xstar
  lambdas <- diag(seq(1, cond.no, ,d))
  # matrix with arbitrary orientation. The seed number decides the orientation.
  if (!exists("glob_umat")) {
    set.seed(1) # change this seed to change the orientation of the quadratic function
    glob_umat <<- qr.Q(qr(matrix(runif(d*d),nrow=d,ncol=d)))
    glob_umat <<- glob_umat[,sample(seq(1,d))]
    #  glob_umat <<- diag(1, d) # to generate a quadratic function whose principal axes are aligned with coordinates
  }
  H <- glob_umat%*%lambdas%*%t(glob_umat)  
  y <- 0.5*t(aa - xstar)%*%H%*%(aa - xstar)
  return(y)
}

#### Tunnel function #######"
tunnel <- function(xx, b=0.5){
  glob_xstar <- rep(0,length(xx))
  d <- length(xx)
  xx <- xx - glob_xstar
  in_tunnel <- TRUE
  if (d > 1) {
    for(i in 2:d){
      if ((xx[i] < -b) | (xx[i] > b)){
        in_tunnel <- FALSE
      }
    }
    if (in_tunnel) {
      y <- -exp(-norm(as.matrix(xx), type="F")/10)
    } else{
      y <- t(xx)%*%xx
    }
  } else{
    y <- xx^2
  }
  return(y)
}

##### function with hierarchical sensitivities and local optima #########
quad_wave <- function(xx){
  cond.no <- 1.e1
  d <- length(xx)
  lambdas <- diag(seq(1, cond.no, ,d))
  glob_umat <<- diag(1, d) # to generate a quadratic function whose 
      # principal axes are aligned with coordinates. Could be made into a rotating matrix.
  H <- glob_umat%*%lambdas%*%t(glob_umat)  
  y <- 1+0.5*t(xx)%*%H%*%(xx)-0.2*sum(cos(4*pi*xx[c(1,2)]))
  return(y)
}


##### function with hierarchical sensitivities and local optima #########
L1norm <- function(xx){
     y <- sum(abs(xx))
  return(y)
}

rosen <- function(xx)
{
  ##########################################################################
  #
  # ROSENBROCK FUNCTION
  # global optimum at (1,...,1)
  # copied and slightly changed from
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  ##########################################################################
  d <- length(xx)
  xi <- xx[1:(d-1)]
  xnext <- xx[2:d]
  y <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
  return(y)
}
