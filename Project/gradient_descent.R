#############################################
# A family of gradient-based descent algorithms
# 
# Rodolphe Le Riche, CNRS LIMOS, july 2021
# 
#
# The algorithms differ in how the search direction is calculated
#   direction_type = (gradient, momentum, NAG) , where NAG= Nesterov Accelerated Gradient
# and they differ on whether a single step is taken or a line search performed
#   linesearch_type = (none, armijo)
#
#   All algorithms have a search direction and a step size
#       x_{t+1} <- x_t + stepSize*direction
#   where direction is a vector of norm == 1 and stepSize a scalar
#   The "direction" is the normalized "step", direction = step/l2norm(step)
#   "step" : vector which depends on "direction_type" (gradient, momentum or NAG)
#     Ex: if direction_type == gradient then step = -gradf
#     and other formula for other direction_type. 
#   Calculation of stepSize. It depends on "linesearch_type":
#   When linesearch_type == "armijo" : stepSize comes from a line search, cf line_searches.R file
#   When linesearch_type == "none" (no line search) : stepSize = stepFactor*l2norm(step)
#   Expl for "gradient" & "none" :  stepSize = stepFactor*normGrad
#
# INPUTS
#
# pbFormulation$... fields related to problem formulation
#   pbFormulation$fun : function to optimize
#   pbFormulation$d : problem dimension
#   pbFormulation$LB : d dimensional vector of lower bounds for the variables
#   pbFormulation$UB : d dimensional vector of upper bounds for the variables
# 
# algoParam$...  optimization algorithm settings
#   algoParam$xinit : d dimensional vector for the starting point. No default.
#   algoParam$direction_type : character string designating the algorithm that 
#     decides on the direction. choices are : "gradient", "momentum", "NAG". 
#     Default is "NAG". 
#   algoParam$linesearch_type : character string for the line search. 
#     Choices are "none", "armijo". Default is "armijo".
#   algoParam$stepFactor : step factor when there is no line search, 
#     stepfactor ~ 1/Lipschitz_constant : the steeper the function, 
#     the smaller stepfactor. Default is 0.1
#   algoParam$beta : momentum term for direction_type == "momentum" or "NAG" 
#     here we keep it constant. Pradeep Ravikumar and Aarti Singh recommend
#     beta <- (iter-2)/(iter+1). Here default is 0.9 
#   algoParam$budget : maximum number of function evaluations. Default is 1000.
#    The current true number of calls to the function is larger because of the 
#    finite differences scheme: nbFun = (d+1)*iter+nbFunLS.
#   algoParam$minGradNorm : stop when gradient norm / sqrt(d) is smaller than it. 
#     Default is 1.e-6
#   algoParam$minStepSize : stop if stepSize < stopStepSize. Default is 1.e-11 
#
# printlevel : controls how much is stored and printed
#                 =0 store overall best only, no plot
#                 =1 store best history, silent (no plot, no output)
#                 =2 store best history and plot it, some messages
#                 =3 store best history and all points, silent, no plot
#                 =4 store all points, best history, plot it + extra plots in particular when d==2
# 
# OUTPUTS
# 
# All results are fields of the res list: 
#   res$xbest : best point found 
#   res$fbest : objectif function of best point
#   res$nbFun : total number of calls to the objective function
#   res$stopCode : code for stopping criterion
#       = 1 for budget exhausted
#       = 2 too small gradient norm
#       = 3 too small step
#   res$rBest : history of best so far points
#   res$rec : history of all computed points
#   Both res$rec and res$rBest are lists with fields $X for the points, 
#     $F for the associated function, $Time for the number of fun calls 
#     when it was discovered.
# 
# Implementation notes: 
# * nbFun, the variable that counts the number of calls to the objective function
#   counts 1 for a finite difference, which is not correct (it should be d+1), 
#   but it would be correct for analytical gradients
# * Delete the global variables if you change dimension and optimize quadratic function. 
#   This is because the Hessian of the quadratic function depends on
#   the global variable "glob_umat" and is not recalculated if glob_mat exists.
#############################################

gradient_descent <- function(pbFormulation,algoParam,printlevel=1){
  ### process input parameters
  # pb formulation
  if (is.null(pbFormulation$fun)) stop("need a function to optimize, pbFormulation$fun")
  else {fun<-pbFormulation$fun}
  if (is.null(pbFormulation$d)) stop("need a problem dimension, pbFormulation$d")
  else {d<-pbFormulation$d}
  if (is.null(pbFormulation$LB)) stop("need variables lower bounds, pbFormulation$LB")
  else if (length(pbFormulation$LB) != d) stop("length LB not equal dimension")
  else {LB<-pbFormulation$LB}
  if (is.null(pbFormulation$UB)) stop("need variables upper bounds, pbFormulation$UB")
  else if (length(pbFormulation$UB) != d) stop("length UB not equal dimension")
  else {UB<-pbFormulation$UB}
  # algorithm settings
  if (is.null(algoParam$xinit)) stop("need to provide an initial point algoParam$xinit")
  else if (length(algoParam$xinit) != d) stop("length xinit not equal to dimension")
  if (is.null(algoParam$direction_type)) algoParam$direction_type <- "NAG" 
  else direction_type <- algoParam$direction_type
  if (is.null(algoParam$linesearch_type)) algoParam$linesearch_type <- "armijo"
  else linesearch_type <- algoParam$linesearch_type
  if (is.null(algoParam$stepFactor)) stepFactor <- 0.1
  else stepFactor<-algoParam$stepFactor 
  if (is.null(algoParam$beta)) beta <- 0.9 
  else beta<- algoParam$beta  
  if (is.null(algoParam$budget)) budget <- 1000 
  else budget <- algoParam$budget 
  if (is.null(algoParam$minGradNorm)) minGradNorm <- 1.e-6 
  else minGradNorm <- algoParam$minGradNorm 
  if (is.null(algoParam$minStepSize)) minStepSize <- 1.e-11
  else minStepSize <- algoParam$minStepSize
  if ((printlevel >=2) & (printlevel != 3)) silent<-FALSE else silent<-TRUE
    
  ### initializations
  x <- algoParam$xinit
  eval <- f.gradf(x=x,f=fun,h=1.e-8) #eval$fofx is the function value
  # eval$gradf the associated gradient
  normGrad <- l2norm(eval$gradf)
  iter <- 1
  nbFun <- 1 # counts the number of calls to the objective function. 
  #       nbFun includes the calls done during line search and 1 call to 
  #       the gradient counts for 1 call to the objective function
  nbFunLS <- 0 # nb of calls to fun done during the line search
  # stopping criteria:
  stopGradNorm<-FALSE
  stopBudget<-FALSE
  stopStepSize<-FALSE
  # recordings of best so far
  fbest <- eval$fofx
  xbest <- x
  if (printlevel >= 1){
    rBest <- list()
    rBest$X <- matrix(x,nrow=1)
    rBest$F <- fbest
    rBest$Time <- c(nbFun)    
  }
  # below are recordings of the whole history. Memory consuming. 
  # If memory issues, set printlevel<2
  if (printlevel >= 3){
    rec <- list()
    rec$X <- matrix(x,nrow=1)
    rec$F <- eval$fofx
    rec$Time <- nbFun
  }
  
  
  ### start the search
  if (!silent) {
    startmsg <- paste0("Start ",direction_type," (direction) + ",linesearch_type," (line search) descent\n")
    cat(startmsg)
  }
  while (isFALSE(stopBudget) & isFALSE(stopGradNorm) & isFALSE(stopStepSize) ){
    
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
      nbFun <- nbFun+2 # +2 , +1 for the obj function, +1 for the gradient
      if (printlevel >=3 ) {rec <- updateRec(arec=rec,x=xnew,f=eval$fofx,t=nbFun)}
    } 
    else if (linesearch_type=="armijo"){
      lsres <- BacktrackLineSearch(x=x,fofx=eval$fofx,gradf=eval$gradf,
                                   direction=direction,f=fun,LB=LB,UB=UB,
                                   printlevel=printlevel,nbFun=nbFun)
      xnew <- lsres$xnew
      nbFun <- nbFun+lsres$nFcalls
      nbFunLS <- nbFunLS + lsres$nFcalls
      # re-evaluate for the gradient and increment nbFun of 1 only 
      # as the point has already been calculated in linesearch
      eval <- f.gradf(x=xnew,f=fun,h=1.e-8)
      nbFun <- nbFun+1
      normGrad <- l2norm(eval$gradf)
      if (printlevel>=3) {
        rec <- updateRec(arec=rec,x=lsres$rec$X,f=lsres$rec$F,t=lsres$rec$Time)
      }
    }
    else {
      stop("unknown algo_type")
    }
    # make the step
    iter <- iter+1
    previous_step <- step
    xprevious <- x
    x <- xnew
    # stopping conditions
    if (nbFun >= budget){stopBudget<-TRUE}
    if ((normGrad/sqrt(d)) <= minGradNorm){stopGradNorm<-TRUE}
    if (l2norm(x-xprevious) <= minStepSize){stopStepSize<-TRUE}
    
    # bookkeeping
    if (eval$fofx < fbest) {
      fbest <- eval$fofx
      xbest <- x
      if (printlevel >=1) {
        rBest <- updateRec(arec=rBest,x=xnew,f=fbest,t=nbFun)
      }
    }
    if (!silent){
      cat(" iteration : ",iter,"\r")
    } 
  } # end while of main loop
  
  # gather final results
  res <- list()
  # lrBest <- dim(rBest$X)[1]
  res$xbest <-xbest 
  res$fbest <-fbest
  res$nbFun <- nbFun
  if (stopBudget) res$stopCode<-1
  if (stopGradNorm) res$stopCode<-2
  if (stopStepSize) res$stopCode<-3
  if (printlevel >=1) res$rBest <- rBest
  if (printlevel >= 3) {res$rec <- rec} 
  #
  if (!silent){
    cat("Search exited after ",iter," iterations, ",nbFun," fct evaluations\n")
    stopMsg <- ""
    if (stopBudget) {stopMsg <- paste(stopMsg,"max budget reached") }
    if (stopGradNorm) {stopMsg <-paste(stopMsg,"small gradient norm") }
    if (stopStepSize) {stopMsg <-paste(stopMsg,"small step size") }
    cat("Active stopping criteria :",stopMsg,"\n")
    cat('best x:',res$xbest,'\n')
    cat('best f:',res$fbest,'\n')
  }
  
  ### Visualization
  if (!silent) {
    plot(x = rBest$Time,y=log(1+rBest$F),type = "l",xlab = "nb. fct. eval",ylab="log(1+f)",col="red",xlim=c(1,nbFun))
    if (printlevel >= 4){
      lines(x = rec$Time,y = log(1+rec$F),col="blue")
    }
    if (d==2) { 
      if (printlevel >=2){
        # png(filename="./contour.png") # save the contour in the current directory
        plot_contour(LB=LB,UB=UB,f=fun)
        if (printlevel >=4){
          # all search points on top
          points(rec$X[,1], rec$X[,2], pch=20, col="blue", cex=0.5)        
        }
        # put a selection of best-so-far on top in red (should resemble iterates)
        points(rBest$X[,1], rBest$X[,2], pch=19, col="red") 
        iplot <- inc.geom.seq(from=1,to=length(rBest$Time),coef=1.4) # selection, otherwise unreadable
        text(x=rBest$X[iplot,1], y=rBest$X[iplot,2],labels=rBest$Time[iplot], pos=3, cex=1.0) # (un)comment for labeling (or not) nb of calls to f when points created
        # dev.off()
      }

    } # end if d==2
  } # end if !silent
  
  return(res)
}

