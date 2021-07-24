#############################################
# Very simple restart mechanism for gradient based searches
#   The total budget is equally split between the local searches
#   The initial point is randomly chosen
# 
# INPUTS
#   the same as for the local search (see gradient_descent.R) +
#   algoParam$nb_restarts : number of restarted local searches
# 
# OUTPUTS
#   
restarted_descent <- function(pbFormulation,algoParam,printlevel=1) {
  if (is.null(algoParam$nb_restarts)) {algoParam$nb_restarts<-5}
  if ((printlevel >=2) & (printlevel != 3)) silent<-FALSE else silent<-TRUE
  # 
  total_budget <- algoParam$budget
  algoParam$budget <- total_budget/algoParam$nb_restarts
  res <- list() # all the results of the local searches as a list
  total_res <- list() # a list for cumulated and processed results 
  cumTime <- 0
  total_res$fbest<-.Machine$double.xmax
  #
  if (!silent) cat("**** Start restarted descent ****\n")
  for (i in 1:algoParam$nb_restarts){

    # generate random initial point
    algoParam$xinit <- runif(n=pbFormulation$d,min=pbFormulation$LB,max=pbFormulation$UB)
    # proceed with local search
    res<-gradient_descent(pbFormulation=pbFormulation,algoParam=algoParam,printlevel=printlevel) 
    # process and cumulate the results
    if (res$fbest<total_res$fbest){
      total_res$fbest<-res$fbest
      total_res$xbest<-res$xbest
    }
    if (printlevel>=1){
      res$rBest$Time <- res$rBest$Time + cumTime
      total_res$rBest <- updateRec(arec=total_res$rBest,x=res$rBest$X,
                                     f = res$rBest$F,t = res$rBest$Time)
    }
    if (printlevel>=4){
      res$rec$Time <- res$rec$Time + cumTime
      total_res$rec <- updateRec(arec=total_res$rec,x=res$rec$X,
                                     f = res$rec$F,t = res$rec$Time)
    }
    cumTime <- cumTime+res$nbFun
    
    if (!silent) cat("**** restarted search no.",i," ended\n")
  } # end loop on restarted local searches
  if (!silent) {
    cat("**** End restarted descent, reporting:\n")
    cat('**** best x:',total_res$xbest,'\n')
    cat('**** best f:',total_res$fbest,'\n')
  }
  # visualization
  if (!silent) {
    plot(x = total_res$rBest$Time,y=log(1+total_res$rBest$F),type = "l",
         xlab = "nb. fct. eval",ylab="log(1+f)",col="red",xlim=c(1,cumTime))
    if (printlevel >= 4){
      lines(x = total_res$rec$Time,y = log(1+total_res$rec$F),col="blue")
    }
  }
  #
  return(total_res)
}