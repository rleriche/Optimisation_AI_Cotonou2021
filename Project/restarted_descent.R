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
  total_budget <- algoParam$budget
  algoParam$budget <- total_budget/algoParam$nb_restarts
  res <- list() # all the results of the local searches as a list
  total_res <- list() # a list for cumulated and processed results 
  cumTime <- 0
  total_res$fbest<-.Machine$double.xmax
  #
  for (i in 1:algoParam$nb_restarts){

    # generate random initial point
    algoParam$xinit <- runif(n=pbFormulation$d,min=pbFormulation$LB,max=pbFormulation$UB)
    res[[i]]<-gradient_descent(pbFormulation=pbFormulation,algoParam=algoParam,printlevel=printlevel) 
    # process and cumulate the results
    cumTime <- cumTime+res[[i]]$nbFun
    if (res[[i]]$fbest<total_res$fbest){
      total_res$fbest<-res[[i]]$fbest
      total_res$xbest<-res[[i]]$xbest
    }
    if (printlevel>=1){
      res[[i]]$recBest$Time <- res[[i]]$recBest$Time + cumTime
      total_res$recBest <- updateRec(rec=total_res$recBest,x=res[[i]]$recBest$X,
                                     f = res[[i]]$recBest$F,t = res[[i]]$recBest$Time)
    }

  } # end loop on restarted local searches
  
  return(total_res)
}