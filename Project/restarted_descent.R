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
  if (printlevel >= 3) locPrintLev <-3
  else if (printlevel >= 1) locPrintLev <- 1
  else locPrintLev <- 0
  # 
  total_budget <- algoParam$budget
  algoParam$budget <- total_budget/algoParam$nb_restarts
  res <- list() # all the results of the local searches as a list
  total_res <- list() # a list for cumulated and processed results 
  cumTime <- 0
  total_res$fbest<-.Machine$double.xmax
  #
  if (!silent) cat("**** START RESTARTED DESCENT\n")
  for (i in 1:algoParam$nb_restarts){

    # generate random initial point
    algoParam$xinit <- runif(n=pbFormulation$d,min=pbFormulation$LB,max=pbFormulation$UB)
    # proceed with local search
    res<-gradient_descent(pbFormulation=pbFormulation,algoParam=algoParam,printlevel=locPrintLev) 
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
    
    if (!silent) {
      cat("**** search no.",i," done:\n")
      cat("     started at x = ",algoParam$xinit,"\n")
      cat("     converged to x = ",res$xbest,"\n")
      cat("     where f = ",res$fbest,"\n")
    }
  } # end loop on restarted local searches
  if (!silent) {
    cat("**** END RESTARTED DESCENT, overall best:\n")
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
    
    if (d==2) { 
      if (printlevel >=2){
        # png(filename="./contour.png") # save the contour in the current directory
        plot_contour(LB=pbFormulation$LB,UB=pbFormulation$UB,f=pbFormulation$fun)
        if (printlevel >=4){
          # all search points on top
          points(total_res$rec$X[,1], total_res$rec$X[,2], pch=20, col="blue", cex=0.5)        
        }
        # put a selection of best-so-far on top in red (should resemble iterates)
        points(total_res$rBest$X[,1], total_res$rBest$X[,2], pch=19, col="red") 
        iplot <- inc.geom.seq(from=1,to=length(total_res$rBest$Time),coef=1.4) # selection, otherwise unreadable
        # text(x=total_res$rBest$X[iplot,1], y=total_res$rBest$X[iplot,2],labels=total_res$rBest$Time[iplot], pos=3, cex=1.0) # (un)comment for labeling (or not) nb of calls to f when points created
        # dev.off()
      }
      
    } # end if d==2
  }
  #
  return(total_res)
}