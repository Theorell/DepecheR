#' Find the smallest optimal sampe size and the best penalty offset for a certain dataset
#'
#'
#' This function is used before the dClust function to identify the smallest sample size that gives rise to the most stable clustering. It is primarily used when datasets are very large (i.e. >100 000 observations), where there is a computational gain in clustering based on a subset of the cells and then assigning all other cells to the correct cluster center. It can however be generally useful when it is not entirely clear which sample size that is needed to perform the optimization of the penalty term, something that this function also works as a wrapper for.
#' @importFrom graphics box
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
#' @param sampleSize The number of observations that are included in each bootstrap subsampling of the data. NB! The algorithm uses resampling, so the same event can be used twice. This is the central argument to this function, that it optimizes over.
#' @param initCenters Number of starting points for clusters. The higher the number, the greater the precision of the clustering, but the computing time is also increased with the number of starting points. Default is 30.
#' @param maxIter The maximal number of iterations that are performed to reach the minimal improvement. 
#' @param minImprovement This is connected to the evaluation of the performance of the algorithm. The more iterations that are run, or the larger the samples are, the smaller will the improvement be, and this sets the threshold when the iterations stop. 
#' @param penaltyOptOnly If this is set to true, no optimization on the sample size is performed. This might be useful in situations when the whole dataset will be used for clustering, such as when the dataset is small, or when the number of variables is very large (>100), so that increasing the sample size too far might lead to computational overload. 
#' @param prePenaltyOptSampleSize The size of the initial samle that will be used to define the penalty offset that will be used in all subsequent optimizations of the sample size. 
#' @param penaltyOffset These values are the ones that are evaluated and the ones that decide the penalization. The number of suggested default values are empirically defined and might not be optimal for a specific dataset, but the algorithm will warn if the most optimal values are on the borders of the range. Note that when this offset is 0, there is no penalization, which means that the algorithm runs normal K-means clustering.
#' @seealso \code{\link{dClust}}, \code{\link{dClustPredict}}, \code{\link{dOptPenalty}}
#' @return If penaltyOptOnly is set to FALSE, which is default, one graph showing the stability of the models with different sample sizes, one showing the distance with different penalties before the optimization of the sample size, and a third showing the same, but after the sample size has been optimized. In addition, a list is returned with the following content:
#' \describe{
#'     \item{sampleSizeOpt}{A dataframe, in which each row represents one sample size, and in which the last row is thus the chosen, optimal sample size. It has the following columns:
#'     \describe{
#'               \item{SampleSize}{This column shows the sample size of each boot strap subsampling in the optimization procedure.}
#'               \item{Lowest distance}{This vector shows the optimal stability, expressed as the lowest distance between the bootstrap subsampling runs at each of the boot strap subsamling sizes.}
#'               \item{Improvement}{Here, the improvement, expressed as a fraction between 0 and 1 is shown. When the improvement is less than minImprovement, the algorithm automatically stops.}
#'              } 
#'     }
#'     \item{penaltyOpt.df}{A dataframe with one row with all the information about which settings that were used to generate the optimal clustering with the optimal sample size. The "withOrigoClust" information tells the user if the solution with or without a cluster in origo gives the most optimal solution. If yes, this origo population is generally small and could be viewed as not fitting in the model.}
#'     \item{meanOptimDf}{A dataframe with the information about the results with all tested penalty values}
#' }
#' @return If instead penaltyOptOnly is set to TRUE, one graph showing the distance with different penalties as well as the two last items in the list described above is returned. 
#' @examples
#' #Generate a dataframe with bimodally distributed data with a million rows.
#' x <- generateBimodalData(observations=1000000)
#'
#' #Scale this datamframe
#' x_scaled <- dScale(x[,2:ncol(x)])
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the function to identify at what sample size the cluster stability plateaus
#' x_optim <- dOpt(x_scaled)
#' @export dOpt
dOpt <- function(inDataFrameScaled, sampleSize=1000*c(2^0, 2^1, 2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10), initCenters=30, maxIter=100, minImprovement=0.01, penaltyOptOnly=FALSE, prePenaltyOptSampleSize=8000, penaltyOffset=c(0,2,4,8,16,32,64,128)){

  #First, the optimal penaltyOffset is identified with a reasonable sample size
  prePenaltyOpt <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, maxIter=maxIter, minImprovement=minImprovement, bootstrapObservations=prePenaltyOptSampleSize, penaltyOffset=penaltyOffset, makeGraph=TRUE, graphName="Pre-optimization of penalty term.pdf")
  if(penaltyOptOnly==TRUE){
    return(prePenaltyOpt)
  } else {

    bestPenaltyOffsetPos <- which(penaltyOffset==prePenaltyOpt[[1]][1,1])
    print("Now, the pre-optimization of the penalty terms is done and the sample size optimization will start.")
  
	  lowestDist <- vector()
	  for(i in 1:length(sampleSize)){
	    ptm <- proc.time()
	    dOptPenaltyResult <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, sampleSize[i], penaltyOffset=c(penaltyOffset[bestPenaltyOffsetPos-2], penaltyOffset[bestPenaltyOffsetPos-1], penaltyOffset[bestPenaltyOffsetPos], penaltyOffset[bestPenaltyOffsetPos+1], penaltyOffset[bestPenaltyOffsetPos+2]), makeGraph=FALSE, disableWarnings=TRUE)
	    #graphName=paste("Distance over penalty values for sample size ", sampleSize[i], ".pdf", sep="")
	    timing <- proc.time() - ptm
	    lowestDist[i] <- min(dOptPenaltyResult[[2]][,1:2])
		
	    if(i<=3){
	      print(paste("Cycle", i, "completed. Jumping to next sample size."))
	    } else {
	      if(abs(lowestDist[i-1]-lowestDist[i])<=minImprovement){
	        print(paste("Cycle", i, "optimal.", "Done."))
	        break
	      } else if(timing[3]>30){
	        ANSWER <- readline(paste("The improvement between the two latest cycles was ", lowestDist[i], " . ", "The next cycle will take approximately ",  timing[3]*2, " seconds. Do you wish to run it? Print no if you do not, and yes if you do.", sep=""))

	        if (substr(ANSWER, 1, 1) == "n"){
	          cat("Ok. Lets stop here then.")
	          break
	        }
			      
			    else {
			      cat("Ok.")
			    }
			      
  			} else{
				print(paste("Cycle ", i, " completed. Jumping to next sample size.", sep=""))
  			}
	
      }
		
  	}

  	pdf("Distance as a function of sample sizes.pdf")
    par(mar=c(5, 4, 4, 6) + 0.1)
    #Plot the lowest distances as a function of the sample size.
    plot(sampleSize[1:length(lowestDist)], lowestDist, pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="",
         type="b",col="black", main="Distance between bootstraps as a as a function of sample size")
    axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("Distance between bootstraps",side=2,line=2.5)
    box()
    
    ## Draw the penalty axis
    axis(1,pretty(range(sampleSize[1:length(lowestDist)]), n=10))
    mtext("Sample size",side=1,col="black",line=2.5)
    
    ## Add Legend
    legend("topleft",legend="Distance (low is good)",
           text.col="black",pch=c(16,15),col="black")
    
    dev.off() 

    #Create a vector of improvements between sample sizes
    improvement <- lowestDist
    #Here, the first position, that cannot be related to any smaller fractions, is turned to one.
    improvement[1] <- 1
    for(i in 2:length(lowestDist)){
      improvement[i] <- lowestDist[i-1]-lowestDist[i]
    }
    
    sampleSizeOpt <- data.frame("SampleSize"=sampleSize[1:length(lowestDist)], "Lowest distance"=lowestDist, "Improvement"=improvement)
    
    #Now, as a final step, an optimiation is performed with the optimal sample size
    optPenalty <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, maxIter=maxIter, bootstrapObservations=sampleSizeOpt$SampleSize[length(lowestDist)], penaltyOffset=penaltyOffset, makeGraph=TRUE, graphName="Optimization of penalty term with optimal sample size.pdf")
    
    
    dClustOptObject <- c(sampleSizeOpt, optPenalty)
    return(dClustOptObject)
  }
}