#' Fint the smallest optimal sampe size
#'
#'
#' This function is used before the dClust function to identify the smallest sample size that gives rise to the most stable clustering. It is primarily used when datasets are very large (i.e. >100 000 observations), where there is a computational gain in clustering based on a subset of the cells and then assigning all other cells to the correct cluster center.
#' @importFrom graphics box
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
#' @param sampleSize The number of observations that are included in each bootstrap subsampling of the data. NB! The algorithm uses resampling, so the same event can be used twice. This is the central argument to this function, that it optimizes over.
#' @param initCenters Number of starting points for clusters. This essentially means that it is the highest possible number of clusters that can be defined. The higher the number, the greater the precision, but the computing time is also increased with the number of starting points. Default is 30
#' @param iterations As it sounds, the number of bootstrap reiterations that are performed.
#' @param minimalImprovement Expressed as a fraction between 0 and 1, this value tells how much better a doubling of the sample size needs to be for the algorithm to continue to another step. Default is 0.001, corresponding to an improvement of one percent.
#' @param penaltyOptSampleSize The size of the initial samle that will be used to define the penalty offset that will be used in all subsequent optimizations of the sample size. 
#' @param penaltyOffset These values are the ones that are evaluated and the ones that decide the penalization. The number of suggested default values are empirically defined and might not be optimal for a specific dataset, but the algorithm will warn if the most optimal values are on the borders of the range. Note that when this offset is 0, there is no penalization, which means that the algorithm runs normal K-means clustering.
#' @param makeOptimGraph If a graph should be created showing the distance between bootstraps under different penalties during the initial pre-optimization.
#' @seealso \code{\link{dClust}}, \code{\link{pClustPredict}}
#' @return A graph showing the stability of the models with different sample sizes. Optionally also the optimization graphs for each step, see makeOptimGraphs above. In addition, a data frame is returned. This one has three columns:
#' \describe{
#'     \item{SampleSize}{This column shows the sample size of each boot strap subsampling in the optimization procedure.}
#'     \item{Lowest distance}{This vector shows the optimal stability, expressed as the lowest distance between the bootstrap subsampling runs at each of the boot strap subsamling sizes.}
#'     \item{Improvement}{Here, the improvement, expressed as a fraction between 0 and 1 is shown. When the improvement is less than 1 percent, the algorithm automatically stops.}
#'}
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
#' x_optim <- dOptSampleSize(x_scaled)
#' @export dOptSampleSize
dOptSampleSize <- function(inDataFrameScaled, sampleSize=1000*c(2^0, 2^1, 2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10), initCenters=30, iterations=50, minimalImprovement=0.01, penaltyOptSampleSize=16000, penaltyOffset=c(0,2,4,8,16,32,64,128), makeOptimGraph=TRUE){
	
  #Here, it is checked  with a rule of thumb if it is worth running the function
  if (nrow(inDataFrameScaled)<100000){
    ANSWER <- readline("The number of rows in the total inDataFrame does not exceed 100000. It is not worth running this function. Do you want to continue anyway? Answer with yes or no.")
    
    if (substr(ANSWER, 1, 1) == "n"){
      stop("Ok. Let us stop here then. This is not an error.")
    } else {
      cat("Ok. Let us go on anyway then.")
    }
    
  }
  
  #First, the optimal penaltyOffset is identified with a reasonable sample size
  bestPenaltyOffset <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, iterations=iterations, bootstrapObservations=penaltyOptSampleSize, penaltyOffset=penaltyOffset, makeGraph=makeOptimGraph)[[1]][1,1]
  
  print("Now, the pre-optimization of the penalty terms is done and the sample size optimization will start.")
  
	lowestDist <- vector()
	for(i in 1:length(sampleSize)){
	  ptm <- proc.time()
		dOptPenaltyResult <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, iterations=iterations, sampleSize[i], penaltyOffset=c(penaltyOffset[which(penaltyOffset==bestPenaltyOffset)], penaltyOffset[which(penaltyOffset==bestPenaltyOffset)-1], penaltyOffset[which(penaltyOffset==bestPenaltyOffset)+1]), makeGraph=FALSE, disableWarnings=TRUE)
		#graphName=paste("Distance over penalty values for sample size ", sampleSize[i], ".pdf", sep="")
		timing <- proc.time() - ptm
		lowestDist[i] <- min(dOptPenaltyResult[[2]][,1:2])
		
		if(i<=3){
			print(paste("Cycle", i, "completed. Jumping to next sample size."))
		}
		if(i>3){
			if(abs(lowestDist[i-1]-lowestDist[i])<=minimalImprovement){
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
			      
			  
			}
			  
			  else{
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
    
    result <- data.frame("SampleSize"=sampleSize[1:length(lowestDist)], "Lowest distance"=lowestDist, "Improvement"=improvement)
    return(result)

}