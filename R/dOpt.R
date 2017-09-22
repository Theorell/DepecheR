#' Find the smallest optimal sampe size and the best penalty offset for a certain dataset
#'
#'
#' This function is used before the dClust function to identify the smallest sample size that gives rise to the most stable clustering. It is primarily used when datasets are very large (i.e. >100 000 observations), where there is a computational gain in clustering based on a subset of the cells and then assigning all other cells to the correct cluster center. It can however be generally useful when it is not entirely clear which sample size that is needed to perform the optimization of the penalty term, something that this function also works as a wrapper for.
#' @importFrom graphics box
#' @importFrom ggplot2 ggplot aes geom_line ggtitle xlab ylab ylim ggsave
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
#' @param sampleSizes The number of observations that are included in each bootstrap subsampling of the data. NB! The algorithm uses resampling, so the same event can be used twice. This is the central argument to this function, that it optimizes over.
#' @param initCenters Number of starting points for clusters. The higher the number, the greater the precision of the clustering, but the computing time is also increased with the number of starting points. Default is 30.
#' @param maxIter The maximal number of iterations that are performed to reach the minimal improvement. 
#' @param minImprovement This is connected to the evaluation of the performance of the algorithm. The more iterations that are run, or the larger the samples are, the smaller will the improvement be, and this sets the threshold when the iterations stop. 
#' @param penaltyOptOnly If this is set to true, no optimization on the sample size is performed. This might be useful in situations when the whole dataset will be used for clustering, such as when the dataset is small, or when the number of variables is very large (>100), so that increasing the sample size too far might lead to computational overload. If this is set to true, only one sample size should be specified in the sampleSizes argument.
#' @param penalties These values are the ones that are evaluated and the ones that decide the penalization. The number of suggested default values are empirically defined and might not be optimal for a specific dataset, but the algorithm will warn if the most optimal values are on the borders of the range. Note that when this offset is 0, there is no penalization, which means that the algorithm runs normal K-means clustering.
#' @seealso \code{\link{dClust}}, \code{\link{dAllocate}}, \code{\link{dOptAndClust}}
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
dOpt <- function(inDataFrameScaled, sampleSizes=1000*c(2^1, 2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10), initCenters=30, maxIter=100, minImprovement=0.01, penaltyOptOnly=FALSE, penalties=c(0,2,4,8,16,32,64,128)){

  
  #First, the optimal penalties is identified with a reasonable sample size
  if(penaltyOptOnly==TRUE){
    if(length(sampleSizes)>1){
      warning(paste("penaltyOptOnly is set to true, which means that no sample size optimization is performed. Still, multiple sample sizes have been added. Thus, only the first sample size will be used, in this case ", sampleSizes[1], ". If you do not like this, change the sampleSize to a number of your preference.", sep=""))
    }
    penaltyOptOnly <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, maxIter=maxIter, minImprovement=minImprovement, bootstrapObservations=sampleSizes[1], penalties=penalties, makeGraph=TRUE, graphName="Optimization of penalty term.pdf")
    return(penaltyOptOnly)
  } else {
    
    dOptPenaltyResultList <- list()
	  lowestDist <- vector()
	  for(i in 1:length(sampleSizes)){
	    ptm <- proc.time()
	    dOptPenaltyResult <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, maxIter=maxIter, bootstrapObservations=sampleSizes[i], penalties=penalties, makeGraph=FALSE, disableWarnings=TRUE)
	    timing <- proc.time() - ptm
	    dOptPenaltyResultList[[i]] <- dOptPenaltyResult
	    lowestDist[i] <- min(dOptPenaltyResult[[2]][,1:2])
		
	    if(i<2){
	      print(paste("Cycle", i, "completed. Jumping to next sample size."))
	    } else {
	      if(abs(lowestDist[i-1]-lowestDist[i])<=minImprovement){
	        print(paste("Cycle", i, "optimal."))
	        break
	      } 
	      #else if(timing[3]>30){
	        #ANSWER <- readline(paste("The improvement between the two latest cycles was ", lowestDist[i], " . ", "The next cycle could take",  timing[3]*2, " seconds. Do you wish to run it? Print no if you do not, and yes if you do.", sep=""))

	        #if (substr(ANSWER, 1, 1) == "n"){
	       #   cat("Ok. Lets stop here then.")
	       #   break
	       # }
			      
			   # else if (substr(ANSWER, 1, 1) == "y"){
			   #   cat("Ok. Lets go on.")
			   # }
			      
  		#	} 
	    else{
				print(paste("Cycle ", i, " completed. Jumping to next sample size.", sep=""))
  			}
	
      }
		
	  }
	  
	  #Now, the curve of distances with different penalties and different sample sizes are plotted together. 
	  #First, the optimal solution is retrieved from the last cycle.
	  dOptPenaltyOptSampleSize <- dOptPenaltyResultList[[length(dOptPenaltyResultList)]]
	  
	  #Then, all values for all sample sizes are made into a long dataframe
	  dOptPenaltyResults <- do.call("rbind", lapply(dOptPenaltyResultList, "[[", 2))
	  
	  #Here, the solution with or without origo cluster is retrieved, depending on the optimal solution
	  
	  if(dOptPenaltyOptSampleSize[[1]][1,2]=="yes"){
	    Distances <- dOptPenaltyResults[,1]
	  }
	  if(dOptPenaltyOptSampleSize[[1]][1,2]=="no"){
	    Distances <- dOptPenaltyResults[,2]
	  }
	  
	  #Now, a vector of penalties is created
	  Penalties <- rep(penalties, times=length(dOptPenaltyResultList))
	  
	  #Here, a vector of sample sizes is created instead
	  SampleSizes <- as.factor(rep(sampleSizes[1:length(dOptPenaltyResultList)], each=length(penalties)))
	  
	  #Now combine these three
	  plottingObject <- data.frame("Sample_sizes"=SampleSizes, Penalties, Distances)
	  
	  ggplot(data=plottingObject,
	         aes(x=sqrt(Penalties), y=Distances, colour=SampleSizes)) +
	    geom_line() + ylim(c(0,1)) +
	    ggtitle("Distance as a function of penalties for different sample sizes") +
	    xlab("Square root of the penalties") + ylab("Distance  between bootstrap subsamplings, [0 to 1]")
	    ggsave("Distance as a function of penalties for different sample sizes.pdf")
	    
    #Create a vector of improvements between sample sizes
    improvement <- lowestDist
    #Here, the first position, that cannot be related to any smaller fractions, is turned to one.
    improvement[1] <- 1
    for(i in 2:length(lowestDist)){
      improvement[i] <- lowestDist[i-1]-lowestDist[i]
    }
    
    sampleSizeOpt <- data.frame("SampleSize"=sampleSizes[1:length(lowestDist)], "Lowest distance"=lowestDist, "Improvement"=improvement)
    
    #Now, as a final step, an optimiation is performed with the optimal sample size
    optPenalty <- dOptPenaltyOptSampleSize
    
    
    dClustOptObject <- c(sampleSizeOpt, optPenalty)
    return(dClustOptObject)
  }
}