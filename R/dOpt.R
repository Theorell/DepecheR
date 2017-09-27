#' Find the smallest optimal sampe size and the best penalty offset for a certain dataset
#'
#'
#' This function is used before the dClust function to identify the smallest sample size that gives rise to the most stable clustering. It is primarily used when datasets are very large (i.e. >100 000 observations), where there is a computational gain in clustering based on a subset of the cells and then assigning all other cells to the correct cluster center. It can however be generally useful when it is not entirely clear which sample size that is needed to perform the optimization of the penalty term, something that this function also works as a wrapper for.
#' @importFrom graphics box
#' @importFrom ggplot2 ggplot aes geom_line ggtitle xlab ylab ylim ggsave
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
#' @param initSampleSize The initial sample size. If penaltyOptOnly==TRUE, this value will define the sample size for the penalty optimization. 
#' @param sampleSizePowerIncrement With what value the sample size will increase. The value is taken to the power of two, i e, the standard value of 0.5 translates to a start value of x*2^0, and the following value of 2^(0.5), and so on. x 
#' @param maxSampleSize At what value the algorithm stops. Defaults to 100 000. 
#' @param initCenters Number of starting points for clusters. The higher the number, the greater the precision of the clustering, but the computing time is also increased with the number of starting points. Default is 30.
#' @param maxIter The maximal number of iterations that are performed to reach the minimal improvement. 
#' @param minImprovement This is connected to the evaluation of the performance of the algorithm. The more iterations that are run, or the larger the samples are, the smaller will the improvement be, and this sets the threshold when the iterations stop. 
#' @param penaltyOptOnly If this is set to true, no optimization on the sample size is performed. This might be useful in situations when the whole dataset will be used for clustering, such as when the dataset is small, or when the number of variables is very large (>100), so that increasing the sample size too far might lead to computational overload. If this is set to true, the sample size will be the value specified in the sampleSizeIncrement argument.
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
dOpt <- function(inDataFrameScaled, initSampleSize=2000, sampleSizePowerIncrement=0.5, maxSampleSize=100000, initCenters=30, maxIter=100, minImprovement=0.01, penaltyOptOnly=FALSE, penalties=c(0,2,4,8,16,32,64,128)){

  
  #First, the optimal penalties is identified with a reasonable sample size
  if(penaltyOptOnly==TRUE){
    penaltyOptOnly <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, maxIter=maxIter, minImprovement=minImprovement, bootstrapObservations=sampleSizeIncrement, penalties=penalties, makeGraph=TRUE, graphName="Optimization of penalty term.pdf")
    return(penaltyOptOnly)
  } else {
    
    dOptPenaltyResultList <- list()
	  lowestDist <- vector()
	  i <- 1
	  sampleSizes <- initSampleSize*2^0
	  #This loop continues to run until two runs produce the lowest distance at the same penalty, and these distances diverge less than 0.01.
	  while(i<=3 && sampleSizes[i]<=maxSampleSize || (dOptPenaltyResultList[[i-1]][[1]][1,1]!=dOptPenaltyResultList[[i-2]][[1]][1,1]) || (abs(lowestDist[i-2]-lowestDist[i-1])>minImprovement)){
	    if(i>1){
	      print(paste("Cycle", i-1, "completed. Jumping to next sample size."))
	    }
	    
	    dOptPenaltyResult <- dOptPenalty(inDataFrameScaled, initCenters=initCenters, maxIter=maxIter, bootstrapObservations=sampleSizes[i], penalties=penalties, makeGraph=FALSE, disableWarnings=TRUE)

	    dOptPenaltyResultList[[i]] <- dOptPenaltyResult
	    lowestDist[i] <- min(dOptPenaltyResult[[2]][,1:2])

	    
	    i <- i+1
	    sampleSizes[i] <- initSampleSize*2^(sampleSizePowerIncrement*(i-1))
	    
	  }
	  
	 print(paste("Cycle", i-1, "optimal."))
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
	  SampleSizes <- as.factor(rep(round(sampleSizes[1:length(sampleSizes)-1]), each=length(penalties)))
	  
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
    
    sampleSizeOpt <- data.frame("SampleSize"=sampleSizes[1:(length(sampleSizes)-1)], "Lowest distance"=lowestDist, "Improvement"=improvement)
    
    #Now, as a final step, an optimiation is performed with the optimal sample size
    optPenalty <- dOptPenaltyOptSampleSize
    
    
    dClustOptObject <- c(sampleSizeOpt, optPenalty)
    return(dClustOptObject)
  }
}