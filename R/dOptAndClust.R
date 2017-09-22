#' Wrapper function for dOpt and dClust
#' 
#' 
#' This function is a user-friendly wrapper integrating the dOpt and dClust functions. It only requires a dataset and an id vector. It starts by doing all necessary optimizations, both on the smallest sample size that is needed to perform the most stable clustering, and to identify the optimal penalty. It then performs clustering based on the values identified in the optimization step. 
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
#' @param ids A vector of the same length as rows in the inDataFrameScaled. It is used to generate the final analysis, where a table of the percentage of observations for each individual and each cluster is created.
#' @param sampleSizes The number of observations that are included in each bootstrap subsampling of the data in the optimization function. NB! The algorithm uses resampling, so the same event can be used twice. 
#' @param penalties These values are evaluated in dOpt. The number of suggested default values are empirically defined and might not be optimal for a specific dataset, but the algorithm will warn if the most optimal values are on the borders of the range. Note that when this offset is 0, there is no penalization, which means that the algorithm runs normal K-means clustering.
#' @param initCenters Number of starting points for clusters in both dOpt and dClust. The higher the number, the greater the precision of the clustering, but the computing time is also increased with the number of starting points. Default is 30.
#' @param maxIter The maximal number of iterations that are performed to reach the minimal improvement in the optimization steps. 
#' @param minImprovement This is connected to the evaluation of the performance of dOpt. The more iterations that are run, or the larger the samples are, the smaller will the improvement be, and this sets the threshold when the iterations stop. 
#' @seealso \code{\link{dAllocate}}, \code{\link{dOpt}}, \code{\link{dClust}}
#' @return A nested list with two component lists, of which the first is inherited from dOpt and the second from dClust. See these functions for details. 
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateBimodalData(samplings=2)
#'
#' #Scale this datamframe
#' x_scaled <- dScale(x[,2:ncol(x)])
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the dOptPenalty function to get good starting points
#' x_opt_clust <- dOptAndClust(x_scaled, ids=x[,1])
#'
#' #Look at the result
#' str(x_opt_clust)
#' @export dOptAndClust
dOptAndClust <- function(inDataFrameScaled, ids, sampleSizes=1000*c(2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10), penalties=c(0,2,4,8,16,32,64,128), initCenters=30, maxIter=100, minImprovement=0.01){
  
  if(missing(ids)){
    stop("Vector of ids is missing. Save youself some time and put it in before running again, as the function will otherwise throw an error at the end.")
  }
  
  if(exists("ids")==FALSE){
    stop("Ids is a non-existent object. Save youself some time and put in an existing one before running again, as the function will otherwise throw an error at the end.")
  }
  
  dOptObject <- dOpt(inDataFrameScaled, sampleSizes=sampleSizes, initCenters=initCenters, maxIter=maxIter, minImprovement=minImprovement, penaltyOptOnly=FALSE, penalties=penalties)

  dClustObject <- dClust(inDataFrameScaled, dOptObject=dOptObject, ids)
  
  dOptAndClustObject <- list(dOptObject, dClustObject)
  return(dOptAndClustObject)

}