#' Performing optimization and penalized K-means clustering
#' 
#'
#' This function is a user-friendly wrapper integrating the dOpt and dClust functions. It only requires a dataset and an id vector. It starts by doing all necessary optimizations, both on the smallest sample size that is needed to perform the most stable clustering, and to identify the optimal penalty. It then performs clustering based on the values identified in the optimization step. 
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
#' @param k Number of initial cluster centers. The higher the number, the greater the precision of the clustering, but the computing time is also increased with the number of starting points. Default is 30. If penalties=0, k-means clustering with k clusters will be performed.
#' @param penalties This argument decides whether a single penalty will be used for clustering, or if multiple penalties will be evaluated to identify the optimal one. The suggested default values are empirically defined and might not be optimal for a specific dataset, but the algorithm will warn if the most optimal values are on the borders of the range. Note that when this offset is 0, there is no penalization, which means that the algorithm runs normal K-means clustering.
#' @param minCRIImprovement This is the stop criterion for the penalty optimization algorithm: the more iterations that are run, the smaller will the improvement of the corrected Rand index be, and this sets the threshold when the inner iterations stop. 
#' @param maxIter The maximal number of iterations that are performed in the penalty optimization. 
#' @param sampleSizes If the full dataset should not be used for the penalty optimization and the clustering, this term specifies one or multiple values that will be evaluated regarding their ability to correctly classify all events. Only meaningful when datasets are very large, i e >1 000 000 data points, as each cycle takes considerable time. "default" tries to adjust to this: if an extensive dataset is added, a number of smaller samples will be tested instead of running the full dataset.
#' @param maxCRI This is the stop criterion for the iterative optimization of the sample size: the maximum corrected Rand index that is acceptable. Defaults to 0.01, or 1 percent.
#' @param withOrigoClust In the event that no optimization of either penalties or sample sizes is performed, this argument specifies if a solution with a cluster in origo should be included or not. No default.
#' @param ids Optionally, a vector of the same length as rows in the inDataFrameScaled can be included. If so, it is used to generate a final analysis, where a table of the fraction of observations for each individual and each cluster is created.
#' @seealso \code{\link{dAllocate}}, \code{\link{dOpt}}, \code{\link{dClust}}
#' @return A nested list with varying components depending on the setup above:
#'  \describe{
#'    \item{clusterVector}{A vector with the same length as number of rows in the inDataFrameUsed, where the cluster identity of each observation is noted.}
#'    \item{clusterCenters}{A matrix containing information about where the centers are in all the variables that contributed to creating the cluster with the given penalty term.}
#'    \item{clusterCentersWZeroVariables}{A matrix similar to clusterCenters, but where all variables are included, regardless of if they contributed to creating the clusters or not. To be used in dAllocate.}
#'    \item{penaltyOptList}{This is only included when multiple penalties are run. A list of two dataframes:
#'     \describe{
#'               \item{penaltyOpt.df}{A dataframe with one row with all the information about which settings that were used to generate the optimal clustering with the optimal sample size. The "withOrigoClust" information tells the user if the solution with or without a cluster in origo gives the most optimal solution. If yes, this origo population is generally small and could be viewed as not fitting in the model.}
#'               \item{meanOptimDf}{A dataframe with the information about the results with all tested penalty values}
#'              }
#'     }
#'     \item{sampleSizeOptList}{This is only included if multiple sample sizes are run. It is a dataframe, in which each row represents one sample size, and in which the last row is thus the chosen, optimal sample size. It has the following columns:
#'     \describe{
#'               \item{SampleSize}{This column shows the sample size of each boot strap subsampling in the optimization procedure.}
#'               \item{Lowest distance}{This vector shows the optimal stability, expressed as the lowest distance between the bootstrap subsampling runs at each of the boot strap subsamling sizes.}
#'               \item{Improvement}{Here, the improvement, expressed as a fraction between 0 and 1 is shown. When the improvement is less than minCRIImprovement, the algorithm automatically stops.}
#'              } 
#'     }
#'     \item{idClusterFractions}{If a valid ids vector is included, this dataframe is returned that contains the what fraction of each id that is present in each cluster. Calculated on a per id basis.}
#' }
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
#' #First, just run with the standard settings
#' xClustObject <- dClust(x_scaled, ids=x[,1])
#'
#' #Look at the result
#' str(xClustObject)
#' 
#' #Then try changing to a few different sample sizes
#' xClustObject <- dClust(x_scaled, ids=x[,1], sampleSizes=c(100, 500, 1000))
#'
#' #Look at the result
#' str(xClustObject)
#' 
#' @useDynLib DepecheR
#' @export dClust
dClust <- function(inDataFrameScaled, ids, k=30, penalties=c(0,2,4,8,16,32,64,128), minCRIImprovement=0.01, sampleSizes="default", maxCRI=0.01, maxIter=100, withOrigoClust){

  if(sampleSizes=="default" && nrow(inDataFrameScaled)>1000000){
    sampleSizes <- c(10000, 50000, 100000)
  } else if(sampleSizes=="default" && nrow(inDataFrameScaled)<=1000000){
    sampleSizes <- nrow(inDataFrameScaled)
  }
  
  #First the two simplest cases, namely no optimization and all cells, and a pre-decided sample size, respectively
  if(length(penalties)==1 && length(sampleSizes)==1){
   dClustResult <-  dClustCoFunction(inDataFrameScaled, sampleSize=sampleSizes, penalty=penalties, k=k, withOrigoClust=withOrigoClust)

  }

  #Now, the cases where penalty optimization is performed, but the full dataset is used.
  if(length(penalties)>1 && sampleSizes[1]==nrow(inDataFrameScaled)){
    dOptPenaltyObject <- dOptPenalty(inDataFrameScaled, k=k, maxIter=maxIter, minCRIImprovement=minCRIImprovement, bootstrapObservations=nrow(inDataFrameScaled), penalties=penalties, makeGraph=TRUE, graphName="Optimization of penalties.pdf", disableWarnings=TRUE, returnClusterCenters=FALSE)
    dClustResult <- dClustCoFunction(inDataFrameScaled, sampleSize=nrow(inDataFrameScaled), dOptPenaltyObject=dOptPenaltyObject)
    newClustResultRow <- length(dClustResult)+1
    dClustResult[[newClustResultRow]] <- dOptPenaltyObject

  }  

  #Now, the cases where sample size(s) is not equal to the full dataset and penalty optimization should be performed.
  if(length(penalties)>1 && sampleSizes[1]!=nrow(inDataFrameScaled)){
      dClustResult <- dOptSubset(inDataFrameScaled=inDataFrameScaled, sampleSizes=sampleSizes, k=k, maxIter=maxIter, maxCRI=maxCRI, minCRIImprovement=minCRIImprovement, penalties=penalties)
  }   
  
  ######################################
  
  #Provided that a viable Id vector is added, a table with the percentage of cells in each cluster for each individual is created
  
  if(missing(ids)==FALSE && length(ids)==nrow(inDataFrameScaled)){

    clusterTable <- table(dClustResult$clusterVector, ids)

    countTable <- table(ids)
    
    clusterFractionsForAllIds <- clusterTable
    
    for(i in 1:length(countTable)){
      x <- clusterTable[,i]/countTable[i]
      clusterFractionsForAllIds[,i] <- x
    }

    nextClustResultPosition <- length(dClustResult)+1
    dClustResult[[nextClustResultPosition]] <- as.data.frame.matrix(clusterFractionsForAllIds)
    names(dClustResult)[[length(dClustResult)]] <- "idClusterFractions"
    
  }

  #Here, a heatmap over the cluster centers is saved. Only true if the number of clusters exceeds one.
  reducedClusterCentersColRow <- dClustResult[[2]]
  if(nrow(reducedClusterCentersColRow)>1){
    pdf("Cluster centers.pdf")
    heatmap.2(as.matrix(reducedClusterCentersColRow), col=colorRampPalette(c("blue", "white", "red"))(100), trace="none")
    dev.off()    
  }
  
  return(dClustResult)

}