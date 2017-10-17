# Function to run penalized K means
#
#
# This function is the core user function of the Depeche package. It clusters the data with a penalized version of K-means.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
#' @importFrom gplots heatmap.2
#' @importFrom dplyr sample_n
# @param inDataFrameUsed A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
# @param dOptPenaltyObject This object contains information about the optimal penalty offset, solution with or without a cluster in origo, and the number of initial cluster centers that were used to find this optimal information.
# @param sampleSize By default inherited from dOptPenaltyObject. Number of observations that shoult be included in the initial clustering step. Three possible values. Either inherited, "All" or a user-specified number. Defaults to inheriting from dClustObject. If a dClustObject is not substituted, all rows in inDataFrameUsed are added by default. If another number, a sample is created from inDataFrameUsed. This is extra useful when clustering very large datasets. Replacement is set to TRUE.
# @param penalty By default inherited from dOptPenaltyObject. The parameter that controls the level of penalization. 
# @param withOrigoClust By default inherited from dOptPenaltyObject. This parameter controls if the generated result should contain a cluster in origo or not. 
# @param k By default inherited from dOptPenaltyObject. Number of starting points for clusters. This essentially means that it is the highest possible number of clusters that can be defined. The higher the number, the greater the precision, but the computing time is also increased with the number of starting points. Default is 30.
# @seealso \code{\link{dAllocate}}, \code{\link{dOpt}}, \code{\link{dOptAndClust}}
# @return A list with three components:
# \describe{
#     \item{clusterVector}{A vector with the same length as number of rows in the inDataFrameUsed, where the cluster identity of each observation is noted.}
#     \item{clusterCenters}{A matrix containing information about where the centers are in all the variables that contributed to creating the cluster with the given penalty term.}
#     \item{clusterPercentagesForAllIds}{A matrix showing the percentage of observations for each id in each cluster.}
# }
# @examples
# #Generate a default size dataframe with bimodally distributed data
# x <- generateBimodalData(samplings=2)
#
# #Scale this datamframe
# x_scaled <- dScale(x[,2:ncol(x)])
#
# #Set a reasonable working directory, e.g.
# setwd("~/Desktop")
#
# #Run the dOptPenalty function to get good starting points
# x_optim <- dOpt(x_scaled)
#
# #Then run the actual function
# x_dClust <- dClust(x_scaled, dOptPenaltyObject=x_optim, ids=x[,1])
#
# #And finally look at your great result
# str(x_dClust)
# @export dClustCoFunction
dClustCoFunction <- function(inDataFrameScaled, sampleSize, dOptPenaltyObject, penalty, withOrigoClust, k=30){

  if(missing(dOptPenaltyObject)==FALSE){
    penalty <- dOptPenaltyObject[[1]][1,1]
    withOrigoClust <- dOptPenaltyObject[[1]][1,2]
    k <- dOptPenaltyObject[[1]][1,3]
  }
  
  if(sampleSize==nrow(inDataFrameScaled)){
    inDataFrameUsed <- inDataFrameScaled
  } else {
    inDataFrameUsed <- sample_n(inDataFrameScaled, size=sampleSize)
  }
  
  penaltyConstant <- ((sampleSize*sqrt(ncol(inDataFrameUsed)))/1450)

  penaltyForRightSize <- penalty*penaltyConstant 

  dataMat<-data.matrix(inDataFrameUsed, rownames.force = NA)

  #Here the number of iterations is chosen. Very many are not needed, but a few will make the clustering even better than if just one was chosen.
  n_cores <- detectCores() - 1
  if(n_cores>=7){
    if(n_cores<=21){
      iterations <- n_cores
    } else {
      iterations <- 21
    }
  } else {
    iterations <- 7
  }
    
#This is the central function of the whole package.
	
  cl <-  parallel::makeCluster(iterations, type = "SOCK")
  registerDoSNOW(cl)
  return_all <- foreach(i=1:iterations) %dopar% sparse_k_means(dataMat,k,penaltyForRightSize,1, i)
  parallel::stopCluster(cl)	
  
  #Alternative deprecated parallelization.
	#if(Sys.info()['sysname']!="Windows"){
	#  cl <- makeCluster(n_cores, type="FORK")
	#  return_all <-parLapply(cl,0:iterations,function(x) sparse_k_means(dataMat,k,penaltyForRightSize,1,x))
	#  stopCluster(cl)
	#} else {
	#  cl <- makeCluster(n_cores, type="PSOCK")
	#  return_all <-parLapply(cl,0:iterations,function(x) sparse_k_means(dataMat,k,penaltyForRightSize,1,x))
	#  stopCluster(cl)
	#}
	
  #Here, the best iteration is retrieved
  logMaxLik <- as.vector(do.call("rbind", lapply(return_all, "[[", 5)))
  minimumN <- max(logMaxLik)
  returnLowest <- return_all[[which(abs(logMaxLik)==minimumN)[1]]]

  #And here, the optimal results, given if an origo cluster should be included or not, are retrieved further
  if(withOrigoClust=="yes"&& length(unique(returnLowest$i))<k){
    clusterVector <- returnLowest$i
    clusterCenters <- returnLowest$c
  } else {
    clusterVector <- returnLowest$o
    clusterCenters <- returnLowest$v
  }
  
  retrieveOrigoOrNotList <- retrieveOrigoOrNot(withOrigoClust=withOrigoClust, clusterVector, clusterCenters, colnamesClusterCenters=colnames(inDataFrameScaled), k=k)

  clusterVectorEquidistant <- retrieveOrigoOrNotList[[1]]
  reducedClusterCentersColRow <- retrieveOrigoOrNotList[[2]]
  reducedClusterCentersRow <- retrieveOrigoOrNotList[[3]]

  if(sampleSize!=nrow(inDataFrameScaled)){
    allocationData <- data.matrix(inDataFrameScaled, rownames.force = NA)
    noZero <- ifelse(withOrigoClust=="yes", 0, 1)
    clusterVectorAllocated <- allocate_points(allocationData, reducedClusterCentersRow, noZero)[[1]]
    if(noZero==0 && length(unique(returnLowest$i))<k){
      clusterVectorEquidistant <-  turnVectorEquidistant(clusterVectorAllocated, startValue=0)
    } else{
      clusterVectorEquidistant <- turnVectorEquidistant(clusterVectorAllocated)
    }
  }

  
  #Here, the results are combined
  dClustResult <- list(clusterVectorEquidistant, as.data.frame.matrix(reducedClusterCentersColRow), reducedClusterCentersRow)
  names(dClustResult) <- c("clusterVector", "clusterCenters", "clusterCentersWZeroVariables")
  return(dClustResult)

}

