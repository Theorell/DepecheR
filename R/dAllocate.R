#' Allocation of observations to pre-established cluster centers.
#'
#'
#' Here, observations of a dataset are allocated to a set of preestablished cluster centers. This is intended to be used for the test set in train-test dataset situations. It is called "predict" as most similar functions of other clustering algorithms have this term.
#' @param inDataFrameScaled A dataframe with the data that that the cluster centers will be allocated to. The data in this dataframe should be scaled in the same way as the dataframe used to generate the clusters in the first place.
#' @param clusterCenters This is a matrix that needs to be inherited from a dClust run. It contains the information about which clusters and variables that have been sparsed away and where the cluster centers are located for the remaining clusters and variables.
#' @param ids A vector of the same length as rows in the inDataFrameScaled. If included, it is used to generate a table of what fraction of observations for each individual that is present in each cluster.
#' @seealso \code{\link{dClust}}
#' @return A list with two components:
#' \describe{
#'     \item{realloClusterVector}{A vector with the same length as number of rows in the inDataFrameScaled, where the cluster identity of each observation is noted.}
#'     \item{realloClusterPercentagesForAllIds}{A matrix showing the percentage of observations for each id in each cluster.}
#' }
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateBimodalData()
#'
#' #Scale this datamframe
#' x_scaled <- dScale(x[,2:ncol(x)])
#'
#' #Divide this scaled dataframe in two parts
#' x_scaled_train <- x_scaled[1:5000,]
#' x_scaled_test <- x_scaled[5001:10000,]
#'
#' #Create two completely meaningless ids vectors
#' id_vector_train <- c(rep("Train 1", 2500), rep("Train 2", 2500))
#' id_vector_test <- c(rep("Test 3", 2500), rep("Test 4", 2500))
#' 
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the dOptAndClust function for the train set
#' x_dClust_train <- dClust(x_scaled_train, maxIter=20, sampleSizes=1000, ids=id_vector_train)
#'
#' #Retrieve the clustering info
#' clusterVector <- x_dClust_train[[1]]
#' clusterCenters <- x_dClust_train[[2]]
#' 
#' #This is followed by running the actual function in question
#' x_dClust_test <- dAllocate(x_scaled_test, 
#' clusterCenters=clusterCenters, ids=id_vector_test)
#'
#' #And finally plot this to see how great the overlap was:
#' xmatrix <- t(cbind(x_dClust_train$idClusterFractions, 
#' x_dClust_test$realloIdClusterFractions))
#' library(gplots)
#' barplot2(xmatrix, beside = TRUE, legend = rownames(xmatrix))
#' title(main = "Difference between train and test set")
#' title(xlab = "Clusters")
#' title(ylab = "Fraction")
#' @export dAllocate
dAllocate <- function(inDataFrameScaled, clusterCenters, ids){

  #If some variables have been excluded as they did not contribute to construction of any cluster, they are removed from the inData here
  inDataFrameReduced <- inDataFrameScaled[,colnames(clusterCenters)]
  
  dataMat <- data.matrix(inDataFrameReduced)
  centersMat <- data.matrix(clusterCenters)
  
  clusterReallocationResult <- allocate_points(dataMat,centersMat,1)[[1]]

  #Here, the individual numbers are changed to accomodate the difference between the inclusion or exclusion of an origo cluster
    newNumbers <- rownames(clusterCenters)
    clusterReallocationResult <- turnVectorEquidistant(clusterReallocationResult, newNumbers=newNumbers)

#A table with the percentage of cells in each cluster for each individual is created
  if(missing(ids)==FALSE && length(ids)==nrow(inDataFrameScaled)){
    clusterTable <- table(clusterReallocationResult, ids)

    countTable <- table(ids)

    idClusterFractions <- clusterTable

    for(i in 1:length(countTable)){
    	x <- clusterTable[,i]/countTable[i]
    	idClusterFractions[,i] <- x
    }
    clusterReallocationResult <- list(clusterReallocationResult, as.data.frame.matrix(idClusterFractions))
    
    names(clusterReallocationResult) <- c("realloClusterVector", "realloIdClusterFractions")
  }

return(clusterReallocationResult)
  
}
