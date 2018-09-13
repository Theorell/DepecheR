# Allocation of observations to pre-established cluster centers.
#
#
# Here, observations of a dataset are allocated to a set of preestablished cluster centers. This is intended to be used for the test set in train-test dataset situations. 
# @param inDataFrameScaled A dataframe or matrix with the data that that the cluster centers will be allocated to. This data does not need to be scaled. Rather, if 
# @param clusterCenters A matrix that needs to be inherited from a depeche run. It contains the information about which clusters and variables that have been sparsed away and where the cluster centers are located for the remaining clusters and variables.
# @param ids A vector of the same length as rows in the inDataFrameScaled. If included, it is used to generate a table of what fraction of observations for each individual that is present in each cluster.
# @seealso \code{\link{depeche}}
# @return A list with two components:
# \describe{
#     \item{realloClusterVector}{A vector with the same length as number of rows in the inDataFrameScaled, where the cluster identity of each observation is noted.}
#     \item{realloClusterPercentagesForAllIds}{A matrix showing the percentage of observations for each id in each cluster.}
# }
# @examples
# #Retrieve some example data
# data(testData)
# 
# #Now arbitrarily (for the sake of the example) divide the data into a 
# #training- and a test set.
# testDataSample <- sample(1:nrow(testData), size=48500)
# testDataTrain <- testData[testDataSample,]
# testDataTest <- testData[-testDataSample,]
# 
# #Set a reasonable working directory, e.g.
# setwd("~/Desktop")
#
# #Run the depeche function for the train set
# x_depeche_train <- depeche(testDataTrain[,2:15], maxIter=20, sampleSizes=1000,
# ids=testDataTrain$ids)
#
# #Allocate the test dataset to the centers of the train dataset
# x_depeche_test <- dAllocate(testDataTest[,2:15], 
# clusterCenters=x_depeche_train$clusterCenters, ids=testDataTest$ids)
#
# #And finally plot this to see how great the overlap was:
# xmatrix <- t(cbind(x_depeche_train$idClusterFractions, 
# x_depeche_test$realloIdClusterFractions))
# library(gplots)
# barplot2(xmatrix, beside = TRUE, legend = rownames(xmatrix))
# title(main = "Difference between train and test set")
# title(xlab = "Clusters")
# title(ylab = "Fraction")
# @export dAllocate
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
