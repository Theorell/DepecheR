#' Allocation of observations to pre-established cluster centers.
#'
#'
#' Here, observations of a dataset are allocated to a set of preestablished cluster centers. This is intended to be used for the test set in train-test dataset situations. 
#' @importFrom moments kurtosis
#' @param inDataFrame A dataframe or matrix with the data that that the cluster centers will be allocated to. This data should be scaled in the same way as the data for the original depeche was scaled  when it entered the algorithm, i.e. in the normal case, not at all.
#' @param clusterCenters A matrix that needs to be inherited from a depeche run. It contains the information about which clusters and variables that have been sparsed away and where the cluster centers are located for the remaining clusters and variables.
#' @param log2Off If the automatic detection for high kurtosis, and followingly, the log2 transformation, should be turned off.
#' @seealso \code{\link{depeche}}
#' @return A vector with the same length as number of rows in the inDataFrame, where the cluster identity of each observation is noted.
#'
#' @examples
#' #Retrieve some example data
#' data(testData)
#' 
#' #Now arbitrarily (for the sake of the example) divide the data into a 
#' #training- and a test set.
#' testDataSample <- sample(1:nrow(testData), size=48500)
#' testDataTrain <- testData[testDataSample,]
#' testDataTest <- testData[-testDataSample,]
#' 
#' #Run the depeche function for the train set
#' \dontrun{
#' 
#' x_depeche_train <- depeche(testDataTrain[,2:15], maxIter=20, sampleSize=1000)
#'
#' #Allocate the test dataset to the centers of the train dataset
#' x_depeche_test <- dAllocate(testDataTest[,2:15], 
#' clusterCenters=x_depeche_train$clusterCenters)
#'
#' #And finally plot the two groups to see how great the overlap was:
#' trainTablePerId <- apply(as.matrix(table(testDataTrain$ids, 
#' x_depeche_train$clusterVector)), 1, function(x) x/sum(x))
#' trainTableCollapsed <- apply(trainTablePerId, 1, sum)
#' trainTableFraction <- trainTableCollapsed/sum(trainTableCollapsed)
#' testTablePerId <- apply(as.matrix(table(testDataTest$ids, x_depeche_test)), 1, function(x) x/sum(x))
#' testTableCollapsed <- apply(testTablePerId, 1, sum)
#' testTableFraction <- testTableCollapsed/sum(testTableCollapsed)
#' xmatrix <- t(cbind(trainTableFraction, testTableFraction))
#' library(gplots)
#' barplot2(xmatrix, beside = TRUE, legend = rownames(xmatrix))
#' title(main = "Difference between train and test set")
#' title(xlab = "Clusters")
#' title(ylab = "Fraction")
#' }
#' @export dAllocate
dAllocate <- function(inDataFrame, clusterCenters, log2Off=FALSE){

  if(class(inDataFrame)=="matrix"){
    inDataFrame <- as.data.frame.matrix(inDataFrame)
  }
  
  if(log2Off==FALSE && kurtosis(as.vector(as.matrix(inDataFrame)))>100){
    kurtosisValue1 <- kurtosis(as.vector(as.matrix(inDataFrame)))
    #Here, the log transformation is performed. In cases where the lowest value is 0, everything is simple. In other cases, a slightly more complicated formula is needed
    if(min(inDataFrame)>=0){
      inDataFrame <- log2(inDataFrame+1)
    } else {
      #First, the data needs to be reasonably log transformed to not too extreme values, but still without loosing resolution.
      inDataMatrixLog <- log2(apply(inDataFrame, 2, function(x) x-min(x))+1)
      #Then, the extreme negative values will be replaced by 0, as they give rise to artefacts.
      inDataMatrixLog[which(is.nan(inDataMatrixLog))] <- 0
      inDataFrame <- as.data.frame(inDataMatrixLog)
      rm(inDataMatrixLog)
    }
    
    kurtosisValue2 <- kurtosis(as.vector(as.matrix(inDataFrame)))
    print(paste("The data was found to be heavily tailed (kurtosis", kurtosisValue1, "). Therefore, it was log2-transformed, leading to a new kurtosis value of", kurtosisValue2, "."))
  }
  
  #If some variables have been excluded as they did not contribute to construction of any cluster, they are removed from the inData here
  inDataFrameReduced <- inDataFrame[,colnames(clusterCenters)]
  
  dataMat <- data.matrix(inDataFrameReduced)
  centersMat <- data.matrix(clusterCenters)
  
  clusterReallocationResult <- allocate_points(dataMat,centersMat,1)[[1]]

  #Here, the individual numbers are changed to accomodate the difference between the inclusion or exclusion of an origo cluster
    newNumbers <- rownames(clusterCenters)
    clusterReallocationResult <- turnVectorEquidistant(clusterReallocationResult, newNumbers=newNumbers)

return(clusterReallocationResult)
  
}
