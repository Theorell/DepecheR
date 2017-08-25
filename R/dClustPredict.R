#' Allocation of observations to pre-established cluster centers.
#'
#'
#' Here, observations of a dataset are allocated to a set of preestablished cluster centers. This is intended to be used for the test set in train-test dataset situations. It is called "predict" as most similar functions of other clustering algorithms have this term.
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the quantileScale function. It should naturally be scaled together with the data used to genreate the cluster centers.
#' @param penalizedClusterCenters This is a matrix that needs to be inherited from a dClust run. It contains the information about which clusters and variables that have been sparsed away and where the cluster centers are located for the remaining clusters and variables.
#' @param withOrWithoutZeroClust This parameter controls if the generated result should contain a cluster in origo or not. This information is given by dClustOpt, again.
#' @param ids A vector of the same length as rows in the inDataFrameScaled. It is used to generate the final analysis, where a table of the percentage of observations for each individual and each cluster is created.
#' @seealso \code{\link{dClustOpt}}, \code{\link{dClust}}
#' @return A list with two components:
#' \describe{
#'     \item{realloClusterVector}{A vector with the same length as number of rows in the inDataFrameScaled, where the cluster identity of each observation is noted.}
#'     \item{realloClusterPercentagesForAllIds}{A matrix showing the percentage of observations for each id in each cluster.}
#' }
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData()
#'
#' #Scale this datamframe
#' x_scaled <- quantileScale(x[,2:ncol(x)])
#'
#' #Divide this scaled dataframe in two parts
#' x_scaled_train <- x_scaled[1:5000,]
#' x_scaled_test <- x_scaled[5001:10000,]
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the dClustOpt function to get good starting points
#' x_optim <- dClustOpt(x_scaled_train, iterations=5, bootstrapObservations=1000)
#'
#' #Create two completely meaningless ids vectors
#' id_vector_train <- c(rep("Train 1", 2500), rep("Train 2", 2500))
#' id_vector_test <- c(rep("Test 3", 2500), rep("Test 4", 2500))
#'
#' #Then run the dClust function for the train set
#' x_dClust_train <- dClust(x_scaled_train, regVec=x_optim[[1]][["optimalRegularizationValue"]], 
#' withOrWithoutZeroClust=x_optim[[1]][["withOrWithoutZeroClust"]], iterations=2, ids=id_vector_train)
#'
#' #This is followed by running the actual function in question
#' x_dClust_test <- dClustPredict(x_scaled_test, 
#' penalizedClusterCenters=x_dClust_train$penalizedClusterCenters, 
#' withOrWithoutZeroClust=x_optim[[1]][["withOrWithoutZeroClust"]], ids=id_vector_test)
#'
#' #And finally plot this to see how great the overlap was:
#' xmatrix <- as.matrix(rbind(x_dClust_train$clusterPercentagesForAllIds, 
#' x_dClust_test$realloClusterPercentagesForAllIds))
#' library(gplots)
#' barplot2(xmatrix, beside = TRUE, legend = rownames(xmatrix))
#' title(main = "Difference between train and test set")
#' title(xlab = "Clusters")
#' title(ylab = "Percentage")
#' @export dClustPredict
dClustPredict <- function(inDataFrameScaled, penalizedClusterCenters, withOrWithoutZeroClust, ids){

myMat<-data.matrix(inDataFrameScaled, rownames.force = NA)

	if(withOrWithoutZeroClust=="stabWZero"){
	newInds <- allocate_points(myMat,penalizedClusterCenters,0)
	}
	if(withOrWithoutZeroClust=="stabWOZero"){
	newInds <- allocate_points(myMat,penalizedClusterCenters,1)
	}

#A table with the percentage of cells in each cluster for each individual is created

clusterTable <- table(newInds[[1]], ids)

countTable <- table(ids)

clusterPercentagesForAllIds <- clusterTable

for(i in 1:length(countTable)){
	x <- 100*clusterTable[,i]/countTable[i]
	clusterPercentagesForAllIds[,i] <- x
}

clusterReallocationResult <- list(newInds$i, as.data.frame.matrix(t(clusterPercentagesForAllIds)))

names(clusterReallocationResult) <- c("realloClusterVector", "realloClusterPercentagesForAllIds")

return(clusterReallocationResult)
}
