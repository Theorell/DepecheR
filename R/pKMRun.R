#' Function to run penalized K means
#'
#'
#' This function is the core user function of the Depeche package. It clusters the data with a penalized version of K-means.
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the quantileScale function.
#' @param regVec The parameter that controls the level of penalization. Preferrably, it should be inherited from a pKMOptim run, as the algorithm will then generate the most stable result.
#' @param withOrWithoutZeroClust This parameter controls if the generated result should contain a cluster in origo or not. This information is given by pKMOptim, again.
#' @param kVec Number of starting points for clusters. This essentially means that it is the highest possible number of clusters that can be defined. The higher the number, the greater the precision, but the computing time is also increased with the number of starting points. Default is 30.
#' @param iterations As it sounds, this controls how many iterations that are performed, among which the most stable is chosen. If pKMOptim has been performed before, this number should not need to be extensive. Default is 10.
#' @param ids A vector of the same length as rows in the inDataFrameScaled. It is used to generate the final analysis, where a table of the percentage of observations for each individual and each cluster is created.
#' @seealso \code{\link{pKMOptim}}, \code{\link{pKMPredict}}
#' @return A list with three components:
#' \describe{
#'     \item{clusterVector}{A vector with the same length as number of rows in the inDataFrameScaled, where the cluster identity of each observation is noted.}
#'     \item{penalizedClusterCenters}{A matrix containing information about which variables and clusters that have been penalizedd out, shown by zero, and for the non-penalized variables and clusters, where the cluster center is located. This is done with one row per cluster.}
#'     \item{clusterPercentagesForAllIds}{A matrix showing the percentage of observations for each id in each cluster.}
#' }
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData(samplings=2)
#'
#' #Scale this datamframe
#' x_scaled <- quantileScale(x[,2:ncol(x)])
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the pKMOptim function to get good starting points
#' x_optim <- pKMOptim(x_scaled, iterations=50, bootstrapObservations=1000)
#'
#' #Then run the actual function
#' x_pKM <- pKMRun(x_scaled, regVec=x_optim[[1]][["optimalRegularizationValue"]], withOrWithoutZeroClust=x_optim[[1]][["withOrWithoutZeroClust"]], iterations=1, ids=x[,1])
#'
#' #And finally look at your great result
#' str(x_pKM)
#' @export pKMRun
pKMRun <- function(inDataFrameScaled, regVec, withOrWithoutZeroClust, kVec=30, iterations=10, ids){

  if(missing(ids)){
    stop("Vector of ids is missing. Save youself some time and put it in before running again, as the function will otherwise throw an error at the end.")
  }

  if(exists("ids")==FALSE){
    stop("Ids is a non-existent object. Save youself some time and put in an existing one before running again, as the function will otherwise throw an error at the end.")
  }

  dataMat<-data.matrix(inDataFrameScaled, rownames.force = NA)

#This is the most central function of the whole package.
	no_cores <- detectCores() - 1
	cl <- makeCluster(no_cores, type="FORK")
	return_all <-parLapply(cl,0:iterations,function(x) sparse_k_means(dataMat,kVec,regVec,1,x))
	stopCluster(cl)

#Here, the lowest iteration is retrieved

	for(i in 1:length(return_all)){
		minimumN <- min(do.call("rbind", lapply(return_all, "[[", 5)))
		if(return_all[[i]][5]==minimumN){
			returnLowest <- return_all[[i]]

		}
	}



	if(withOrWithoutZeroClust=="stabWZero"){
			clusterVector <- returnLowest$i
			penalizedClusterCenters <- returnLowest$c

	}

	if(withOrWithoutZeroClust=="stabWOZero"){
			clusterVector <- returnLowest$o
			penalizedClusterCenters <- returnLowest$v
	}


	colnames(penalizedClusterCenters) <- colnames(inDataFrameScaled)
	rownames(penalizedClusterCenters) <- c(seq(0,kVec-1))

#A table with the percentage of cells in each cluster for each individual is created

clusterTable <- table(clusterVector, ids)
	if(nrow(clusterTable)==1){
		print("Warning. The number of clusters is one, i.e. no separation of observations is present with the current settings. Try lowering the regVec.")
	}

countTable <- table(ids)

clusterPercentagesForAllIds <- clusterTable

for(i in 1:length(countTable)){
	x <- 100*clusterTable[,i]/countTable[i]
	clusterPercentagesForAllIds[,i] <- x
}


pKMResult <- list(clusterVector, penalizedClusterCenters, as.data.frame.matrix(t(clusterPercentagesForAllIds)))

names(pKMResult) <- c("clusterVector", "penalizedClusterCenters", "clusterPercentagesForAllIds")


	return(pKMResult)

}

