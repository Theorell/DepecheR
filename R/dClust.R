#' Function to run penalized K means
#'
#'
#' This function is the core user function of the Depeche package. It clusters the data with a penalized version of K-means.
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @importFrom gplots heatmap.2
#' @importFrom dplyr sample_n
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the quantileScale function.
#' @param sampleSize Number of observations that shoult be included in the initial clustering step. Defaults to all rows in inDataFrameScaled. If another number, a sample is created from inDataFrameScaled. This is extra useful when clustering very large datasets. Replacement is set to TRUE.
#' @param penaltyOffset The parameter that controls the level of penalization. Preferrably, it should be inherited from a dClustOpt run, as the algorithm will then generate the most stable result.
#' @param withOrWithoutZeroClust This parameter controls if the generated result should contain a cluster in origo or not. This information is given by dClustOpt, again.
#' @param initCenters Number of starting points for clusters. This essentially means that it is the highest possible number of clusters that can be defined. The higher the number, the greater the precision, but the computing time is also increased with the number of starting points. Default is 30.
#' @param iterations As it sounds, this controls how many iterations that are performed, among which the most stable is chosen. If dClustOpt has been performed before, this number should not need to be extensive. Default is 10.
#' @param ids A vector of the same length as rows in the inDataFrameScaled. It is used to generate the final analysis, where a table of the percentage of observations for each individual and each cluster is created.
#' @seealso \code{\link{dClustOpt}}, \code{\link{dClustPredict}}
#' @return A list with three components:
#' \describe{
#'     \item{clusterVector}{A vector with the same length as number of rows in the inDataFrameScaled, where the cluster identity of each observation is noted.}
#'     \item{clusterCenters}{A matrix containing information about where the centers are in all the variables that contributed to creating the cluster with the given penalty term.}
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
#' #Run the dClustOpt function to get good starting points
#' x_optim <- dClustOpt(x_scaled, iterations=10, bootstrapObservations=1000)
#'
#' #Then run the actual function
#' x_dClust <- dClust(x_scaled, penaltyOffset=x_optim[[1]][["bestPenaltyOffset"]], 
#' withOrigoClust=x_optim[[1]][["withOrigoClust"]], iterations=1, ids=x[,1])
#'
#' #And finally look at your great result
#' str(x_dClust)
#' @export dClust
dClust <- function(inDataFrameScaled, sampleSize, penaltyOffset, withOrigoClust, initCenters=30, iterations=10, ids){

  if(missing(ids)){
    stop("Vector of ids is missing. Save youself some time and put it in before running again, as the function will otherwise throw an error at the end.")
  }

  if(exists("ids")==FALSE){
    stop("Ids is a non-existent object. Save youself some time and put in an existing one before running again, as the function will otherwise throw an error at the end.")
  }

  if(missing(sampleSize)){
    sampleSize <- nrow(inDataFrameScaled)
  }
  
  if(sampleSize==nrow(inDataFrameScaled)){
    inDataFrameUsed <- inDataFrameScaled
  } else {
    inDataFrameUsed <- sample_n(inDataFrameScaled, sampleSize, replace=TRUE)
  }

  #k <- ((sampleSize*sqrt(ncol(inDataFrameUsed)))/1450) #This is for the old algorithm.
  k <- ((sampleSize*sqrt(ncol(inDataFrameUsed)))/1450)
  penalty <- penaltyOffset*k
  
  dataMat<-data.matrix(inDataFrameUsed, rownames.force = NA)
  
#This is the most central function of the whole package.
	
	if(Sys.info()['sysname']!="Windows"){
	  no_cores <- detectCores() - 1
	  cl <- makeCluster(no_cores, type="FORK")
	  return_all <-parLapply(cl,0:iterations,function(x) sparse_k_means(dataMat,initCenters,penalty,1,x))
	  stopCluster(cl)
	} else {
	  no_cores <- detectCores() - 1
	  cl <- makeCluster(no_cores, type="PSOCK")
	  return_all <-parLapply(cl,0:iterations,function(x) sparse_k_means(dataMat,initCenters,penalty,1,x))
	  stopCluster(cl)
	}
	
  #Here, the best iteration is retrieved
  logDistance <- as.vector(do.call("rbind", lapply(return_all, "[[", 5)))
  minimumN <- max(logDistance)
  returnLowest <- return_all[[which(abs(logDistance)==minimumN)[1]]]


	if(withOrigoClust=="yes"){
			clusterVector <- returnLowest$i
			#Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is included, this cluster gets the denomination 0.
			clusterVectorEquidistant <- turnVectorEquidistant(clusterVector, startValue=0)
			penalizedClusterCenters <- returnLowest$c			  
			colnames(penalizedClusterCenters) <- colnames(inDataFrameUsed)
			#Remove all rows that do not contain any information.
			
			reducedPenalizedClusterCenters <- penalizedClusterCenters[which(rowSums(penalizedClusterCenters)!=0),which(colSums(penalizedClusterCenters)!=0)]
			
			#Add the zero cluster back
			reducedPenalizedClusterCentersOrigo <- penalizedClusterCenters[1,which(colSums(penalizedClusterCenters)!=0)]
			reducedPenalizedClusterCenters <- rbind(reducedPenalizedClusterCentersOrigo, reducedPenalizedClusterCenters)
	}

	if(withOrigoClust=="no"){
			clusterVector <- returnLowest$o
			#Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is not included, the first cluster gets the denomination 1.
			clusterVectorEquidistant <- turnVectorEquidistant(clusterVector)			
			penalizedClusterCenters <- returnLowest$v		
			colnames(penalizedClusterCenters) <- colnames(inDataFrameUsed)
			#Remove all rows that do not contain any information
			reducedPenalizedClusterCenters <- penalizedClusterCenters[which(rowSums(penalizedClusterCenters)!=0),which(colSums(penalizedClusterCenters)!=0)]

	}

	#Make the row names the same as the cluster names in the clusterVectorEquidistant
  rownames(reducedPenalizedClusterCenters) <- sort(unique(clusterVectorEquidistant))
  
	
	#If the number of observations used for the clustering is not equal to the number of rows in the inDataFrameScaled, here a prediction is made for all the other events. 
	if(sampleSize!=nrow(inDataFrameScaled)){
	  myMat<-data.matrix(inDataFrameScaled, rownames.force = NA)
	  
	  if(withOrigoClust=="yes"){
	    clusterVectorEquidistant <- unlist(allocate_points(myMat,reducedPenalizedClusterCenters,0))
	  }
	  if(withOrigoClust=="no"){
	    clusterVectorEquidistant <- unlist(allocate_points(myMat,reducedPenalizedClusterCenters,1))
	  }
	}
		
  #A table with the percentage of cells in each cluster for each individual is created

  clusterTable <- table(clusterVectorEquidistant, ids)
  	if(nrow(clusterTable)==1){
  		print("Warning. The number of clusters is one, i.e. no separation of observations is present with the current settings. Try lowering the penalty.")
  	}

  countTable <- table(ids)

  clusterFractionsForAllIds <- clusterTable

  for(i in 1:length(countTable)){
  	x <- clusterTable[,i]/countTable[i]
  	clusterFractionsForAllIds[,i] <- x
  }


  dClustResult <- list(clusterVectorEquidistant, as.data.frame.matrix(reducedPenalizedClusterCenters), as.data.frame.matrix(t(clusterFractionsForAllIds)))

  names(dClustResult) <- c("clusterVector", "clusterCenters", "idClusterFractions")

  #Here, a heatmap over the cluster centers is saved
  pdf("Cluster centers.pdf")
  heatmap.2(reducedPenalizedClusterCenters, col=colorRampPalette(c("blue", "white", "red"))(100), trace="none")
  dev.off()

	return(dClustResult)

}

