#' Function to run penalized K means
#'
#'
#' This function is the core user function of the Depeche package. It clusters the data with a penalized version of K-means.
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @importFrom gplots heatmap.2
#' @importFrom dplyr sample_n
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
#' @param dOptObject This object contains information about optimal sample size, penalty offset, solution with or without a cluster in origo, and the number of initial cluster centers that were used to find this optimal information.
#' @param ids A vector of the same length as rows in the inDataFrameScaled. It is used to generate the final analysis, where a table of the percentage of observations for each individual and each cluster is created.
#' @param sampleSize By default inherited from dOptObject. Number of observations that shoult be included in the initial clustering step. Three possible values. Either inherited, "All" or a user-specified number. Defaults to inheriting from dClustObject. If a dClustObject is not substituted, all rows in inDataFrameScaled are added by default. If another number, a sample is created from inDataFrameScaled. This is extra useful when clustering very large datasets. Replacement is set to TRUE.
#' @param penalty By default inherited from dOptObject. The parameter that controls the level of penalization. 
#' @param withOrigoClust By default inherited from dOptObject. This parameter controls if the generated result should contain a cluster in origo or not. 
#' @param initCenters By default inherited from dOptObject. Number of starting points for clusters. This essentially means that it is the highest possible number of clusters that can be defined. The higher the number, the greater the precision, but the computing time is also increased with the number of starting points. Default is 30.
#' @seealso \code{\link{dAllocate}}, \code{\link{dOpt}}, \code{\link{dOptAndClust}}
#' @return A list with three components:
#' \describe{
#'     \item{clusterVector}{A vector with the same length as number of rows in the inDataFrameScaled, where the cluster identity of each observation is noted.}
#'     \item{clusterCenters}{A matrix containing information about where the centers are in all the variables that contributed to creating the cluster with the given penalty term.}
#'     \item{clusterPercentagesForAllIds}{A matrix showing the percentage of observations for each id in each cluster.}
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
#' #Run the dOptPenalty function to get good starting points
#' x_optim <- dOpt(x_scaled)
#'
#' #Then run the actual function
#' x_dClust <- dClust(x_scaled, dOptObject=x_optim, ids=x[,1])
#'
#' #And finally look at your great result
#' str(x_dClust)
#' @export dClust
dClust <- function(inDataFrameScaled, dOptObject, ids, sampleSize, penalty, withOrigoClust, initCenters=30){

  
  if(missing(ids)){
    stop("Vector of ids is missing. Save youself some time and put it in before running again, as the function will otherwise throw an error at the end.")
  }

  if(exists("ids")==FALSE){
    stop("Ids is a non-existent object. Save youself some time and put in an existing one before running again, as the function will otherwise throw an error at the end.")
  }

  if(missing(dOptObject)==FALSE){
    if(missing(sampleSize)==TRUE){
      sampleSize <- dOptObject[[1]][length(dOptObject[[1]])]
    } else if(sampleSize!="All"){
      sampleSize <- nrow(inDataFrameScaled)
      inDataFrameUsed <- inDataFrameScaled
    }
    penalty <- dOptObject[[4]][1,1]
    withOrigoClust <- dOptObject[[4]][1,2]
    initCenters <- dOptObject[[4]][1,3]
  }
  
  if(sampleSize!=nrow(inDataFrameScaled)){
    inDataFrameUsed <- sample_n(inDataFrameScaled, sampleSize, replace=TRUE)
  }

  k <- ((sampleSize*sqrt(ncol(inDataFrameUsed)))/1450)

  penaltyForRightSize <- penalty*k 

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
    
#This is the most central function of the whole package.
	
	if(Sys.info()['sysname']!="Windows"){
	  cl <- makeCluster(n_cores, type="FORK")
	  return_all <-parLapply(cl,0:iterations,function(x) sparse_k_means(dataMat,initCenters,penaltyForRightSize,1,x))
	  stopCluster(cl)
	} else {
	  cl <- makeCluster(n_cores, type="PSOCK")
	  return_all <-parLapply(cl,0:iterations,function(x) sparse_k_means(dataMat,initCenters,penaltyForRightSize,1,x))
	  stopCluster(cl)
	}
	
  #Here, the best iteration is retrieved
  logMaxLik <- as.vector(do.call("rbind", lapply(return_all, "[[", 5)))
  minimumN <- max(logMaxLik)
  returnLowest <- return_all[[which(abs(logMaxLik)==minimumN)[1]]]

	if(withOrigoClust=="yes"){
			clusterVector <- returnLowest$i

			#Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is included, this cluster gets the denomination 0.
			clusterVectorEquidistant <- turnVectorEquidistant(clusterVector, startValue=0)
			clusterCenters <- returnLowest$c			  
			colnames(clusterCenters) <- colnames(inDataFrameUsed)

			#Remove all columns and rows that do not contain any information.
			reducedClusterCentersColRow <- clusterCenters[which(rowSums(clusterCenters)!=0),which(colSums(clusterCenters)!=0)]
			
			#Here, only the rows that do not contain any information is removed. To be used in dClustPredict. 
			reducedClusterCentersRow <- clusterCenters[which(rowSums(clusterCenters)!=0),]
			
			#Add the zero cluster back. This is not done when there is an origo cluster.
			if(nrow(clusterCenters)>nrow(reducedClusterCentersRow)){
			  reducedClusterCentersColRowOrigo <- clusterCenters[1,which(colSums(clusterCenters)!=0)]
			  reducedClusterCentersColRow <- rbind(reducedClusterCentersColRowOrigo, reducedClusterCentersColRow)
			  reducedClusterCentersRow <- rbind(clusterCenters[1,], reducedClusterCentersRow)
			  }

	}

	if(withOrigoClust=="no"){
			clusterVector <- returnLowest$o
			#Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is not included, the first cluster gets the denomination 1.
			clusterVectorEquidistant <- turnVectorEquidistant(clusterVector)			
			clusterCenters <- returnLowest$v		
			colnames(clusterCenters) <- colnames(inDataFrameUsed)
			#Remove all rows that do not contain any information
			reducedClusterCentersColRow <- clusterCenters[which(rowSums(clusterCenters)!=0),which(colSums(clusterCenters)!=0)]

			#Here, only the rows that do not contain any information is removed. To be used in dClustPredict. 
			reducedClusterCentersRow <- clusterCenters[which(rowSums(clusterCenters)!=0),]
			
	}

	#Make the row names the same as the cluster names in the clusterVectorEquidistant
  rownames(reducedClusterCentersColRow) <- sort(unique(clusterVectorEquidistant))
  rownames(reducedClusterCentersRow) <- sort(unique(clusterVectorEquidistant))
	
	#If the number of observations used for the clustering is not equal to the number of rows in the inDataFrameScaled, here a prediction is made for all the other events. 
	if(sampleSize!=nrow(inDataFrameScaled)){
	  myMat<-data.matrix(inDataFrameScaled, rownames.force = NA)
	  
	  if(withOrigoClust=="yes"){
	    clusterVectorEquidistant <- unlist(allocate_points(myMat,reducedClusterCentersRow,0))
	  }
	  if(withOrigoClust=="no"){
	    clusterVectorEquidistant <- unlist(allocate_points(myMat,reducedClusterCentersRow,1))
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


  dClustResult <- list(clusterVectorEquidistant, as.data.frame.matrix(reducedClusterCentersColRow), reducedClusterCentersRow, as.data.frame.matrix(t(clusterFractionsForAllIds)))

  names(dClustResult) <- c("clusterVector", "clusterCenters", "clusterCentersWZeroVariables", "idClusterFractions")

  #Here, a heatmap over the cluster centers is saved
  pdf("Cluster centers.pdf")
  heatmap.2(reducedClusterCentersColRow, col=colorRampPalette(c("blue", "white", "red"))(100), trace="none")
  dev.off()

	return(dClustResult)

}

