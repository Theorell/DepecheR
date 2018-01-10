#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
#' @importFrom gplots heatmap.2
#' @importFrom dplyr sample_n
#' @export dClustAllData
dClustAllData <- function(inDataFrameScaled, penalty, firstClusterNumber=1, k=20){

  penaltyForRightSize <- penalty*((nrow(inDataFrameScaled)*sqrt(ncol(inDataFrameScaled)))/1450)

  dataMat<-data.matrix(inDataFrameScaled)
  
  #Here the number of iterations is chosen. Very many are not needed, but a few will make the clustering even better than if just one was chosen.
  n_cores <- detectCores() - 1
  #if(n_cores>=7){
  #  if(n_cores<=21){
 #     iterations <- n_cores
 #   } else {
 #     iterations <- 21
 #   }
 # } else {
 #   iterations <- 7
 # }
  
  #This is the central function of the whole package.
  
  cl <-  parallel::makeCluster(n_cores, type = "SOCK")
  registerDoSNOW(cl)
  return_all <- foreach(i=1:21) %dopar% sparse_k_means(dataMat,k,penaltyForRightSize,1, i)
  parallel::stopCluster(cl)	

  #Here, the best iteration is retrieved
  logMaxLik <- as.vector(do.call("rbind", lapply(return_all, "[[", 5)))
  minimumN <- max(logMaxLik)
  returnLowest <- return_all[[which(abs(logMaxLik)==minimumN)[1]]]
  
  #And here, the optimal results, given if an origo cluster should be included or not, are retrieved further
  clusterVector <- returnLowest$i
  clusterCenters <- returnLowest$c

  #And here, the optimal results are made more dense by removing empty rows and columns, etc.
  #Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is not included, the first cluster gets the denomination 1.
  clusterVectorEquidistant <- turnVectorEquidistant(clusterVector, startValue=firstClusterNumber)			
  colnames(clusterCenters) <- colnames(inDataFrameScaled)
    
  #Remove all rows and columns that do not contain any information
  reducedClusterCenters <- clusterCenters[which(rowSums(clusterCenters)!=0),which(colSums(clusterCenters)!=0)]
    
  #In the specific case that only one row is left, due to a high penalty, the data needs to be converted back to a matrix from a vector. The same is done if the number of informative variables is just one.
  if(length(which(rowSums(clusterCenters)!=0))==1){
    reducedClusterCenters <- t(reducedClusterCenters)
  } else if(length(which(colSums(clusterCenters)!=0)==1)){
    reducedClusterCenters <- as.matrix(reducedClusterCenters)
  }
    
  #Make the row names the same as the cluster names in the clusterVectorEquidistant
  rownames(reducedClusterCenters) <- rep(firstClusterNumber:(firstClusterNumber+(nrow(reducedClusterCenters))-1))
    
    
  #Here, the results are combined
  dClustResult <- list(clusterVectorEquidistant, reducedClusterCenters)
  names(dClustResult) <- c("clusterVector", "clusterCenters")
  return(dClustResult)
  
}

