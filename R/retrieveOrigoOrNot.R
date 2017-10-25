#This function is a support function for dOptSubset, ddOptPenaltyObj och dClustCoFunction that retrieves the right dataset depending on if an origo cluster should be returned or not

retrieveOrigoOrNot <- function(withOrigoClust, clusterVector, clusterCenters, colnamesClusterCenters, k, firstClusterNumber){
  if(withOrigoClust=="yes" && length(unique(clusterVector))<k){
    
    #Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is included, this cluster gets the denomination 0.
    clusterVectorEquidistant <- turnVectorEquidistant(clusterVector, startValue=firstClusterNumber)
    colnames(clusterCenters) <- colnamesClusterCenters
    
    #Remove all columns and rows that do not contain any information.
    reducedClusterCenters <- clusterCenters[which(rowSums(clusterCenters)!=0),which(colSums(clusterCenters)!=0)]
    
    #In the specific case that only one row is left, due to a high penalty, the data needs to be converted back to a matrix from a vector
    if(class(reducedClusterCenters)=="numeric"){
      reducedClusterCenters <- t(reducedClusterCenters)
    }
 
    #Add the origo cluster back. This is not done when there is no sparsity.
    
    reducedClusterCentersOrigo <- clusterCenters[1,which(colSums(clusterCenters)!=0)]
    reducedClusterCenters <- rbind(rep(0, times=ncol(reducedClusterCenters)), reducedClusterCenters)

    #Make the row names the same as the cluster names in the clusterVectorEquidistant
    
    rownames(reducedClusterCenters) <- rep(firstClusterNumber:(firstClusterNumber+nrow(reducedClusterCenters)-1))

  } 
  if(withOrigoClust=="no" || (withOrigoClust=="yes" && length(unique(clusterVector))==k)){
    #Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is not included, the first cluster gets the denomination 1.
    clusterVectorEquidistant <- turnVectorEquidistant(clusterVector, startValue=firstClusterNumber)			
    colnames(clusterCenters) <- colnamesClusterCenters
    #Remove all rows that do not contain any information
    reducedClusterCenters <- clusterCenters[which(rowSums(clusterCenters)!=0),which(colSums(clusterCenters)!=0)]
    
    #In the specific case that only one row is left, due to a high penalty, the data needs to be converted back to a matrix from a vector
    if(class(reducedClusterCenters)=="numeric"){
      reducedClusterCenters <- t(reducedClusterCenters)
    }
    
    #Make the row names the same as the cluster names in the clusterVectorEquidistant
    
    rownames(reducedClusterCenters) <- rep(firstClusterNumber:(firstClusterNumber+(nrow(reducedClusterCenters))-1))

    
  }

  return(list(clusterVectorEquidistant, reducedClusterCenters))
}
