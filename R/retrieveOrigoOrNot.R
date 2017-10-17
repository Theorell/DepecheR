#This function is a support function for dOptSubset, ddOptPenaltyObj och dClustCoFunction that retrieves the right dataset depending on if an origo cluster should be returned or not

retrieveOrigoOrNot <- function(withOrigoClust, clusterVector, clusterCenters, colnamesClusterCenters, k){
  if(withOrigoClust=="yes" && length(unique(clusterVector))<k){
    
    #Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is included, this cluster gets the denomination 0.
    clusterVectorEquidistant <- turnVectorEquidistant(clusterVector, startValue=0)
    colnames(clusterCenters) <- colnamesClusterCenters
    
    #Remove all columns and rows that do not contain any information.
    reducedClusterCentersColRow <- clusterCenters[which(rowSums(clusterCenters)!=0),which(colSums(clusterCenters)!=0)]
    
    #Here, only the rows that do not contain any information is removed. To be used in dClustPredict. 
    reducedClusterCentersRow <- clusterCenters[which(rowSums(clusterCenters)!=0),]
    
    #In the specific case that only one row is left, due to a high penalty, the data needs to be converted back to a matrix from a vector
    if(class(reducedClusterCentersColRow)=="numeric"){
      reducedClusterCentersColRow <- t(reducedClusterCentersColRow)
    }
    if(class(reducedClusterCentersRow)=="numeric"){
      reducedClusterCentersRow <- t(reducedClusterCentersRow)
    }
    
    #Add the origo cluster back. This is not done when there is no sparsity.
    
    reducedClusterCentersColRowOrigo <- clusterCenters[1,which(colSums(clusterCenters)!=0)]
    reducedClusterCentersColRow <- rbind(rep(0, times=ncol(reducedClusterCentersColRow)), reducedClusterCentersColRow)
    reducedClusterCentersRow <- rbind(rep(0, times=ncol(reducedClusterCentersRow)), reducedClusterCentersRow)
    
    #Make the row names the same as the cluster names in the clusterVectorEquidistant
    
    rownames(reducedClusterCentersColRow) <- rep(0:(nrow(reducedClusterCentersColRow)-1))
    rownames(reducedClusterCentersRow) <- rep(0:(nrow(reducedClusterCentersColRow)-1))
    
  } 
  if(withOrigoClust=="no" || (withOrigoClust=="yes" && length(unique(clusterVector))==k)){
    #Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is not included, the first cluster gets the denomination 1.
    clusterVectorEquidistant <- turnVectorEquidistant(clusterVector)			
    colnames(clusterCenters) <- colnamesClusterCenters
    #Remove all rows that do not contain any information
    reducedClusterCentersColRow <- clusterCenters[which(rowSums(clusterCenters)!=0),which(colSums(clusterCenters)!=0)]
    
    #Here, only the rows that do not contain any information is removed. To be used in dClustPredict. 
    reducedClusterCentersRow <- clusterCenters[which(rowSums(clusterCenters)!=0),]
    
    #In the specific case that only one row is left, due to a high penalty, the data needs to be converted back to a matrix from a vector
    if(class(reducedClusterCentersColRow)=="numeric"){
      reducedClusterCentersColRow <- t(reducedClusterCentersColRow)
    }
    if(class(reducedClusterCentersRow)=="numeric"){
      reducedClusterCentersRow <- t(reducedClusterCentersRow)
    }
    
    #Make the row names the same as the cluster names in the clusterVectorEquidistant
    
    rownames(reducedClusterCentersColRow) <- rep(1:(nrow(reducedClusterCentersColRow)))
    rownames(reducedClusterCentersRow) <- rep(1:(nrow(reducedClusterCentersColRow)))
    
    
  }

  return(list(clusterVectorEquidistant, reducedClusterCentersColRow, reducedClusterCentersRow))
}
