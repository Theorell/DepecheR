removeEmptyVariablesAndAllocatePoints <- function(selectionDataSet, clusterCenters){
  #First, all variables that do not contribute to defining a single cluster is removed
  clusterCentersNoEmptyVariables <- clusterCenters[,which(colSums(clusterCenters)!=0)]
  #After this, the same variables are removed from the selectionDataSet
  selectionDataSetNoEmptyVariables <- selectionDataSet[,which(colSums(clusterCenters)!=0)]
  allocationResult <- allocate_points(selectionDataSetNoEmptyVariables, clusterCentersNoEmptyVariables, no_zero=1)[[1]]
  return(allocationResult)
}