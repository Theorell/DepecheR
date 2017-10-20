removeEmptyVariablesAndAllocatePoints <- function(fMeasureDataSet, clusterCenters, noZero){
  #First, all variables that do not contribute to defining a single cluster is removed
  clusterCentersNoEmptyVariables <- clusterCenters[,which(colSums(clusterCenters)!=0)]
  #After this, the same variables are removed from the fMeasureDataSet
  fMeasureDataSetNoEmptyVariables <- fMeasureDataSet[,which(colSums(clusterCenters)!=0)]
  allocationResult <- allocate_points(fMeasureDataSetNoEmptyVariables, clusterCentersNoEmptyVariables, no_zero=noZero)[[1]]
  return(allocationResult)
}