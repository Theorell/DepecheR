removeEmptyVariablesAndAllocatePoints <- function(selectionDataSet, 
                                                  clusterCenters) {
    # First, all variables that do not
    # contribute to defining a single cluster
    # is removed. A specific case, namely
    # that only one variable contains
    # meaningful information, is taken into
    # account
    clusterCentersNoEmptyVariables <- 
        clusterCenters[, which(colSums(clusterCenters) != 0)]
    if (length(which(colSums(clusterCenters) != 0) == 1)) {
        clusterCentersNoEmptyVariables <- as.matrix(clusterCentersNoEmptyVariables)
    }
    # After this, the same variables are
    # removed from the selectionDataSet
    selectionDataSetNoEmptyVariables <- 
        selectionDataSet[, which(colSums(clusterCenters) != 0)]
    if (length(which(colSums(clusterCenters) != 0) == 1)) {
        selectionDataSetNoEmptyVariables <- 
            as.matrix(selectionDataSetNoEmptyVariables)
    }
    allocationResult <- 
        allocate_points(selectionDataSetNoEmptyVariables, 
                        clusterCentersNoEmptyVariables, no_zero = 1)[[1]]
    return(allocationResult)
}
