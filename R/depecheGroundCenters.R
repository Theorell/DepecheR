# This function is part of hte depeche function
# Here, a cluster center matrix is created that relates its values to the
# original indata variables, which increases the interpretability vastly.
# "Unique" flags: 
# reducedClusterCenters: the selected cluster centers, where only variables
# that have not been completely sparsed out are kept. 
depecheGroundCenters <- function(reducedClusterCenters, logCenterSd){
    
    # First, the cluster centers are multiplied by the standard deviation of
    # the data
    clusterCentersMultSd <- reducedClusterCenters * logCenterSd[[3]]
    
    # Then, the center term is added to all
    # variables separately
    if (is.logical(logCenterSd[[2]]) == FALSE) {
        
        groundClusterCenters <- 
            do.call("cbind",lapply(seq_len(ncol(clusterCentersMultSd)), 
                                   function(i){
                clusterCentersMultSd[, i] + logCenterSd[[2]][i]
        }))
        groundClusterCenters[which(reducedClusterCenters == 0)] <- 0
    } else {
        groundClusterCenters <- clusterCentersMultSd
    }
    dimnames(groundClusterCenters) <- dimnames(reducedClusterCenters)
    groundClusterCenters
}