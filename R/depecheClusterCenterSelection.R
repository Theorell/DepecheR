#This function is used internally in the depeche function.
#allSolutions: this is an object that should be derived from the 
#depechePenaltyOpt function output, which saves all the cluster matrices that
#are generated in the optimization process. 
#selectionDataSet: the name speaks for itself
#k: same meaning as in depeche, i.e. the number of startpoint cluster centers.
depecheClusterCenterSelection <- function(allSolutions, selectionDataSet, k,
                                          nCores){
    # Now, all clusterCenters are used to
    # allocate the selectionDataSet.
    allocationResultList <-
        lapply(allSolutions, function(x) {
            dAllocate(inDataFrame = selectionDataSet, depModel = 
                                   list("clusterCenters" = x, 
                               "logCenterSd" = FALSE))
        })
    
    # Here, the corrected Rand index with
    # each allocationResult as the first
    # vector vector and all the others as
    # individual second vectors is identified
    cl <- makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)
    meanARIList <-
        foreach(i = seq_len(length(allocationResultList))) %dopar%
        mean(vapply(allocationResultList,
                    FUN.VALUE = 0.5, rand_index,
                    inds2 = allocationResultList[[i]], k = k
        ))
    stopCluster(cl)
    meanARIVector <- unlist(meanARIList)
    
    # Now the solution being the most similar
    # to all the others is retrieved
    optimalClusterCenters <- 
        unlist(allSolutions[[which(meanARIVector == 
                                       max(meanARIVector))[1]]])
    colnames(optimalClusterCenters) <- colnames(selectionDataSet)
    # Remove all rows and columns that do not
    # contain any information
    reducedClusterCenters <-
        optimalClusterCenters[
            which(rowSums(optimalClusterCenters) != 0),
            which(colSums(optimalClusterCenters) != 0)
            ]
    
    # In the specific case that only one row
    # is left, due to a high penalty, the
    # data needs to be converted back to a
    # matrix from a vector. The same is done
    # if the number of informative variables
    # is just one.
    if (length(which(rowSums(optimalClusterCenters) != 0)) == 1) {
        reducedClusterCenters <- t(reducedClusterCenters)
    } else if (length(which(colSums(optimalClusterCenters) != 0) == 1)) {
        reducedClusterCenters <- as.matrix(reducedClusterCenters)
    }
    reducedClusterCenters
}