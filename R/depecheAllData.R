#This function is used by depecheCoFunction. 
#It is only used for datasets with less than 10 000 total observations. 
#"Unique" parameters
#inDataFrameScaled: the scaled version of inDataFrame. See initial part of
#depeche function for details. 
#firstClusterNumber: this defines if the first number for the cluster 
#definitions.
#For information on the other parameters, see depeche. 
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom moments kurtosis
#' @importFrom foreach foreach %dopar%
#' @importFrom gplots heatmap.2
#' @importFrom dplyr sample_n
depecheAllData <- function(inDataFrameScaled, 
    penalty, firstClusterNumber = 1, k, nCores) {
    penaltyForRightSize <- 
        penalty * ((nrow(inDataFrameScaled) * 
                        sqrt(ncol(inDataFrameScaled)))/1450)
    
    dataMat <- data.matrix(inDataFrameScaled)
    
    if( nCores=="default"){
        nCores <- floor(detectCores()*0.875) 
        if(nCores>10){ nCores <- 10}
    }
    
    cl <- makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)
    return_all <- foreach(i = seq_len(21), .packages = "DepecheR") %dopar% 
        sparse_k_means(dataMat, round(k * 3), penaltyForRightSize, 1, i)
    stopCluster(cl)
    
    # Here, the best iteration is retrieved,
    # if this is not one with only 1 cluster
    
    logMaxLik <- as.vector(do.call("rbind", lapply(return_all, "[[", 5)))
    nClust <- vapply(return_all, FUN.VALUE = 1, 
                     function(x) sum(rowSums(x[[3]] != 0) != 0))
    logMaxLikNotOne <- logMaxLik[which(nClust > 1)]
    maxN <- max(logMaxLikNotOne)
    returnLowest <- return_all[[which(logMaxLik == maxN)[1]]]
    
    # And here, the optimal results, given if
    # an origo cluster should be included or
    # not, are retrieved further
    clusterVector <- returnLowest$i
    clusterCenters <- returnLowest$c
    
    # And here, the optimal results are made
    # more dense by removing empty rows and
    # columns, etc.  Here, the numbers of the
    # removed clusters are removed as well,
    # and only the remaining clusters are
    # retained and numbered equidistantly.
    originalNumbers <- sort(unique(clusterVector))
    newNumbers <- seq(firstClusterNumber, (length(originalNumbers) + 
                                               (firstClusterNumber - 1)))
    clustVecEquidist <- clusterVector
    for (i in seq_along(originalNumbers)) {
        clustVecEquidist[clustVecEquidist == originalNumbers[i]] <- 
            newNumbers[i]
    }

    colnames(clusterCenters) <- colnames(inDataFrameScaled)
    
    # Remove all rows and columns that do not
    # contain any information
    reducedClusterCenters <- clusterCenters[which(rowSums(clusterCenters) != 0),
                                            which(colSums(clusterCenters) != 0)]
    
    # In the specific case that only one row
    # is left, due to a high penalty, the
    # data needs to be converted back to a
    # matrix from a vector. The same is done
    # if the number of informative variables
    # is just one.
    if (length(which(rowSums(clusterCenters) != 0)) == 1) {
        reducedClusterCenters <- t(reducedClusterCenters)
    } else if (length(which(colSums(clusterCenters) != 0) == 1)) {
        reducedClusterCenters <- as.matrix(reducedClusterCenters)
    }
    
    # Make the row names the same as the
    # cluster names in the
    # clustVecEquidist
    rownames(reducedClusterCenters) <- 
        seq(firstClusterNumber, (firstClusterNumber + 
                                     (nrow(reducedClusterCenters)) - 1))
    
    
    # Here, the results are combined
    dClustResult <- list(clustVecEquidist, reducedClusterCenters)
    names(dClustResult) <- c("clusterVector", "clusterCenters")
    return(dClustResult)
}
