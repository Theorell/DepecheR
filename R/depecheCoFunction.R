#This function is used by depeche
#The reason this level of complexity is needed is to handle both single-run
#depeche and dual-run depeche. In the latter case, this function is called
#multiple times. 
#"Unique" parameters
#inDataFrameScaled: the scaled version of inDataFrame. See initial part of
#depeche function for details. 
#firstClusterNumber: this defines if the first number for the cluster 
#definitions.
#createDirectory: should a directory be created? In practice, FALSE for single
#depeche and TRUE for the second part of dual depeche. 
#For information on the other parameters, see depeche. 
#' @importFrom gplots heatmap.2
#' @importFrom graphics box
#' @importFrom ggplot2 ggplot aes geom_line ggtitle xlab ylab ylim ggsave
#' @importFrom dplyr sample_n
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %do% %dopar%
#' @useDynLib DepecheR
depecheCoFunction <- function(inDataFrameScaled, firstClusterNumber = 1, 
                              directoryName, penalties, sampleSize, 
                              selectionSampleSize, k, minARIImprovement, 
                              optimARI, maxIter, newNumbers, 
                              createDirectory = FALSE, nCores, 
                              createOutput, logCenterSd) {
    
    if (createDirectory == TRUE) { dir.create(directoryName) }
    
    # First, if the dataset is very, very
    # big, a subset of it is used to subset
    # from. Otherwise the system memory
    # needed to just perform the boot
    # strapping becomes so consuming, that
    # the process stalls.
    if (nrow(inDataFrameScaled) > 1e+06) {
        sampleRows <- sample(seq_len(nrow(inDataFrameScaled)), 1e+06)
        inDataFrameUsed <- inDataFrameScaled[sampleRows,]
    } else {
        inDataFrameUsed <- inDataFrameScaled
    }
    
    # Here, the sampleSize is set in cases it
    # is 'default'.
    if (sampleSize == "default") {
        if (nrow(inDataFrameUsed) <= 10000) {
            sampleSize <- nrow(inDataFrameScaled)
        } else {
            sampleSize <- 10000
        }
    }
    
    dOptPenaltyResult <- dOptPenalty(inDataFrameUsed, k = k, maxIter = maxIter,
                                     sampleSize = sampleSize, 
                                     penalties = penalties, 
                                     createOutput = createOutput, 
                                     minARIImprovement = minARIImprovement, 
                                     optimARI = optimARI, nCores=nCores, 
                                     createDirectory=createDirectory, 
                                     directoryName=directoryName)
    
    # Now over to creating the final solution
    # Here, the selectionDataSet is created
    if (selectionSampleSize == "default") {
        selectionSampleSize <- sampleSize
    }
    if (nrow(inDataFrameUsed) <= selectionSampleSize) {
        selectionDataSet <- inDataFrameUsed
    } else {
        selectionDataSet <- 
            inDataFrameUsed[sample(seq_len(nrow(inDataFrameUsed)), 
            selectionSampleSize, replace = TRUE), ]
    }
    
    # If the dataset is small, a new set of
    # seven clusterings are performed (on all
    # the data or on a subsample, depending
    # on the sample size), and the maximum
    # likelihood solution is returned as the
    # result
    if (nrow(inDataFrameUsed) < 10000) {
        penalty <- dOptPenaltyResult[[1]][1, 1]
        depecheAllDataResult <- 
            depecheAllData(inDataFrameUsed, penalty = penalty, k = k, 
                           firstClusterNumber = firstClusterNumber,
                           nCores=nCores)
        clusterVectorEquidistant <- depecheAllDataResult[[1]]
        reducedClusterCenters <- depecheAllDataResult[[2]]
    } else {
        # Now, the best run amongst all the runs
        # with the largest sample size is
        # defined, by identifying the solution
        # that gives the highest mean f-measure
        # for all the others.  First, all
        # solutions are retrieved
        allSolutions <- dOptPenaltyResult[[3]]
        
        # Now, all clusterCenters are used to
        # allocate the selectionDataSet.
        allocationResultList <- list()
        selectionDataSetMatrix <- data.matrix(selectionDataSet)
        allocationResultList <- 
            foreach(i = seq_len(length(allSolutions))) %dopar% 
            dAllocate(inDataMatrix = selectionDataSetMatrix, 
                      clusterCenters = allSolutions[[i]], log2Off=TRUE, 
                      noZeroNum==FALSE)
        
        # Here, the corrected Rand index with
        # each allocationResult as the first
        # vector vector and all the others as
        # individual second vectors is identified
        if( nCores=="default"){
            nCores <- floor(detectCores()*0.875) 
            if(nCores>10){
                nCores <- 10
            }
        }
        cl <- makeCluster(nCores, type = "SOCK")
        registerDoSNOW(cl)
        meanARIList <- 
            foreach(i = seq_len(length(allocationResultList))) %dopar% 
            mean(vapply(allocationResultList, FUN.VALUE = 0.5, rand_index, 
                inds2 = allocationResultList[[i]], k = k))
        stopCluster(cl)
        meanARIVector <- unlist(meanARIList)
        
        # Now the solution being the most similar
        # to all the others is retrieved
        optimalClusterCenters <- unlist(allSolutions[[which(meanARIVector == 
            max(meanARIVector))[1]]])
        colnames(optimalClusterCenters) <- colnames(inDataFrameUsed)
        # Remove all rows and columns that do not
        # contain any information
        reducedClusterCenters <- 
            optimalClusterCenters[which(rowSums(optimalClusterCenters) != 0), 
                                  which(colSums(optimalClusterCenters) != 0)]
        
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
        
        # Make the row names nice
        rownames(reducedClusterCenters) <- 
            seq(firstClusterNumber,(firstClusterNumber + 
                                        (nrow(reducedClusterCenters)) - 1))
        
        # And here, the optimal solution is
        # created with the full dataset
        clusterVectorEquidistant <- dAllocate(inDataFrameScaled, 
            reducedClusterCenters, log2Off=TRUE, noZeroNum==FALSE)
    }
    
    # Here, the optPenalty information is
    # retrieved from the optimal sample size
    # run.
    optPenalty <- list(dOptPenaltyResult[[1]], dOptPenaltyResult[[2]])
    
    # Here, a cluster center matrix is
    # created that relates its values to the
    # original, indata variables, which
    # increases the interpretability vastly.
    # First, the cluster centers are
    # multiplied by the standar deviation of
    # the data
    clusterCentersMultSd <- reducedClusterCenters * logCenterSd[[3]]
    
    # Then, the center term is added to all
    # variables separately
    if (is.logical(logCenterSd[[2]]) == FALSE) {
        correctClusterCentersList <- list()
        for (i in seq_len(ncol(clusterCentersMultSd))) {
            focusMean <- logCenterSd[[2]][which(names(logCenterSd[[2]]) %in% 
                colnames(clusterCentersMultSd)[i])]
            correctClusterCentersList[[i]] <- clusterCentersMultSd[,i] 
            + focusMean
        }
        correctClusterCenters <- do.call("cbind", correctClusterCentersList)
        colnames(correctClusterCenters) <- colnames(reducedClusterCenters)
    } else {
        correctClusterCenters <- clusterCentersMultSd
        colnames(correctClusterCenters) <- colnames(reducedClusterCenters)
    }
    
    # Here, a sparsity matrix is generated,
    # as the sparsed out variables will seem
    # like they are situated in the middle,
    # when they are in fact pushed to the
    # center of the variable in question.
    sparsityMatrix <- reducedClusterCenters
    sparsityMatrix[sparsityMatrix != 0] <- 1
    
    # Here, a heatmap over the cluster
    # centers is saved. Only true if the
    # number of clusters exceeds one.
    if (ncol(reducedClusterCenters) > 500) {
        message("As the number of variables in the result exceeds 500, 
                it is not meaningful to produce a cluster center heatmap, 
                so it is omitted")
    } else if (nrow(reducedClusterCenters) > 1 && 
               ncol(reducedClusterCenters) > 1) {
        graphicClusterCenters <- reducedClusterCenters
        # Here we scale each center value to the
        # range between the lowest and highest
        # permille of the observations in the
        # inDataScaled for that variable
        for (i in seq_len(ncol(graphicClusterCenters))) {
            scaledFocus <- 
                inDataFrameUsed[, colnames(inDataFrameUsed) == 
                                    colnames(graphicClusterCenters)[i]]
            graphicClusterCenters[, i] <- dScale(reducedClusterCenters[, i], 
                                                 scaledFocus, 
                                                 robustVarScale = FALSE, 
                                                 center = FALSE, 
                                                 truncate = TRUE)
        }
        
        graphicClusterCenters[sparsityMatrix == 0] <- NA
        
        colorLadder <- dColorVector(seq_len(11), 
            colorScale = c("#0D0887FF", "#6A00A8FF", "#900DA4FF", "#B12A90FF", 
                           "#CC4678FF", "#E16462FF", "#F1844BFF", "#FCA636FF", 
                           "#FCCE25FF"))
        
        if (createOutput == TRUE) {
            #First, the name is defined, depending on two criteria: if the data
            #was log transformed, and if a directory should be created or not. 
            logOrNoLogName <- ifelse (logCenterSd[1] == FALSE, "Cluster",
                                          "Log2_transformed_cluster")
            clusterCenterName <- ifelse (createDirectory == TRUE, 
                                         file.path(directoryName, 
                                                   paste0(logOrNoLogName,
                                                   "_centers.pdf")), 
                                         paste0(logOrNoLogName, "_centers.pdf"))
            pdf(clusterCenterName)
            heatmap.2(as.matrix(graphicClusterCenters), Rowv = FALSE, 
                      Colv = FALSE, dendrogram = "none", scale = "none", 
                      col = colorLadder, breaks = seq(0, 1, length.out = 12), 
                      trace = "none", keysize = 1.5, density.info = "none", 
                      key.xlab = 
                          "0=no expression, 1=high expression\ngrey=penalized", 
                      na.color = "#A2A2A2")
            dev.off()
        }
    }
    
    # And now, the result that should be
    # returned is compiled
    depecheResult <- list(clusterVectorEquidistant, correctClusterCenters, 
                          sparsityMatrix, optPenalty)
    names(depecheResult) <- c("clusterVector", "clusterCenters", 
                              "sparsityMatrix", "penaltyOptList")
    if(logCenterSd[1] == TRUE){
        names(depecheResult)[2] <- "log2ClusterCenters"
        }

    return(depecheResult)
}
