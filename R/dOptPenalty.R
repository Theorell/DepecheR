#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom Rcpp evalCpp
#' @importFrom graphics box

dOptPenalty <- function(inDataFrameScaled, 
    k, maxIter, minARIImprovement, sampleSize, 
    penalties, makeGraph, disableWarnings = FALSE, 
    optimARI, nCores=nCores) {
    
    # The constant k is empirically
    # identified by running a large number of
    # penalty values for a few datasets.
    penaltyConstant <- ((sampleSize * sqrt(ncol(inDataFrameScaled)))/1450)
    realPenalties <- penalties * penaltyConstant
    roundPenalties <- round(penalties, digits = 1)
    
    if( nCores=="default"){
        nCores <- floor(detectCores()*0.875) 
        if(nCores>10){
            nCores <- 10
        }
    }
    
    dataMat <- data.matrix(inDataFrameScaled, 
        rownames.force = NA)
    
    # This command is reiterated the number
    # of times that is needed to reach a
    # minimal improvement of the distance.
    iter <- 1
    std <- 1
    # distanceBetweenMinAndBestPrevious=-1
    iterTimesNCores <- 1
    allClusterCentersPenaltySorted <- list()
    cl <- makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)
    
    interestingPenalties <- realPenalties
    usedPositions <- seq_len(length(realPenalties))
    
    while (iterTimesNCores < 20 || (iterTimesNCores < 
        maxIter && (std >= minARIImprovement || 
        distanceBetweenMaxAndSecond > 0))) {
        ptm <- proc.time()
        optimList <- foreach(i = seq_len(nCores), 
            .packages = "DepecheR") %dopar% 
            grid_search(dataMat, k, interestingPenalties, 
                1, sampleSize, i)
        
        # Before any further analyses are
        # performed, any penalty that can result
        # in a trivial solution are practically
        # eliminated.
        optimListNonTrivial <- optimList
        for (i in seq_len(length(optimListNonTrivial))) {
            optimListNonTrivial[[i]]$d[which(optimList[[i]]$n == 
                1)] <- 0
        }
        
        # Now, the new list is combined with the
        # older, if there are any
        if (iter == 1) {
            optimListFull <- optimListNonTrivial
        } else {
            # First, the excluded values are included
            # as NAs for each of the iterations in
            # the chunk
            for (i in seq_len(length(optimListNonTrivial))) {
                d <- rep("NA", length(optimListFull[[1]][[1]]))
                d[usedPositions] <- optimListNonTrivial[[i]]$d
                optimListNonTrivial[[i]]$d <- d
                z <- rep("NA", length(optimListFull[[1]][[1]]))
                z[usedPositions] <- optimListNonTrivial[[i]]$z
                optimListNonTrivial[[i]]$z <- z
                n <- rep("NA", length(optimListFull[[1]][[1]]))
                n[usedPositions] <- optimListNonTrivial[[i]]$n
                optimListNonTrivial[[i]]$n <- n
                m <- rep("NA", length(optimListFull[[1]][[1]]))
                m[usedPositions] <- optimListNonTrivial[[i]]$m
                optimListNonTrivial[[i]]$m <- m
                # The same is done for all the matrices,
                # to make the lists the right length
                mockCenterMatrix <- matrix("NA", 
                  nrow = nrow(optimListFull[[1]][[5]][[1]][[1]]), 
                  ncol = ncol(optimListFull[[1]][[5]][[1]][[1]]))
                mockCenterMatrixList <- list(mockCenterMatrix, 
                  mockCenterMatrix, mockCenterMatrix, 
                  mockCenterMatrix)
                mockCenterList <- rep(list(mockCenterMatrixList), 
                  times = length(optimListFull[[1]][[5]]))
                mockCenterList[usedPositions] <- optimListNonTrivial[[i]][[5]]
                optimListNonTrivial[[i]][[5]] <- mockCenterList
            }
            
            optimListFull <- c(optimListFull, 
                optimListNonTrivial)
        }
        
        # Here, the average and standard
        # deviation of the error of the mean (or
        # something like that) is retrieved
        meanOptimList <- list()
        stdOptimList <- list()
        for (i in seq_len(4)) {
            x <- do.call("rbind", lapply(optimListFull, 
                "[[", i))
            xNumeric <- suppressWarnings(apply(x, 
                2, as.numeric))
            meanOptimList[[i]] <- apply(xNumeric, 
                2, mean, na.rm = TRUE)
            stdOptimList[[i]] <- apply(xNumeric, 
                2, sd, na.rm = TRUE)
        }
        
        stdOptimDf <- (as.data.frame(stdOptimList))/(sqrt(iter * 
            nCores))
        meanOptimDf <- data.frame(meanOptimList[[1]], 
            meanOptimList[[3]])
        
        # Turn these into vectors
        meanOptimVector <- meanOptimDf[, 
            1]
        stdOptimVector <- stdOptimDf[, 1]
        # Return the position of the mazimum
        # value
        maxPos <- which.max(meanOptimVector)[1]
        
        
        if (length(interestingPenalties) > 
            1) {
            # Add the standard deviation of this
            # position to its mean
            meanMinus2StdMax <- meanOptimVector[maxPos] - 
                (2 * stdOptimVector[maxPos])
            
            # Now add the standard deviation of each
            # of the non-minimal values to the mean
            
            meanPlus2StdAll <- meanOptimVector + 
                (2 * stdOptimVector)
            
            # Identify the highest value among these
            maxMeanPlus2StdAllNonMax <- max(meanPlus2StdAll[-maxPos])
            
            # Now, the distance between
            # minMeanMinusSdAllNonMin and the lowest
            # value. If they overlap, the iteration
            # has not made it totally clear which
            # point is optimal.
            distanceBetweenMaxAndSecond <- meanMinus2StdMax - 
                maxMeanPlus2StdAllNonMax
            
            # And now, the interesting positions are
            # defined.  These are the ones that
            # either overlap with uncertainity with
            # the optimal solution, or that are very
            # similar to it. In the first iterations,
            # the critera are more inclusive.
            if (iter == 1) {
                usedPositions <- which(realPenalties %in% 
                  interestingPenalties & 
                  (meanPlus2StdAll >= (meanMinus2StdMax - 
                    (2 * stdOptimVector[maxPos])) | 
                    meanOptimDf[, 1] >= max(meanOptimDf[, 
                      1]) - ((1 - optimARI) * 
                      4)))
            } else if (iter == 2) {
                usedPositions <- which(realPenalties %in% 
                  interestingPenalties & 
                  (meanPlus2StdAll >= (meanMinus2StdMax - 
                    stdOptimVector[maxPos]) | 
                    meanOptimDf[, 1] >= max(meanOptimDf[, 
                      1]) - ((1 - optimARI) * 
                      3)))
            } else {
                usedPositions <- which(realPenalties %in% 
                  interestingPenalties & 
                  (meanPlus2StdAll >= meanMinus2StdMax | 
                    meanOptimDf[, 1] >= max(meanOptimDf[, 
                      1]) - ((1 - optimARI) * 
                      2)))
            }
            
            # Here, all penalties and solutions that
            # do not overlap with the optimal
            # solution are excluded from further
            # optimiations, to reduce calculation
            # time.
            interestingPenalties <- realPenalties[usedPositions]
        } else {
            distanceBetweenMaxAndSecond <- 1
        }
        
        # Finally, another criterion on the gain
        # of adding more rows is included
        std <- stdOptimVector[maxPos]
        iterTimesNCores <- iter * nCores
        
        # Here, the cluster center information
        # for each run is saved: First all
        # cluster center information is saved in
        # one place.
        funval <- optimListNonTrivial[[1]][5]
        allClusterCenters <- vapply(optimListNonTrivial, 
            FUN.VALUE = funval, "[", 5)
        
        # Then each penalty and the solutions are
        # reorganized and, if iter>1, integrated
        # with previous runs.
        for (i in seq_len(length(allClusterCenters[[1]]))) {
            funval <- allClusterCenters[[1]][i]
            tempPenaltyList <- vapply(allClusterCenters, 
                FUN.VALUE = funval, "[", 
                i)
            funval1 <- allClusterCenters[[1]][2]
            funval2 <- allClusterCenters[[1]][2]
            tempClusterCenters <- c(vapply(tempPenaltyList, 
                FUN.VALUE = funval1, "[", 
                1), vapply(tempPenaltyList, 
                FUN.VALUE = funval2, "[", 
                2))
            
            if (iter == 1) {
                allClusterCentersPenaltySorted[[i]] <- tempClusterCenters
            } else {
                allClusterCentersPenaltySorted[[i]] <- c(allClusterCentersPenaltySorted[[i]], 
                  tempClusterCenters)
            }
        }
        
        fullTime <- proc.time() - ptm
        print(paste("Set ", iter, " with ", 
            nCores, " iterations completed in ", 
            round(fullTime[3]), " seconds.", 
            sep = ""))
        iter <- iter + 1
    }
    
    stopCluster(cl)
    
    print(paste("The optimization was iterated ", 
        (iter - 1) * nCores, " times.", 
        sep = ""))
    
    if (iter * nCores >= maxIter && std > 
        minARIImprovement) {
        warning("The maximum number of iterations was reached before stable optimal solution was found")
    }
    
    rownames(meanOptimDf) <- roundPenalties
    colnames(meanOptimDf) <- c("ARI", "nClust")
    
    # Here, the optimal penalty is selected.
    # This is defined as the lowest penalty
    # that yields an ARI that is not lower
    # than 0.01 less than the best ARI.
    optimalPenalties <- roundPenalties[which(meanOptimDf[, 
        1] >= max(meanOptimDf[, 1]) - (1 - 
        optimARI))]
    
    # The best penalty is defined as the
    # median penalty or the optimal penalty
    # with the even number among the two most
    # centrally placed, in the case of an
    # even number of optimal solutions.
    penaltyOpt.df <- data.frame(bestPenalty = optimalPenalties[round(mean(c(seq_len(length(optimalPenalties)))))], 
        k)
    
    lowestPenalty <- roundPenalties[1]
    highestPenalty <- roundPenalties[length(roundPenalties)]
    
    if (length(penalties) > 1) {
        if (penaltyOpt.df$bestPenalty == 
            lowestPenalty) {
            warning("The lowest penalty was the most optimal in the range. This might be either due to using a to small samle size, or because the penalty range is suboptimal. ")
        }
        if (penaltyOpt.df$bestPenalty == 
            highestPenalty) {
            warning("The highest penalty was the most optimal in the range. This might be either due to using a to small samle size, or because the penalty range is suboptimal. ")
        }
    }
    
    # Here, the optimization is plotted if
    # wanted.
    if (makeGraph == TRUE) {
        pdf("Distance_as_a_function_of_penalty_values.pdf")
        par(mar = c(5, 4, 4, 6) + 0.1)
        # Plot the data
        plot(log10(roundPenalties), meanOptimDf[, 
            1], pch = 16, axes = FALSE, ylim = c(0, 
            1), xlab = "", ylab = "", type = "b", 
            col = "black", main = "Distance between bootstraps as a function of penalties values")
        axis(2, ylim = c(0, 1), col = "black", 
            las = 1)  ## las=1 makes horizontal labels
        mtext("Adjusted rand index (ARI) between data resamplings", 
            side = 2, line = 2.5)
        graphics::box()
        
        # Draw the penalty axis
        axis(1, pretty(range(log10(roundPenalties)), 
            n = 10))
        mtext("Log10 of penalty values", 
            side = 1, col = "black", line = 2.5)
        
        # Add Legend
        legend("topleft", legend = "ARI (high is good)", 
            text.col = "black", pch = c(16, 
                15), col = "black")
        
        dev.off()
    }
    
    # Here, the list of solutions with the
    # best penalty is exported
    bestClusterCenters <- allClusterCentersPenaltySorted[[which(roundPenalties == 
        penaltyOpt.df$bestPenalty)]]
    
    penaltyOptList <- list(penaltyOpt.df, 
        meanOptimDf, bestClusterCenters)
    
    
    return(penaltyOptList)
}
