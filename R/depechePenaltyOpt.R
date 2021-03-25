# This function is used by depecheCoFunction, and thus indirectly by depeche.
# Here, the optimal penalty is defined.
# "Unique" parameters
# inDataFrameUsed: the scaled version of inDataFrame. See initial part of
# depeche function for details.
# disableWarnings: Should warnings be disabled?
# For information on the other parameters, see depeche.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom Rcpp evalCpp
#' @importFrom graphics box

depechePenaltyOpt <- function(inDataFrameUsed, k, maxIter, minARIImprovement,
                        sampleSize, penalties, createOutput,
                        disableWarnings = FALSE, optimARI, nCores,
                        plotDir) {

    # The constant 1450 was empirically
    # identified by running a large number of
    # penalty values for a few test datasets.
    penaltyConstant <- ((sampleSize * sqrt(ncol(inDataFrameUsed))) / 1450)
    realPenalties <- penalties * penaltyConstant
    roundPenalties <- round(penalties, digits = 1)

    dataMat <- data.matrix(inDataFrameUsed, rownames.force = NA)

    # This command is reiterated the number
    # of times that is needed to reach a
    # minimal improvement of the distance.
    iter <- 1
    lastStd <- 1
    stdChange <- 1
    # distanceBetweenMinAndBestPrevious=-1
    iterTimesNCores <- 1
    allClusterCentersPenaltySorted <- list()
    cl <- makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)

    usedPenalties <- realPenalties
    usedRoundPenalties <- roundPenalties
    allPositions <- seq_along(realPenalties)
    usedPositions <- allPositions

    # That the while loop should never terminate before 20 iterations is
    # arbitrary. 
    while ((iterTimesNCores < 20) || (iterTimesNCores < maxIter && 
           stdChange >= minARIImprovement)) {
        ptm <- proc.time()
        optimList <- foreach(
            i = seq_len(nCores),
            .packages = "DepecheR"
        ) %dopar%
            grid_search(dataMat, k, usedPenalties, 1, sampleSize, i)[
                c(2,4,5)] # This exclusion is due to a previous version keeping
        # a cluster in origo, which has since been deprecated.
        
        # Here, we name the output, and further remove half of the data which
        # is solely there due to the old origo solution thing. We also turn the
        # matrices into vectors. 
        optimList <- lapply(optimList, function(x){
            names(x) <- c("ARI", "Clusters", "ClustCenters")
            x[c(1,2)] <- lapply(x[c(1,2)], as.vector)
            x <- lapply(x, 'names<-', usedRoundPenalties)
            x[[3]] <- lapply(x[[3]], function(y){
                y[c(1,2)]
            })
            x
        })

        # Before any further analyses are
        # performed, any penalty that can result
        # in a trivial solution are practically
        # eliminated.
        optimListNonTrivial <- lapply(optimList, function(x){
            x$ARI[which(x$Clusters == 1)] <- 0
            x
        })
        
        # Now, the new list is combined with the older, if there are any
        if (iter == 1) {
            optimListFull <- optimListNonTrivial
        } else {
            optimListFull <- c(optimListFull, optimListNonTrivial)
        }

        # Here, the average and standard deviation of the ARI for each penalty 
        # is retrieved
        ARIList <- lapply(optimListFull, "[[", 1)
        
        ARIVecList <- lapply(as.character(roundPenalties), function(x){
            unlist(lapply(ARIList, function(y){
                if(any(names(y) %in% x)){
                    return(y[which(names(y) == x)])
                }
            }))
        })
        
        meanARIVec <- unlist(lapply(ARIVecList, mean))
        stdOptVec <- unlist(lapply(ARIVecList, sd))
        
        # Return the position of the maximum value
        maxPos <- which.max(meanARIVec)[1]
        
        # Add the standard deviation of this
        # position to its mean
        meanMinus2StdMax <- meanARIVec[maxPos] -
            (2 * stdOptVec[maxPos])
        
        # Now add the standard deviation of each
        # of the non-minimal values to the mean
        
        meanPlus2StdAll <- meanARIVec +
            (2 * stdOptVec)
        
        # And now, the interesting positions are
        # defined. These are the ones that
        # either overlap with uncertainity with
        # the optimal solution, or that are very
        # similar to it. In the first iterations,
        # the critera are more inclusive, to avoid
        # discarding solutions due to coincidence.
        localIter <- ifelse(iter < 3, iter, 3)
        k1 <- switch(EXPR = localIter, 2, 1, 0)
        k2 <- switch(EXPR = localIter, 4, 3, 2)
        usedPositions <-
            which(realPenalties %in% usedPenalties &
                      (meanPlus2StdAll >=
                           (meanMinus2StdMax -
                                (k1 * stdOptVec[maxPos])) |
                           meanARIVec >= max(meanARIVec) -
                           ((1 - optimARI) * k2)))
        
        # Here, all penalties and solutions that
        # do not overlap with the optimal
        # solution are excluded from further
        # optimiations, to reduce calculation
        # time.
        usedPenalties <- realPenalties[usedPositions]
        usedRoundPenalties <- roundPenalties[usedPositions]

        # Finally, another criterion on the gain
        # of adding more iterations is included. Here, we also kill the
        #optimization if we are running only one penalty. 
        if (length(usedPenalties) > 1) {
            if(lastStd != 0){
                stdChange <- ((lastStd - stdOptVec[maxPos])/lastStd)/
                    sqrt(iter * nCores)
            } else {
                stdChange <- 0
            }
        } else {
            stdChange <- 0
        }
        
        lastStd <- stdOptVec[maxPos]
        
        iterTimesNCores <- iter * nCores


        fullTime <- proc.time() - ptm
        message(
            "Set ", iter, " with ", nCores, " iterations completed in ",
            round(fullTime[3]), " seconds."
        )
        iter <- iter + 1
    }

    stopCluster(cl)

    message("The optimization was iterated ", (iter - 1) * nCores, " times.")

    if (iter * nCores >= maxIter && stdChange >
        minARIImprovement) {
        warning("The maximum number of iterations was reached before stable 
                optimal solution was found")
    }


    # Here, the optimal penalty is selected.
    # This is defined as the lowest penalty
    # that yields an ARI that is not lower
    # than 0.01 less than the best ARI.
    optimalPenalties <-
        roundPenalties[which(meanARIVec >= max(meanARIVec) - (1 - optimARI))]
    
    # The best penalty is defined as the
    # median penalty or the optimal penalty
    # with the even number among the two most
    # centrally placed, in the case of an
    # even number of optimal solutions.
    bestPenalty <- optimalPenalties[round(mean(c(seq_along(optimalPenalties))))]
    lowestPenalty <- roundPenalties[1]
    highestPenalty <- roundPenalties[length(roundPenalties)]

    if (length(penalties) > 1) {
        if (bestPenalty == lowestPenalty) {
            warning("The lowest penalty was the most optimal in the range. ",
                    "This might either be due a suboptimal set of penalties,",
                    " or to the use of a too small sample size")
        }
        if (bestPenalty == highestPenalty) {
            warning("The highest penalty was the most optimal in the range. ",
                    "This might either be due a suboptimal set of penalties,",
                    " or to the use of a too small sample size")
        }
    }

    # Here, the optimization is plotted if wanted.
    
    if (createOutput) {
        pdf(file.path(plotDir, "ARI_as_a_function_of_penalties.pdf"))
        par(mar = c(5, 4, 4, 6) + 0.1)
        # Plot the data
        plot(log10(roundPenalties), meanARIVec,
            pch = 16, axes = FALSE,
            ylim = c(0, 1), xlab = "", ylab = "", type = "b", col = "black",
            main = "Difference between samples as a function of penalties"
        )
        axis(2, ylim = c(0, 1), col = "black", las = 1)
        mtext("Adjusted rand index (ARI) between data resamplings",
            side = 2,
            line = 2.5
        )
        graphics::box()

        # Draw the penalty axis
        axis(1, pretty(range(log10(roundPenalties)), n = 10))
        mtext("Log10 of penalty values", side = 1, col = "black", line = 2.5)

        # Add Legend
        legend("topleft",
            legend = "ARI (high is good)", text.col = "black",
            pch = c(16, 15), col = "black"
        )

        dev.off()
    }

    # Now, we create the final output
    #The average number of clusters: 
    nClustList <- lapply(optimListFull, "[[", 2)
    meanClustVec <- unlist(lapply(as.character(roundPenalties), function(x){
        mean(unlist(lapply(nClustList, function(y){
            if(any(names(y) %in% x)){
                return(y[which(names(y) == x)])
            }
        })), na.rm = TRUE)
    }))
    
    meanOptimDf <- data.frame("ARI" = meanARIVec, "nClust" = meanClustVec)
    row.names(meanOptimDf) <- roundPenalties
    # Here, the cluster center information
    # for each run is saved.
    allClustCentRaw <- lapply(optimListFull, "[[", 3)

    bestClusterCenters <- unlist(lapply(allClustCentRaw, "[[", 
                                 as.character(bestPenalty)), recursive = FALSE)

    penaltyOptList <- list("bestPenalty" = bestPenalty, 
                           meanOptimDf, bestClusterCenters)


    return(penaltyOptList)
}
