dSplsdaPreCalculations <- function(clusterVector, 
    idsVector, groupVector, pairingVector, groupName1, groupName2) {
    # Here, the statistical evaluation is
    # performed.  First, the data is divided
    # into each group.
    clusterVectorGroup1 <- clusterVector[groupVector == unique(groupVector)[1]]
    clusterVectorGroup2 <- clusterVector[groupVector == unique(groupVector)[2]]
    idsVectorGroup1 <- 
        as.character(idsVector[groupVector == unique(groupVector)[1]])
    idsVectorGroup2 <- 
        as.character(idsVector[groupVector == unique(groupVector)[2]])
    
    # Now, a table with the percentage of
    # cells in each cluster for each
    # individual is created for both groups
    
    clusterTable1 <- 
        as.matrix(as.data.frame.matrix(table(clusterVectorGroup1, 
                                             idsVectorGroup1)))
    clusterTable2 <- 
        as.matrix(as.data.frame.matrix(table(clusterVectorGroup2, 
                                             idsVectorGroup2)))
    
    # In the very unlikely event that there
    # is not a single observation for one
    # cluster from one of the groups, this
    # cluster is substituted to that table
    # with a row of zeros.
    if (nrow(clusterTable1) < length(unique(clusterVector))) {
        zeroMat <- matrix(data = 0, nrow = length(unique(clusterVector)) - 
                              nrow(clusterTable1), ncol = ncol(clusterTable1))
        
        colnames(zeroMat) <- colnames(clusterTable1)
        allRowNames <- as.character(sort(unique(clusterVector)))
        row.names(zeroMat) <- 
            allRowNames[-which(allRowNames %in% row.names(clusterTable1))]
        # Here, rows are added to the cluster
        # table to make the number of rows the
        # same as the unique values of the
        # cluster vector.
        clusterTable1big <- rbind(clusterTable1, zeroMat)
        
        # The rows of the table are re-sorted
        clusterTable1bigResorted <- 
            clusterTable1big[order(as.numeric(row.names(clusterTable1big))), ]
        clusterTable1 <- clusterTable1bigResorted
    }
    
    # And the same procedure is done for the
    # second group
    if (nrow(clusterTable2) < length(unique(clusterVector))) {
        zeroMat <- matrix(data = 0, nrow = length(unique(clusterVector)) - 
            nrow(clusterTable2), ncol = ncol(clusterTable2))
        
        colnames(zeroMat) <- colnames(clusterTable2)
        allRowNames <- as.character(sort(unique(clusterVector)))
        row.names(zeroMat) <- allRowNames[-which(allRowNames %in% 
            row.names(clusterTable2))]
        # Here, rows are added to the cluster
        # table to make the number of rows the
        # same as the unique values of the
        # cluster vector.
        clusterTable2big <- rbind(clusterTable2, zeroMat)
        
        # The rows of the table are re-sorted
        clusterTable2bigResorted <- 
            clusterTable2big[order(as.numeric(row.names(clusterTable2big))), ]
        clusterTable2 <- clusterTable2bigResorted
    }
    
    
    countTable1 <- table(idsVectorGroup1)
    countTable2 <- table(idsVectorGroup2)
    
    clusterFractionsForAllIds1 <- clusterTable1
    clusterFractionsForAllIds2 <- clusterTable2
    
    
    for (i in seq_len(length(countTable1))) {
        x <- clusterTable1[, i]/countTable1[i]
        clusterFractionsForAllIds1[, i] <- x
    }
    
    for (i in seq_len(length(countTable2))) {
        x <- clusterTable2[, i]/countTable2[i]
        clusterFractionsForAllIds2[, i] <- x
    }
    
    if (missing(pairingVector) == FALSE) {
        # Here, a comparable pairing vector pair
        # is created if a multilevel sPLS-DA
        # should be performed.
        
        pairingVectorGroup1 <- 
            as.character(pairingVector[groupVector == unique(groupVector)[1]])
        pairingVectorGroup2 <- 
            as.character(pairingVector[groupVector == unique(groupVector)[2]])
        
        pairingShortGroup1 <- clusterFractionsForAllIds1[1, ]
        
        for (i in seq_len(ncol(clusterFractionsForAllIds1))) {
            pairingShortGroup1[i] <- 
                pairingVectorGroup1[
                    which(colnames(clusterFractionsForAllIds1)[i] 
                          == idsVectorGroup1)[1]]
        }
        
        pairingShortGroup2 <- clusterFractionsForAllIds2[1, ]
        
        for (i in seq_len(ncol(clusterFractionsForAllIds2))) {
            pairingShortGroup2[i] <- 
                pairingVectorGroup2[
                    which(colnames(clusterFractionsForAllIds2)[i] 
                          == idsVectorGroup2)[1]]
        }
        
        pairingAll <- c(pairingShortGroup1, pairingShortGroup2)
    } else {
        pairingAll <- NULL
    }
    
    # A group vector is created with the same
    # length as the number of columns in the
    # tables
    groupId <- as.factor(c(rep(groupName1, ncol(clusterFractionsForAllIds1)), 
                           rep(groupName2, ncol(clusterFractionsForAllIds2))))
    
    # These two tables are combined to one
    clusterFractionsForAllIds <- as.matrix(cbind(clusterFractionsForAllIds1, 
                                                 clusterFractionsForAllIds2))
    
    # Here, the number of possible clusters
    # to be saved in the sPLS-DA is chosen.
    # The number of tested clusters can be up
    # to five, but if the number of clusters
    # is low, the lowest tested variant will
    # be 1 and that will be tested only once.
    testKeepAlternatives <- vector()
    for (i in seq_len(5)) {
        if (i == 1) {
            testKeepAlternatives[i] <- nrow(clusterFractionsForAllIds)
        } else {
            if (testKeepAlternatives[i - 1] > 1) {
                testKeepAlternatives[i] <- round(testKeepAlternatives[i - 1]/2)
            } else {
                break
            }
        }
    }
    
    return(list(clusterFractionsForAllIds, groupId, testKeepAlternatives, 
                pairingAll))
}
