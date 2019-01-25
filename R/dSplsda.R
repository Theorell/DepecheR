#' Sparse partial least squares discriminant analysis with paired and unpaired data
#'
#'
#' This function is used to compare groups of individuals from whom comparable cytometry or other complex data has been generated. It is superior to just running a Wilcoxon analysis in that it does not consider each cluster individually, but instead uses a sparse partial least squares discriminant analysis to first identify which vector thourgh the multidimensional data cloud, created by the cluster-donor matrix, that optimally separates the groups, and as it is a sparse algorithm, applies a penalty to exclude the clusters that are orthogonal, or almost orthogonal to the discriminant vector, i.e. that do not contribute to separating the groups.
#' @importFrom mixOmics splsda tune.splsda
#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual scale_x_continuous theme element_blank element_rect ggsave
#' @param xYData A dataframe or matrix with two columns. Each row contains information about the x and y positition in the field for that observation.
#' @param idsVector Vector with the same length as xYData containing information about the id of each observation.
#' @param groupVector Vector with the same length as xYData containing information about the group identity of each observation.
#' @param clusterVector Vector with the same length as xYData containing information about the cluster identity of each observation.
#' @param displayVector Optionally, if the dataset is very large (>100 000 observations) and hence the SNE calculation becomes impossible to perform for the full dataset, this vector can be included. It should contain the set of rows from the data used for statistics, that has been used to generate the xYData.
#' @param testSampleRows Optionally, if a train-test setup is wanted, the rows specified in this vector are used to divide the dataset into a training set, used to generate the analysis, and a test set, where the outcome is predicted based on the outcome of the training set. All rows that are not labeled as test rows are assumed to be train rows.
#' @param paired Defaults to FALSE, i.e. no assumption of pairing is made and Wilcoxon rank sum-test is performed. If true, the software will by default pair the first id in the first group with the firs id in the second group and so forth, so make sure the order is correct!
#' @param name The main name for the graph and the analysis.
#' @param densContour If density contours should be created for the plot(s) or not. Defaults to TRUE. a
#' @param groupName1 The name for the first group
#' @param groupName2 The name for the second group
#' @param thresholdMisclassRate This threshold corresponds to the usefulness of the model in separating the groups: a misclassification rate of the default 0.05 means that 5 percent of the individuals are on the wrong side of the theoretical robust middle line between the groups along the sPLS-DA axis, defined as the middle point between the 3:rd quartile of the lower group and the 1:st quartile of the higher group.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved as a png, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @param createOutput For testing purposes. Defaults to TRUE. If FALSE, no output is generated.
#' @seealso \code{\link{dColorPlot}}, \code{\link{dDensityPlot}}, \code{\link{dResidualPlot}}
#' @return This function returns the full result of the sPLS-DA. It also returns a SNE based plot showing which events that belong to a cluster dominated by the first or the second group defined by the sparse partial least squares loadings of the clusters.
#' @examples
#' 
#' # Load some data
#' data(testData)
#' 
#' \dontrun{
#' # Run Barnes Hut tSNE on this. For more rapid example execution, a SNE of the
#' # data is inluded
#' # library(Rtsne)
#' # testDataSNE <- Rtsne(testData[,2:15], pca=FALSE)
#' data(testDataSNE)
#' 
#' # Run the clustering function. For more rapid example execution,
#' # a depeche clustering of the data is inluded
#' # testDataDepeche <- depeche(testData[,2:15])
#' data(testDataDepeche)
#' 
#' 
#' # Run the function. This time without pairing.
#' sPLSDAObject <- dSplsda(
#'   xYData = testDataSNE$Y, idsVector = testData$ids,
#'   groupVector = testData$label, clusterVector = testDataDepeche$clusterVector
#' )
#' 
#' # Here, pairing is used. NB!! This artificial example is only present to show how to
#' # use the function. In reality, pairing should only be used in situations 
#' # where true paired data is present! The only reason this works although this is
#' # non-paired data is that the number of donors is identical. As it is, the algorithm
#' # internally converts the idsVector so that the first individual in group1 is 
#' # associated with the first individual in group2. This can lead to erratic 
#' # problems, so make sure that either a valid id vector, with the same id occuring two
#' # times for each individual is provided, or that the individuals occur in the
#' # exact same order in both groups. 
#' 
#' sPLSDAObject <- dSplsda(
#'    xYData = testDataSNE$Y, idsVector = testData$ids,
#'    groupVector = testData$label, clusterVector = testDataDepeche$clusterVector, 
#'    paired = TRUE, name = 'd_sPLSDAPlot_paired', groupName1 = 'Stimulation 1', 
#'    groupName2 = 'Stimulation 2')
#' 
#' # Here is an example of how the display vector can be used.
#' subsetVector <- sample(1:nrow(testData), size = 10000)
#' 
#' # Now, the SNE for this displayVector could be created
#' # testDataSubset <- testData[subsetVector, 2:15]
#' # testDataSNESubset <- Rtsne(testDataDisplay, pca=FALSE)$Y
#' # But we will just subset the testDataSNE immediately
#' testDataSNESubset <- testDataSNE$Y[subsetVector, ]
#' 
#' # And now, this new SNE can be used for display, although all
#' # the data is used for the sPLS-DA calculations
#' sPLSDAObject <- dSplsda(
#'   xYData = testDataSNESubset, idsVector = testData$ids,
#'   groupVector = testData$label, clusterVector = testDataDepeche$clusterVector,
#'   displayVector = subsetVector
#' )
#' 
#' # Finally, an example of a train-test set situation, where a random half the dataset
#' # is used for training and the second half is used for testing. It is naturally more
#' # biologically interesting to use two independent datasets for training and testing
#' # in the real world.
#' testDataRows <- sample(1:nrow(testData), size = 48500)
#' sPLSDAObject <- dSplsda(
#'   xYData = testDataSNE$Y, idsVector = testData$ids,
#'   groupVector = testData$label, clusterVector = testDataDepeche$clusterVector,
#'   testSampleRows = testDataRows
#' )
#' }
#' @export dSplsda
dSplsda <- function(xYData, idsVector, groupVector, 
    clusterVector, displayVector, testSampleRows, 
    paired = FALSE, densContour = TRUE, name = "default", 
    groupName1 = unique(groupVector)[1], 
    groupName2 = unique(groupVector)[2], 
    thresholdMisclassRate = 0.05, title = FALSE, 
    createDirectory = FALSE, directoryName = "dSplsda", 
    bandColor = "black", dotSize = 500/sqrt(nrow(xYData)), 
    createOutput = TRUE) {
    if (createDirectory == TRUE) {
        dir.create(directoryName)
        workingDirectory <- getwd()
        setwd(paste(workingDirectory, directoryName, 
            sep = "/"))
    }
    
    if (length(unique(groupVector)) != 2) {
        stop("More or less than two groups are present. This is currently not supported.")
    }
    
    if (length(unique(idsVector)) < 8) {
        warning("NB! The number of unique ids is smaller than 8, so statistical comparison is not suitable. Use dResidualPlot instead to view differences.")
    }
    
    if (name == "default") {
        name <- paste0(groupName1, "_vs_", 
            groupName2)
    }
    
    if (paired == TRUE) {
        # As the algorithm does not like repeated
        # ids, if the id vector contains the same
        # values for the first group and the
        # second, a new id vector is introduced
        # here
        if (identical(unique(idsVector[groupVector == 
            unique(groupVector)[1]]), unique(idsVector[groupVector == 
            unique(groupVector)[2]]))) {
            pairingVector <- idsVector
            idsVector <- paste0(idsVector, 
                groupVector)
        } else if (length(unique(idsVector[groupVector == 
            unique(groupVector)[1]])) == 
            length(unique(idsVector[groupVector == 
                unique(groupVector)[2]]))) {
            pairingVector <- c(idsVector[groupVector == 
                unique(groupVector)[1]], 
                idsVector[groupVector == 
                  unique(groupVector)[1]])
        } else {
            print("Pairing cannot be performed, as the first and second datasets contain different number of individual Ids")
        }
    }
    
    if (any(is(xYData) == "matrix")) {
        xYData <- as.data.frame(xYData)
    }
    
    if (missing(testSampleRows) == FALSE) {
        clusterVectorTrain <- clusterVector[-testSampleRows]
        idsVectorTrain <- idsVector[-testSampleRows]
        groupVectorTrain <- groupVector[-testSampleRows]
        clusterVectorTest <- clusterVector[testSampleRows]
        idsVectorTest <- idsVector[testSampleRows]
        groupVectorTest <- groupVector[testSampleRows]
        
        if (paired == TRUE) {
            pairingVectorTrain <- pairingVector[-testSampleRows]
            pairingVectorTest <- pairingVector[testSampleRows]
        }
    } else {
        clusterVectorTrain <- clusterVector
        idsVectorTrain <- idsVector
        groupVectorTrain <- groupVector
        if (paired == TRUE) {
            pairingVectorTrain <- pairingVector
        }
    }
    
    if (paired == FALSE) {
        dSplsdaInData <- dSplsdaPreCalculations(clusterVectorTrain, 
            idsVectorTrain, groupVectorTrain, 
            groupName1 = groupName1, groupName2 = groupName2)
    } else {
        dSplsdaInData <- dSplsdaPreCalculations(clusterVectorTrain, 
            idsVectorTrain, groupVectorTrain, 
            groupName1 = groupName1, groupName2 = groupName2, 
            pairingVector = pairingVectorTrain)
    }
    
    # Here, the number of clusters that
    # should be kept in the sPLS-DA is chosen
    nVarSPLSDA <- tune.splsda(X = t(dSplsdaInData[[1]]), 
        Y = dSplsdaInData[[2]], ncomp = 1, 
        logratio = "none", test.keepX = dSplsdaInData[[3]], 
        validation = "loo", dist = "mahalanobis.dist", 
        multilevel = dSplsdaInData[[4]])
    
    # And here the sPLS-DA is performed.
    # Scaling is performed internally in the
    # algorithm. The number of components is
    # set to one, as this situation is by far
    # the most interpretable.
    sPLSDAObject <- splsda(X = t(dSplsdaInData[[1]]), 
        Y = dSplsdaInData[[2]], ncomp = 1, 
        keepX = nVarSPLSDA$choice.keepX, 
        multilevel = dSplsdaInData[[4]])
    
    # Retrieve the x variates for plotting
    sPLSDA_vector <- data.frame(sPLSDAObject$variates$X)
    Group <- dSplsdaInData[[2]]
    densityHist <- cbind(sPLSDA_vector, Group)
    
    colnames(densityHist) <- c("sPLSDA_vector", 
        "Group")
    
    # Density plots with semi-transparent
    # fill
    ggplot(densityHist, aes(x = sPLSDA_vector, 
        fill = Group)) + geom_density(adjust = 0.2, 
        alpha = 0.4) + scale_fill_manual(values = c("red", 
        "blue")) + scale_x_continuous(limits = c(min(densityHist$sPLSDA_vector) - 
        abs(max(densityHist$sPLSDA_vector) - 
            min(densityHist$sPLSDA_vector)) * 
            0.3, max(densityHist$sPLSDA_vector) + 
        abs(max(densityHist$sPLSDA_vector) - 
            min(densityHist$sPLSDA_vector)) * 
            0.3)) + theme(line = element_blank(), 
        panel.background = element_rect(fill = "white"))
    if (createOutput == TRUE) {
        ggsave("Individuals_distributed_along_sPLS-DA_vector.pdf", 
            dpi = 300)
    }
    
    # Retrieve the sparse loadings
    sPLSDALoadings <- sPLSDAObject$loadings$X
    
    # Here, the maximum values for the
    # plotting are defined. If not added by
    # the user, they are related to the
    # overlap between the populations.
    # Furthermore, the data is re-scaled in
    # here, to be more visually
    # understandable.
    group1SplsDa <- densityHist$sPLSDA_vector[which(densityHist$Group == 
        groupName1)]
    group2SplsDa <- densityHist$sPLSDA_vector[which(densityHist$Group == 
        groupName2)]
    
    # Now, the border between the populations
    # is defined
    
    if (min(group1SplsDa) > min(group2SplsDa)) {
        lowGroup <- group2SplsDa
        highGroup <- group1SplsDa
    } else {
        lowGroup <- group1SplsDa
        highGroup <- group2SplsDa
    }
    groupBorder <- (quantile(lowGroup, probs = 0.75) + 
        quantile(highGroup, probs = 0.25))/2
    lowErrors <- length(which(highGroup < 
        groupBorder))
    highErrors <- length(which(lowGroup > 
        groupBorder))
    misclassRate <- (lowErrors + highErrors)/sum(length(highGroup), 
        length(lowGroup))
    if (max(lowGroup) < min(highGroup)) {
        print("The separation of the datasets was perfect, with no overlap between the groups")
        lowestPlottedOverlap <- 0
        absSPLSDALoadings <- abs(sPLSDALoadings)
    } else if (misclassRate <= thresholdMisclassRate) {
        print(paste("The separation of the datasets was very clear, with the misclassification rate being ", 
            round(100 * misclassRate), " percent", 
            sep = ""))
        lowestPlottedOverlap <- misclassRate
        absSPLSDALoadings <- abs(sPLSDALoadings)
    } else {
        print(paste("The separation of the datasets was not clear, with the misclassification rate being ", 
            round(100 * misclassRate), " percent", 
            sep = ""))
        scalingValue <- thresholdMisclassRate/misclassRate
        absSPLSDALoadings <- abs(sPLSDALoadings * 
            scalingValue)
    }
    
    
    # Now, the sign is changed of the
    # splsda-loadings, to make the function
    # generate a similar visual output to the
    # dWilcox function.  Now the median for
    # each group and cluster is calculated
    group1Data <- dSplsdaInData[[1]][, which(as.character(dSplsdaInData[[2]]) == 
        groupName1)]
    group2Data <- dSplsdaInData[[1]][, which(as.character(dSplsdaInData[[2]]) == 
        groupName2)]
    
    median1 <- apply(group1Data, 1, median)
    median2 <- apply(group2Data, 1, median)
    correctSPLSDALoadings <- absSPLSDALoadings
    for (i in seq_len(nrow(group1Data))) {
        if (median1[i] >= median2[i]) {
            correctSPLSDALoadings[i] <- absSPLSDALoadings[i]
        } else {
            correctSPLSDALoadings[i] <- -absSPLSDALoadings[i]
        }
    }
    
    # Here, a vector with the same length as
    # the cluster vector is generated, but
    # where the cluster info has been
    # substituted with the statistic.  If a
    # displayVector has been included, it is
    # used here, to subset the clusterVector
    if (missing(displayVector) == FALSE) {
        statisticVector <- clusterVector[displayVector]
        clusterVectorUsed <- clusterVector[displayVector]
    } else {
        statisticVector <- clusterVector
        clusterVectorUsed <- clusterVector
    }
    
    for (i in seq_len(nrow(correctSPLSDALoadings))) {
        statisticVector[clusterVectorUsed == 
            rownames(correctSPLSDALoadings)[i]] <- correctSPLSDALoadings[i]
    }
    
    # Here the data that will be used for
    # plotting is scaled.
    colnames(xYData) <- c("V1", "V2")
    
    brks <- seq(-1, 1, length.out = 10)
    
    # assign each value to a bin
    grps <- cut(statisticVector, breaks = brks, 
        include.lowest = TRUE)
    colors <- colorRampPalette(c("#FF0000", 
        "white", "#0000FF"))(9)
    xYData$col <- colors[grps]
    
    # Create the density matrix for xYData.
    if (is.logical(densContour)) {
        if (densContour == TRUE) {
            densContour <- dContours(xYData)
        }
    }
    if (length(densContour) > 1) {
        xlim <- c(min(densContour[[1]]), 
            max(densContour[[1]]))
        ylim <- c(min(densContour[[2]]), 
            max(densContour[[2]]))
    } else {
        minX <- min(xYData[, 1])
        maxX <- max(xYData[, 1])
        minY <- min(xYData[, 2])
        maxY <- max(xYData[, 2])
        xlim <- c(minX - abs(minX * 0.05), 
            maxX + abs(maxX * 0.05))
        ylim <- c(minY - abs(minY * 0.05), 
            maxY + abs(maxY * 0.05))
    }
    
    png(paste(name, "_sPLSDA_result.png", 
        sep = ""), width = 2500, height = 2500, 
        units = "px", bg = "transparent")
    if (createOutput == TRUE) {
        if (title == TRUE) {
            plot(V2 ~ V1, data = xYData, 
                main = name, pch = 20, cex = dotSize, 
                cex.main = 5, col = col, 
                xlim = xlim, ylim = ylim, 
                axes = FALSE, xaxs = "i", 
                yaxs = "i")
        }
        
        if (title == FALSE) {
            plot(V2 ~ V1, data = xYData, 
                main = "", pch = 20, cex = dotSize, 
                cex.main = 5, col = col, 
                xlim = xlim, ylim = ylim, 
                axes = FALSE, xaxs = "i", 
                yaxs = "i")
        }
        if (length(densContour) > 1) {
            par(fig = c(0, 1, 0, 1), mar = c(6, 
                4.5, 4.5, 2.5), new = TRUE)
            contour(x = densContour$x, y = densContour$y, 
                z = densContour$z, xlim = xlim, 
                ylim = ylim, nlevels = 10, 
                col = bandColor, lwd = 8, 
                drawlabels = FALSE, axes = FALSE, 
                xaxs = "i", yaxs = "i")
        }
    }
    dev.off()
    # Create a color legend with text
    
    yname <- "Misclass-corrected sPLS-DA loadings"
    topText <- paste(groupName1, " is more abundant", 
        sep = "")
    bottomText <- paste(groupName2, " is more abundant", 
        sep = "")
    legendTitle <- paste("Color scale for", 
        name, "sPLS-DA analysis.pdf", sep = " ")
    
    if (createOutput == TRUE) {
        pdf(legendTitle)
        par(fig = c(0.35, 0.65, 0, 1), xpd = NA)
        z <- matrix(seq_len(9), nrow = 1)
        x <- 1
        y <- seq(-1, 1, len = 9)
        image(x, y, z, col = colors, axes = FALSE, 
            xlab = "", ylab = yname)
        axis(2)
        text(1, 1 * 1.2, labels = topText, 
            cex = 1.1)
        text(1, -1 * 1.2, labels = bottomText, 
            cex = 1.1)
        box()
        dev.off()
    }
    
    # Return data from the sPLS-DA that was
    # needed for the generation of the graphs
    if (createOutput == TRUE) {
        write.csv(sPLSDALoadings, paste0(name, 
            "_sPLSDALoadings.csv"))
        write.csv(sPLSDA_vector, paste0(name, 
            "_sPLSDA_vector.csv"))
    }
    
    
    if (createDirectory == TRUE) {
        setwd(workingDirectory)
    }
    
    # Now, prediction is performed, if the
    # setup is train-test.
    if (missing(testSampleRows) == FALSE) {
        if (paired == FALSE) {
            dSplsdaInDataTest <- dSplsdaPreCalculations(clusterVectorTest, 
                idsVectorTest, groupVectorTest, 
                groupName1 = groupName1, 
                groupName2 = groupName2)
        } else {
            dSplsdaInDataTest <- dSplsdaPreCalculations(clusterVectorTest, 
                idsVectorTest, groupVectorTest, 
                groupName1 = groupName1, 
                groupName2 = groupName2, 
                pairingVector = pairingVectorTest)
        }
        # And here the sPLS-DA is performed.
        # Scaling is performed internally in the
        # algorithm. The number of components is
        # set to one, as this situation is by far
        # the most interpretable.
        sPLSDAPredictObject <- predict(object = sPLSDAObject, 
            newdata = t(dSplsdaInDataTest[[1]]), 
            multilevel = dSplsdaInDataTest[[4]])
        
        # Retrieve the x variates for plotting
        sPLSDAX <- data.frame(sPLSDAPredictObject$variates)
        densityHist <- cbind(sPLSDAX, dSplsdaInDataTest[[2]])
        
        colnames(densityHist) <- c("sPLSDA_vector", 
            "Group")
        
        # Density plots with semi-transparent
        # fill
        ggplot(densityHist, aes(x = sPLSDA_vector, 
            fill = Group)) + geom_density(adjust = 0.2, 
            alpha = 0.4) + scale_fill_manual(values = c("red", 
            "blue")) + scale_x_continuous(limits = c(min(densityHist$sPLSDA_vector) - 
            abs(max(densityHist$sPLSDA_vector) - 
                min(densityHist$sPLSDA_vector)) * 
                0.3, max(densityHist$sPLSDA_vector) + 
            abs(max(densityHist$sPLSDA_vector) - 
                min(densityHist$sPLSDA_vector)) * 
                0.3)) + theme(line = element_blank(), 
            panel.background = element_rect(fill = "white"))
        if (createOutput == TRUE) {
            ggsave("Predicted_individuals_distributed_along_sPLS-DA_vector.pdf", 
                dpi = 300)
        }
    }
    
    print(paste0("Files were saved at ", 
        getwd()))
    
    if (missing(testSampleRows) == TRUE) {
        return(sPLSDAObject)
    } else {
        return(list(sPLSDAObject, sPLSDAPredictObject))
    }
}
