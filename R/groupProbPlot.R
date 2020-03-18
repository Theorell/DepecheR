#' Define and plot group probabilities
#'
#'
#' This function defines and plots the single-observation probability for
#' belonging to either of two groups. It builds on the same idea as has been put
#' forward in the Sconify package:
#' -Burns TJ (2019). Sconify: A toolkit for performing KNN-based statistics for
#' flow and mass cytometry data. R package version 1.4.0 and
#' -Hart GT, Tran TM, Theorell J, Schlums H, Arora G, Rajagopalan S, et al.
#' Adaptive NK cells in people exposed to Plasmodium falciparum correlate
#' with protection from malaria. J Exp Med. 2019 Jun 3;216(6):1280–90.
#' First, the k nearest neighbors are defined for each individual cell, and the
#' cell in question thereafter gets a group probability assigned to it, which
#' is calculated by defining the percentage of neighbors belonging to each
#' respective groups. In other words, if 20 out of 100 neighbors belong to group
#' A and 80 belong to group B, and the value for the cell will be 20% for group
#' A or 80% for group B, depending on the point of view. This will also be
#' accordingly reflected in the color scale on the resulting plot.
#'
#' @param xYData A dataframe or matrix with two columns. Each row contains
#' information about the x and y positition in the field for that observation.
#' @param groupVector Vector with the same length as xYData containing
#' information about the group identity of each observation.
#' @param dataTrans The data cloud in which the nearest neighbors for the
#' events should be identified.
#' @param kNeighK The number of nearest neighbors.
#' @param kMeansK The number of clusters in the initial step of the algorithm.
#' A higher number leads to shorter runtime, but potentially lower accuracy.
#' @param plotName The main name for the graph and the analysis.
#' @param densContour If density contours should be created for the plot(s) or
#' not. Defaults to TRUE. a
#' @param groupName1 The name for the first group
#' @param groupName2 The name for the second group
#' @param title If there should be a title displayed on the plotting field. As
#' the plotting field is saved as a png, this title cannot be removed as an
#' object afterwards, as it is saved as coloured pixels. To simplify usage for
#' publication, the default is FALSE, as the files are still named, eventhough
#' no title appears on the plot.
#' @param plotDir If different from the current directory. If specified and
#' non-existent, the function creates it. If "." is specified, the plots will be
#' saved at the current directory.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots
#' smaller the more observations that are included.
#' @param returnProb Should a probability vector be returned? Mutually exclusive
#' with returnProbColVec. 
#' @param returnProbColVec Should the color vector be returned as part of the
#' output? Mutually exclusive with returnProb. 
#' @param createOutput For testing purposes. Defaults to TRUE. If FALSE, no
#' output is generated.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom FNN knnx.index
#' @importFrom stats kmeans
#' @return A graph showing the probability as a color scale from blue over white
#' to red for each event to belong to one group or the other, with a separate 
#' color scale. Optionally also the color vector, if returnProbColVec is TRUE.
#' @examples
#' data(testData)
#' data(testDataSNE)
#' dataTrans <-
#'   testData[, c("SYK", "CD16", "CD57", "EAT.2", "CD8", "NKG2C", "CD2", "CD56")]
#' \dontrun{
#' groupProbPlot(xYData = testDataSNE$Y, groupVector = testData$label, 
#' dataTrans)
#' }
#'
#' @export groupProbPlot

groupProbPlot <- function(xYData, groupVector, dataTrans,
                          kNeighK = max(100, round(nrow(dataTrans) / 10000)),
                          kMeansK = round(nrow(dataTrans) / 1000),
                          densContour = TRUE,
                          groupName1 = unique(groupVector)[1],
                          groupName2 = unique(groupVector)[2],
                          plotName = "default", title = FALSE,
                          bandColor = "black", plotDir = ".",
                          dotSize = 400 / sqrt(nrow(xYData)), 
                          returnProb = FALSE,
                          returnProbColVec = FALSE, createOutput = TRUE) {
    if (plotDir != ".") {
        dir.create(plotDir)
    }

    if (length(unique(groupVector)) != 2) {
        stop("More or less than two groups are present. Please correct this.")
    }

    if (is.matrix(xYData)) {
        xYData <- as.data.frame(xYData)
    }

    if (plotName == "default") {
        plotName <- paste0(groupName1, "_vs_", groupName2)
    }

    rowNumbers <- seq_len(nrow(dataTrans))
    rowsGroup1 <- rowNumbers[which(groupVector == unique(groupVector)[1])]
    rowsGroup2 <- rowNumbers[which(groupVector == unique(groupVector)[2])]
    groupValVec <- rep(100, length.out = nrow(dataTrans))
    groupValVec[rowsGroup2] <- -100

    # Here, a new, balanced group neighbor vector is introduced
    if (length(rowsGroup1) != length(rowsGroup2)) {
        groupLengths <- c(length(rowsGroup1), length(rowsGroup2))
        if (min(groupLengths) < 10000) {
            if (min(groupLengths) < 5000) {
                stop(paste0("At least one of the groups has less than 5000",
                            " events, which will lead to a too low resolution ",
                            "for this analysis to be meaningful. Use ",
                            "dResidualPlot instead."))
            }
            message(paste0("The groups are not matched in size and the ",
                           "smaller group has only ", min(groupLengths), 
                           " rows, leading to a reduction in neighbor",
                    " resolution to this level for both groups"))
        }
        rowsGroup1 <- rowsGroup1[sample(
            seq(1, length(rowsGroup1)),
            min(groupLengths)
        )]
        rowsGroup2 <- rowsGroup2[sample(
            seq(1, length(rowsGroup2)),
            min(groupLengths)
        )]
    }

    groupNeighRowVec <- c(rowsGroup1, rowsGroup2)


    # Each of the raw channels are scaled using the standard
    # robustVarScale method in the dScale function
    dataScaled <- dScale(dataTrans, center = FALSE)

    # Now, the cells are clustered according to this analysis
    kMeansResult <- kmeans(dataScaled, kMeansK, iter.max = max(200, kMeansK))
    kMeansCenters <- kMeansResult$centers
    kMeansClusters <- kMeansResult$cluster
    print("Done with k-means")

    dataTransList <- split(dataTrans, kMeansClusters)
    rowNumbersList <- split(rowNumbers, kMeansClusters)
    
    if(kMeansK > 11){
        distCenters <- knnx.index(kMeansCenters, kMeansCenters, 11)  
    } else {
        distCenters <- t(matrix(seq_len(kMeansK), ncol = kMeansK, nrow = kMeansK))
    }

    print("Now the first bit is done, and the iterative part takes off")
    nCores <- detectCores() - 1

    cl <- parallel::makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)

    allClusters <- unique(kMeansClusters)
    firstCluster <- 1
    resultList <- list()
    x <- 1
    while (firstCluster <= length(allClusters)) {
        timeBefore <- Sys.time()
        if (((firstCluster + nCores) - 1) < length(allClusters)) {
            clusterRange <- seq(firstCluster, ((firstCluster + nCores) - 1))
        } else {
            clusterRange <- seq(firstCluster, length(allClusters))
        }

        # And here, all the clusters that are close to the nonIndata clusters 
        #are selected
        if (length(clusterRange) > 1) {

            # Now, the datasets are constructed from this cluster range
            dataTransListFocus <- dataTransList[clusterRange]

            distCentClustRange <- distCenters[clusterRange, ]
            rowNumbersNeighListFocus <-
                lapply(seq_len(nrow(distCentClustRange)), function(y)
                    unlist(rowNumbersList[allClusters %in% distCentClustRange[y,]]))
            # This list is now focused on the group neighbors, which will only
            # have an effect in the cases where the number of events differ 
            #between groups.
            groupRowNumNeighListFocus <-
                lapply(rowNumbersNeighListFocus, function(y)
                    return(y[y %in% groupNeighRowVec]))

            dataTransNeighListFocus <- lapply(groupRowNumNeighListFocus, 
                                              function(y)
                                                return(dataTrans[y, ]))

            dataReturnListFocus <- lapply(groupRowNumNeighListFocus, 
                                          function(y)
                                            return(groupValVec[y]))
            i <- 1
            resultFocus <-
                foreach(
                    i = seq_along(dataTransListFocus),
                    .packages = "DepecheR"
                ) %dopar% microClust(
                    dataCenter = dataTransListFocus[[i]],
                    dataNeigh = dataTransNeighListFocus[[i]],
                    dataReturn = dataReturnListFocus[[i]], method = "mean",
                    k = kNeighK, trim = 0
                )

            resultList[[x]] <- unlist(resultFocus)
        } else {
            dataTransFocus <- dataTransList[[clusterRange]]

            rowNumbersNeighFocus <-
                unlist(rowNumbersList[allClusters %in% 
                                        distCenters[clusterRange, ]])
            # This list is now focused on the group neighbors, which will only
            # have an effect in the cases where the number of events differ
            #between groups.
            groupRowNumNeighFocus <-
                rowNumbersNeighFocus[rowNumbersNeighFocus %in% groupNeighRowVec]

            dataTransNeighFocus <- dataTrans[groupRowNumNeighFocus, ]

            dataReturnFocus <- groupValVec[groupRowNumNeighFocus]

            resultList[[x]] <- microClust(
                dataCenter = dataTransFocus,
                dataNeigh = dataTransNeighFocus,
                dataReturn = dataReturnFocus,
                method = "mean", k = kNeighK
            )
        }


        timeAfter <- as.numeric(Sys.time() - timeBefore)
        print(paste0("Clusters ", clusterRange[1], " to ", 
                     clusterRange[length(clusterRange)], 
                     " smoothed in ", timeAfter, " ", 
                     attributes(timeAfter)$units, ". Now, ", 
                     length(allClusters) - clusterRange[length(clusterRange)], 
                     " clusters are left."))
        firstCluster <- (clusterRange[length(clusterRange)] + 1)
        x <- x + 1
    }

    parallel::stopCluster(cl)
    fullResult <- unlist(resultList)
    fullResultOrdered <- fullResult[order(unlist(rowNumbersList))]


    # Here the data that will be used for
    # plotting is scaled.
    colnames(xYData) <- c("V1", "V2")

    # Make a color vector with the same
    # length as the data
    residual.df <- as.data.frame(fullResultOrdered)

    # make a breaks vector to define each bin
    # for the colors
    brks <- with(residual.df, seq(-100, 100,
        length.out = 10
    ))

    # assign each value to a bin
    grps <- with(residual.df, cut(residual.df[, 1],
        breaks = brks,
        include.lowest = TRUE
    ))
    colors <- colorRampPalette(c("#FF0000", "white", "#0000FF"))(9)
    xYData$col <- colors[grps]

    dPlotCoFunction(
        colorVariable = xYData$col, plotName =
            paste0(plotName, "_probabilities"),
        xYData = xYData, title = title,
        densContour = densContour, bandColor = bandColor,
        dotSize = dotSize, plotDir = plotDir,
        createOutput = createOutput
    )

    # Create a color legend with text

    yname <- "Probability for group identity, percent"
    topText <- paste0(groupName1, " more probable")
    bottomText <- paste0(groupName2, " more probable")
    z <- matrix(seq_len(9), nrow = 1)
    x <- 1
    y <- seq(-100, 100, len = 9)

    if (createOutput) {
        pdf(file.path(plotDir, paste0(plotName, "_probability_scale.pdf")))
        par(fig = c(0.35, 0.65, 0, 1), xpd = NA)
        image(x, y, z, col = colors, axes = FALSE, xlab = "", ylab = yname)
        text(0.4, -100, labels = 100)
        text(0.4, -50, labels = 75)
        text(x = 0.4, y = 0, labels = 50)
        text(0.4, 50, labels = 75)
        text(0.4, 100, labels = 100)
        text(1, 120, labels = topText, cex = 1.1)
        text(1, -120, labels = bottomText, cex = 1.1)
        box()
        dev.off()
    }
    if(returnProb){
        return(fullResultOrdered)
    } else if (returnProbColVec){
        return(xYData$col)
    }
}