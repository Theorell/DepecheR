#' Define and plot group probabilities
#'
#'
#' This function defines and plots the single-observation probability for
#' belonging to either of two groups. It uses the \code{\link{neighSmooth}}
#' function with the special case that the values are binary: For each set of
#' k nearest neighbors, cell x is assigned a probability to belong to one group
#' or the other based on the percentage of the neighbors belonging to each
#' group. In other words, if 20 out of 100 neighbors belong to group
#' A and 80 belong to group B, and the value for the cell will be 20% for group
#' A or 80% for group B, depending on the point of view. This will also be
#' accordingly reflected in the color scale on the resulting plot.
#'
#' @param xYData A dataframe or matrix with two columns. Each row contains
#' information about the x and y positition in the field for that observation.
#' @param groupVector Vector with the same length as xYData containing
#' information about the group identity of each observation.
#' @param euclidSpaceData The data cloud in which the nearest neighbors for the
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
#' euclidSpaceData <-
#'     testData[, c(
#'         "SYK", "CD16", "CD57", "EAT.2",
#'         "CD8", "NKG2C", "CD2", "CD56"
#'     )]
#' \dontrun{
#' groupProbPlot(
#'     xYData = testDataSNE$Y, groupVector = testData$label,
#'     euclidSpaceData
#' )
#' }
#'
#' @export groupProbPlot

groupProbPlot <- function(xYData, groupVector, euclidSpaceData,
                          kNeighK = max(
                              100, round(nrow(euclidSpaceData) / 10000)
                          ),
                          kMeansK = round(nrow(euclidSpaceData) / 1000),
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

    rowNumbers <- seq_along(groupVector)
    rowsGroup1 <- rowNumbers[which(groupVector == unique(groupVector)[1])]
    rowsGroup2 <- rowNumbers[which(groupVector == unique(groupVector)[2])]
    groupValVec <- rep(100, length.out = nrow(euclidSpaceData))
    groupValVec[rowsGroup2] <- -100

    # Here, a new, balanced group neighbor vector is introduced
    if (length(rowsGroup1) != length(rowsGroup2)) {
        groupLengths <- c(length(rowsGroup1), length(rowsGroup2))
        if (min(groupLengths) < 10000) {
            if (min(groupLengths) < 5000) {
                stop(paste0(
                    "At least one of the groups has less than 5000",
                    " events, which will lead to a too low resolution ",
                    "for this analysis to be meaningful. Use ",
                    "dResidualPlot instead."
                ))
            }
            message(paste0(
                "The groups are not matched in size and the ",
                "smaller group has only ", min(groupLengths),
                " rows, leading to a reduction in neighbor",
                " resolution to this level for both groups"
            ))
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

    # And here, we use the neighbor smoothing algorithm to identify the
    # probabilities
    fullResult <- neighSmooth(
        focusData = groupValVec,
        euclidSpaceData = euclidSpaceData,
        neighRows = c(rowsGroup1, rowsGroup2),
        kNeighK = kNeighK, kMeansK = kMeansK
    )

    #Now, we make a vital change to the scale: the values never go below
    #50%, as we are projecting two scale on each other, thus showing that
    #a value of -49 also means +51 and vice versa. Therefore the values
    #need to be changed.
    fullTempResult <- fullResult+(50-(50*fullResult)/100)
    fullResult <- vapply(fullTempResult, function(x) if(x<50){x-100}else{x}, 100)

    # Here the data that will be used for
    # plotting is scaled.
    colnames(xYData) <- c("V1", "V2")

    # Make a color vector with the same
    # length as the data
    residual.df <- as.data.frame(fullResult)

    # make a breaks vector to define each bin
    # for the colors
    brks <- with(residual.df, c(-100, -90, -80, -70, -60, 0, 60, 70, 80, 90, 100))

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
    if (returnProb) {
        return(fullResult)
    } else if (returnProbColVec) {
        return(xYData$col)
    }
}
