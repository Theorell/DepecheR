#' Showing the residuals when subtracting the values from one group from 
#' another on a SNE plot
#'
#'
#' This function is used to visually compare groups of individuals from whom 
#' comparable cytometry or other complex data has been generated, but where the
#' number of individuals does not permit any statistical comparisons.
#' @param xYData A dataframe or matrix with two columns. Each row contains 
#' information about the x and y positition in the field for that observation.
#' @param groupVector Vector with the same length as xYData containing 
#' information about the group identity of each observation.
#' @param clusterVector Vector with the same length as xYData containing 
#' information about the cluster identity of each observation.
#' @param densContour If density contours should be created for the plot(s) or
#' not. Defaults to TRUE.
#' @param plotName The main name for the graph and the analysis.
#' @param groupName1 The name for the first group
#' @param groupName2 The name for the second group
#' @param maxAbsPlottingValues If multiple plots should be compared, it might 
#' be useful to define a similar color scale for all plots, so that the same 
#' color always means the same value. Such a value can be added here. It 
#' defaults to the maximum Wilcoxon statistic that is generated in the analysis.
#' @param title If there should be a title displayed on the plotting field. As
#' the plotting field is saved as a png, this title cannot be removed as an 
#' object afterwards, as it is saved as coloured pixels. To simplify usage for 
#' publication, the default is FALSE, as the files are still named, eventhough 
#' no title appears on the plot.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots 
#' smaller the more observations that are included.
#' @param plotDir If different from the current directory. If specified and 
#' non-existent, the function creates it. If "." is specified, the plots will be
#' saved at the current directory.
#' @param createOutput For testing purposes. Defaults to TRUE. If FALSE, no 
#' plots are generated.
#' @seealso \code{\link{dColorPlot}}, \code{\link{dDensityPlot}}, 
#' \code{\link{dWilcox}}
#' @return A sne based plot showing which events that belong to a cluster 
#' dominated by the first or the second group.
#' @examples
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
#' # And finally run the function
#' dResidualPlot(
#'   xYData = testDataSNE$Y, groupVector = testData[, 16],
#'   clusterVector = testDataDepeche$clusterVector
#' )
#' }
#' @export dResidualPlot
dResidualPlot <- function(xYData, groupVector, clusterVector, 
                          densContour = TRUE, 
                          groupName1 = unique(groupVector)[1], 
                          groupName2 = unique(groupVector)[2], 
                          plotName = "default", title = FALSE, maxAbsPlottingValues,
                          bandColor = "black", plotDir = ".", 
                          dotSize = 400/sqrt(nrow(xYData)), 
                          createOutput = TRUE) {
    
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
    
    # Here, the residuals are identified.  A
    # table with the percentage of cells in
    # each cluster for each group is created
    # in analogy with XXX pKMRun.
    
    clusterTable <- table(clusterVector, groupVector)
    
    countTable <- table(groupVector)
    
    clusterPercentagesForGroups <- clusterTable
    
    for (i in seq_len(length(countTable))) {
        x <- 100 * clusterTable[, i]/countTable[i]
        clusterPercentagesForGroups[, i] <- x
    }
    
    # Now the residual matrix is constructed
    residualVector <- as.vector(clusterPercentagesForGroups[, 1] - 
                                    clusterPercentagesForGroups[, 2])
    
    names(residualVector) <- as.numeric(row.names(clusterPercentagesForGroups))
    
    # Here, the maximum values for the
    # plotting is defined if not added by the
    # user.
    if (missing(maxAbsPlottingValues)) {
        maxAbsPlottingValues <- max(abs(clusterPercentagesForGroups))
    }
    
    # It might be that the
    # maxAbsPlottingValues defined by the
    # user are lower than the most extreme
    # values in the residual data. Then, the
    # data is truncated to fit the range.
    residualVector[residualVector > maxAbsPlottingValues] <- 
        maxAbsPlottingValues
    residualVector[residualVector < -maxAbsPlottingValues] <- 
        -maxAbsPlottingValues
    
    # Here, a vector with the same length as
    # the cluster vector is generated, but
    # where the cluster info has been
    # substituted with the statistic.
    residualVectorLong <- clusterVector
    for (i in seq_len(length(residualVector))) {
        residualVectorLong[clusterVector == 
            names(residualVector)[i]] <- residualVector[i]
    }
    
    # Here the data that will be used for
    # plotting is scaled.
    colnames(xYData) <- c("V1", "V2")
    
    # Make a color vector with the same
    # length as the data
    residual.df <- as.data.frame(residualVectorLong)
    
    # make a breaks vector to define each bin
    # for the colors
    brks <- with(residual.df, seq(-maxAbsPlottingValues, maxAbsPlottingValues, 
                                  length.out = 10))
    
    # assign each value to a bin
    grps <- with(residual.df, cut(residual.df[, 1], breaks = brks, 
                                  include.lowest = TRUE))
    colors <- colorRampPalette(c("#FF0000", "white", "#0000FF"))(9)
    xYData$col <- colors[grps]
    
    dPlotCoFunction(colorVariable = xYData$col, plotName = 
                        paste0(plotName, "_residuals"), 
                    xYData = xYData, title = title, 
                    densContour = densContour, bandColor = bandColor, 
                    dotSize = dotSize, plotDir = plotDir, 
                    createOutput = createOutput)
    
    # Create a color legend with text
    
    yname <- "Residual values"
    topText <- paste0(groupName1, " is more abundant")
    bottomText <- paste0(groupName2, " is more abundant")

    if (createOutput) {
        pdf(file.path(plotDir, paste0(plotName, "_residual_scale.pdf")))
        par(fig = c(0.35, 0.65, 0, 1), xpd = NA)
        z <- matrix(seq_len(9), nrow = 1)
        x <- 1
        y <- seq(-maxAbsPlottingValues, maxAbsPlottingValues, len = 9)
        image(x, y, z, col = colors, axes = FALSE, xlab = "", ylab = yname)
        axis(2)
        text(1, maxAbsPlottingValues * 1.2, labels = topText, cex = 1.1)
        text(1, -maxAbsPlottingValues * 1.2, labels = bottomText, cex = 1.1)
        box()
        dev.off()
    }
}
