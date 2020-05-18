# This function is used by dDensityPlot.
# It calculates the density at each point of a 2-dimensional field of
# data, creates a color density gradient based on this, and resorts the events
# so that the lowest density areas will be plotted first, and so on.
# The multiple color option allows for multiple separate densities to be 
# plotted simulataneously, for example if multiple clusters should be displayed
# together with density.
# For information on the parameters, see dDensityPlot

dDensityPlotCoFunction <- function(xYData, idsVector, uniqueIds, color,
                                   colorList, plotName, densContour, bandColor,
                                   dotSize, title, plotDir,
                                   createOutput = FALSE) {
    if (missing(idsVector)) {
        df <- xYData
        colnames(df) <- c("V1", "V2")
        ## Use densCols() output to get density at
        ## each point. The colors here are only
        ## supporting the coming order of the rows
        ## further down the script.
        df$col <- densCols(df$V1, df$V2,
            colramp =
                colorRampPalette(c("black", "grey", color))
        )
        df$dens <- col2rgb(df$col)[1, ] + 1L
    } else {
        # Divide the dataframe according to which
        # color annotation the event has
        dfList <- list()
        for (i in seq_along(uniqueIds)) {
            df <- xYData[idsVector == uniqueIds[i], ]
            colnames(df) <- c("V1", "V2")
            ## Use densCols() output to get density at
            ## each point. The colors here are only
            ## supporting the coming order of the rows
            ## further down the script.
            df$col <- densCols(df$V1, df$V2,
                colramp =
                    colorRampPalette(c(
                        "black",
                        "grey", color[i]
                    ))
            )
            df$dens <- col2rgb(df$col)[1, ] + 1L
            dfList[[i]] <- df
        }

        df <- do.call("rbind", dfList)
    }

    # Plot it, reordering rows so that
    # densest points are plotted on top
    newXyData <- df[order(df$dens), ]
    dPlotCoFunction(
        colorVariable = newXyData$col, plotName = plotName,
        xYData = newXyData[, seq_len(2)], title = title,
        densContour = densContour, bandColor = bandColor,
        dotSize = dotSize, plotDir = plotDir,
        createOutput = createOutput
    )
}
