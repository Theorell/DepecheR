#This function is used by dDensityPlot. 
#It calculates the density at each point of a 2-dimensional field of
#data, creates a color density gradient based on this, and resorts the events
#so that the lowest density areas will be plotted first, and so on. 
#The multiple color option allows for multiple separate densities to be plotted
#simulataneously, for example if multiple clusters should be displayed together
#with density. 
#"Unique" parameters
#multipleColors: if multiple densities should be calculated in the sample plot. 
#cols: The set of colors to use for the density distribution. 
#colorList if multipleColors==TRUE, this provides the list of colors. See cols.
#name: see commonName. If multiple ids are plotted individually, the name is 
#composed of the unique id and the commonName. 
#For information on the other parameters, see dDensityPlot

dDensityPlotCoFunction <- function(xYData, multipleColors = FALSE, cols, 
                                   colorList, name, densContour, bandColor, 
                                   dotSize, title, createDirectory, 
                                   directoryName, createOutput = FALSE) {

    if (multipleColors == FALSE) {
        V1 <- xYData[, 1]
        V2 <- xYData[, 2]
        df <- data.frame(V1, V2)
        
        ## Use densCols() output to get density at
        ## each point. The colors here are only
        ## supporting the coming order of the rows
        ## further down the script.
        x <- densCols(V1, V2, colramp = colorRampPalette(c("black", "white")))
        df$dens <- col2rgb(x)[1, ] + 1L
        df$col <- cols[df$dens]
    }
    
    if (multipleColors == TRUE) {
        # Divide the dataframe according to which
        # color annotation the event has
        colors <- colorList[[length(colorList) - 1]]
        color <- colorList[[length(colorList)]]
        dfList <- list()
        for (i in seq_along(colors)) {
            V1 <- xYData[color == colors[i],1]
            V2 <- xYData[color == colors[i],2]
            df <- data.frame(V1, V2)
            cols <- colorList[[i]]
            ## Use densCols() output to get density at
            ## each point. The colors here are only
            ## supporting the coming order of the rows
            ## further down the script.
            x <- densCols(V1, V2, colramp = colorRampPalette(c("black",
                                                               "white")))
            df$dens <- col2rgb(x)[1, ] + 1L
            df$col <- cols[df$dens]
            dfList[[i]] <- df
        }
        
        df <- do.call("rbind", dfList)
    }
    
    # Plot it, reordering rows so that
    # densest points are plotted on top
    newXyData <- df[order(df$dens), ]
    dPlotCoFunction(colorVariable = newXyData$col, name = name, 
                    xYData = newXyData[,seq_len(2)], title = title, 
                    densContour = densContour, bandColor = bandColor, 
                    dotSize = dotSize, 
                    createDirectory = createDirectory, 
                    directoryName = directoryName, 
                    createOutput = createOutput)
}
