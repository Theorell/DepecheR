dDensityPlotCoFunction <- function(xYData, 
    multipleColors = FALSE, cols, colorList, 
    name, densContour, bandColor, dotSize, 
    title, createDirectory, directoryName,
    createPlot = FALSE) {
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
    
    if (multipleColors == FALSE) {
        x1 <- xYData[, 1]
        x2 <- xYData[, 2]
        df <- data.frame(x1, x2)
        
        ## Use densCols() output to get density at
        ## each point. The colors here are only
        ## supporting the coming order of the rows
        ## further down the script.
        x <- densCols(x1, x2, colramp = colorRampPalette(c("black", 
            "white")))
        df$dens <- col2rgb(x)[1, ] + 1L
        df$col <- cols[df$dens]
    }
    
    if (multipleColors == TRUE) {
        # Divide the dataframe according to which
        # color annotation the event has
        colors <- colorList[[length(colorList) - 
            1]]
        color <- colorList[[length(colorList)]]
        dfList <- list()
        for (i in seq_along(colors)) {
            x1 <- xYData[color == colors[i],1]
            x2 <- xYData[color == colors[i],2]
            df <- data.frame(x1, x2)
            cols <- colorList[[i]]
            ## Use densCols() output to get density at
            ## each point. The colors here are only
            ## supporting the coming order of the rows
            ## further down the script.
            x <- densCols(x1, x2, colramp = colorRampPalette(c("black",
                                                               "white")))
            df$dens <- col2rgb(x)[1, ] + 1L
            df$col <- cols[df$dens]
            dfList[[i]] <- df
        }
        
        df <- do.call("rbind", dfList)
    }
    
    if (createDirectory == TRUE) {
        png(file.path(directoryName, paste0(name, ".png")), width = 2500, 
                      height = 2500, units = "px", bg = "transparent")
    } else {
        png(paste0(name, ".png"), width = 2500, 
            height = 2500, units = "px", bg = "transparent")
    }
    
    
    # Plot it, reordering rows so that
    # densest points are plotted on top
    if (createPlot == TRUE) {
        if (title == TRUE) {
            plot(x2 ~ x1, data = df[order(df$dens), 
                ], main = name, pch = 20, 
                cex = dotSize, cex.main = 5, 
                col = col, xlim = xlim, ylim = ylim, 
                axes = FALSE, xaxs = "i", 
                yaxs = "i")
        }
        if (title == FALSE) {
            plot(x2 ~ x1, data = df[order(df$dens), 
                ], main = NULL, pch = 20, 
                cex = dotSize, col = col, 
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
}
