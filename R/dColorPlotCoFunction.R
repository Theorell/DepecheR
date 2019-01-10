dColorPlotCoFunction <- function(colorVariable, name, xYData, 
    title = FALSE, densContour, bandColor, dotSize, createPlot = TRUE) {
    colnames(xYData) <- c("V1", "V2")
    
    
    if (length(densContour) > 1) {
        xlim <- c(min(densContour[[1]]), max(densContour[[1]]))
        ylim <- c(min(densContour[[2]]), max(densContour[[2]]))
    } else {
        minX <- min(xYData[, 1])
        maxX <- max(xYData[, 1])
        minY <- min(xYData[, 2])
        maxY <- max(xYData[, 2])
        xlim <- c(minX - abs(minX * 0.05), maxX + abs(maxX * 
            0.05))
        ylim <- c(minY - abs(minY * 0.05), maxY + abs(maxY * 
            0.05))
    }
    
    png(paste(name, ".png", sep = ""), width = 2500, height = 2500, 
        units = "px", bg = "transparent")
    
    # Plot it
    if (createPlot == TRUE) {
        if (title == TRUE) {
            plot(V2 ~ V1, data = xYData, main = name, pch = 20, 
                cex = dotSize, cex.main = 5, col = colorVariable, 
                xlim = xlim, ylim = ylim, axes = FALSE, xaxs = "i", 
                yaxs = "i")
        }
        if (title == FALSE) {
            plot(V2 ~ V1, data = xYData, main = NULL, pch = 20, 
                cex = dotSize, cex.main = 5, col = colorVariable, 
                xlim = xlim, ylim = ylim, axes = FALSE, xaxs = "i", 
                yaxs = "i")
        }
        if (length(densContour) > 1) {
            par(fig = c(0, 1, 0, 1), mar = c(6, 4.5, 4.5, 2.5), 
                new = TRUE)
            contour(x = densContour$x, y = densContour$y, z = densContour$z, 
                xlim = xlim, ylim = ylim, nlevels = 10, col = bandColor, 
                lwd = 8, drawlabels = FALSE, axes = FALSE, xaxs = "i", 
                yaxs = "i")
        }
    }
    dev.off()
}
