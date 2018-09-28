dDensityPlotCoFunction <- function(xYDataScaled, multipleColors=FALSE, cols, colorList, name, densContour, bandColor, dotSize, title, createPlot=TRUE){
  
  if(multipleColors==FALSE){
    
    x1 <- xYDataScaled[,1]
    x2 <- xYDataScaled[,2]
    df <- data.frame(x1,x2)
    
    ## Use densCols() output to get density at each point. The colors here are only supporting the coming order of the rows further down the script.
    x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(x)[1,] + 1L
    df$col <- cols[df$dens]
    
  }
  
  if(multipleColors==TRUE){
    
    #Divide the dataframe according to which color annotation the event has
    colors <- colorList[[length(colorList)-1]]
    color <- colorList[[length(colorList)]]
    dfList <- list()
    for(i in 1:length(colors)){
      
      x1 <- xYDataScaled[color==colors[i],1]
      x2 <- xYDataScaled[color==colors[i],2]
      df <- data.frame(x1,x2)
      cols <- colorList[[i]]
      ## Use densCols() output to get density at each point. The colors here are only supporting the coming order of the rows further down the script.
      x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
      df$dens <- col2rgb(x)[1,] + 1L
      df$col <- cols[df$dens]
      dfList[[i]] <- df
    }
    
    df <- do.call("rbind", dfList)
  }
  
  png(paste(name, ".png", sep=""), width = 2500, height = 2500, units = "px", bg="transparent")
  # Plot it, reordering rows so that densest points are plotted on top
  if(createPlot==TRUE){
    if(title==TRUE){
      plot(x2~x1, data=df[order(df$dens),], main=name, pch=20, cex=dotSize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
    }
    if(title==FALSE){
      plot(x2~x1, data=df[order(df$dens),], main=NULL, pch=20, cex=dotSize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
    }
    if(length(densContour)>1){
      par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
      contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")
    } 
  }
  dev.off()
  
}