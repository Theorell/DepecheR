dColorPlotCoFunction <- function(colorVariable, name, xYDataFraction, title=FALSE, densContour, bandColor, dotSize){
  
  colnames(xYDataFraction) <- c("V1", "V2")
  
  png(paste(name, ".png", sep=""), width = 2500, height = 2500, units = "px", bg="transparent")
  # Plot it
  if(title==TRUE){
    plot(V2~V1, data=xYDataFraction, main=name, pch=20, cex=dotSize, cex.main=5, col=colorVariable, xlim=c(-0.05, 10.5), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }
  if(title==FALSE){
    plot(V2~V1, data=xYDataFraction, main=NULL, pch=20, cex=dotSize, cex.main=5, col=colorVariable, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }
  if(length(densContour)>1){
    par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
    contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")
  } 
  dev.off()
}
