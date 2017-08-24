#' Display third variable as color on 2D plot
#'
#'
#' Function to overlay one variable for a set of observations on a field created by two other variables known for the same observations. The plot is constructed primarily for displaying variables on 2D-stochastic neighbour embedding fields, but can be used for any sets of (two or) three variables known for the same observations. As the number of datapoints is often very high, the files would, if saved as pdf of another vector based file type become extremely big. For this reason, the plots are saved as jpeg and no axes or anything alike are added, to simplify usage in publications.
#' @param colorData A vector or a dataframe of numeric observations that will be displayed as color on the plot.
#' @param names The name(s) for the plots. The default alternative, "default" returns the column names of the colorData object in the case this is a dataframe and otherwise returns the somewhat generic name "testVariable". It can be substitutet with a string (in the case colorData is a vector) or vector of strings, as long as it has the same length as the number of columns in colorData.
#' @param xYData These variables create the field on which the colorData will be displayed. It needs to be a dataframe with two columns and the same number of rows as the colorData object.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved a jpeg, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param densContour An object to create the density contours for the plot. If not present, it will be generated with the xYData. Useful when only a subfraction of a dataset is plotted, and a superimposition of the distribution of the whole dataset is of interest.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @param drawColorPalette If a separate plot with the color palette used for the plots should be printed and saved.
#' @param addLegend If this is set to true, a separate legend plot is produced. This is most useful when the color data contains specific info about separate ids, such as clusters.
#' @param idsVector If a legend is added, this argument controls the naming in the legend.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @seealso \code{\link{depecheDensity}}, \code{\link{depecheResidual}}, \code{\link{depecheWilcox}}, \code{\link{colorVector}}
#' @return Plots showing the colorData displayed as color on the field created by xYData.
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData(samplings=3, observations=3000)
#'
#' #Scale the data (not actually necessary in this artificial 
#' #example due to the nature of the generated data)
#' x_scaled <- quantileScale(x=x[2:ncol(x)])
#' 
#' #Run Barnes Hut tSNE on this. 
#' library(Rtsne.multicore)
#' xSNE <- Rtsne.multicore(x_scaled, pca=FALSE)
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the function for all the variables
#' depecheColor(colorData=x_scaled, xYData=as.data.frame(xSNE$Y), drawColorPalette=TRUE)
#'
#' #Create a color vector and display it on the SNE field.
#' xColor <- colorVector(x[,1])
#' depecheColor(colorData=xColor, xYData=as.data.frame(xSNE$Y), names="separate samplings", addLegend=TRUE, idsVector=x[,1])
#' 
#' @export depecheColor
depecheColor <- function(colorData, names="default", xYData,  title=FALSE, densContour, bandColor="black", dotSize=30000/nrow(xYData), drawColorPalette=FALSE, addLegend=FALSE, idsVector, createDirectory=TRUE, directoryName="Variables displayed as color on SNE field"){

  if(class(colorData)!="numeric" && class(colorData)!="data.frame" && class(colorData)!="character"){
    stop("colorData needs to be either a numeric, vector, a character vector of colors or a dataframe. Change the class and try again.")
  }

  if(class(xYData)!="data.frame"){
    stop("xYData needs to be a dataframe. Change the class and try again.")
  }

  if(createDirectory==TRUE){
    dir.create(directoryName)
    workingDirectory <- getwd()
    setwd(paste(workingDirectory, directoryName, sep="/"))

  }
  if(names=="default" && class(colorData)=="numeric"){
    names <- "testVariable"
  }
  if(names=="default" && class(colorData)=="character"){
    names <- "Ids"
  }

  if(names=="default" && class(colorData)=="data.frame"){
    names <- colnames(colorData)
  }

  #If there is no matrix present to construct the contour lines, create the density matrix for all_data to make them.
  if(missing("densContour")){
    densContour <- densityContours(xYData)
  }
  
  if(class(colorData)!="character" && drawColorPalette==TRUE){
    pdf("palette.pdf")
    palette(rev(rich.colors(100, plot=TRUE)))
    dev.off()
  }

  if(class(colorData)!="character" && drawColorPalette==FALSE){
    palette(rev(rich.colors(100, plot=FALSE)))
  }

  xYDataFraction <- quantileScale(xYData, robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE)
  
  if(class(colorData)=="numeric"){
    depecheColorCoFunction(colorVariable=colorData, name=names, xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize)
  }
  if(class(colorData)=="data.frame"){
    mapply(depecheColorCoFunction, colorData, names, MoreArgs=list(xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize))
  }
  if(class(colorData)=="character"){
    depecheColorCoFunctionSetCols(colorVariable=colorData, name=names, xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize)
  }

  if(addLegend==TRUE){
    pdf(paste("Legend for ", names, ".pdf", sep=""))
    plot.new()
    legend("center",legend = unique(idsVector), col=unique(colorData), cex=5, pch=19)
    dev.off()
    }
  
  
  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }

}

depecheColorCoFunction <- function(colorVariable, name, xYDataFraction, title=FALSE, densContour, bandColor, dotSize){

  colorVariableTruncated <- truncateData(colorVariable)
  colorVariablePercent <- quantileScale(colorVariableTruncated, robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE, multiplicationFactor=100)
  
  colnames(xYDataFraction) <- c("V1", "V2")
    
  png(paste(name, ".png", sep=""), width = 2500, height = 2500, units = "px", bg="transparent")
  # Plot it
  if(title==TRUE){
    plot(V2~V1, data=xYDataFraction, main=name, pch=20, cex=dotSize, cex.main=5, col=(1 + 0.98*(102-colorVariablePercent)), xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }
  if(title==FALSE){
    plot(V2~V1, data=xYDataFraction, main=NULL, pch=20, cex=dotSize, cex.main=5, col=(1 + 0.98*(102-colorVariablePercent)), xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }
  
  
  par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
  contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")
  
  dev.off()
  
}


depecheColorCoFunctionSetCols <- function(colorVariable, name, xYDataFraction, title=FALSE, densContour, bandColor, dotSize){

  colnames(xYDataFraction) <- c("V1", "V2")
  
  png(paste(name, ".png", sep=""), width = 2500, height = 2500, units = "px", bg="transparent")
  # Plot it
  if(title==TRUE){
    plot(V2~V1, data=xYDataFraction, main=name, pch=20, cex=dotSize, cex.main=5, col=colorVariable, xlim=c(-0.05, 10.5), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }
  if(title==FALSE){
    plot(V2~V1, data=xYDataFraction, main=NULL, pch=20, cex=dotSize, cex.main=5, col=colorVariable, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }
  
  par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
  contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")
  
  dev.off()
  
  }
