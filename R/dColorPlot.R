#' Display third variable as color on 2D plot
#'
#'
#' Function to overlay one variable for a set of observations on a field created by two other variables known for the same observations. The plot is constructed primarily for displaying variables on 2D-stochastic neighbour embedding fields, but can be used for any sets of (two or) three variables known for the same observations. As the number of datapoints is often very high, the files would, if saved as pdf of another vector based file type become extremely big. For this reason, the plots are saved as jpeg and no axes or anything alike are added, to simplify usage in publications.
#' @importFrom gplots rich.colors
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
#' @param colorData A vector or a dataframe of numeric observations that will be displayed as color on the plot.
#' @param xYData These variables create the field on which the colorData will be displayed. It needs to be a dataframe with two columns and the same number of rows as the colorData object.
#' @param names The name(s) for the plots. The default alternative, "default" returns the column names of the colorData object in the case this is a dataframe and otherwise returns the somewhat generic name "testVariable". It can be substitutet with a string (in the case colorData is a vector) or vector of strings, as long as it has the same length as the number of columns in colorData.
#' @param densContour An object to create the density contours for the plot. Three possible values: 
#' \describe{
#'               \item{densContour}{A densContour object generated previously with dContours}
#'               \item{TRUE}{a densContour object will be generated internally}
#'               \item{FALSE}{No density contours will be displayed.}
#'              }
#' If not present, it will be generated with the xYData. Useful when only a subfraction of a dataset is plotted, and a superimposition of the distribution of the whole dataset is of interest.
#' @param addLegend If this is set to true, a separate legend plot is produced. This is most useful when the color data contains specific info about separate ids, such as clusters. Default is FALSE.
#' @param idsVector If a legend is added, this argument controls the naming in the legend.
#' @param drawColorPalette If a separate plot with the color palette used for the plots should be printed and saved.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved a jpeg, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @param multiCore If the algorithm should be performed on multiple cores. This increases the speed if the dataset is medium-large (>100000 rows) and has at least 5 columns. Default is true, as it only affects datasets with more than one column.
#' @seealso \code{\link{dDensityPlot}}, \code{\link{dResidualPlot}}, \code{\link{dWilcoxPlot}}, \code{\link{dColorVector}}
#' @return Plots showing the colorData displayed as color on the field created by xYData.
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateBimodalData(samplings=3, observations=3000)
#'
#' #Scale the data 
#' x_scaled <- dScale(x=x[2:ncol(x)])
#' 
#' #Run Barnes Hut tSNE on this. 
#' library(Rtsne.multicore)
#' xSNE <- Rtsne.multicore(x_scaled, pca=FALSE)
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the function for all the variables
#' dColorPlot(colorData=x_scaled, xYData=as.data.frame(xSNE$Y), drawColorPalette=TRUE)
#'
#' #Create a color vector and display it on the SNE field.
#' xColor <- dColorVector(x[,1], colorScale="plasma")
#' dColorPlot(colorData=xColor, xYData=as.data.frame(xSNE$Y), names="separate samplings", addLegend=TRUE, idsVector=x[,1])
#' 
#' @export dColorPlot
dColorPlot <- function(colorData, xYData,  names="default", densContour=TRUE, addLegend=FALSE, idsVector, drawColorPalette=FALSE, title=FALSE, createDirectory=TRUE, directoryName="Variables displayed as color on SNE field", bandColor="black", dotSize=400/sqrt(nrow(xYData)), multiCore=TRUE){

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

  #If there is no matrix present to construct the contour lines and these are wanted, create the density matrix for xYData to make them.
  if(is.logical(densContour)==TRUE){
    if(densContour==TRUE){
      densContour <- dContours(xYData)      
    }
  }
  
  if(drawColorPalette==TRUE){
    pdf("palette.pdf")
    palette(rev(rich.colors(100, plot=TRUE)))
    dev.off()
  }


  xYDataFraction <- dScale(xYData, scale=c(0,1), robustVarScale=FALSE, center=FALSE)
  
  if(class(colorData)=="numeric"){
    colorDataTruncated <- truncateData(colorData)
    colorDataPercent <- dScale(colorDataTruncated, scale=c(0,1), robustVarScale=FALSE, center=FALSE, multiplicationFactor=100)
    colorVector <- dColorVector(round(colorDataPercent), colorScale="rich.colors", order=c(1:100))
    dColorPlotCoFunction(colorVariable=colorVector, name=names, xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize, drawColorPalette=drawColorPalette)
  }
  if(class(colorData)=="data.frame"){
    colorDataTruncated <- apply(colorData, 2, truncateData)
    colorDataPercent <- apply(colorDataTruncated, 2, dScale, scale=c(0,1), robustVarScale=FALSE, center=FALSE, multiplicationFactor=100)
    colorVectors <- apply(round(colorDataPercent), 2, dColorVector, colorScale="rich.colors", order=c(1:100))
    if(multiCore==TRUE){
      no_cores <- detectCores() - 1
      cl = makeCluster(no_cores, type = "SOCK")
      registerDoSNOW(cl)
      foreach(i=1:ncol(colorVectors), .inorder=FALSE) %dopar% dColorPlotCoFunction(colorVariable=colorVectors[,i], name=names[i], xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize)
      stopCluster(cl)
    } else {
       mapply(dColorPlotCoFunction, as.data.frame.matrix(colorVectors, stringsAsFactors =FALSE), names, MoreArgs=list(xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize))
    }
  }
  if(class(colorData)=="character"){
    dColorPlotCoFunction(colorVariable=colorData, name=names, xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize)
  }

  if(addLegend==TRUE){

    #Some preparations for the legend
    #Create a dataframe from the ids and the color vectors
    colorIdsDataFrame <- data.frame(unique(colorData), unique(idsVector), stringsAsFactors = FALSE)
    colorIdsDataFrame <- colorIdsDataFrame[order(colorIdsDataFrame[,2]),]
    pdf(paste("Legend for ", names, ".pdf", sep=""))
    plot.new()
    legend("center",legend = colorIdsDataFrame[,2], col=colorIdsDataFrame[,1], cex=15/length(unique(idsVector)), pch=19)
    dev.off()
  }
  
  
  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }

}

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
