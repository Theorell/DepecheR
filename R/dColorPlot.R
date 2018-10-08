#' Display third variable as color on a 2D plot
#'
#'
#' Function to overlay one variable for a set of observations on a field created by two other variables known for the same observations. The plot is constructed primarily for displaying variables on 2D-stochastic neighbour embedding fields, but can be used for any sets of (two or) three variables known for the same observations. As the number of datapoints is often very high, the files would, if saved as pdf of another vector based file type become extremely big. For this reason, the plots are saved as jpeg and no axes or anything alike are added, to simplify usage in publications.
#' @importFrom gplots rich.colors
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
#' @param colorData A vector, matrix or dataframe of numeric observations that will be displayed as color on the plot.
#' @param controlData Optional. A numeric/integer vector or dataframe of values that could be used to define the range of the colorData. If no control data is present, the function defaults to using the colorData as control data.
#' @param xYData These variables create the field on which the colorData will be displayed. It needs to be a matrix or dataframe with two columns and the same number of rows as the colorData object.
#' @param names The name(s) for the plots. The default alternative, "default" returns the column names of the colorData object in the case this is a dataframe and otherwise returns the somewhat generic name "testVariable". It can be substituted with a string (in the case colorData is a vector) or vector of strings, as long as it has the same length as the number of columns in colorData.
#' @param densContour If density contours should be created for the plot(s) or not. Defaults to TRUE. a
#' @param addLegend If this is set to true, a separate legend plot is produced. This is most useful when the color data contains specific info about separate ids, such as clusters. Default is FALSE.
#' @param idsVector If a legend is added, this argument controls the naming in the legend.
#' @param drawColorPalette If a separate plot with the color palette used for the plots should be printed and saved.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved a jpeg, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @param truncate If truncation of the most extreme values should be performed for the visualizations. Three possible values: TRUE, FALSE, and a vector with two values indicating the low and high threshold quantiles for truncation. 
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @param multiCore If the algorithm should be performed on multiple cores. This increases the speed if the dataset is medium-large (>100000 rows) and has at least 5 columns. Default is TRUE when the rows exceed 100000 rows and FALSE otherwise.
#' @param createPlot For testing purposes. Defaults to TRUE. If FALSE, no plots are generated.
#' @seealso \code{\link{dDensityPlot}}, \code{\link{dResidualPlot}}, \code{\link{dWilcox}}, \code{\link{dColorVector}}
#' @return Plots showing the colorData displayed as color on the field created by xYData.
#' @examples
#' #Load some data
#' data(testData)
#' 
#' #Run Barnes Hut tSNE on this. For more rapid example execution, a pre-run SNE is inluded
#' #library(Rtsne)
#' #testDataSNE <- Rtsne(testData[,2:15], pca=FALSE)
#' data(testDataSNE)
#'
#' #Run the function for two of the variables
#' dColorPlot(colorData=testData[2:3], xYData=testDataSNE$Y, drawColorPalette=TRUE)
#'
#' #Create a color vector and display it on the SNE field. For this purpose,
#' #four individual donors are extracted from the test dataset
#' testDataSubset <- rbind(testData[1:2000,], testData[95001:97000,])
#' testDataSNESubset <- rbind(testDataSNE$Y[1:2000,], testDataSNE$Y[95001:97000,])
#' testColor <- dColorVector(testDataSubset$ids, colorScale="plasma")
#' dColorPlot(colorData=testColor, xYData=testDataSNESubset, 
#'            names="separate samplings", addLegend=TRUE, idsVector=testDataSubset$ids)
#' 
#' 
#' @export dColorPlot
dColorPlot <- function(colorData, controlData, xYData,  names="default", densContour=TRUE, addLegend=FALSE, idsVector, drawColorPalette=FALSE, title=FALSE, createDirectory=TRUE, directoryName="Variables displayed as color on SNE field", truncate=TRUE, bandColor="black", dotSize=500/sqrt(nrow(xYData)), multiCore="default", createPlot=TRUE){

  if(class(colorData)=="matrix"){
    colorData <- as.data.frame(colorData)
  }
  
  if(class(colorData)!="numeric" && class(colorData)!="data.frame" && class(colorData)!="character"){
    stop("ColorData needs to be either a numeric vector, a character vector of colors or a matrix or dataframe of numbers.")
  }

  if(class(xYData)=="matrix"){
    xYData <- as.data.frame(xYData)
  }
  if(ncol(xYData)!=2){
    stop("xYData needs to contain two vectors to be displayed.")
  }

  if(createDirectory==TRUE){
    workingDirectory <- getwd()
    dir.create(directoryName)
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
  
  if(missing(controlData)==TRUE){
    controlData <- colorData
  }

  #Create the density matrix for xYData.
  if(logial(densContour)){
    if(densContour==TRUE){
      densContour <- dContours(xYData)
    }
  }  
  
  if(drawColorPalette==TRUE && createPlot==TRUE){
    pdf("palette.pdf")
    palette(rev(rich.colors(100, plot=TRUE)))
    dev.off()
  }


  xYDataFraction <- dScale(xYData, scale=c(0,1), robustVarScale=FALSE, center=FALSE)
  
  if(class(colorData)=="numeric"){
    colorDataPercent <- dScale(colorData, control=controlData, scale=c(0,1), robustVarScale=FALSE, center=FALSE, multiplicationFactor=100, truncate=truncate)
    colorVector <- dColorVector(round(colorDataPercent), colorScale="rich_colors", order=c(0:100))
    dColorPlotCoFunction(colorVariable=colorVector, name=names, xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize, createPlot=createPlot)
  }
  if(class(colorData)=="data.frame"){
    colorDataPercent <- dScale(x=colorData, control=controlData, scale=c(0,1), robustVarScale=FALSE, center=FALSE, multiplicationFactor=100, truncate=truncate)
    colorVectors <- apply(round(colorDataPercent), 2, dColorVector, colorScale="rich_colors", order=c(0:100))
    if(multiCore=="default"){
      if(nrow(colorData)>100000){
        multiCore <- TRUExYDataFraction <- dScale(xYData, scale=c(0,1), robustVarScale=FALSE, center=FALSE)

      } else {
        multiCore <- FALSE
      }
    }
    if(multiCore==TRUE){
      no_cores <- detectCores() - 1
      cl = makeCluster(no_cores, type = "SOCK")
      registerDoSNOW(cl)
      foreach(i=1:ncol(colorVectors), .inorder=FALSE) %dopar% dColorPlotCoFunction(colorVariable=colorVectors[,i], name=names[i], xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize, createPlot=createPlot)
      stopCluster(cl)
    } else {
       mapply(dColorPlotCoFunction, as.data.frame.matrix(colorVectors, stringsAsFactors =FALSE), names, MoreArgs=list(xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize, createPlot=createPlot))
    }
  }
  if(class(colorData)=="character"){
    dColorPlotCoFunction(colorVariable=colorData, name=names, xYDataFraction=xYDataFraction, title=title, densContour=densContour, bandColor=bandColor, dotSize=dotSize, createPlot=createPlot)
  }

  if(addLegend==TRUE && createPlot==TRUE){

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

  print(paste0("Files were saved at ", getwd()))
}