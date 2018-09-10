#' Display density on 2D plot
#'
#'
#' Function to show density for a set of observations on a field created by two variables. The plot is constructed primarily for displaying density of 2D-stochastic neighbour embedding fields, but can be used for any sets of two known for the same observations. As the number of datapoints is often very high, the files would, if saved as pdf of another vector based file type become big. For this reason, the plots are saved as jpeg and no axes or anything alike are added, to simplify usage in publications.
#' @param xYData A dataframe with two columns. Each row contains information about the x and y positition in the field for that observation.
#' @param color This gives the specific color of the plot. It has three potential return values:
#' \describe{
#'     \item{A specific color, e.g. "red" or "#FF0000"}{This will be the color in the densest part of the plot(s)}
#'     \item{"rainbowCols"}{Creates a specific set of colors similar to the rainbow.}
#'     \item{A vector of colors, the same length as nrow(xYData)}{Here, a separate density estimate will be performed for each color subcompartment of the xYData and these will be plotted together.}
#' }
#' @param commonName A name that is common to all density plots created. It can be the groups name, e.g. "Malaria patients" or "Clusters". If only one plot is created, the name is still taken from here.
#' @param plotEachIdSeparately Separates the xYData into subframes specified in the idsVector and plots all ids together but with separate colors. Especially useful when all clusters in a simultaneous analysis should be shown. Colors are inherited from color. Defaults to FALSE.
#' @param idsVector Needed in two situations:
#' \describe{
#'     \item{When "plotEachIdSeparately"=TRUE}{Provides information about which rows that belong to which id and the names for the individual plots.}
#'     \item{When "color" is a vector of colors}{The ids are used to create the legend.}
#' }

#' @param densContour An object to create the density contours for the plot. Three possible values: 
#' \describe{
#'               \item{densContour}{A densContour object generated previously with dContours}
#'               \item{TRUE}{a densContour object will be generated internally}
#'               \item{FALSE}{No density contours will be displayed.}
#'              }
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved as a png, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @param scalingControl A dataframe with two columns. This argument is only necessary if the xYData should be scaled to another range than its internal maximum and minimum values.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.

#' @seealso \code{\link{dColorPlot}}, \code{\link{dResidualPlot}}, \code{\link{dWilcox}}
#' @return Plots showing the densities of the specific xYData (subset) displayed as color on the field created by the same xYData (subset).
#' @examples
#' #Generate a dataframe with bimodally distributed data and a few separate subsamplings
#' x <- generateBimodalData(samplings=5, observations=500)
#'
#' #Run Barnes Hut tSNE on this. 
#' library(Rtsne.multicore)
#' xSNE <- Rtsne.multicore(x[,2:ncol(x)], pca=FALSE)
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Plot all ids together and use a fixed color
#' dDensityPlot(xYData=as.data.frame(xSNE$Y), commonName="All_samplings", 
#' color="blue", createDirectory=FALSE)
#'
#' #Now plot each id separately using a predefined colorscale separating each cluster
#' xColor <- dColorVector(x[,1], colorScale="plasma")
#' dDensityPlot(xYData=as.data.frame(xSNE$Y), color=xColor, plotEachIdSeparately=TRUE, 
#' idsVector=x[,1], commonName="sampling")
#' 
#' #Now all clusters are plotted together using the same predefined colorscale
#' dDensityPlot(xYData=as.data.frame(xSNE$Y), color=xColor, idsVector=x[,1],
#' commonName="all samplings")
#'
#' @export dDensityPlot
dDensityPlot <- function(xYData, color=c("blue", "rainbowCols", "a colorVector"), commonName, plotEachIdSeparately=FALSE, idsVector, densContour=TRUE, title=FALSE, createDirectory=TRUE, directoryName=paste("Density plots for ", commonName, "s", sep=""), scalingControl,  bandColor="black", dotSize=400/sqrt(nrow(xYData))){

  if(createDirectory==TRUE){
    dir.create(directoryName)
    workingDirectory <- getwd()
    setwd(paste(workingDirectory, directoryName, sep="/"))
  }

	if(missing(scalingControl)){
		scalingControl <- xYData
	}

  xYDataScaled <- dScale(xYData, scalingControl, scale=c(0,1), robustVarScale=FALSE, center=FALSE)

  #If there is no matrix present to construct the contour lines and these are wanted, create the density matrix for xYData to make them.
  if(is.logical(densContour)==TRUE){
    if(densContour==TRUE){
      densContour <- dContours(xYData)
    }
  }

  if(length(color)==1){
    #Here, the colors are defined
    ifelse(color=="rainbowCols", cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                                             "#FCFF00", "#FF9400", "#FF3100"))(256), cols <- colorRampPalette(c("black", "grey", color))(256))

    if(plotEachIdSeparately==FALSE){
      dDensityPlotCoFunction(xYDataScaled=xYDataScaled, cols=cols, name=commonName, densContour=densContour, bandColor=bandColor, dotSize=dotSize, title=title)
    }   

    if(plotEachIdSeparately==TRUE){
      uniqueIds <- unique(idsVector)
      
      for (i in 1:length(uniqueIds)){
        dDensityPlotCoFunction(xYDataScaled=xYDataScaled[idsVector==uniqueIds[i],], cols=cols, name=paste(commonName, uniqueIds[i], "density", sep = " "), densContour=densContour, bandColor=bandColor, dotSize=dotSize, title=title)
      }
      
    }
      
  }
  
  if(length(color)>1){
    #Here, the colors are defined
    colors <- sapply(unique(color), as.character)
    colorList <- list()
    for(i in 1:length(colors)){
      colorList[[i]] <- colorRampPalette(c("black", "grey", colors[i]))(256)
    }
    colorList[[length(colors)+1]] <- colors
    colorList[[length(colors)+2]] <- color
  
    if(plotEachIdSeparately==FALSE && length(color)>1){
      dDensityPlotCoFunction(xYDataScaled=xYDataScaled, multipleColors=TRUE, colorList=colorList, name=commonName, densContour=densContour, bandColor=bandColor, dotSize=dotSize, title=title)

      #Some preparations for the legend
      #Create a dataframe from the ids and the color vectors
      colorIdsDataFrame <- data.frame(unique(color), unique(idsVector), stringsAsFactors = FALSE)
      colorIdsDataFrame <- colorIdsDataFrame[order(colorIdsDataFrame[,2]),]
      pdf(paste("Legend for ", commonName, ".pdf", sep=""))
      plot.new()
      legend("center",legend = colorIdsDataFrame[,2], col=colorIdsDataFrame[,1], cex=15/length(unique(idsVector)), pch=19)
      dev.off()
    }

    if(plotEachIdSeparately==TRUE){
      uniqueIds <- unique(idsVector)
      
      for (i in 1:length(uniqueIds)){
        dDensityPlotCoFunction(xYDataScaled=xYDataScaled[idsVector==uniqueIds[i],], cols=colorList[[i]], name=paste(commonName, uniqueIds[i], "density", sep = " "), densContour=densContour, bandColor=bandColor, dotSize=dotSize, title=title)
      }
      
    }
        
  }

  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }

}

