#' Display density on 2D plot
#'
#'
#' Function to show density for a set of observations on a field created by two variables. The plot is constructed primarily for displaying density of 2D-stochastic neighbour embedding fields, but can be used for any sets of two known for the same observations. As the number of datapoints is often very high, the files would, if saved as pdf of another vector based file type become big. For this reason, the plots are saved as jpeg and no axes or anything alike are added, to simplify usage in publications.
#' @param xYData A dataframe with two columns. Each row contains information about the x and y positition in the field for that observation.
#' @param scalingControl A dataframe with two columns. This argument is only necessary if the xYData should be scaled to another range than its internal maximum and minimum values.
#' @param plotEachIdSeparately If the xYData is made up of data from multiple, individual subjects, these can be plotted separately. Defaults to FALSE.
#' @param idsVector If "plotEachIdSeparately" is TRUE, then this argument is needed to provide information about which rows that belong to which id.
#' @param commonName A name that is common to all density plots created. It can be the groups name, e.g. "Malaria patients" or "Clusters". If only one plot is created, the name is still taken from here.
#' @param color This gives the specific color of the plot. It has four potential return values. Either it is a specific color ("red" or "#FF0000"), it is "rainbowCols", which creates a specific set of colors similar to the rainbow. The third alternative is "idSpecific", which will hten give each plot a different color depending on their id. The fourth alternative is that a vector of colors (e.g. either "red" or "#FF0000") with the same length as the number of rows in xYData is provided. Here, a separate density estimate will be performed for each color subcompartment of the xYData and these will be plotted together.
#' @param densContour An object to create the density contours for the plot. If not present, it will be generated with the xYData. Useful when only a subfraction of a dataset is plotted, and a superimposition of the distribution of the whole dataset is of interest.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved as a png, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @seealso \code{\link{depecheColor}}, \code{\link{depecheResidual}}, \code{\link{depecheWilcox}}
#' @return Plots showing the densities of the specific xYData (subset) displayed as color on the field created by the same xYData (subset).
#' @examples
#' #Generate a dataframe with bimodally distributed data and a few separate subsamplings
#' x <- generateFlowCytometryData(samplings=5, observations=2000)
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
#' #Plot all ids together and use rainbowColors
#' depecheDensity(xYData=as.data.frame(xSNE$Y), commonName="All_samplings", 
#' color="rainbowCols", createDirectory=FALSE)
#'
#' #Now plot each id separately
#' depecheDensity(xYData=as.data.frame(xSNE$Y), plotEachIdSeparately=TRUE, 
#' idsVector=x[,1], commonName="sampling", color="idSpecific")
#'
#' @export depecheDensity
depecheDensity <- function(xYData, scalingControl, plotEachIdSeparately=FALSE, idsVector, commonName, color=c("blue", "rainbowCols", "idSpecific", "multipleColors"), densContour, bandColor="black", dotSize=400/sqrt(nrow(xYData)), title=FALSE, createDirectory=TRUE, directoryName=paste("Density plots for ", commonName, "s", sep="")){

    if(plotEachIdSeparately==FALSE && color=="idSpecific"){
    stop("A different color can only be created for each object if multiple objects are plotted. Submit an idsVector and turn plotEachIdSeparately to TRUE, ant try again.")
  }

  if(createDirectory==TRUE){
    dir.create(directoryName)
    workingDirectory <- getwd()
    setwd(paste(workingDirectory, directoryName, sep="/"))
  }

	if(missing(scalingControl)){
		scalingControl <- xYData
	}

  xYDataScaled <- quantileScale(xYData, scalingControl, robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE)

  #If there is no matrix present to construct the contour lines, create the density matrix for all_data to make them.
  if(missing("densContour")){
    densContour <- densityContours(xYData)
  }

  if(color!="idSpecific" && length(color)==1){
    #Here, the colors are defined
    ifelse(color=="rainbowCols", cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                                             "#FCFF00", "#FF9400", "#FF3100"))(256), cols <- colorRampPalette(c("black", "grey", color))(256))

    if(plotEachIdSeparately==FALSE){

      depecheDensityCoFunction(xYDataScaled=xYDataScaled, cols=cols, name=commonName, densContour=densContour, bandColor=bandColor, dotSize=dotSize, title=title)
    }

    if(plotEachIdSeparately==TRUE){
      number <- unique(idsVector)

      for (i in 1:length(number)){

        depecheDensityCoFunction(xYDataScaled=xYDataScaled[idsVector==number[i],], cols=cols, name=paste(commonName, number[i], "density", sep = " "), densContour=densContour, bandColor=bandColor, dotSize=dotSize, title=title)
      }
    }
  }
  
  if(plotEachIdSeparately==TRUE && color=="idSpecific"){
    equidistantIds <- turnVectorEquidistant(idsVector)
    normIds <- quantileScale(equidistantIds, robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE, multiplicationFactor=100)
    normNumber <- unique(normIds)
    palette(rev(rich.colors(100, plot=FALSE)))
    for (i in 1:length(normNumber)){
      color <- (1 + 0.98*(102-normNumber[i]))
      cols <- colorRampPalette(c("black", "grey", color))(256)
      depecheDensityCoFunction(xYDataScaled=xYDataScaled[normIds==normNumber[i],], cols=cols, name=paste(commonName, unique(idsVector)[i], "density", sep = " "), densContour=densContour, bandColor=bandColor, dotSize=dotSize, title=title)
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

      depecheDensityCoFunction(xYDataScaled=xYDataScaled, multipleColors=TRUE, colorList=colorList, name=commonName, densContour=densContour, bandColor=bandColor, dotSize=dotSize, title=title)


  }

  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }

}

depecheDensityCoFunction <- function(xYDataScaled, multipleColors=FALSE, cols, colorList, name, densContour, bandColor, dotSize, title){

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
  if(title==TRUE){
    plot(x2~x1, data=df[order(df$dens),], main=name, pch=20, cex=dotSize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }
if(title==FALSE){
  plot(x2~x1, data=df[order(df$dens),], main=NULL, pch=20, cex=dotSize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
}


par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")

dev.off()

}


