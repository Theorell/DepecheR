#' Display third variable as color on 2D plot
#'
#'
#' Function to overlay one variable for a set of observations on a field created by two other variables known for the same observations. The plot is constructed primarily for displaying variables on 2D-stochastic neighbour embedding fields, but can be used for any sets of (two or) three variables known for the same observations. As the number of datapoints is often very high, the files would, if saved as pdf of another vector based file type become extremely big. For this reason, the plots are saved as jpeg and no axes or anything alike are added, to simplify usage in publications.
#' @importFrom ggplot2 ggplot aes geom_point geom_density_2d xlim ylim ggtitle theme element_blank element_text ggplot_gtable ggplot_build element_rect ggsave
#' @importFrom grid grid.draw
#' @param colorData A vector or a dataframe of numeric observations that will be displayed as color on the plot.
#' @param names The name(s) for the plots. The default alternative, "default" returns the column names of the colorData object in the case this is a dataframe and otherwise returns the somewhat generic name "testVariable". It can be substitutet with a string (in the case colorData is a vector) or vector of strings, as long as it has the same length as the number of columns in colorData.
#' @param xYData These variables create the field on which the colorData will be displayed. It needs to be a dataframe with two columns and the same number of rows as the colorData object.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved a jpeg, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @param drawColorPalette If a separate plot with the color palette used for the plots should be printed and saved.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @seealso \code{\link{depecheDensity}}, \code{\link{depecheResidual}}, \code{\link{depecheWilcox}}
#' @return Plots showing the colorData displayed as color on the field created by xYData.
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData()
#'
#' #Scale the data (not actually necessary in this artificial example due to the nature of the generated data)
#' x_scaled <- quantileScale(x[2:ncol(x)])
#'
#' #Run Barnes Hut tSNE on this. NB! This takes quite a while (a minute or so) as the algorithm is slow.
#' library(Rtsne)
#' xSNE <- Rtsne(x_scaled, pca=FALSE)
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the function for all the variables
#' depecheColor(colorData=x_scaled, xYData=as.data.frame(xSNE$Y), drawColorPalette=TRUE)
#'
#' @export depecheColor
depecheColor <- function(colorData, names="default", xYData,  title=FALSE, dotSize=20000/nrow(colorData), drawColorPalette=FALSE, createDirectory=TRUE, directoryName="Variables displayed as color on SNE field"){

  if(class(colorData)!="numeric" && class(colorData)!="data.frame"){
    stop("colorData needs to be either a numeric vector or a dataframe. Change the class and try again.")
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

  if(names=="default" && class(colorData)=="data.frame"){
    names <- colnames(colorData)
  }

  if(drawColorPalette==TRUE){
    pdf("palette.pdf")
    palette(rev(rich.colors(100, plot=TRUE)))
    dev.off()
  }

  if(drawColorPalette==FALSE){
    palette(rev(rich.colors(100, plot=FALSE)))
  }

  xYDataPercent <- minMaxScale(xYData, multiplicationFactor=100)

  if(class(colorData)=="numeric"){
    depecheColorCoFunction(colorVariable=colorData, name=names, xYDataPercent=xYDataPercent, title=title, dotSize=dotSize)
  }
  if(class(colorData)=="data.frame"){
    mapply(depecheColorCoFunction, colorData, names, MoreArgs=list(xYDataPercent=xYDataPercent, title=title, dotSize=dotSize))
  }

  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }

}

depecheColorCoFunction <- function(colorVariable, name, xYDataPercent, title=FALSE, dotSize=20000/nrow(colorVariable)){
	fname = paste(name,'.png', sep="")

	colorVariablePercent <- minMaxScale(colorVariable, multiplicationFactor=100)

colnames(xYDataPercent) <- c("V1", "V2")

allData <- cbind(colorVariablePercent, xYDataPercent)

if(title==TRUE){
  p <- ggplot(allData,aes(x=V1,y=V2)) +
    geom_point(col=(1 + 0.98*(102-allData[,1])), size=dotSize, alpha=0.7) 	+ geom_density_2d(col="black", size=I(0.6)) + xlim(-5, 105) + ylim(-5, 105) + ggtitle(name) +
    theme (rect = element_rect(fill = "transparent"),
            line = element_blank(),
            title = element_text(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.background = element_rect(fill = "transparent"))
}

if(title==FALSE){
  p <- ggplot(allData,aes(x=V1,y=V2)) +
    geom_point(col=(1 + 0.98*(102-allData[,1])), size=dotSize, alpha=0.7) 	+ geom_density_2d(col="black", size=I(0.6)) + xlim(-5, 105) + ylim(-5, 105) +
    theme (rect = element_rect(fill = "transparent"),
            line = element_blank(),
            title = element_text(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.background = element_rect(fill = "transparent"))
}

gt <- ggplot_gtable(ggplot_build(p))
ge <- subset(gt$layout, name == "panel")
grid.draw(gt[ge$t:ge$b, ge$l:ge$r])
ggsave(filename = fname, dpi=300,  bg = "transparent")
}
