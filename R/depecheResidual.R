#' Showing the residuals when subtracting the values from one group from another on a sne plot
#'
#'
#' This function is used to compare groups of individuals from whom comparable cytometry or other complex data has been generated where the number of individuals does not allow any statistical comparisons.
#' @param groupClusterData This dataframe should contain at least two columns  with information about the cluster (named "cluster") and the group identity (named "group") for each observation.
#' @param xYData A dataframe with two columns. Each row contains information about the x and y positition in the field for that observation. It needs to have the same number of rows as groupClusterData.
#' @param densContour An object to create the density contours for the plot. If not present, it will be generated with the xYData. Useful when only a subfraction of a dataset is plotted, and a superimposition of the distribution of the whole dataset is of interest.
#' @param name The main name for the graph and the analysis.
#' @param groupName1 The name for the first group
#' @param groupName2 The name for the second group
#' @param maxAbsPlottingValues If multiple plots should be compared, it might be useful to define a similar color scale for all plots, so that the same color always means the same value. Such a value can be added here. It defaults to the maximum Wilcoxon statistic that is generated in the analysis.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved as a png, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @seealso \code{\link{depecheColor}}, \code{\link{depecheDensity}}, \code{\link{depecheWilcox}}
#' @return A sne based plot showing which events that belong to a cluster dominated by the first or the second group.
#' @examples
#' #Generate a dataframe with bimodally distributed data and 2 subsamplings.
#' x <- generateFlowCytometryData(samplings=2, ncols=7, observations=1000)
#'
#' #Scale the data (not actually necessary in this artificial 
#' #example due to the nature of the generated data)
#' x_scaled <- quantileScale(x[2:ncol(x)])
#'
#' #Create the optimized number of clusters for this dataset
#' x_optim <- pKMOptim(x_scaled, iterations=50, bootstrapObservations=1000)
#' x_pKM <- pKMRun(x_scaled, regVec=x_optim[[1]][["optimalRegularizationValue"]], 
#' withOrWithoutZeroClust=x_optim[[1]][["withOrWithoutZeroClust"]], iterations=1, ids=x[,1])
#'
#' #Run Barnes Hut tSNE on this. 
#' library(Rtsne.multicore)
#' xSNE <- Rtsne.multicore(x_scaled, pca=FALSE)
#'
#' #Create the groupClusterData object
#' groupClusterData <- as.data.frame(cbind(x[,1], x_pKM$clusterVector))
#' colnames(groupClusterData) <- c("group", "cluster")
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #And finally run the function
#' depecheResidual(groupClusterData=groupClusterData, xYData=as.data.frame(xSNE$Y))
#' @export depecheResidual
depecheResidual <- function(groupClusterData, xYData, densContour, name="depecheResidual", groupName1="Group 1", groupName2="Group 2", title=FALSE,  maxAbsPlottingValues, bandColor="black", createDirectory=FALSE, directoryName="depecheResidual", dotSize=400/sqrt(nrow(xYData))){

  if(createDirectory==TRUE){
    dir.create(directoryName)
    workingDirectory <- getwd()
    setwd(paste(workingDirectory, directoryName, sep="/"))

  }

  if(length(unique(groupClusterData$group))!=2){
    stop("More or less than two groups are present. Please correct this.")
  }

  #Here, the residuals are identified.
  #A table with the percentage of cells in each cluster for each group is created in analogy with XXX pKMRun.

  clusterTable <- table(groupClusterData$cluster, groupClusterData$group)

  countTable <- table(groupClusterData$group)

  clusterPercentagesForGroups <- clusterTable

  for(i in 1:length(countTable)){
     x <- 100*clusterTable[,i]/countTable[i]
     clusterPercentagesForGroups[,i] <- x
  }

  #Now the residual matrix is constructed
  residualVector <- as.vector(clusterPercentagesForGroups[,1]-clusterPercentagesForGroups[,2])

  names(residualVector) <- as.numeric(row.names(clusterPercentagesForGroups))

  #Here, the maximum values for the plotting is defined if not added by the user.
  if(missing(maxAbsPlottingValues)){
    maxAbsPlottingValues <- max(abs(clusterPercentagesForGroups))
  }

  #It might be that the maxAbsPlottingValues defined by the user are lower than the most extreme values in the residual data. Then, the data is truncated to fit the range.
  residualVector[residualVector > maxAbsPlottingValues] <- maxAbsPlottingValues
  residualVector[residualVector < -maxAbsPlottingValues] <- -maxAbsPlottingValues

  #Here, a vector with the same length as the cluster vector is generated, but where the cluster info has been substituted with the statistic.
  residualVectorLong <- groupClusterData$cluster
  for(i in 1:length(residualVector)){
    residualVectorLong[groupClusterData$cluster==names(residualVector)[i]] <- residualVector[i]
  }

  #Here the data that will be used for plotting is scaled.
  xYDataScaled <- minMaxScale(xYData)
  colnames(xYDataScaled) <- c("V1", "V2")

  #Make a color vector with the same length as the data
  residual.df <- as.data.frame(residualVectorLong)

  #make a breaks vector to define each bin for the colors
  brks <- with(residual.df, seq(-maxAbsPlottingValues, maxAbsPlottingValues, length.out = 22))

  #assign each value to a bin
  grps <- with(residual.df, cut(residual.df[,1], breaks = brks, include.lowest = TRUE))
  colors <- colorRampPalette(c("#FF0000",  "white","#0000FF"))(21)
  xYDataScaled$col <- rev(colors)[grps]

  #If there is no matrix present to construct the contour lines, create the density matrix from xYData to make them.
  if(missing("densContour")){
    densContour <- densityContours(xYData)
  }

  if(title==TRUE){
  	png(paste(name,'.png', sep=""), width = 2500, height = 2500, units = "px", bg="transparent")
  plot(V2~V1, data=xYDataScaled, main=name, pch=20, cex=dotSize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
  contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")

  dev.off()
  }

  if(title==FALSE){
  	png(paste(name,'.png', sep=""), width = 2500, height = 2500, units = "px", bg="transparent")
  plot(V2~V1, data=xYDataScaled, main="", pch=20, cex=dotSize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
  contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")

  dev.off()
  }

#Create a color legend with text

	yname <- "Residual values"
	topText <- paste(groupName1, " is more abundant", sep="")
	bottomText <- paste(groupName2, " is more abundant", sep="")
	legendTitle <- paste("Color scale for", name, "residual analysis.pdf", sep=" ")

  pdf(legendTitle)
  par(fig=c(0.35,0.65,0,1), xpd=NA)
  z=matrix(1:21,nrow=1)
  x=1
  y=seq(-maxAbsPlottingValues,maxAbsPlottingValues,len=21)
  image(x,y,z,col=rev(colors),axes=FALSE,xlab="",ylab=yname)
  axis(2)
  text(1,maxAbsPlottingValues*1.1, labels=topText, cex=1.1)
  text(1,-maxAbsPlottingValues*1.1, labels=bottomText, cex=1.1)

  box()
  dev.off()

  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }

}
