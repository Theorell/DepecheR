#' Wilcoxon rank-sum or signed rank test comparison of subject groups in a dClust result
#'
#'
#' This function is used to compare groups of individuals from whom comparable cytometry or other complex data has been generated.
#' @param xYData A dataframe with two columns. Each row contains information about the x and y positition in the field for that observation.
#' @param idsVector Vector with the same length as xYData containing information about the id of each observation.
#' @param groupVector Vector with the same length as xYData containing information about the group identity of each observation.
#' @param clusterVector Vector with the same length as xYData containing information about the cluster identity of each observation.
#' @param displayVector Optionally, if the dataset is very large and the SNE calculation hence becomes impossible to perform for the full dataset, this vector can be included. It should contain the set of rows from the data used for statistics, that has been used to generate the xYData. 
#' @param paired If the data is paired, so that Wilcoxon signed rank test instead of Wilcoxon rank-sum test/Mann_Whitney test can be used. Defaults to false, i.e. no assumption of pairing is made and Wilcoxon rank sum-test.
#' @param multipleCorrMethod Which method that should be used for adjustment ofmultiple comparisons. Defaults to Benjamini-Hochberg, but all other methods available in \code{\link{p.adjust}} can be used.
#' @param densContour An object to create the density contours for the plot. Three possible values: 
#' \describe{
#'               \item{densContour}{A densContour object generated previously with dContours}
#'               \item{TRUE}{a densContour object will be generated internally}
#'               \item{FALSE}{No density contours will be displayed.}
#'              }
#' @param name The main name for the graph and the analysis.
#' @param groupName1 The name for the first group
#' @param groupName2 The name for the second group
#' @param maxAbsPlottingValues If multiple plots should be compared, it might be useful to define a similar color scale for all plots, so that the same color always means the same statistical value. Such a value can be added here. It defaults to the maximum Wilcoxon statistic that is generated in the analysis.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved as a png, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @seealso \code{\link{dColorPlot}}, \code{\link{dDensityPlot}}, \code{\link{dResidualPlot}}
#' @return This function always returns a dataframe showing the Wilcoxon statistic and the p-value for each cluster, with an included adjustment for multiple comparisons (see above). It also returns a sne based plot showing which events that belong to a cluster dominated by the first or the second group.
#' @examples
#' #Generate a dataframe with bimodally distributed data and 20 subsamplings.
#' xindividuals <- generateBimodalData(samplings=40, ncols=7, observations=500)
#'
#' #Now add three columns that will separate the first ten from 
#' #the second ten individuals and merge the datasets
#' xgroups <- generateBimodalData(samplings=2, ncols=3, observations=10000)
#' colnames(xgroups)[2:4] <- c("X8", "X9", "X10")
#' x <- cbind(xindividuals[,1], xgroups[,1], 
#' xindividuals[,2:ncol(xindividuals)], xgroups[,2:ncol(xgroups)])
#' 
#' colnames(x)[1:2] <- c("ids", "group")
#'
#' #Scale the data
#' x_scaled <- dScale(x[3:ncol(x)])
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#' 
#' #Optimize and run the clustering function.
#' xClustObject <- dClust(x_scaled)
#' clusterVector <- xClustObject[[1]]
#'
#' #Run Barnes Hut tSNE on this. 
#' library(Rtsne.multicore)
#' xSNE <- Rtsne.multicore(x_scaled, pca=FALSE)
#'
#' #Run the function
#' dWilcox(xYData=as.data.frame(xSNE$Y), idsVector=x$ids, groupVector=x$group, clusterVector=clusterVector)
#' @export dWilcox
dWilcox <- function(xYData, idsVector, groupVector, clusterVector, paired=FALSE, multipleCorrMethod="hochberg", densContour=TRUE, name="dWilcox", groupName1=unique(groupVector)[1], groupName2=unique(groupVector)[2], title=FALSE, maxAbsPlottingValues, createDirectory=FALSE, directoryName="dWilcox", bandColor="black", dotSize=400/sqrt(nrow(xYData))){

  if(createDirectory==TRUE){
    dir.create(directoryName)
    workingDirectory <- getwd()
    setwd(paste(workingDirectory, directoryName, sep="/"))

  }

  if(length(unique(groupVector))!=2){
    stop("More or less than two groups are present. Please correct this.")
  }

  if(length(unique(idsVector))<8){
    warning("NB! The number of unique ids is smaller than 8, so statistical comparison is not suitable. Use dResidualPlot instead to view differences.")
  }

  #Here, the statistical evaluation is performed. First, the data is divided into each group.
  clusterVectorGroup1 <- clusterVector[groupVector==unique(groupVector)[1]]
  clusterVectorGroup2 <- clusterVector[groupVector==unique(groupVector)[2]]
  idsVectorGroup1 <- as.character(idsVector[groupVector==unique(groupVector)[1]])
  idsVectorGroup2 <- as.character(idsVector[groupVector==unique(groupVector)[2]])

  #Now, a table with the percentage of cells in each cluster for each individual is created for both groups, in analogy with XXX pKMRun.

  clusterTable1 <- table(clusterVectorGroup1, idsVectorGroup1)
  clusterTable2 <- table(clusterVectorGroup2, idsVectorGroup2)

  countTable1 <- table(idsVectorGroup1)
  countTable2 <- table(idsVectorGroup2)

  clusterFractionsForAllIds1 <- clusterTable1
  clusterFractionsForAllIds2 <- clusterTable2

  for(i in 1:length(countTable1)){
     x <- clusterTable1[,i]/countTable1[i]
    clusterFractionsForAllIds1[,i] <- x
  }

  for(i in 1:length(countTable2)){
    x <- clusterTable2[,i]/countTable2[i]
    clusterFractionsForAllIds2[,i] <- x
  }

  #And here the statistical test is performed for each cluster individually
  statisticList <- mapply(wilcox.test, as.data.frame.matrix(t(clusterFractionsForAllIds1)), as.data.frame.matrix(t(clusterFractionsForAllIds2)), MoreArgs=list(alternative="two.sided", paired=paired, exact=FALSE), SIMPLIFY=FALSE)

  #Now, the statistics and the p-values are retrieved
  statistic <- unlist(lapply(statisticList, `[[`, 1))
  p_values <- unlist(lapply(statisticList, `[[`, 3))

  #Here, adjustments for multiple comparisons are performed
  p_adjusted <- p.adjust(p_values, method=multipleCorrMethod)

  #Now the median for each group and cluster is calculated
  median1 <- 100*apply(clusterFractionsForAllIds1, 1, median)
  median2 <- 100*apply(clusterFractionsForAllIds2, 1, median)
  
  #Combine the four
  result <- data.frame(as.numeric(names(p_values)), median1, median2, statistic, p_values, p_adjusted)
  row.names(result) <- c(1:nrow(result))
  colnames(result) <- c("Cluster", paste("Median percentage for", groupName1, sep=" "), paste("Median percentage for", groupName2, sep=" "), "Wilcoxon_statistic", "p-value", paste(multipleCorrMethod, "corrected p-value", sep=" "))

  #Here, a vector with the same length as the cluster vector is generated, but where the cluster info has been substituted with the statistic.
  #If a displayVector has been included, it is used here, to subset the clusterVector
  if(is.null(displayVector)==FALSE){
    statisticVector <- clusterVector[displayVector]
    clusterVectorUsed <- clusterVector[displayVector]
  } else {
    statisticVector <- clusterVector
    clusterVectorUsed <- clusterVector
  }
  for(i in 1:nrow(result)){
    statisticVector[clusterVectorUsed==result$Cluster[i]] <- result$Wilcoxon_statistic[i]
  }

  #Here the data that will be used for plotting are scaled.
  xYDataScaled <- dScale(xYData, scale=c(0,1), robustVarScale=FALSE, center=FALSE)
  colnames(xYDataScaled) <- c("V1", "V2")

  #Make a color vector with the same length as the data
  statistic.df <- as.data.frame(statisticVector)

  #Make a breaks vector to define each bin for the colors. 
  #To do this, the Wilcoxon statistic that corresponds to a p-value of 1 needs to be obtained. This is done by running a perfectly overlapping Wilcoxon.
  
  nullWilcox <- wilcox.test(rep(c(1,2), length=ncol(clusterFractionsForAllIds1)), rep(c(2,1), length=ncol(clusterFractionsForAllIds2)), paired=paired, exact=FALSE)

  centerWilcoxStatistic <- nullWilcox$statistic
  
  #Here, the maximum values for the plotting are defined. If not added by the user, they are obtained from the data.
  if(missing(maxAbsPlottingValues)){
    maxValue <- max(statistic)-centerWilcoxStatistic
    minValue <- centerWilcoxStatistic-min(statistic)
    maxAbsPlottingValues <- max(c(maxValue, minValue))
  }
    
  #And here, the lower border 
  
  brks <- with(statistic.df, seq(centerWilcoxStatistic-maxAbsPlottingValues, centerWilcoxStatistic+maxAbsPlottingValues, length.out = 22))

  #assign each value to a bin
  grps <- with(statistic.df, cut(statistic.df[,1], breaks = brks, include.lowest = TRUE))
  colors <- colorRampPalette(c("#FF0000",  "white","#0000FF"))(21)
  xYDataScaled$col <- rev(colors)[grps]

  #If there is no matrix present to construct the contour lines and these are wanted, create the density matrix for xYData to make them.
  if(densContour==TRUE){
    densContour <- dContours(xYData)
  }

  if(title==TRUE){
  	png(paste(name,'.png', sep=""), width = 2500, height = 2500, units = "px", bg="transparent")
  plot(V2~V1, data=xYDataScaled, main=name, pch=20, cex=dotSize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }

  if(title==FALSE){
  	png(paste(name,'.png', sep=""), width = 2500, height = 2500, units = "px", bg="transparent")
  plot(V2~V1, data=xYDataScaled, main="", pch=20, cex=dotSize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
  }
  if(length(densContour)>1){
    par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
    contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")
  } 
  dev.off()
#Create a color legend with text

	yname <- "Wilcoxon values"
	topText <- paste(groupName1, " is more abundant", sep="")
	bottomText <- paste(groupName2, " is more abundant", sep="")
	legendTitle <- paste("Color scale for", name, "analysis.pdf", sep=" ")

  pdf(legendTitle)
  par(fig=c(0.35,0.65,0,1), xpd=NA)
  z=matrix(1:21,nrow=1)
  x=1
  y=seq(centerWilcoxStatistic-maxAbsPlottingValues, centerWilcoxStatistic+maxAbsPlottingValues,len=21)
  image(x,y,z,col=rev(colors),axes=FALSE,xlab="",ylab=yname)
  axis(2)
  text(1,centerWilcoxStatistic+maxAbsPlottingValues*1.1, labels=topText, cex=1.1)
  text(1,centerWilcoxStatistic-maxAbsPlottingValues*1.1, labels=bottomText, cex=1.1)

  box()
  dev.off()

  write.csv(result, "dWilcoxResult.csv", row.names=FALSE)

  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }

  return(result)

}
