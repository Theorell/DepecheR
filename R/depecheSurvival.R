#' Wilcoxon rank-sum or signed rank test comparison of subject groups in pKMRun result
#'
#'
#' This function is used to compare groups of individuals from whom comparable cytometry or other complex data has been generated.
#' @importFrom survival coxph
#' @param idGroupClusterData This dataframe should contain three columns with information about the cluster identity (named "cluster"), the id (named "ids")and the group (named "group")identity for each observation.
#' @param surv
#' @param xYData A dataframe with two columns. Each row contains information about the x and y positition in the field for that observation. It needs to have the same number of rows as idGroupClusterData.
#' @param densContour An object to create the density contours for the plot. If not present, it will be generated with the xYData. Useful when only a subfraction of a dataset is plotted, and a superimposition of the distribution of the whole dataset is of interest.
#' @param name The main name for the graph and the analysis.
#' @param maxAbsPlottingValues If multiple plots should be compared, it might be useful to define a similar color scale for all plots, so that the same color always means the same statistical value. Such a value can be added here. It defaults to the maximum Wilcoxon statistic that is generated in the analysis.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved as a png, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @return This function always returns a dataframe showing the Wilcoxon statistic and the p-value for each cluster, with an included adjustment for multiple comparisons (see above). It also returns a sne based plot showing which events that belong to a cluster dominated by the first or the second group.
#' @examples
#' #Generate a dataframe with bimodally distributed data and 20 subsamplings.
#' x <- generateFlowCytometryData(samplings=100, ncols=7, observations=200)
#'
#' #Split it into the samplings and order these according to the value in one of the columns
#' xSplitted <- split(x, x[,1])
#'
#' names(xSplitted) <- seq(1:length(xSplitted))
#'
#' orderingList <- list()
#' for(i in 2:4){
#'   y <- sapply(xSplitted, `[[`, i)
#'   yMean <- apply(y, 2, mean)
#'   yMeanRound <- sapply(yMean, round)
#'   orderingList[[i]] <- yMeanRound
#' }
#' orderingDf <- as.data.frame(t(do.call(rbind, orderingList)))
#' orderingDf$splitNumber <- seq(1:nrow(orderingDf))
#' orderingDfOrdered <- orderingDf[order(orderingDf[,1]),]
#'
#' #Simulate some survival data
#' library(survsim)
#' surv.data <- simple.surv.sim(n=100, foltime=3600, dist.ev=c('llogistic'), 
#' anc.ev=c(0.69978200185280),beta0.ev=c(5.84298525742252),,anc.cens=1.17783687569519, 
#' beta0.cens=7.39773677281100,z=list(c("unif", 0.8, 1.2)))
#'
#' #Combine the relevant parts from these survival objects
#' survival <- data.frame(surv.data$stop, surv.data$status)
#' colnames(survival) <- c("survivalTime", "status")
#'
#' #Now sort the survival data in the same order as the means of 
#' #the first column in the synthetic flow data.
#' survivalOrdered <- survival[order(survival$survivalTime),]
#' survivalOrdered$Identity <- orderingDfOrdered$splitNumber
#' survivalMatched <- survivalOrdered[order(survivalOrdered$Identity),]
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
#' #Create the idClusterData object
#' idClusterData <- data.frame(x[,1], x_pKM$clusterVector)
#' colnames(idClusterData) <- c("id", "cluster")
#'
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #And finally run the function
#' depecheSurvival(idGroupClusterData=idGroupClusterData, survival=survival, 
#' xYData=as.data.frame(xSNE$Y))
#' @export depecheSurvival
depecheSurvival <- function(idClusterData, survival, xYData, densContour, name=depecheSurvival, title=FALSE, dotsize=400/sqrt(length(data)), highestAbsValues=max(abs(data)), bandColor="black", createDirectory=TRUE, directoryName="depeche Survival"){

  if(createDirectory==TRUE){
    dir.create(directoryName)
    workingDirectory <- getwd()
    setwd(paste(workingDirectory, directoryName, sep="/"))

  }

  #Here, the residuals are identified.
  #A table with the percentage of cells in each cluster for each group is created in analogy with XXX pKMRun.

  clusterTable <- table(idClusterData$cluster, idClusterData$id)

  countTable <- table(idClusterData$id)

  clusterPercentagesForIds <- clusterTable

  for(i in 1:length(countTable)){
    x <- 100*clusterTable[,i]/countTable[i]
    clusterPercentagesForIds[,i] <- x
  }

  clusterPercentagesForIdsDfScaled <- as.data.frame.matrix(scale(t(clusterPercentagesForIds)))

  #Now identify the cluster that contributes the least by itself to separating the groups and turn all its values to 0. This is done not to lose information on the last cluster, as the analysis cannot be performed with completely matching data, which is the case when each row becomes a 100 as in this case with the percentages of clusters.

  coxPhList <- list()
  for(i in 1:ncol(clusterPercentagesForIdsDfScaled)){
    coxPhList[[i]] <-  coxph(Surv(survival$survivalTime, survival$status) ~ clusterPercentagesForIdsDfScaled[,i], data=clusterPercentagesForIdsDfScaled)

  }
  #Retrieve the Wald test statistic results
  waldTestsSingle <- unlist(lapply(coxPhList, `[[`, 15))

  clusterPercentagesForIdsWorstExcluded <-clusterPercentagesForIdsDfScaled[,-which(waldTestsSingle==min(waldTestsSingle))]
  clusterPercentagesForIdsWorstExcluded[,which(waldTestsSingle==min(waldTestsSingle))] <- 0

  #A survival analysis is performed with the reduced dataset
  coxphData <- coxph(Surv(survival$survivalTime, survival$status) ~., data=clusterPercentagesForIdsWorstExcluded)

  #Extract the coefficients
  coxphCoefficients <- as.data.frame(summary(coxphData)$coefficients)

  coxPhWaldValues <- coxphCoefficients$z

  #Add a row of zeros for the cluster that did not contribute to the clustering
  coxPhWaldValuesFull <- c(coxPhWaldValues, 0)
  names(coxPhWaldValuesFull) <- c(numeric(row.names(coxphCoefficients)), as.numeric(colnames(clusterPercentagesForIdsDfScaled)[which(waldTestsSingle==min(waldTestsSingle))]))

  #Here the data that will be used for plotting are scaled.
  xYDataScaled <- minMaxScale(xYData)
  colnames(xYDataScaled) <- c("V1", "V2")

  #Make a color vector with the same length as the data

  #Make a "breaks" vector to define each bin for the colors
  brks <- with(wald.df, seq(-highestAbsValues, highestAbsValues, length.out = 22))

  #assign each value to a bin
  grps <- with(wald.df, cut(wald.df[,1], breaks = brks, include.lowest = TRUE))

  colors <- colorRampPalette(c("#FF0000",  "white","#0000FF"))(21)

  xYDataScaled$col <- rev(colors)[grps]

  #If there is no matrix present to construct the contour lines, create the density matrix from xYData to make them.
  if(missing("densContour")){
    densContour <- densityContours(xYData)
  }


if(title==TRUE){
	jpeg(paste(name,'.jpeg', sep=""), width = 2500, height = 2500, units = "px")
plot(V2~V1, data=xYDataScaled, main=name, pch=20, cex=dotsize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")

dev.off()
}

if(title==FALSE){
	jpeg(paste(name,'.jpeg', sep=""), width = 2500, height = 2500, units = "px")
plot(V2~V1, data=xYDataScaled, main="", pch=20, cex=dotsize, cex.main=5, col=col, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), axes=FALSE, xaxs="i", yaxs="i")
par(fig=c(0,1,0,1), mar=c(6,4.5,4.5,2.5), new=TRUE)
contour(x=densContour$x, y=densContour$y, z=densContour$z, xlim=c(-0.05, 1.05), ylim=c(-0.05, 1.05), nlevels=10, col=bandColor, lwd=8, drawlabels = FALSE, axes=FALSE, xaxs="i", yaxs="i")

dev.off()
}

#Create a color legend with text
	yname <- "Wald values"
	topText <- "Higher risk"
	bottomText <- "Lower risk"
	legendTitle <- paste("Color scale for", name, ".pdf", sep="")

pdf(legendTitle)
par(fig=c(0.35,0.65,0,1), xpd=NA)
z=matrix(1:21,nrow=1)
x=1
y=seq(-highestAbsValues,highestAbsValues,len=21)
image(x,y,z,col=rev(colors),axes=FALSE,xlab="",ylab=yname)
axis(2)
text(1,highestAbsValues*1.13, labels=topText, cex=1.1)
text(1,highestAbsValues*-1.13, labels=bottomText, cex=1.1)
box()
dev.off()

}


