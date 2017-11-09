#' Sparse partial least squares discriminant analysis with paired and unpaired data
#'
#'
#' This function is used to compare groups of individuals from whom comparable cytometry or other complex data has been generated. It is superior to just running a Wilcoxon analysis in that it does not consider each cluster individually, but instead uses a sparse partial least squares discriminant analysis to first identify which vector thourgh the multidimensional data cloud, created by the cluster-donor matrix, that optimally separates the groups, and as it is a sparse algorithm, applies a penalty to exclude the clusters that are orthogonal, or almost orthogonal to the discriminant vector, i.e. that do not contribute to separating the groups.  
#' @importFrom mixOmics splsda tune.splsda
#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual scale_x_continuous theme element_blank element_rect ggsave
#' @param xYData A dataframe with two columns. Each row contains information about the x and y positition in the field for that observation.
#' @param idsVector Vector with the same length as xYData containing information about the id of each observation.
#' @param groupVector Vector with the same length as xYData containing information about the group identity of each observation.
#' @param clusterVector Vector with the same length as xYData containing information about the cluster identity of each observation.
#' @param pairingVector If this vector is present, a multilevel spls-da will be performed, that considers the within-donor variation between different stimuli. Defaults to NULL.
#' @param displayVector Optionally, if the dataset is very large and the SNE calculation hence becomes impossible to perform for the full dataset, this vector can be included. It should contain the set of rows from the data used for statistics, that has been used to generate the xYData. 
#' @param testSampleRows Optionally, if a train-test setup is wanted, the rows specified in this vector are used to divide the dataset into a training set, used to generate the analysis, and a test set, where the outcome is predicted based on the outcome of the training set. All rows that are not labeled as test rows are assumed to be train rows. 
#' @param densContour An object to create the density contours for the plot. If not present, it will be generated with the xYData. Useful when only a subfraction of a dataset is plotted, and a superimposition of the distribution of the whole dataset is of interest.
#' @param name The main name for the graph and the analysis.
#' @param densContour An object to create the density contours for the plot. Three possible values: 
#' \describe{
#'               \item{densContour}{A densContour object generated previously with dContours}
#'               \item{TRUE}{a densContour object will be generated internally}
#'               \item{FALSE}{No density contours will be displayed.}
#'              }
#' @param groupName1 The name for the first group
#' @param groupName2 The name for the second group
#' @param maxAbsPlottingValues If multiple plots should be compared, it might be useful to define a similar color scale for all plots, so that the same color always means the same statistical value. Such a value can be added here. It defaults to the maximum statistic that is generated in the analysis.
#' @param title If there should be a title displayed on the plotting field. As the plotting field is saved as a png, this title cannot be removed as an object afterwards, as it is saved as coloured pixels. To simplify usage for publication, the default is FALSE, as the files are still named, eventhough no title appears on the plot.
#' @param createDirectory If a directory (i.e. folder) should be created. Defaults to TRUE.
#' @param directoryName The name of the created directory, if it should be created.
#' @param bandColor The color of the contour bands. Defaults to black.
#' @param dotSize Simply the size of the dots. The default makes the dots smaller the more observations that are included.
#' @seealso \code{\link{dColorPlot}}, \code{\link{dDensityPlot}}, \code{\link{dResidualPlot}}
#' @return This function returns the full result of the sPLS-DA. It also returns a sne based plot showing which events that belong to a cluster dominated by the first or the second group defined by the sparse partial least squares loadings of the clusters.
#' @examples
#' #Generate a dataframe with bimodally distributed data and 20 subsamplings.
#' xindividuals <- generateBimodalData(samplings=40, dataCols=7, observations=500)
#'
#' #Now add three columns that will separate the first ten from 
#' #the second ten individuals and merge the datasets
#' xgroups <- generateBimodalData(samplings=2, dataCols=3, observations=10000)
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
#' xClustObject <- dClust(x_scaled, sampleSizes=100, selectionSampleSize=1000, maxIter=20)
#' clusterVector <- xClustObject[[1]]
#'
#' #Run Barnes Hut tSNE on this. 
#' library(Rtsne.multicore)
#' xSNE <- Rtsne.multicore(x_scaled, pca=FALSE)
#'
#' #Run the function. This time without pairing.
#' sPLSDAObject <- dSplsda(xYData=as.data.frame(xSNE$Y), idsVector=x$ids, groupVector=x$group, clusterVector=clusterVector)
#' 
#' #Here, pairing is used
#' #First, an artificial pairing vector, making the first donor amongst the first ten connected to the first donor among the second ten.
#' pairingVector <- c(rep(1:20, each=500), rep(1:20, each=500))
#' 
#' #Then the actual multilevel sPLS-DA is run. 
#' sPLSDAObject <- dSplsda(xYData=as.data.frame(xSNE$Y), idsVector=x$ids, groupVector=x$group, clusterVector=clusterVector, pairingVector=pairingVector, name="d_sPLSDAPlot_paired", groupName1="Stimulation 1", groupName2="Stimulation 2")
#' @export dSplsda
dSplsda <- function(xYData, idsVector, groupVector, clusterVector, pairingVector, displayVector, testSampleRows, densContour=TRUE, name="dSplsda", groupName1=unique(groupVector)[1], groupName2=unique(groupVector)[2], title=FALSE, maxAbsPlottingValues, createDirectory=FALSE, directoryName="dSplsda", bandColor="black", dotSize=400/sqrt(nrow(xYData))){

  if(createDirectory==TRUE){
    dir.create(directoryName)
    workingDirectory <- getwd()
    setwd(paste(workingDirectory, directoryName, sep="/"))
  }

  if(length(unique(groupVector))!=2){
    stop("More or less than two groups are present. This is currently not supported.")
  }

  if(length(unique(idsVector))<8){
    warning("NB! The number of unique ids is smaller than 8, so statistical comparison is not suitable. Use dResidualPlot instead to view differences.")
  }

  if(missing(testSampleRows)==FALSE){
    
    clusterVectorTrain <- clusterVector[-testSampleRows]
    idsVectorTrain <- idsVector[-testSampleRows]
    groupVectorTrain <- groupVector[-testSampleRows]
    clusterVectorTest <- clusterVector[testSampleRows]
    idsVectorTest <- idsVector[testSampleRows]
    groupVectorTest <- groupVector[testSampleRows]
    
    if(missing(pairingVector)==FALSE){
      pairingVectorTrain <- pairingVector[-testSampleRows]
      pairingVectorTest <- pairingVector[testSampleRows]
    }
  } else{
    clusterVectorTrain <- clusterVector
    idsVectorTrain <- idsVector
    groupVectorTrain <- groupVector
    if(missing(pairingVector)==FALSE){
      pairingVectorTrain <- pairingVector
    }
  }
  
  if(missing(pairingVector)){
    dSplsdaInData <- dSplsdaPreCalculations(clusterVectorTrain, idsVectorTrain, groupVectorTrain, groupName1=groupName1, groupName2=groupName2)
  } else {
    dSplsdaInData <- dSplsdaPreCalculations(clusterVectorTrain, idsVectorTrain, groupVectorTrain, groupName1=groupName1, groupName2=groupName2, pairingVector=pairingVectorTrain)
  }
  
  #Here, the number of clusters that should be kept in the sPLS-DA is chosen
  nVarSPLSDA = tune.splsda(X=t(dSplsdaInData[[1]]), Y=dSplsdaInData[[2]], ncomp=1, logratio = "none",
                     test.keepX = dSplsdaInData[[3]], validation = 'loo', dist = "mahalanobis.dist", multilevel = dSplsdaInData[[4]])
  
  #And here the sPLS-DA is performed. Scaling is performed internally in the algorithm. The number of components is set to one, as this situation is by far the most interpretable.
  sPLSDAObject <-  splsda(X=t(dSplsdaInData[[1]]),Y=dSplsdaInData[[2]],ncomp=1, keepX=nVarSPLSDA$choice.keepX, multilevel=dSplsdaInData[[4]])
  
  #Retrieve the x variates for plotting
  sPLSDAX <- data.frame(sPLSDAObject$variates$X)
  densityHist <- cbind(sPLSDAX, dSplsdaInData[[2]])
  
  colnames(densityHist) <- c("sPLSDA_vector","Group")
  
  # Density plots with semi-transparent fill
  ggplot(densityHist, aes(x=sPLSDA_vector, fill=Group)) + geom_density(adjust=0.2, alpha=.4)+scale_fill_manual(values = c("blue", "red")) + scale_x_continuous(limits = c(min(densityHist$sPLSDA_vector)-abs(max(densityHist$sPLSDA_vector)-min(densityHist$sPLSDA_vector))*0.3, max(densityHist$sPLSDA_vector)+abs(max(densityHist$sPLSDA_vector)-min(densityHist$sPLSDA_vector))*0.3)) +
    theme (line = element_blank(),
           panel.background = element_rect(fill = "white"))
  ggsave("Individuals_distributed_along_sPLS-DA_vector.pdf", dpi=300)
  
  #Retrieve the sparse loadings
  sPLSDALoadings <- sPLSDAObject$loadings$X
  
  #Here, a vector with the same length as the cluster vector is generated, but where the cluster info has been substituted with the statistic.
  #If a displayVector has been included, it is used here, to subset the clusterVector
  if(missing(displayVector)==FALSE){
    statisticVector <- clusterVector[displayVector]
    clusterVectorUsed <- clusterVector[displayVector]
  } else {
    statisticVector <- clusterVector
    clusterVectorUsed <- clusterVector
  }
    
  for(i in 1:nrow(sPLSDALoadings)){
    statisticVector[clusterVectorUsed==rownames(sPLSDALoadings)[i]] <- sPLSDALoadings[i]
  }

  #Here, the maximum values for the plotting are defined. If not added by the user, they are obtained from the data.
  if(missing(maxAbsPlottingValues)){
    maxAbsPlottingValues <- max(abs(sPLSDALoadings))
  }

  #Here the data that will be used for plotting is scaled.
  xYDataScaled <- dScale(xYData, scale=c(0,1), robustVarScale=FALSE, center=FALSE)
  colnames(xYDataScaled) <- c("V1", "V2")

  #Make a color vector with the same length as the data
  statistic.df <- as.data.frame(statisticVector)

  #make a breaks vector to define each bin for the colors
  brks <- with(statistic.df, seq(-maxAbsPlottingValues, maxAbsPlottingValues, length.out = 22))

  #assign each value to a bin
  grps <- with(statistic.df, cut(statistic.df[,1], breaks = brks, include.lowest = TRUE))
  colors <- colorRampPalette(c("#FF0000",  "white","#0000FF"))(21)
  xYDataScaled$col <- rev(colors)[grps]

  #If there is no matrix present to construct the contour lines and these are wanted, create the density matrix for xYData to make them.
  if(is.logical(densContour)==TRUE){
    if(densContour==TRUE){
      densContour <- dContours(xYData)
    }
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

	yname <- "sPLS-DA values"
	topText <- paste(groupName1, " is more abundant", sep="")
	bottomText <- paste(groupName2, " is more abundant", sep="")
	legendTitle <- paste("Color scale for", name, "analysis.pdf", sep=" ")

  pdf(legendTitle)
  par(fig=c(0.35,0.65,0,1), xpd=NA)
  z=matrix(1:21,nrow=1)
  x=1
  y=seq(0,maxAbsPlottingValues,len=21)
  image(x,y,z,col=rev(colors),axes=FALSE,xlab="",ylab=yname)
  axis(2)
  text(1,maxAbsPlottingValues*1.1, labels=topText, cex=1.1)
  text(1,maxAbsPlottingValues-maxAbsPlottingValues*1.1, labels=bottomText, cex=1.1)

  box()
  dev.off()

  
  #Return data from the sPLS-DA that was needed for the generation of the graphs
  
  write.csv(sPLSDALoadings, "sPLSDALoadings.csv")
  
  write.csv(sPLSDAX, "sPLSDAVariatesX.csv")
  
  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }

  #Now, prediction is performed, if the setup is train-test.
  if(missing(testSampleRows)==FALSE){
    if(missing(pairingVector)){
      dSplsdaInDataTest <- dSplsdaPreCalculations(clusterVectorTest, idsVectorTest, groupVectorTest, groupName1=groupName1, groupName2=groupName2)
    } else {
      dSplsdaInDataTest <- dSplsdaPreCalculations(clusterVectorTest, idsVectorTest, groupVectorTest, groupName1=groupName1, groupName2=groupName2, pairingVector=pairingVectorTest)
    }
    #And here the sPLS-DA is performed. Scaling is performed internally in the algorithm. The number of components is set to one, as this situation is by far the most interpretable.
    sPLSDAPredictObject <-  predict(object=sPLSDAObject,newdata=t(dSplsdaInDataTest[[1]]), multilevel=dSplsdaInDataTest[[4]])
    
    #Retrieve the x variates for plotting
    sPLSDAX <- data.frame(sPLSDAPredictObject$variates)
    densityHist <- cbind(sPLSDAX, dSplsdaInDataTest[[2]])
    
    colnames(densityHist) <- c("sPLSDA_vector","Group")
    
    # Density plots with semi-transparent fill
    ggplot(densityHist, aes(x=sPLSDA_vector, fill=Group)) + geom_density(adjust=0.2, alpha=.4)+scale_fill_manual(values = c("blue", "red")) + scale_x_continuous(limits = c(min(densityHist$sPLSDA_vector)-abs(max(densityHist$sPLSDA_vector)-min(densityHist$sPLSDA_vector))*0.3, max(densityHist$sPLSDA_vector)+abs(max(densityHist$sPLSDA_vector)-min(densityHist$sPLSDA_vector))*0.3)) +
      theme (line = element_blank(),
             panel.background = element_rect(fill = "white"))
    ggsave("Predicted_individuals_distributed_along_sPLS-DA_vector.pdf", dpi=300)
  }
  if(missing(testSampleRows)==TRUE){
    return(sPLSDAObject)
  } else {
    return(list(sPLSDAObject, sPLSDAPredictObject))
  }

}
