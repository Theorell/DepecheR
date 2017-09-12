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
#' @param densContour An object to create the density contours for the plot. If not present, it will be generated with the xYData. Useful when only a subfraction of a dataset is plotted, and a superimposition of the distribution of the whole dataset is of interest.
#' @param name The main name for the graph and the analysis.
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
#' xindividuals <- generateFlowCytometryData(samplings=40, ncols=7, observations=500)
#'
#' #Now add three columns that will separate the first ten from 
#' #the second ten individuals and merge the datasets
#' xgroups <- generateFlowCytometryData(samplings=2, ncols=3, observations=10000)
#' colnames(xgroups)[2:4] <- c("X8", "X9", "X10")
#' x <- cbind(xindividuals[,1], xgroups[,1], 
#' xindividuals[,2:ncol(xindividuals)], xgroups[,2:ncol(xgroups)])
#' 
#' colnames(x)[1:2] <- c("ids", "group")
#'
#' #Scale the data
#' x_scaled <- quantileScale(x[3:ncol(x)])
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#' 
#' #Create the optimized number of clusters for this dataset
#' x_optim <- dClustOpt(x_scaled, iterations=50, bootstrapObservations=1000)
#' x_pKM <- dClust(x_scaled, penaltyOffset=x_optim[[1]][["bestPenaltyOffset"]], 
#' withOrigoClust=x_optim[[1]][["withOrigoClust"]], iterations=1, ids=x[,1])
#'
#' #Run Barnes Hut tSNE on this. 
#' library(Rtsne.multicore)
#' xSNE <- Rtsne.multicore(x_scaled, pca=FALSE)
#'
#' #Run the function. This time without pairing.
#' sPLSDAObject <- dsPLSDAPlot(xYData=as.data.frame(xSNE$Y), idsVector=x$ids, groupVector=x$group, clusterVector=x_pKM$clusterVector)
#' 
#' #Here, pairing is used
#' #First, an artificial pairing vector, making the first donor amongst the first ten connected to the first donor among the second ten.
#' pairingVector <- c(rep(1:20, each=500), rep(1:20, each=500))
#' 
#' #Then the actual multilevel sPLS-DA is run. 
#' sPLSDAObject <- dsPLSDAPlot(xYData=as.data.frame(xSNE$Y), idsVector=x$ids, groupVector=x$group, clusterVector=x_pKM$clusterVector, pairingVector=pairingVector, name="d_sPLSDAPlot_paired", groupName1="Stimulation 1", groupName2="Stimulation 2")
#' @export dsPLSDAPlot
dsPLSDAPlot <- function(xYData, idsVector, groupVector, clusterVector, pairingVector=NULL, densContour, name="dsPLSDAPlot", groupName1=unique(groupVector)[1], groupName2=unique(groupVector)[2], title=FALSE, maxAbsPlottingValues, createDirectory=FALSE, directoryName="dsPLSDAPlot", bandColor="black", dotSize=400/sqrt(nrow(xYData))){

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

  #In the very unlikely event that there is not a single observation for one cluster from one of the groups, this cluster is substituted to that table with a row of zeros.
  if(length(clusterTable1)<length(unique(clusterVector))){
    clusterTable1big <- clusterTable1
    zeroVector <- rep(0, times=ncol(clusterTable1))
    
    #Here, rows are added to the cluster table to make the number of rows the same as the unique values of the cluster vector.
    for(i in length(setdiff(unique(clusterVector), as.numeric(row.names(clusterTable1))))){
      clusterTable1big <- rbind(clusterTable1big, zeroVector)
      #Row names are added to the new rows, that correspond to the clusters that were missing in the table
      row.names(clusterTable1big)[nrow(clusterTable1big)] <- setdiff(unique(clusterVector), as.numeric(row.names(clusterTable1)))[i]
      
    }
    
    #The rows of the table are re-sorted
    clusterTable1bigResorted <- clusterTable1big[order(as.numeric(row.names(clusterTable1big))),]
    
    clusterTable1 <- clusterTable1bigResorted  
    }
  #And the same procedure is done for the second group 
  if(length(clusterTable2)<length(unique(clusterVector))){
    clusterTable2big <- clusterTable2
    zeroVector <- rep(0, times=ncol(clusterTable2))
    
    #Here, rows are added to the cluster table to make the number of rows the same as the unique values of the cluster vector.
    for(i in length(setdiff(unique(clusterVector), as.numeric(row.names(clusterTable2))))){
      clusterTable2big <- rbind(clusterTable2big, zeroVector)
      #Row names are added to the new rows, that correspond to the clusters that were missing in the table
      row.names(clusterTable2big)[nrow(clusterTable2big)] <- setdiff(unique(clusterVector), as.numeric(row.names(clusterTable2)))[i]
      
    }
  
    #The rows of the table are re-sorted
    clusterTable2bigResorted <- clusterTable2big[order(as.numeric(row.names(clusterTable2big))),]
    
    clusterTable2 <- clusterTable2bigResorted
  }    

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
  
  if(is.null(pairingVector)==FALSE){
    #Here, a comparable pairing vector pair is created if a multilevel sPLS-DA should be performed.
  
    pairingVectorGroup1 <- as.character(pairingVector[groupVector==unique(groupVector)[1]])
    pairingVectorGroup2 <- as.character(pairingVector[groupVector==unique(groupVector)[2]])
    
    pairingShortGroup1 <- clusterFractionsForAllIds1[1,]
    
    for(i in 1:ncol(clusterFractionsForAllIds1)){
      pairingShortGroup1[i] <- pairingVectorGroup1[which(as.numeric(colnames(clusterFractionsForAllIds1)[i])==idsVectorGroup1)[1]]
    }
    
    pairingShortGroup2 <- clusterFractionsForAllIds2[1,]
    
    for(i in 1:ncol(clusterFractionsForAllIds2)){
      pairingShortGroup2[i] <- pairingVectorGroup2[which(as.numeric(colnames(clusterFractionsForAllIds2)[i])==idsVectorGroup2)[1]]
    }
    
    pairingAll <- c(pairingShortGroup1, pairingShortGroup2)
  } else{
    pairingAll <- NULL
  }
  
  #A group vector is created with the same length as the number of columns in the tables
  groupId <- as.factor(c(rep(groupName1, ncol(clusterFractionsForAllIds1)), rep(groupName2, ncol(clusterFractionsForAllIds2))))
  
  #These two tables are combined to one
  clusterFractionsForAllIds <- as.matrix(cbind(clusterFractionsForAllIds1, clusterFractionsForAllIds2))
  
  #Here, the number of possible clusters to be saved in the sPLS-DA is chosen. The number of tested clusters can be up to five, but if the number of clusters is low, the lowest tested variant will be 1 and that will be tested only once.
  testKeepAlternatives <- vector()
  for(i in 1:5){
    if(i==1){
      testKeepAlternatives[i] <- nrow(clusterFractionsForAllIds)
    } else {
      if(testKeepAlternatives[i-1]>1){
        testKeepAlternatives[i] <- round(testKeepAlternatives[i-1]/2)
      } else {break}
    }
  }

  
  #Here, the number of clusters that should be kept in the sPLS-DA is chosen
  nVarSPLSDA = tune.splsda(X=t(clusterFractionsForAllIds), Y=groupId, ncomp=1, logratio = "none",
                     test.keepX = testKeepAlternatives, validation = 'loo', dist = "mahalanobis.dist", multilevel = pairingAll)
  
  #And here the sPLS-DA is performed. Scaling is performed internally in the algorithm. The number of components is set to one, as this situation is by far the most interpretable.
  sPLSDAObject <-  splsda(X=t(clusterFractionsForAllIds),Y=groupId,ncomp=1, keepX=nVarSPLSDA$choice.keepX, multilevel=pairingAll)
  
  #Retrieve the x variates for plotting
  sPLSDAX <- data.frame(sPLSDAObject$variates$X)
  densityHist <- cbind(sPLSDAX, groupId)
  
  colnames(densityHist) <- c("sPLSDA_vector","Group")
  
  # Density plots with semi-transparent fill
  ggplot(densityHist, aes(x=sPLSDA_vector, fill=Group)) + geom_density(alpha=.4)+scale_fill_manual(values = c("blue", "red")) + scale_x_continuous(limits = c(min(densityHist$sPLSDA_vector)-abs(max(densityHist$sPLSDA_vector)-min(densityHist$sPLSDA_vector))*0.3, max(densityHist$sPLSDA_vector)+abs(max(densityHist$sPLSDA_vector)-min(densityHist$sPLSDA_vector))*0.3)) +
    theme (line = element_blank(),
           panel.background = element_rect(fill = "white"))
  ggsave("Individuals_distributed_along_sPLS-DA_vector.pdf", dpi=300)
  
  #Retrieve the sparse loadings
  sPLSDALoadings <- sPLSDAObject$loadings$X
  
  #if(nVarSPLSDA$choice.keepX==1){
    #ANSWER <- readline("With the optimal sparsity, only one cluster is selected. Do you want to plot non-sparse PLS-DA info instead? The sparse PLS-DA information will still be retrieved in the result of the function. Answer with yes or no.")
    
    #if (substr(ANSWER, 1, 1) == "n"){
      #print("Ok, then we will plot the only non-sparse cluster")
    #} else if(substr(ANSWER, 1, 1) == "y"){
     # cat("Ok, let us plot all clusters as if they were in a PLS-DA")
      #sPLSDALoadings <- sPLSDAObject$mat.c
    #}
 # }

  #Here, a vector with the same length as the cluster vector is generated, but where the cluster info has been substituted with the statistic.
  statisticVector <- clusterVector
  for(i in 1:nrow(sPLSDALoadings)){
    statisticVector[clusterVector==rownames(sPLSDALoadings)[i]] <- sPLSDALoadings[i]
  }

  #Here, the maximum values for the plotting are defined. If not added by the user, they are obtained from the data.
  if(missing(maxAbsPlottingValues)){
    maxAbsPlottingValues <- max(abs(sPLSDALoadings))
  }

  #Here the data that will be used for plotting is scaled.
  xYDataScaled <- quantileScale(xYData, robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE)
  colnames(xYDataScaled) <- c("V1", "V2")

  #Make a color vector with the same length as the data
  statistic.df <- as.data.frame(statisticVector)

  #make a breaks vector to define each bin for the colors
  brks <- with(statistic.df, seq(-maxAbsPlottingValues, maxAbsPlottingValues, length.out = 22))

  #assign each value to a bin
  grps <- with(statistic.df, cut(statistic.df[,1], breaks = brks, include.lowest = TRUE))
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

  return(sPLSDAObject)

}
