#' Create violin plots for all non-penalized variable for all clusters
#'
#' Here, violin plots of a specific cluster and the total population are created for each variable that has not been penalized away in the penalized K-means analysis. As al such plots are generated for each cluster, this function creates a great number of plots in most instances.
#' @importFrom ggplot2 ggplot aes geom_violin scale_color_manual scale_fill_manual theme_classic labs ggsave
#' @importFrom viridis inferno
#' @param clusterCenters A matrix containing information about where the centers are in all the variables that contributed to creating the cluster with the given penalty term.
#' @param clusterVector A vector with information about the cluster identity of all observations. Needs to have the same length as the number of rows in the inDataFrame.
#' @param order The order that the unique features of the cluster vector should appear in. For harmonization with colorVector and all subsequent functions.
#' @param inDataFrame A dataframe that has been used to generate the cluster vector and the clusterCenters. Note that the scaling does not matter in this case, as each variable wil be plotted separately.
#' @param plotAll If all parameters, including the non-contributing, should be plotted for each cluster. Defaults to FALSE.
#' @return One graph is created for each non-penalized variable in each non-penalized cluster, which often means that the function creates a vast number of graphs. The graphs are sorted into subfolders for each cluster.
#' @seealso \code{\link{dDensityPlot}}, \code{\link{dColorPlot}}, \code{\link{colorVector}}
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateBimodalData(samplings=2, ncols=8)
#'
#' #Scale this dataframe
#' x_scaled <- dScale(x[2:ncol(x)])
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Optimize and run the clustering function.
#' xOptAndClustObject <- dOptAndClust(x_scaled,ids=x[,1])
#' xClustObject <- xOptAndClustObject[[2]]
#'
#' #Create the plots of the variables that contribute to creating each cluster
#' dViolins(clusterCenters=xClustObject$clusterCenters, clusterVector=as.numeric(xClustObject$clusterVector), inDataFrame=x[,2:ncol(x)])
#' 
#' #Now, finally, create plots of all clusters, regardless of if they contributed or not
#' dViolins(clusterCenters=xClustObject$clusterCenters, clusterVector=as.numeric(xClustObject$clusterVector), inDataFrame=x[,2:ncol(x)], plotAll=TRUE)
#' @export dViolins
dViolins <- function(clusterCenters, clusterVector, order=unique(clusterVector), inDataFrame, plotAll=FALSE){

  percentClusterVector <- dScale(clusterVector, robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE, multiplicationFactor=100)

  paletteColors <- inferno(length(order))

  #Here, a directory for all the subdirectories for each cluster is made
  directoryName <- "Cluster expressions"
  dir.create(directoryName)
  workingDirectory <- getwd()
  setwd(paste(workingDirectory, directoryName, sep="/"))

  #Here, the columns in the inDataFrame that are not selected as contributing are excluded from further analysis, if plotAll is not TRUE
  if(plotAll==FALSE){
    inDataFocused <- subset(inDataFrame, select=colnames(clusterCenters)) 
  } else {
    inDataFocused <- inDataFrame
  }
  
  
  for(i in 1:length(order)){

    #Here, a specific directory for the graphics are made.
    directoryName <- paste("Cluster", order[i])
    dir.create(directoryName)
    workingDirectoryClusters <- getwd()
    setwd(paste(workingDirectoryClusters, directoryName, sep="/"))


    #This code is an efficient way of giving all rows in the "Clusters" column the same name, except for the rows with the cluster of interest.

    clustIndicesSpecific <- sapply(clusterVector, singleEventClusterNaming, n=order[i])

    #Create a color vector for the visualzation
    clustColorsSpecific <- c(paletteColors[i], "#d3d3d3")

    #Here, the mu variables for the specific cluster is extracted, if not all clusters should be shown. 
    if(plotAll==FALSE){
      oneClustAllMu <- clusterCenters[rownames(clusterCenters)==order[i],]
    } else {
      oneClustAllMu <- rep(1, ncol(inDataFrame))
    }
    
    #Here the variable names is exported
    allVarNames <- colnames(inDataFocused)
    #Then a list is created that contains the objects for the plot creation
    oneClustAllVarList <- mapply(createAllClustOneVarMu, inDataFocused, oneClustAllMu, allVarNames, MoreArgs=list(clust=clustIndicesSpecific, cols=clustColorsSpecific, clustNum=order[i]), SIMPLIFY=FALSE)
    #And then the plots are created
    sapply(oneClustAllVarList, createOneViolin, plotAll=plotAll)

    setwd(workingDirectoryClusters)

  }

  setwd(workingDirectory)

}

#This function is the core of this whole thing, creating the plot
createOneViolin <- function(allClustOneVarOneMuList, plotAll){
  if(plotAll==FALSE){
    if(allClustOneVarOneMuList[[2]]!=0){
      
      plotname = paste("Cluster", "_", allClustOneVarOneMuList[[4]], "_", allClustOneVarOneMuList[[3]], ".pdf", sep="")
      dp <- ggplot(allClustOneVarOneMuList[[1]], aes(x=Cluster, y=var, fill=Cluster)) +
        geom_violin(trim=FALSE)+
        scale_color_manual(values=allClustOneVarOneMuList[[5]])+
        scale_fill_manual(values=allClustOneVarOneMuList[[5]])+
        labs(title=paste("Plot of", allClustOneVarOneMuList[[3]], "for cluster", allClustOneVarOneMuList[[4]]), x="Cluster", y = "Intensity")
      dp + theme_classic()
      ggsave(filename = plotname, dpi=300)
    }
  } else {

      plotname = paste("Cluster", "_", allClustOneVarOneMuList[[4]], "_", allClustOneVarOneMuList[[3]], ".pdf", sep="")
      dp <- ggplot(allClustOneVarOneMuList[[1]], aes(x=Cluster, y=var, fill=Cluster)) +
        geom_violin(trim=FALSE)+
        scale_color_manual(values=allClustOneVarOneMuList[[5]])+
        scale_fill_manual(values=allClustOneVarOneMuList[[5]])+
        labs(title=paste("Plot of", allClustOneVarOneMuList[[3]], "for cluster", allClustOneVarOneMuList[[4]]), x="Cluster", y = "Intensity")
      dp + theme_classic()
      ggsave(filename = plotname, dpi=300)
  }

}

#Function needed to create the right kind of object for the createOneViolinOneClust function
createAllClustOneVarMu <- function(var, oneClustOneMu, varName, clust, cols, clustNum){

	allClustOneVar <- data.frame(clust, var)

	colnames(allClustOneVar) <- c("Cluster", "var")

	allClustOneVarOneMuList <- list(allClustOneVar, oneClustOneMu, varName, clustNum, cols)
	return(allClustOneVarOneMuList)
}


#Function to give the right cluster annotation to each event for each cluster investigation
singleEventClusterNaming <- function(event, n){
	ifelse(event==n, return(n), return("All clusters"))
}

