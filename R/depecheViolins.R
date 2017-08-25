#' Create violin plots for all non-penalized variable for all clusters
#'
#' Here, violin plots of a specific cluster and the total population are created for each variable that has not been penalized away in the penalized K-means analysis. As al such plots are generated for each cluster, this function creates a great number of plots in most instances.
#' @importFrom ggplot2 ggplot aes geom_violin scale_color_manual scale_fill_manual theme_classic labs ggsave
#' @importFrom viridis inferno
#' @param clusterCenters A matrix containing information about where the centers are in all the variables that contributed to creating the cluster with the given penalty term.
#' @param clusterVector A vector with information about the cluster identity of all observations. Needs to have the same length as the number of rows in the inDataFrame.
#' @param order The order that the unique features of the cluster vector should appear in. For harmonization with colorVector and all subsequent functions.
#' @param inDataFrame A dataframe that has been used to generate the cluster vector and the clusterCenters. Note that the scaling does not matter in this case, as each variable wil be plotted separately.
#' @return One graph is created for each non-penalized variable in each non-penalized cluster, which often means that the function creates a vast number of graphs. The graphs are sorted into subfolders for each cluster.
#' @seealso \code{\link{depecheDensity}}, \code{\link{depecheColor}}, \code{\link{colorVector}}
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData(samplings=2, ncols=8)
#'
#' #Scale this dataframe
#' x_scaled <- quantileScale(x[2:ncol(x)])
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the Optim function to get good starting points
#' x_optim <- pKMOptim(x_scaled, iterations=5, bootstrapObservations=1000)
#'
#' #Then run the clustering function
#' x_clustered <- pKMRun(x_scaled, regVec=x_optim[[1]][1,1], 
#' withOrWithoutZeroClust=x_optim[[1]][1,2], iterations=2, ids=x[,1])
#'
#' #And finally create all the clusters
#' depecheViolins(clusterCenters=x_clustered$clusterCenters, clusterVector=as.numeric(x_clustered$clusterVector), inDataFrame=x[,2:ncol(x)])
#' @export depecheViolins
depecheViolins <- function(clusterCenters, clusterVector, order=unique(clusterVector), inDataFrame){

  percentClusterVector <- quantileScale(clusterVector, robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE, multiplicationFactor=100)

  paletteColors <- inferno(length(order))

  #Here, a directory for all the subdirectories for each cluster is made
  directoryName <- "Cluster expressions"
  dir.create(directoryName)
  workingDirectory <- getwd()
  setwd(paste(workingDirectory, directoryName, sep="/"))

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

    #Here, the mu variables for the specific cluster is extracted
    oneClustAllMu <- clusterCenters[rownames(clusterCenters)==order[i],]

    #Here the variable names is exported
    allVarNames <- colnames(inDataFrame)
    #Then a list is created that contin the objects for the prot creation
    oneClustAllVarList <- mapply(createAllClustOneVarMu, inDataFrame, oneClustAllMu, allVarNames, MoreArgs=list(clust=clustIndicesSpecific, cols=clustColorsSpecific, clustNum=order[i]), SIMPLIFY=FALSE)
    #And then the plots are created
    sapply(oneClustAllVarList, createOneViolin)

    setwd(workingDirectoryClusters)

  }

  setwd(workingDirectory)

}

#This function is the core of this whole thing, creating the plot
createOneViolin <- function(allClustOneVarOneMuList){
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

