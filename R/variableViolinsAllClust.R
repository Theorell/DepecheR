#' Create violin plots for all non-penalized variable for all clusters
#'
#' Here, violin plots of a specific cluster and the total population are created for each variable that has not been penalized away in the penalized K-means analysis. As al such plots are generated for each cluster, this function creates a great number of plots in most instances.
#' @importFrom ggplot2 ggplot aes geom_violin scale_color_manual scale_fill_manual theme_classic labs ggsave
#' @importFrom gplots rich.colors
#' @param sparsityMatrix A matrix containing information about
#' @param clusterVector A vector with information about the cluster identity of all observations. Needs to have the same length as the number of rows in the inDataFrame.
#' @param inDataFrame A dataframe that has been used to generate the cluster vector and the sparsityMatrix. Note that the scaling does not matter in this case, as each variable wil be plotted separately.
#' @return One graph is created for each non-penalized variable in each non-penalized cluster, which often means that the function creates a vast number of graphs. The graphs are sorted into subfolders for each cluster.
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData(samplings=2, ncols=8)
#'
#' #Scale this datamframe
#' x_scaled <- quantileScale(x[2:ncol(x)])
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the Optim function to get good starting points
#' x_optim <- Optim(x_scaled, iterations=5, bootstrapObservations=1000)
#'
#' #Then run the clustering function
#' x_ <- Run(x_scaled, regVec=x_optim[[1]][["optimalRegularizationValue"]], withOrWithoutZeroClust=x_optim[[1]][["withOrWithoutZeroClust"]], iterations=2, ids=x[,1])
#'
#' #And finally create all the clusters
#' variableViolinsAllClust(x_$penalizedClusterCenters, as.numeric(x_$clusterVector), x[,2:ncol(x)])
#' @export variableViolinsAllClust
variableViolinsAllClust <- function(sparsityMatrix, clusterVector, inDataFrame){

  number <- sort(unique(clusterVector))

  percentClusterVector <- minMaxScale(clusterVector, multiplicationFactor=100)


  percentNumber <- minMaxScale(number, multiplicationFactor=100)
  paletteColors <- palette(rev(rich.colors(100, plot=FALSE)))[1 + 0.98*(102-percentNumber)]
  dev.off()

  #Here, a directory for all the subdirectories for each cluster is made
  directoryName <- "Cluster expressions"
  dir.create(directoryName)
  workingDirectory <- getwd()
  setwd(paste(workingDirectory, directoryName, sep="/"))

  for(i in 1:length(number)){

    #Here, a specific directory for the graphics are made.
    directoryName <- paste("Cluster", number[i])
    dir.create(directoryName)
    workingDirectoryClusters <- getwd()
    setwd(paste(workingDirectoryClusters, directoryName, sep="/"))


    #This code is an efficient way of giving all rows in the "Clusters" column the same name, except for the rows with the cluster of interest.

    clustIndicesSpecific <- sapply(clusterVector, singleEventClusterNaming, n=number[i])

    #Create a color vector for the visualzation
    clustColorsSpecific <- c(paletteColors[i], "#d3d3d3")

    #Here, the mu variables for the specific cluster is extracted
    oneClustAllMu <- sparsityMatrix[rownames(sparsityMatrix)==number[i],]

    #Here the variable names is exported
    allVarNames <- colnames(inDataFrame)
    #Then a list is created that contin the objects for the prot creation
    oneClustAllVarList <- mapply(createAllClustOneVarMu, inDataFrame, oneClustAllMu, allVarNames, MoreArgs=list(clust=clustIndicesSpecific, cols=clustColorsSpecific, clustNum=number[i]), SIMPLIFY=FALSE)
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

