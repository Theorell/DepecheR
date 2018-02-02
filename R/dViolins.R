#' Create violin plots for all non-penalized variable for all clusters
#'
#' Here, violin plots of a specific cluster and the total population are created for each variable that has not been penalized away in the penalized K-means analysis. As al such plots are generated for each cluster, this function creates a great number of plots in most instances.
#' @importFrom ggplot2 ggplot aes geom_violin scale_color_manual scale_fill_manual theme_classic labs ggsave
#' @importFrom viridis inferno magma plasma viridis
#' @importFrom gplots rich.colors
#' @importFrom grDevices rainbow
#' @param depecheObject A list object generated with the depeche function, containing a cluster vector and a cluster centers matrix.
#' @param inDataFrame The data used to generate the depecheObject
#' @param order The order that the unique features of the cluster vector should appear in. For harmonization with colorVector and all subsequent functions.
#' @param colorScale The color scale. Inherited from the viridis, gplots and grDevices packages (and the package-specific "dark_rainbow"). Seven possible scales are pre-made: inferno, magma, plasma, viridis, rich_colors, rainbow and dark_rainbow. User specified vectors of colors (e.g. c("#FF0033", "#03AF49")) are also accepted.
#' @param plotAll If all parameters, including the non-contributing, should be plotted for each cluster. Defaults to FALSE.
#' @return One graph is created for each non-penalized variable in each non-penalized cluster, which often means that the function creates a vast number of graphs. The graphs are sorted into subfolders for each cluster.
#' @seealso \code{\link{dDensityPlot}}, \code{\link{dColorPlot}}, \code{\link{colorVector}}
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateBimodalData(samplings=2, dataCols=8)
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Optimize and run the clustering function.
#' xDepecheObject <- depeche(x[2:ncol(x)])
#'
#' #Create the plots of the variables that contribute to creating each cluster
#' dViolins(xDepecheObject)
#' 
#' #Now create plots of all clusters, regardless of if they contributed or not
#' dViolins(xDepecheObject, plotAll=TRUE)
#' @export dViolins
dViolins <- function(depecheObject, inDataFrame, order=unique(clusterVector), colorScale="viridis", plotAll=FALSE){

  clusterVector <- depecheObject$clusterVector
  clusterCenters <- depecheObject$clusterCenters
  
  percentClusterVector <- dScale(clusterVector, scale=c(0,1), robustVarScale=FALSE, center=FALSE, multiplicationFactor=100)

  if(length(colorScale)>1){
    orderColors <- colorRampPalette(colorScale)(length(order)) 
  } else {
    if(colorScale=="inferno"){
      paletteColors <- inferno(length(order)) 
    }
    if(colorScale=="viridis"){
      paletteColors <- viridis(length(order)) 
    }
    if(colorScale=="plasma"){
      paletteColors <- plasma(length(order)) 
    }
    if(colorScale=="magma"){
      paletteColors <- magma(length(order)) 
    }
    if(colorScale=="rich_colors"){
      paletteColors <- rich.colors(length(order)) 
    }
    if(colorScale=="rainbow"){
      paletteColors <- rainbow(length(order)) 
    }
    if(colorScale=="dark_rainbow"){
      paletteColors <- colorRampPalette(c("#990000", "#FFCC00", "#336600", "#000066", "#660033"))(length(order)) 
    }
  }

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

