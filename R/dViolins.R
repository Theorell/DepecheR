#' Create violin plots for any variables of choise
#'
#' Here, assymetrical violin plots for each cluster vs all other clusters are
#' plotted for variables either retrieved from a depeche analysis or 
#' user-defined.
#' @param clusterVector Vector with the same length as inDataFrame containing 
#' information about the cluster identity of each observation.
#' @param plotElements This provides information on which features to plot. In 
#' the typical case, this is the essenceElementList from a depeche run. Other
#' input formats are however accepted: if a vector of column names is provided, 
#' then these features will be plotted for all clusters. A custom list of 
#' features specific for each cluster is also accepted. A final alternative is 
#' to return "all" (default), in which case all markers will be plotted for all
#' clusters.If more than a 100 markers are provided, however, this will return 
#' an error. 
#' @param inDataFrame The data used to generate the depecheObject
#' @param colorOrder The order of the cluster colors. Defaults to the order 
#' that the unique values in clusterVector occurs. 
#' @param colorScale The color scale. Options identical to dColorVector.
#' @param directoryName The name of the created directory.
#' @param createOutput For testing purposes. Defaults to TRUE. If FALSE, no 
#' plots are generated.
#' @return One graph is created for each cluster, containing a bean per 
#' specified variable.
#' @seealso \code{\link{dDensityPlot}}, \code{\link{dColorPlot}}, 
#' \code{\link{dColorVector}}, \code{\link{depeche}}
#' @examples
#' # Load some data
#' data(testData)
#' 
#' \dontrun{
#' # Run the clustering function. For more rapid example execution,
#' # a depeche clustering of the data is inluded
#' # testDataDepeche <- depeche(testData[,2:15])
#' data(testDataDepeche)
#' 
#' # Create the plots of the variables that contribute to creating each cluster
#' dViolins(testDataDepeche$clusterVector, testDataDepeche$essenceElementList, 
#' inDataFrame=testData)
#' }
#' 
#' @export dViolins
dViolins <- function(clusterVector, plotElements="all", inDataFrame, 
                     colorOrder = unique(clusterVector), colorScale = "viridis",
                     directoryName = "dViolin_result", createOutput = TRUE) {
    if(createOutput == TRUE){
        dir.create(directoryName)   
    }

    
    allClusterNums <- unique(clusterVector)
    allClusterCols <- dColorVector(allClusterNums, colorOrder = colorOrder, 
                                   colorScale = colorScale)
    
    if(length(plotElements) == 1){
        if(plotElements == "all"){
            plotElements <- colnames(inDataFrame)
            if(length(plotElements)>100){
                stop("All features should be plotted, and these are more than 
                    100 in total, which is unfeasible. Please specify something
                    more reasonable for plotElements and try again.")
            }
        }
        
    }
    
    if(is.list(plotElements)==FALSE){
        localPlotElementList <- 
            list(plotElements)[rep(1,length(allClusterNums))] 
        plotElements <- localPlotElementList
    }
    
    #To reduce the size of the objects involved, all features that should not be
    #plotted at all are excluded. 
    
    allPlotElements <- unique(unlist(plotElements))
    
    inDataFrameReduced <- inDataFrame[,allPlotElements]
    
    #Here the remaining data is truncated. 
    inDataFrameTrunc <- dScale(inDataFrameReduced, scale=TRUE, 
                               robustVarScale = FALSE, center = FALSE, 
                               truncate = TRUE)

    #And now, the inner function is applied once for each cluster
    lapply(seq_along(allClusterNums), function(x) dViolinsCoFunction(
        clusterNum = x, plotElement = plotElements[[x]], 
        clusterVector = clusterVector, inDataFrameTrunc = inDataFrameTrunc,
        clusterCol = allClusterCols[x], directoryName = directoryName, 
        createOutput=createOutput))

    message("Files were saved at ", getwd())
}
