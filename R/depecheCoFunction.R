#' @importFrom dplyr sample_n
#' @importFrom gplots heatmap.2
#' @useDynLib DepecheR
depecheCoFunction <- function(inDataFrameScaled, firstClusterNumber=1, directoryName, penalties=c(0,2,4,8,16,32,64,128), sampleSizes="default", selectionSampleSize="default", k=30, minARIImprovement=0.01, minARI=0.95, maxIter=100, ids, newNumbers, createDirectory=FALSE){

    if(createDirectory==TRUE){
    workingDirectory <- getwd()
    dir.create(directoryName)
    setwd(paste(workingDirectory, directoryName, sep="/"))
    }
  
  #First, if the dataset is very, very big, a subset of it is used to subset from. Otherwise the system memory needed to just perform the boot strapping becomes so consuming, that the process stalls.
  if(nrow(inDataFrameScaled)>1000000){
    inDataFrameUsed <- sample_n(inDataFrameScaled, 1000000)
  } else {
    inDataFrameUsed <- inDataFrameScaled
  }
  
  #Here, the sampleSize is set in cases it is "default".
  if(length(sampleSizes)==1){
    if(sampleSizes=="default"){
      if(nrow(inDataFrameScaled)<=10000){
        sampleSizes <- nrow(inDataFrameScaled)
      } else {
        sampleSizes <- 10000
      }
    }
  }

  depecheResult <- dOptSubset(inDataFrameScaled=inDataFrameUsed, firstClusterNumber=firstClusterNumber, sampleSizes=sampleSizes, k=k, maxIter=maxIter, minARI=minARI, minARIImprovement=minARIImprovement, penalties=penalties, selectionSampleSize=selectionSampleSize)
    
  #
  #Here the data is added back, in the cases where very large datasets are used
  
  if(nrow(inDataFrameScaled)>1000000){
    depecheResult$clusterVector <- dAllocate(inDataFrameScaled, depecheResult$clusterCenters)
  }
  ######################################
  
  #Provided that a viable Id vector is added, a table with the percentage of cells in each cluster for each individual is created
  
  if(missing(ids)==FALSE && length(ids)==nrow(inDataFrameScaled)){
    
    clusterTable <- table(depecheResult$clusterVector, ids)
    
    countTable <- table(ids)
    
    clusterFractionsForAllIds <- clusterTable
    
    for(i in 1:length(countTable)){
      x <- clusterTable[,i]/countTable[i]
      clusterFractionsForAllIds[,i] <- x
    }
    
    nextClustResultPosition <- length(depecheResult)+1
    depecheResult[[nextClustResultPosition]] <- as.data.frame.matrix(clusterFractionsForAllIds)
    names(depecheResult)[[length(depecheResult)]] <- "idClusterFractions"
    
  }
  
  #Here, a heatmap over the cluster centers is saved. Only true if the number of clusters exceeds one.
  reducedClusterCentersColRow <- depecheResult[[2]]
  if(nrow(reducedClusterCentersColRow)>1 && ncol(reducedClusterCentersColRow)>1){
    reducedClusterCentersColRow[reducedClusterCentersColRow==0] <- NA
    colorLadder <- colorRampPalette(c("blue", "white", "red"))(21)
    pdf("Cluster centers.pdf")
    heatmap.2(as.matrix(reducedClusterCentersColRow),Rowv=FALSE, Colv=FALSE, dendrogram="none", scale="none", col=colorLadder, trace="none", na.color="#A2A2A2")
    dev.off()    
  }
  
  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }
  
  return(depecheResult)
  
}