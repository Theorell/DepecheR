#' @importFrom dplyr sample_n
#' @importFrom gplots heatmap.2
#' @useDynLib DepecheR
#' @export dClustCoFunction
dClustCoFunction <- function(inDataFrameScaled, firstClusterNumber=1, penalties=c(0,2,4,8,16,32,64,128), sampleSizes="default", selectionSampleSize="default", k=30, minCRIImprovement=0.01, maxCRI=0.05, maxIter=100, ids, newNumbers, clusterOnAllData){
  
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

  dClustResult <- dOptSubset(inDataFrameScaled=inDataFrameUsed, firstClusterNumber=firstClusterNumber, sampleSizes=sampleSizes, k=k, maxIter=maxIter, maxCRI=maxCRI, minCRIImprovement=minCRIImprovement, penalties=penalties, selectionSampleSize=selectionSampleSize, clusterOnAllData=clusterOnAllData)
    
  #
  #Here the data is added back, in the cases where very large datasets are used
  
  if(nrow(inDataFrameScaled)>1000000){
    dClustResult$clusterVector <- dAllocate(inDataFrameScaled, dClustResult$clusterCenters)
  }
  ######################################
  
  #Provided that a viable Id vector is added, a table with the percentage of cells in each cluster for each individual is created
  
  if(missing(ids)==FALSE && length(ids)==nrow(inDataFrameScaled)){
    
    clusterTable <- table(dClustResult$clusterVector, ids)
    
    countTable <- table(ids)
    
    clusterFractionsForAllIds <- clusterTable
    
    for(i in 1:length(countTable)){
      x <- clusterTable[,i]/countTable[i]
      clusterFractionsForAllIds[,i] <- x
    }
    
    nextClustResultPosition <- length(dClustResult)+1
    dClustResult[[nextClustResultPosition]] <- as.data.frame.matrix(clusterFractionsForAllIds)
    names(dClustResult)[[length(dClustResult)]] <- "idClusterFractions"
    
  }
  
  #Here, a heatmap over the cluster centers is saved. Only true if the number of clusters exceeds one.
  reducedClusterCentersColRow <- dClustResult[[2]]
  if(nrow(reducedClusterCentersColRow)>1){
    pdf("Cluster centers.pdf")
    heatmap.2(as.matrix(reducedClusterCentersColRow), col=colorRampPalette(c("blue", "white", "red"))(100), trace="none")
    dev.off()    
  }
  
  return(dClustResult)
  
}