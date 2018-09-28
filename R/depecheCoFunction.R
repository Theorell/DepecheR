#' @importFrom dplyr sample_n
#' @importFrom gplots heatmap.2
#' @importFrom graphics box
#' @importFrom ggplot2 ggplot aes geom_line ggtitle xlab ylab ylim ggsave
#' @importFrom dplyr sample_n
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %do% %dopar%
#' @useDynLib DepecheR
depecheCoFunction <- function(inDataFrameScaled, firstClusterNumber=1, directoryName, penalties, sampleSize, selectionSampleSize, k, minARIImprovement, minARI, maxIter, ids, newNumbers, createDirectory=FALSE, createOutput){

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
  if(length(sampleSize)==1){
    if(sampleSize=="default"){
      if(nrow(inDataFrameScaled)<=10000){
        sampleSize <- nrow(inDataFrameScaled)
      } else {
        sampleSize <- 10000
      }
    }
  }

  dOptPenaltyResult <- dOptPenalty(inDataFrameScaled, k=k, maxIter=maxIter, sampleSize=sampleSize, penalties=penalties, makeGraph=createOutput, minARI=minARI)
  
  #Now over to creating the final solution
  #Here, the selectionDataSet is created
  if(selectionSampleSize=="default"){
    selectionSampleSize <- sampleSize
  }
  if(nrow(inDataFrameScaled)<=selectionSampleSize){
      selectionDataSet <- inDataFrameScaled
    } else {
      selectionDataSet <- sample_n(inDataFrameScaled, selectionSampleSize, replace=TRUE)
  }
  
  #If the dataset is small, a new set of seven clusterings are performed (on all the data or on a subsample, depending on the sample size), and the maximum likelihood solution is returned as the result
  if(nrow(inDataFrameScaled)<10000){
    penalty <- dOptPenaltyResult[[1]][1,1]
    depecheAllDataResult <- depecheAllData(inDataFrameScaled, penalty=penalty, k=k, firstClusterNumber=firstClusterNumber)
    clusterVectorEquidistant <- depecheAllDataResult[[1]]
    reducedClusterCenters <- depecheAllDataResult[[2]]
  } else {
    #Now, the best run amongst all the runs with the largest sample size is defined, by identifying the solution that gives the highest mean f-measure for all the others.
    
    #First, all solutions are retrieved
    allSolutions <- dOptPenaltyResult[[3]]
    
    #Now, all clusterCenters are used to allocate the selectionDataSet.
    allocationResultList <- list()
    
    selectionDataSetMatrix <- data.matrix(selectionDataSet)
    
    allocationResultList <- foreach(i=1:length(allSolutions)) %do% removeEmptyVariablesAndAllocatePoints(selectionDataSet=selectionDataSetMatrix, clusterCenters=allSolutions[[i]])
    
    
    #Here, the corrected Rand index with each allocationResult as the first vector vector and all the others as individual second vectors is identified
    n_cores <- detectCores() - 1
    cl <-  parallel::makeCluster(n_cores, type = "SOCK")
    registerDoSNOW(cl)
    meanARIList <- foreach(i=1:length(allocationResultList)) %dopar% mean(sapply(allocationResultList, rand_index, inds2=allocationResultList[[i]], k=k))
    parallel::stopCluster(cl)	
    meanARIVector <- unlist(meanARIList)
    
    #Now the solution being the most similar to all the others is retrieved
    optimalClusterCenters <- unlist(allSolutions[[which(meanARIVector==max(meanARIVector))[1]]])
    
    #And here, the optimal solution is created with the full dataset
    optimalClusterVector <- allocate_points(data.matrix(inDataFrameScaled, rownames.force = NA), optimalClusterCenters, no_zero=1)[[1]]
    
    #And here, the optimal results are made more dense by removing empty rows and columns, etc.
    #Here, the numbers of the removed clusters are removed as well, and only the remaining clusters are retained. As the zero-cluster is not included, the first cluster gets the denomination 1.
    clusterVectorEquidistant <- turnVectorEquidistant(optimalClusterVector, startValue=firstClusterNumber)			
    colnames(optimalClusterCenters) <- colnames(inDataFrameScaled)
    
    #Remove all rows and columns that do not contain any information
    reducedClusterCenters <- optimalClusterCenters[which(rowSums(optimalClusterCenters)!=0),which(colSums(optimalClusterCenters)!=0)]
    
    #In the specific case that only one row is left, due to a high penalty, the data needs to be converted back to a matrix from a vector. The same is done if the number of informative variables is just one.
    if(length(which(rowSums(optimalClusterCenters)!=0))==1){
      reducedClusterCenters <- t(reducedClusterCenters)
    } else if(length(which(colSums(optimalClusterCenters)!=0)==1)){
      reducedClusterCenters <- as.matrix(reducedClusterCenters)
    }
    
    
    #Make the row names the same as the cluster names in the clusterVectorEquidistant
    rownames(reducedClusterCenters) <- rep(firstClusterNumber:(firstClusterNumber+(nrow(reducedClusterCenters))-1))
    
  }
  
  #Here, the optPenalty information is retrieved from the optimal sample size run. 
  optPenalty <- list(dOptPenaltyResult[[1]],dOptPenaltyResult[[2]])
  
  depecheResult <- list(clusterVectorEquidistant, reducedClusterCenters, optPenalty)
  names(depecheResult) <- c("clusterVector", "clusterCenters", "penaltyOptList")
  
  
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
    if(createOutput==TRUE){
      pdf("Cluster centers.pdf")
      heatmap.2(as.matrix(reducedClusterCentersColRow),Rowv=FALSE, Colv=FALSE, dendrogram="none", scale="none", col=colorLadder, trace="none", na.color="#A2A2A2")
      dev.off()  
    }
  }
  
  if(createDirectory==TRUE){
    setwd(workingDirectory)
  }
  
  return(depecheResult)
  
}