#' @importFrom graphics box
#' @importFrom ggplot2 ggplot aes geom_line ggtitle xlab ylab ylim ggsave
#' @importFrom dplyr sample_n
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %do% %dopar%
#' @export dOptSubset
dOptSubset <- function(inDataFrameScaled, sampleSizes, fMeasureSampleSize="default", k=30, maxIter=100, maxCRI=0.01, minCRIImprovement=0.01, penalties, withOrigoClust){

    dOptPenaltyResultList <- list()
	  lowestDist <- vector()

	 #This loop continues to run until two runs produce distances that diverge less than 0.01.
	  for(i in 1:length(sampleSizes)){
	    if(i>1){
	      print(paste("Sample size ", sampleSizes[i-1], "completed. Jumping to next sample size."))
	    }
	    
	    dOptPenaltyResult <- dOptPenalty(inDataFrameScaled, k=k, maxIter=maxIter, bootstrapObservations=sampleSizes[i], penalties=penalties, makeGraph=FALSE, disableWarnings=TRUE)
	    dOptPenaltyResultList[[i]] <- dOptPenaltyResult
	    lowestDist[i] <- min(dOptPenaltyResult[[2]][,1:2])
	    print(paste("The lowest corrected Rand index with sample size ", sampleSizes[i], " is ", lowestDist[i], sep=""))
	    
	    if(lowestDist[length(lowestDist)]<=maxCRI){
	      break
	    }
	  }

	 if(length(sampleSizes)>1 && lowestDist[i]<=maxCRI){
	   print(paste("Cycle", i, "optimal. No need to run larger datasets."))	   
	 }  
	 #Here, the optimal solution is retrieved from the last cycle.
	 resultOfOptimalSettings <- dOptPenaltyResultList[[length(dOptPenaltyResultList)]]
	  
	 #Now, the curve of distances with different penalties and different sample sizes are plotted together if more than one sample size is tested. 
	  if(length(sampleSizes)>1){
	    
	    #Here, all values for all sample sizes are made into a long dataframe
	    dOptPenaltyResults <- do.call("rbind", lapply(dOptPenaltyResultList, "[[", 2))
	    
	    #Here, the solution with or without origo cluster is retrieved, depending on the optimal solution
	    
	    if(resultOfOptimalSettings[[1]][1,2]=="yes"){
	      Distances <- dOptPenaltyResults[,1]
	    }
	    if(resultOfOptimalSettings[[1]][1,2]=="no"){
	      Distances <- dOptPenaltyResults[,2]
	    }
	    
	    #Now, a vector of penalties is created
	    Penalties <- rep(penalties, times=length(dOptPenaltyResultList))
	    #Here, a vector of sample sizes is created instead
	    SampleSizes <- as.factor(rep(round(sampleSizes), each=length(penalties)))
	    #Now combine these three
	    plottingObject <- data.frame("Sample_sizes"=SampleSizes, Penalties, Distances)
	    
	    ggplot(data=plottingObject,
	           aes(x=sqrt(Penalties), y=Distances, colour=SampleSizes)) +
	      geom_line() + ylim(c(0,1)) +
	      ggtitle("Distance as a function of penalties for different sample sizes") +
	      xlab("Square root of the penalties") + ylab("Distance  between bootstrap subsamplings, [0 to 1]")
	    ggsave("Distance as a function of penalties for different sample sizes.pdf")
	    
	    #Create a vector of improvements between sample sizes
	    improvement <- lowestDist
	    #Here, the first position, that cannot be related to any smaller fractions, is turned to one.
	    improvement[1] <- 1
	    if(length(lowestDist)>1){
	      for(i in 2:length(lowestDist)){
	        improvement[i] <- lowestDist[i-1]-lowestDist[i]
	      } 
	    }
	    
	    sampleSizeOpt <- data.frame("SampleSize"=sampleSizes[1:(length(sampleSizes))], "Lowest distance"=lowestDist, "Improvement"=improvement)
	    
	  }
	  
    #Now, the best run amongst all the runs with the largest sample size is defined, by identifying the solution that gives the highest mean f-measure for all the others.
    
    #First, all solutions for the largest sample size is retrieved
    allSolutions <- resultOfOptimalSettings[[3]]

    #Here, the fMeasureDataSet is created
    if(is.numeric(fMeasureSampleSize)==TRUE){
      fMeasureDataSet <- data.matrix(sample_n(inDataFrameScaled, fMeasureSampleSize, replace=TRUE))
    }
    if(fMeasureSampleSize=="default"){
      if(nrow(inDataFrameScaled)<=10000){
        fMeasureDataSet <- data.matrix(inDataFrameScaled)
      } else {
        fMeasureDataSet <- data.matrix(sample_n(inDataFrameScaled, 10000, replace=TRUE))
      }
    }
    
    #Now, all clusterCenters are used to allocate the fMeasureDataSet.
    allocationResultList <- list()
    noZero <- ifelse(withOrigoClust=="no", 1, 0)
    
    if(ncol(fMeasureDataSet)<500){
        allocationResultList <- foreach(i=1:length(allSolutions)) %do% removeEmptyVariablesAndAllocatePoints(fMeasureDataSet=fMeasureDataSet, clusterCenters=allSolutions[[i]], noZero=noZero)
    } else {
      cl <-  parallel::makeCluster((detectCores() - 1), type = "SOCK")
      registerDoSNOW(cl)
      allocationResultList <- foreach(i=1:length(allSolutions)) %dopar% removeEmptyVariablesAndAllocatePoints(fMeasureDataSet=fMeasureDataSet, clusterCenters=allSolutions[[i]], noZero=noZero)
      parallel::stopCluster(cl)	
    }

    #Here, the mean F-measure with each allocationResult as the id vector and all the others as cluster vectors is identified
    meanFmeasureList <- list()
    n_cores <- detectCores() - 1
    cl <-  parallel::makeCluster(n_cores, type = "SOCK")
    registerDoSNOW(cl)
    meanFmeasureList <- foreach(i=1:length(allocationResultList)) %dopar% mean(unlist(sapply(allocationResultList, fMeasureCalculation, idVector=allocationResultList[[i]], writeFiles=FALSE)))
    parallel::stopCluster(cl)	

    meanFmeasureVector <- unlist(meanFmeasureList)
    #Now the solution being the most similar to all the others is retrieved
    optimalClusterCenters <- unlist(allSolutions[[which(meanFmeasureVector==max(meanFmeasureVector))[1]]])
    
    #And here, the optimal solution is created with the full dataset
    optimalClusterVector <- allocate_points(data.matrix(inDataFrameScaled, rownames.force = NA), optimalClusterCenters, no_zero=noZero)[[1]]
    
    #Here, the optPenalty information is retrieved from the optimal sample size run. 
    optPenalty <- list(dOptPenaltyResultList[[length(dOptPenaltyResultList)]][[1]],dOptPenaltyResultList[[length(dOptPenaltyResultList)]][[2]])
    
    #And here, the optimal results, given if an origo cluster should be included or not, are retrieved further
    retrieveOrigoOrNotList <- retrieveOrigoOrNot(withOrigoClust = withOrigoClust, clusterVector=optimalClusterVector, clusterCenters=optimalClusterCenters, colnamesClusterCenters=colnames(inDataFrameScaled), k=k)
    clusterVectorEquidistant <- retrieveOrigoOrNotList[[1]]
    reducedClusterCenters <- retrieveOrigoOrNotList[[2]]
    
    #Here, the decision to go with a solution excluding origo is tested
    if(withOrigoClust=="no" && min(resultOfOptimalSettings[[2]][,1])-min(resultOfOptimalSettings[[2]][,2])>0.01){
      warning(paste("The best solution with a cluster in origo was ", 100*min(resultOfOptimalSettings[[2]][,1])-min(resultOfOptimalSettings[[2]][,2]), " percent better than the best solution without an origo cluster. Consider changing \"withOrigoClust\" to \"yes\".", sep=""))
    }
    
    #And now, the final product is produced, whose form depends on the input
    if(length(sampleSizes)>1){
      dOptSubsetObject <- list(clusterVectorEquidistant, reducedClusterCenters, optPenalty, sampleSizeOpt)
      names(dOptSubsetObject) <- c("clusterVector", "clusterCenters", "penaltyOptList", "sampleSizeOptList")
    } else {
      dOptSubsetObject <- list(clusterVectorEquidistant, reducedClusterCenters, optPenalty)
      names(dOptSubsetObject) <- c("clusterVector", "clusterCenters", "penaltyOptList")
    }
    return(dOptSubsetObject)
}