#' @importFrom graphics box
#' @importFrom ggplot2 ggplot aes geom_line ggtitle xlab ylab ylim ggsave
#' @importFrom dplyr sample_n
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %do% %dopar%
#' @export dOptSubset
dOptSubset <- function(inDataFrameScaled, sampleSizes, selectionSampleSize="default", k=30, maxIter=100, minARI=0.95, minARIImprovement=0.01, penalties, firstClusterNumber){

    dOptPenaltyResultList <- list()
	  highestARI <- vector()

	  makeGraph <- ifelse(length(sampleSizes)>1, FALSE, TRUE)
	  
	  #if((clusterOnAllData=="default" && nrow(inDataFrameScaled)<10000) || (is.logical(clusterOnAllData)==TRUE && clusterOnAllData==TRUE)){
	 #   returnClusterCenters <- FALSE
	  #}  else {
	   # returnClusterCenters <- TRUE
	  #}
	  
	  
	 #This loop continues to run until two runs produce distances that diverge less than 0.01.
	  for(i in 1:length(sampleSizes)){
	    if(i>1){
	      print(paste("Sample size ", sampleSizes[i-1], "completed. Jumping to next sample size."))
	    }
	    
	    dOptPenaltyResult <- dOptPenalty(inDataFrameScaled, k=k, maxIter=maxIter, bootstrapObservations=sampleSizes[i], penalties=penalties, makeGraph=makeGraph, disableWarnings=TRUE, minARI=minARI)
	    dOptPenaltyResultList[[i]] <- dOptPenaltyResult

	    highestARI[i] <- max(dOptPenaltyResult[[2]][,1])
	    print(paste("The highest adjusted Rand index with sample size ", sampleSizes[i], " is ", highestARI[i], sep=""))
	    
	    if(highestARI[length(highestARI)]>=minARI){
	      sampleSizes <- sampleSizes[1:i]
	      break
	    }
	  }

	  if(length(sampleSizes)>1 && highestARI[length(highestARI)]>=minARI){
	    print(paste("Cycle", i, "optimal. No need to run larger datasets."))	   
	  }  
	    
	 #Here, the optimal solution is retrieved from the last cycle.
	 resultOfOptimalSettings <- dOptPenaltyResultList[[length(dOptPenaltyResultList)]]
	  
	 #Now, the curve of distances with different penalties and different sample sizes are plotted together if more than one sample size is tested. 
	  if(length(sampleSizes)>1){
	    
	    #Here, all values for all sample sizes are made into a long dataframe
	    dOptPenaltyResults <- do.call("rbind", lapply(dOptPenaltyResultList, "[[", 2))
	    
	    #Here, the solution with or without origo cluster is retrieved, depending on the optimal solution
	    Distances <- dOptPenaltyResults[,1]

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
	    improvement <- highestARI
	    #Here, the first position, that cannot be related to any smaller fractions, is turned to one.
	    improvement[1] <- 1
	    if(length(highestARI)>1){
	      for(i in 2:length(highestARI)){
	        improvement[i] <- highestARI[i]-highestARI[i-1]
	      } 
	    }
	    
	    sampleSizeOpt <- data.frame("SampleSize"=sampleSizes[1:(length(sampleSizes))], "Highest ARI"=highestARI, "Improvement"=improvement)
	    
	  }
	  
	 #Now over to creating the final solution
	 #Here, the selectionDataSet is created
	 if(is.numeric(selectionSampleSize)==TRUE){
	   selectionDataSet <- data.matrix(sample_n(inDataFrameScaled, selectionSampleSize, replace=TRUE))
	 }
	 if(selectionSampleSize=="default"){
	   if(nrow(inDataFrameScaled)<=10000){
	     selectionDataSet <- inDataFrameScaled
	   } else if(sampleSizes[length(highestARI)]<=10000){
	     selectionDataSet <- sample_n(inDataFrameScaled, 10000, replace=TRUE)
	   } else {
	     selectionDataSet <- sample_n(inDataFrameScaled, sampleSizes[length(highestARI)], replace=TRUE)
	   }
	 }
	  
	 #If the dataset is small, a new set of seven clusterings are performed (on all the data or on a subsample, depending on the sample size), and the maximum likelihood solution is returned as the result
	 #This is not practivally removed
	 if(nrow(inDataFrameScaled)<1){
	   penalty <- dOptPenaltyResultList[[length(dOptPenaltyResultList)]][[1]][1,1]
	   dClustAllDataResult <- dClustAllData(selectionDataSet, penalty=penalty, k=k, firstClusterNumber=firstClusterNumber)
	   clusterVectorEquidistant <- dClustAllDataResult[[1]]
	   reducedClusterCenters <- dClustAllDataResult[[2]]
	 } else {
	   #Now, the best run amongst all the runs with the largest sample size is defined, by identifying the solution that gives the highest mean f-measure for all the others.
	   
	   #First, all solutions for the largest sample size is retrieved
	   allSolutions <- resultOfOptimalSettings[[3]]
	   
	   #Now, all clusterCenters are used to allocate the selectionDataSet.
	   allocationResultList <- list()
	   
	   selectionDataSetMatrix <- data.matrix(selectionDataSet)
	   
	   #if(ncol(selectionDataSet)<500){
	     allocationResultList <- foreach(i=1:length(allSolutions)) %do% removeEmptyVariablesAndAllocatePoints(selectionDataSet=selectionDataSetMatrix, clusterCenters=allSolutions[[i]])
	   #} #else {
	     #cl <-  parallel::makeCluster((detectCores() - 1), type = "SOCK")
	     #registerDoSNOW(cl)
	     #allocationResultList <- foreach(i=1:length(allSolutions)) %dopar% removeEmptyVariablesAndAllocatePoints(selectionDataSet=selectionDataSetMatrix, clusterCenters=allSolutions[[i]])
	     #parallel::stopCluster(cl)	
	   #}
	   
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
	 optPenalty <- list(dOptPenaltyResultList[[length(dOptPenaltyResultList)]][[1]],dOptPenaltyResultList[[length(dOptPenaltyResultList)]][[2]])
	 
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