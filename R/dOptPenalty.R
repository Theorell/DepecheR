#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
#' @importFrom Rcpp evalCpp
#' @importFrom graphics box
#' @export dOptPenalty
dOptPenalty <- function(inDataFrameScaled, k=30, maxIter=100, minCRIImprovement=0.01, bootstrapObservations=10000, penalties=c(0,2,4,8,16,32,64,128), makeGraph=TRUE, graphName="Distance as a function of penalty values.pdf", disableWarnings=FALSE, returnClusterCenters=TRUE, withOrigoClust="no"){

  #The constant k is empirically identified by running a large number of penalty values for a few datasets.
  penaltyConstant <- ((bootstrapObservations*sqrt(ncol(inDataFrameScaled)))/1450)
  penalty <- penalties*penaltyConstant

  chunkSize <- detectCores() - 1

  dataMat<-data.matrix(inDataFrameScaled, rownames.force = NA)

  #This command is reiterated the number of times that is needed to reach a minimal improvement of the distance. 
  iter <- 1
  std=1
  distanceBetweenMinAndBestPrevious=-1
  iterTimesChunkSize <- 1
  allClusterCentersPenaltySorted <- list()
  cl <-  parallel::makeCluster(chunkSize, type = "SOCK")
  registerDoSNOW(cl)  
  
  while(iterTimesChunkSize<=6 || (std>=minCRIImprovement && iterTimesChunkSize<=maxIter && distanceBetweenMinAndBestPrevious<0)){

    optimList <- foreach(i=1:chunkSize, .packages="DepecheR") %dopar% grid_search(dataMat,k,penalty,1,bootstrapObservations,i)
    
    #Before any further analyses are performed, any penalty that can result in a trivial solution are practically eliminated.
    optimListNonTrivial <- optimList
    for(i in 1:length(optimListNonTrivial)){	  
	  optimListNonTrivial[[i]]$d[which(optimList[[i]]$n<=2)] <- 1
  	  optimListNonTrivial[[i]]$z[which(optimList[[i]]$m==1)] <- 1	 
    }	
    
    #Now, the new list is combined with the older, if there are any
    if(iter==1){
    	optimListFull <- optimListNonTrivial
    } else {
    	optimListFull <- c(optimListFull, optimListNonTrivial)
    }
 
    #Here, the  average and standard deviation of the error of the mean (or something like that) is retrieved
    meanOptimList <- list()	
    stdOptimList <- list()	
    for(i in 1:4){
      x <- do.call("rbind", lapply(optimListFull, "[[", i))
      meanOptimList[[i]] <- apply(x, 2, mean)
      stdOptimList[[i]] <- apply(x, 2, sd)
    }

    stdOptimDf <- (as.data.frame(stdOptimList))/(sqrt(iter*chunkSize))
    meanOptimDf <- as.data.frame(meanOptimList)

    #Turn these into vectors
    meanOptimVector <- as.vector(t(meanOptimDf[,1:2]))
    stdOptimVector <- as.vector(t(stdOptimDf[,1:2]))
	  #Return the position of the minmum value
    minPos <- which(meanOptimVector==min(meanOptimVector))
	  
    #Add the standard deviation of this position to its mean
    meanPlus2StdMin <- meanOptimVector[minPos]+(2*stdOptimVector[minPos])
    
	  #Return the positions of all values that are not minimum
	  allNonMinPos <- which(meanOptimVector!=min(meanOptimVector))
	  
	  #Now subtract the standard deviation of each of these values from the mean
	  meanMinus2StdAllNonMin <- meanOptimVector[allNonMinPos]-(2*stdOptimVector[allNonMinPos])

	  #Identify the lowest value among these
	  minMeanMinus2StdAllNonMin <- min(meanMinus2StdAllNonMin)
	  
	  #Now, the distance between minMeanMinusSdAllNonMin and . If they overlap, the iteration has not made it totally clear which point is optimal. 
	  distanceBetweenMinAndBestPrevious <- minMeanMinus2StdAllNonMin-meanPlus2StdMin
	  
	  #Finally, another criterion on the gain of adding more rows is included
	  std <- stdOptimVector[minPos]
	  iterTimesChunkSize <- iter*chunkSize
	  
	  #print(std)
	  #print(iterTimesChunkSize)
	  #print(distanceBetweenMinAndBestPrevious)


	  if(returnClusterCenters==TRUE){
	    #Here, the cluster center information for each run is saved, one list for the origoCLust solution and another for the nonOrigoClust:
	    #First all cluster center information is saved in one place.
	    allClusterCenters <- sapply(optimList, "[", 5)
	    
	    #Then each penalty and the solutions with and without origo clusters are reorganized and, if iter>1, integrated with previous runs. 
	    for(i in 1:length(allClusterCenters[[1]])){
	      
	      tempPenaltyList <- sapply(allClusterCenters, "[", i)
	      tempOrigoClusterCenters <- c(sapply(tempPenaltyList, "[", 1), sapply(tempPenaltyList, "[", 2))
	      tempNonOrigoClusterCenters <- c(sapply(tempPenaltyList, "[", 3), sapply(tempPenaltyList, "[", 4))
	      
	      if(iter==1){
	        allClusterCentersPenaltySorted[[i]] <- list(tempOrigoClusterCenters, tempNonOrigoClusterCenters)
	        
	      } else {
	        allClusterCentersPenaltySorted[[i]][[1]] <- c(allClusterCentersPenaltySorted[[i]][[1]], tempOrigoClusterCenters)
	        allClusterCentersPenaltySorted[[i]][[2]] <- c(allClusterCentersPenaltySorted[[i]][[2]], tempNonOrigoClusterCenters)
	      }
	    }
	  }
	  print(paste("Set ", iter, " with ", chunkSize, " iterations completed.", sep=""))
	  iter <- iter+1
  }
	
  parallel::stopCluster(cl)	
  
  print(paste("The optimization was iterated ", (iter-1)*chunkSize, " times.", sep=""))

  if(iter>=maxIter && (std>minCRIImprovement || distanceBetweenMinAndBestPrevious<0)){
    warning("An optimal value was not identified before maxIter was reached")
  }
  
  rownames(meanOptimDf) <- round(penalties, digits=1)
  colnames(meanOptimDf) <- c("distWZero", "distWOZero", "nClustWZero", "nClustWOZero")

  penaltyOpt.df <- as.data.frame(as.numeric(row.names(which(meanOptimDf[,1:2]==min(meanOptimDf[,1:2]), arr.ind=TRUE))))

  colnames(penaltyOpt.df)[1] <- "bestPenalty"

  #Export if the solution with or without zero clusters give the optimal result
  penaltyOpt.df$withOrigoClust <- colnames(meanOptimDf)[which(meanOptimDf[,1:2]==min(meanOptimDf[,1:2]), arr.ind=TRUE, useNames=TRUE)[2]]

  #In the exceptional event that one solution is as good both with and without a origo cluster, this argument is added that always prefers the solution without an origo cluster. It also prefers solutions with a higher penalty in cases where two origo-cluster free solutions give identical results.
  if(nrow(penaltyOpt.df)>1 && length(grep("nClustWOZero", penaltyOpt.df$withOrigoClust))>0){
  
    penaltyOpt.df <- penaltyOpt.df[which(penaltyOpt.df[,2]=="nClustWOZero"),]
  }

  if(nrow(penaltyOpt.df)>1){
  
    penaltyOpt.df <- penaltyOpt.df[min(which(penaltyOpt.df[,1]==max(penaltyOpt.df[,1]))),]
  }

  lowestPenalty <- as.numeric(row.names(meanOptimDf[1,]))
  highestPenalty <- as.numeric(row.names(meanOptimDf[nrow(meanOptimDf),]))

  if(disableWarnings==FALSE){
  
    if(penaltyOpt.df$bestPenalty==lowestPenalty){
      warning("The lowest penalty was the most optimal in the range. It might be a good idea to run with a few lower penalty values to make sure that the most optimal has been found")
    }
    if(penaltyOpt.df$bestPenalty==highestPenalty){
      warning("The highest penalty was the most optimal in the range. It might be a good idea to run with a few higher penalty values to make sure that the most optimal has been found")
    }
  
  }

  #Export the used k, as this needs to be used also when running dClust based on the optimizations.
  penaltyOpt.df$k <- k

  #Here, the optimization is plotted if wanted.
  if(makeGraph==TRUE){
    pdf(graphName)
    par(mar=c(5, 4, 4, 6) + 0.1)
    #Plot the data
    plot(row.names(meanOptimDf), meanOptimDf[[penaltyOpt.df$withOrigoClust]], pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="", type="b",col="black", main="Distance between bootstraps as a function of penalties values")
    axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("Distance between bootstraps",side=2,line=2.5)
    graphics::box()
    
    # Draw the penalty axis
    axis(1,pretty(range(as.numeric(row.names(meanOptimDf))), n=10))
    mtext("Penalty values",side=1,col="black",line=2.5)
    
    # Add Legend
    legend("topleft",legend="Distance (low is good)",
           text.col="black",pch=c(16,15),col="black")
    
    dev.off() 
  }


  #Change the resulting names in withOrigoClust to something more meaningful
  penaltyOpt.df$withOrigoClust <- ifelse(penaltyOpt.df$withOrigoClust=="distWZero", "yes", "no")

  if(returnClusterCenters==TRUE){
    #Here, the list of solutions with the best penalty and with or without origo cluster is exported
    allClusterCentersBestPenalty <- allClusterCentersPenaltySorted[[which(round(penalties, digits=1)==penaltyOpt.df$bestPenalty)]]
    
    if(withOrigoClust=="yes"){
      bestClusterCenters <- allClusterCentersBestPenalty[[1]]
    } else {
      bestClusterCenters <- allClusterCentersBestPenalty[[2]]
    }
  }

  #Return the list, that changes depending on the preference 
  if(returnClusterCenters==TRUE){
    penaltyOptList <- list(penaltyOpt.df, meanOptimDf, bestClusterCenters)
  } else{
    penaltyOptList <- list(penaltyOpt.df, meanOptimDf)
  }

	return(penaltyOptList)
}
