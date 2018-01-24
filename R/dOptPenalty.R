#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
#' @importFrom Rcpp evalCpp
#' @importFrom graphics box
#' @export dOptPenalty
dOptPenalty <- function(inDataFrameScaled, k=30, maxIter=100, minARIImprovement=0.01, bootstrapObservations=10000, penalties=c(0,2,4,8,16,32,64,128), makeGraph=TRUE, graphName="Distance as a function of penalty values.pdf", disableWarnings=FALSE, returnClusterCenters=TRUE, minARI=0.95){

  #The constant k is empirically identified by running a large number of penalty values for a few datasets.
  penaltyConstant <- ((bootstrapObservations*sqrt(ncol(inDataFrameScaled)))/1450)
  penalty <- penalties*penaltyConstant
  roundPenalties <- round(penalties, digits=1)

  chunkSize <- detectCores() - 1

  dataMat<-data.matrix(inDataFrameScaled, rownames.force = NA)

  #This command is reiterated the number of times that is needed to reach a minimal improvement of the distance. 
  iter <- 1
  std=1
  #distanceBetweenMinAndBestPrevious=-1
  iterTimesChunkSize <- 1
  allClusterCentersPenaltySorted <- list()
  cl <-  parallel::makeCluster(chunkSize, type = "SOCK")
  registerDoSNOW(cl)  
  
  while(iterTimesChunkSize< 20 || (std>=minARIImprovement && iterTimesChunkSize<maxIter)){
    ptm <- proc.time()
    optimList <- foreach(i=1:chunkSize, .packages="DepecheR") %dopar% grid_search(dataMat,k,penalty,1,bootstrapObservations,i)
    
    #Before any further analyses are performed, any penalty that can result in a trivial solution are practically eliminated. 
    optimListNonTrivial <- optimList
    for(i in 1:length(optimListNonTrivial)){	  
	    optimListNonTrivial[[i]]$d[which(optimList[[i]]$n==1)] <- 0
	    #Further, solutions with only one dimension and two clusters are eliminated, as they are artifactual and always results in superior ARI.
	    for(j in 1:length(optimListNonTrivial[[i]]$c)){
	      if(optimList[[i]]$n[j]==2){
	        if(length(which(apply(optimList[[i]]$c[[j]][[1]],2,function(x) !all(x==0))))==1 || length(which(apply(optimList[[i]]$c[[j]][[2]],2,function(x) !all(x==0))))==1){
	          optimListNonTrivial[[i]]$d[j] <- 0
	        }
	      }
	    }
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
    meanOptimDf <- data.frame(meanOptimList[[1]], meanOptimList[[3]])

    #Turn these into vectors
    meanOptimVector <- meanOptimDf[,1]
    stdOptimVector <- stdOptimDf[,1]
	  #Return the position of the minmum value
    maxPos <- which(meanOptimVector==max(meanOptimVector))[1]

	  
    if(length(roundPenalties)>1){
      #Add the standard deviation of this position to its mean
      #meanPlus2StdMin <- meanOptimVector[minPos]+(2*stdOptimVector[minPos])
      
      #Now subtract the standard deviation of each of the non-minimal values from the mean
      
      #meanMinus2StdAllNonMin <- meanOptimVector[-minPos]-(2*stdOptimVector[-minPos])
      
      #Identify the lowest value among these
      #minMeanMinus2StdAllNonMin <- min(meanMinus2StdAllNonMin)
      
      #Now, the distance between minMeanMinusSdAllNonMin and the lowest value. If they overlap, the iteration has not made it totally clear which point is optimal. 
      #distanceBetweenMinAndBestPrevious <- minMeanMinus2StdAllNonMin-meanPlus2StdMin
    }
	  
	  #Finally, another criterion on the gain of adding more rows is included
	  std <- stdOptimVector[maxPos]
	  iterTimesChunkSize <- iter*chunkSize

	  if(returnClusterCenters==TRUE){
	    #Here, the cluster center information for each run is saved, one list for the origoCLust solution and another for the nonOrigoClust:
	    #First all cluster center information is saved in one place.
	    allClusterCenters <- sapply(optimList, "[", 5)
	    
	    #Then each penalty and the solutions with and without origo clusters are reorganized and, if iter>1, integrated with previous runs. 
	    for(i in 1:length(allClusterCenters[[1]])){
	      
	      tempPenaltyList <- sapply(allClusterCenters, "[", i)
	      tempClusterCenters <- c(sapply(tempPenaltyList, "[", 1), sapply(tempPenaltyList, "[", 2))

	      if(iter==1){
	        allClusterCentersPenaltySorted[[i]] <- tempClusterCenters
	        
	      } else {
	        allClusterCentersPenaltySorted[[i]] <- c(allClusterCentersPenaltySorted[[i]], tempClusterCenters)
	      }
	    }
	  }
	  fullTime <- proc.time()-ptm
	  print(paste("Set ", iter, " with ", chunkSize, " iterations completed in ", fullTime[3], " seconds.", sep=""))
	  iter <- iter+1
	  	  
  }
	
  parallel::stopCluster(cl)	
  
  print(paste("The optimization was iterated ", (iter-1)*chunkSize, " times.", sep=""))

  if(iter*chunkSize>=maxIter && std>minARIImprovement){
    warning("The maximum number of iterations was reached before stable optimal solution was found")
  }
  
  rownames(meanOptimDf) <- roundPenalties
  colnames(meanOptimDf) <- c("ARI", "nClust")

  #Here, the optimal penalty is selected. This is defined as the lowest penalty that yields an ARI that is not lower than 0.01 less than the best ARI. 
  penaltyOpt.df <- data.frame("bestPenalty"=roundPenalties[which(meanOptimDf[,1]>max(meanOptimDf[,1])-(1-minARI))][1], k)

  lowestPenalty <- roundPenalties[1]
  highestPenalty <- roundPenalties[length(roundPenalties)]

  if(disableWarnings==FALSE){
  
    if(penaltyOpt.df$bestPenalty==lowestPenalty){
      warning("The lowest penalty was the most optimal in the range. It might be a good idea to run with a few lower penalty values to make sure that the most optimal has been found")
    }
    if(penaltyOpt.df$bestPenalty==highestPenalty){
      warning("The highest penalty was the most optimal in the range. It might be a good idea to run with a few higher penalty values to make sure that the most optimal has been found")
    }
  
  }

  #Here, the optimization is plotted if wanted.
  if(makeGraph==TRUE){
    pdf(graphName)
    par(mar=c(5, 4, 4, 6) + 0.1)
    #Plot the data
    plot(log10(roundPenalties), meanOptimDf[,1], pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="", type="b",col="black", main="Distance between bootstraps as a function of penalties values")
    axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("Adjusted rand index (ARI) between data resamplings",side=2,line=2.5)
    graphics::box()
    
    # Draw the penalty axis
    axis(1,pretty(range(log10(roundPenalties)), n=10))
    mtext("Log10 of penalty values",side=1,col="black",line=2.5)
    
    # Add Legend
    legend("topleft",legend="ARI (high is good)",
           text.col="black",pch=c(16,15),col="black")
    
    dev.off() 
  }

  if(returnClusterCenters==TRUE){
    #Here, the list of solutions with the best penalty is exported
    bestClusterCenters <- allClusterCentersPenaltySorted[[which(roundPenalties==penaltyOpt.df$bestPenalty)]]
  }

  #Return the list, that changes depending on the preference 
  if(returnClusterCenters==TRUE){
    penaltyOptList <- list(penaltyOpt.df, meanOptimDf, bestClusterCenters)
  } else{
    penaltyOptList <- list(penaltyOpt.df, meanOptimDf)
  }

	return(penaltyOptList)
}
