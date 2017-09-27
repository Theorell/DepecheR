# Function to find the optimal penalty value for penalized K means
#
#
# This function is used before the dClust function to identify the optimal penalty value for the specific dataset. This value decides the penalization and consequently also the number of clusters that are identified.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
#' @importFrom Rcpp evalCpp
#' @importFrom graphics box
# @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the dScale function.
# @param initCenters Number of starting points for clusters. This essentially means that it is the highest possible number of clusters that can be defined. The higher the number, the greater the precision, but the computing time is also increased with the number of starting points. Default is 30
# @param maxIter The maximal number of iterations that are performed to reach the minimal improvement. 
# @param minImprovement This is connected to the evaluation of the performance of the algorithm. The more iterations that are run, the smaller will the improvement be, and this sets the threshold when the iterations stop. 
# @param bootstrapObservations The number of observations that are included in each bootstrap subsampling of the data. NB! The algorithm uses resampling, so the same event can be used twice.
# @param penalties These values are the ones that are evaluated and the ones that decide the penalization. The number of suggested default values are empirically defined and might not be optimal for a specific dataset, but the algorithm will warn if the most optimal values are on the borders of the range. Note that when this offset is 0, there is no penalization, which means that the algorithm runs normal K-means clustering.
# @param makeGraph If a graph should be created showing the distance between bootstraps under different penalties.
# @param graphName The name of the graph. 
# @param disableWarnings If the lowest or highest penalty is the most optimal, the function gives a warning. This is suppressed with this command. Mostly to simplify use within other functions, such as dOptSampleSize.
# @seealso \code{\link{dClust}}, \code{\link{pClustPredict}}, \code{\link{dOptSampleSize}}
# @return A graph showing the performance of the algorithm under the different penalty values and a list with two components:
# \describe{
#     \item{penaltyOpt.df}{A dataframe with one row with all the information about which settings that were used to generate the optimal clustering. The "withOrigoClust" information tells the user if the solution with or without a cluster in origo gives the most optimal solution. If yes, this origo population is generally small and could be viewed as not fitting in the model.}
#     \item{meanOptimDf}{A dataframe with the information about the results with all tested penalty values}
# }
# @examples
# #Generate a default size dataframe with bimodally distributed data
# x <- generateBimodalData()
#
# #Scale this datamframe
# x_scaled <- dScale(x[,2:ncol(x)])
#
# #Set a reasonable working directory, e.g.
# setwd("~/Desktop")
#
# #Run the function
# x_optim <- dOptPenalty(x_scaled, bootstrapObservations=1000)
# @export dOptPenalty
#' @useDynLib DepecheR
dOptPenalty <- function(inDataFrameScaled, initCenters=30, maxIter=100, minImprovement=0.01, bootstrapObservations=10000, penalties=c(0,2,4,8,16,32,64,128), makeGraph=TRUE, graphName="Distance as a function of penalty values.pdf", disableWarnings=FALSE){

  #The constant k is empirically identified by running a large number of penalty values for a few datasets.
  k <- ((bootstrapObservations*sqrt(ncol(inDataFrameScaled)))/1450)
  penalty <- penalties*k

  chunkSize <- detectCores() - 1

  dataMat<-data.matrix(inDataFrameScaled, rownames.force = NA)

  #This command is reiterated the number of times that is needed to reach a minimal improvement of the distance. 
  iter <- 1
  std=1
  distanceBetweenMinAndBestPrevious=-1
  iterTimesChunkSize <- 1
  
  while((iterTimesChunkSize<=maxIter && std>=minImprovement && distanceBetweenMinAndBestPrevious<0) || iterTimesChunkSize<=14){

    cl <-  parallel::makeCluster(chunkSize, type = "SOCK")
    registerDoSNOW(cl)
    optimList <- foreach(i=1:chunkSize) %dopar% grid_search(dataMat,initCenters,penalty,1,bootstrapObservations)
    parallel::stopCluster(cl)	

    #Alternative deprecated parallelization. 
    #if(Sys.info()['sysname']!="Windows"){
    #  cl <- makeCluster(chunkSize, type="FORK")
    #  } else {
    #  cl <- makeCluster(spec="PSOCK", names=chunkSize)
    #  }
    #optimList <- parLapply(cl,1:chunkSize,function(x) grid_search(dataMat,initCenters,penalty, 1,bootstrapObservations,x))
    
      
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
 
    #Here, the standard deviation of the error of the mean (or something like that) is instead retrieved
    stdOptimList <- list()	
    for(i in 1:length(optimListFull[[1]])){
      x <- do.call("rbind", lapply(optimListFull, "[[", i))
      stdOptimList[[i]] <- apply(x, 2, sd)
    }
    stdOptimDf <- (as.data.frame(stdOptimList))/sqrt(iter*chunkSize) 

    #Now average these runs, to retrieve the best penalty value between the runs. 
    meanOptimList <- list()	
    for(i in 1:length(optimListFull[[1]])){
      x <- do.call("rbind", lapply(optimListFull, "[[", i))
      meanOptimList[[i]] <- apply(x, 2, mean)
    }
    meanOptimDf <- as.data.frame(meanOptimList)

    #Turn these into vectors
    meanOptimVector <- as.vector(t(meanOptimDf[,1:2]))
    stdOptimVector <- as.vector(t(stdOptimDf[,1:2]))
	  #Return the position of the minmum value
    minPos <- which(meanOptimVector==min(meanOptimVector))
	  
    #Add the standard deviation of this position to its mean
    meanPlusStdMin <- meanOptimVector[minPos]+stdOptimVector[minPos]
    
	  #Return the positions of all values that are not minimum
	  allNonMinPos <- which(meanOptimVector!=min(meanOptimVector))
	  
	  #Now subtract the standard deviation of each of these values from the mean
	  meanMinusStdAllNonMin <- meanOptimVector[allNonMinPos]-stdOptimVector[allNonMinPos]

	  #Identify the lowest value among these
	  minMeanMinusStdAllNonMin <- min(meanMinusStdAllNonMin)
	  
	  #Now, the distance between minMeanMinusSdAllNonMin and . If they overlap, the iteration has not made it totally clear which point is optimal. 
	  distanceBetweenMinAndBestPrevious <- minMeanMinusStdAllNonMin-meanPlusStdMin
	  
	  #Finally, another criterion on the gain of adding more rows is included
	  std <- stdOptimVector[minPos]
	  iterTimesChunkSize <- iter*chunkSize

	  iter <- iter+1
  }
	
  print(paste("The optimization was iterated ", iter*chunkSize, " times.", sep=""))
  
  if(iter>maxIter && std>minImprovement){
    warning("The improvement was still larger than minImprovement when maxIter was reached")
  }
  
  #Here, the penalty values for the real dataset, not the bootstrap subsamples, is used.
  realK <- ((nrow(inDataFrameScaled)*sqrt(ncol(inDataFrameScaled)))/1450)
  rownames(meanOptimDf) <- penalties
  colnames(meanOptimDf) <- c("distWZero", "distWOZero", "nClustWZero", "nClustWOZero")

  penaltyOpt.df <- as.data.frame(as.numeric(row.names(which(meanOptimDf[,1:2]==min(meanOptimDf[,1:2]), arr.ind=TRUE))))

  colnames(penaltyOpt.df)[1] <- "bestPenalty"

  #Export if the solution with or without zero clusters give the optimal result
  penaltyOpt.df$withOrigoClust <- colnames(meanOptimDf)[which(meanOptimDf[,1:2]==min(meanOptimDf[,1:2]), arr.ind=TRUE, useNames=TRUE)[2]]


  #In the exceptional event that one solution is as good for both with and without a origo cluster,  this argument is added that always prefers the solution without an origo cluster. It also prefers solutions with a higher penalty in cases where two origo-cluster free solutions give identical results
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
      print("Warning: the lowest penalty was the most optimal in the range. It is suggested to run with a few lower penalty values to make sure that the most optimal has been found")
    }
    if(penaltyOpt.df$bestPenalty==highestPenalty){
      print("Warning: the highest penalty was the most optimal in the range. It is suggested to run with a few higher penalty values to make sure that the most optimal has been found")
    }
  
  }

  #Export the used initCenters, as this needs to be used also when running dClust based on the optimizations.
  penaltyOpt.df$initCenters <- initCenters

  #Here, the optimization is plotted if wanted.
  if(makeGraph==TRUE){
    pdf(graphName)
    par(mar=c(5, 4, 4, 6) + 0.1)
    #Plot the data
    plot(row.names(meanOptimDf), meanOptimDf[[penaltyOpt.df$withOrigoClust]], pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="", type="b",col="black", main="Distance between bootstraps as a function of penalties values")
    axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("Distance between bootstraps",side=2,line=2.5)
    box()
    
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

#Return the list.
	penaltyOptList <- list(penaltyOpt.df, meanOptimDf)
	return(penaltyOptList)
}
