#' Function to find the optimal penalty value for penalized K means
#'
#'
#' This function is used before the dClust function to identify the optimal penalty value for the specific dataset. This value decides the penalization and consequently also the number of clusters that are identified.
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @importFrom Rcpp evalCpp
#' @importFrom graphics box
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the quantileScale function.
#' @param initCenters Number of starting points for clusters. This essentially means that it is the highest possible number of clusters that can be defined. The higher the number, the greater the precision, but the computing time is also increased with the number of starting points. Default is 30
#' @param iterations As it sounds, the number of bootstrap reiterations that are performed.
#' @param bootstrapObservations The number of observations that are included in each bootstrap subsampling of the data. NB! The algorithm uses resampling, so the same event can be used twice.
#' @param penaltyOffset These values are the ones that are evaluated and the ones that decide the penalization. The number of suggested default values are empirically defined and might not be optimal for a specific dataset, but the algorithm will warn if the most optimal values are on the borders of the range. Note that when this offset is 0, there is no penalization, which means that the algorithm runs normal K-means clustering.
#' @param makeGraph If a graph should be created showing the distance between bootstraps under different penalties.
#' @param graphName The name of the graph. 
#' @param disableWarnings If the lowest or highest penalty is the most otpimal, the function gives a warning. This is suppressed with this command. Mostly to simplify use within other functions, such as dOptSampleSize.
#' @seealso \code{\link{dClust}}, \code{\link{pClustPredict}}, \code{\link{dOptSampleSize}}
#' @return A graph showing the performance of the algorithm under the different penalty values and a list with two components:
#' \describe{
#'     \item{regOpt.df}{A dataframe with one row with all the information about which settings that were used to generate the optimal clustering. The "withOrigoClust" information tells the user if the solution with or without a cluster in origo gives the most optimal solution. If yes, this origo population is generally small and could be viewed as not fitting in the model.}
#'     \item{meanOptimDf}{A dataframe with the information about the results with all tested penalty values}
#' }
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData()
#'
#' #Scale this datamframe
#' x_scaled <- quantileScale(x[,2:ncol(x)])
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Run the function
#' x_optim <- dClustOpt(x_scaled, iterations=10, bootstrapObservations=1000)
#' @export dClustOpt
#' @useDynLib DepecheR
dClustOpt <- function(inDataFrameScaled, initCenters=30, iterations=50, bootstrapObservations=10000, penaltyOffset=c(0,2,4,8,16,32,64,128), makeGraph=TRUE, graphName="Distance as a function of penalty values.pdf", disableWarnings=FALSE){

#The constant k is empirically identified by running a large number of penalty values for a few datasets.
k <- ((bootstrapObservations*sqrt(ncol(inDataFrameScaled)))/1450)
penalty <- penaltyOffset*k

dataMat<-data.matrix(inDataFrameScaled, rownames.force = NA)

if(Sys.info()['sysname']!="Windows"){
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores, type="FORK")
  optimList <- parLapply(cl,1:iterations,function(x) grid_search(dataMat,initCenters,penalty,1,bootstrapObservations,x))
  stopCluster(cl)
} else {
  no_cores <- detectCores() - 1
  cl <- makeCluster(no_cores, type="PSOCK")
  optimList <- parLapply(cl=NULL,1:iterations,function(x) grid_search(dataMat,initCenters,penalty,1,bootstrapObservations,x))
  stopCluster(cl=NULL)
}

#Before any further analyses are performed, any penalty that can result in a trivial solution are practically eliminated.
	optimListNonTrivial <- optimList
	for(i in 1:length(optimListNonTrivial)){
	  
	  optimListNonTrivial[[i]]$d[which(optimList[[i]]$n<=2)] <- 1
	 
	  optimListNonTrivial[[i]]$z[which(optimList[[i]]$m==1)] <- 1
	  
	}	
  
	#Now average these runs	
	meanOptimList <- list()
	
	for(i in 1:length(optimListNonTrivial[[1]])){

		x <- do.call("rbind", lapply(optimListNonTrivial, "[[", i))
		meanOptimList[[i]] <- apply(x, 2, mean)

	}
	meanOptimDf <- as.data.frame(meanOptimList)
#Here, the penalty values for the real dataset, not the bootstrap subsamples, is used.
realK <- ((nrow(inDataFrameScaled)*sqrt(ncol(inDataFrameScaled)))/1450)
	rownames(meanOptimDf) <- penaltyOffset
	colnames(meanOptimDf) <- c("distWZero", "distWOZero", "nClustWZero", "nClustWOZero")

	regOpt.df <- as.data.frame(as.numeric(row.names(which(meanOptimDf[,1:2]==min(meanOptimDf[,1:2]), arr.ind=TRUE))))

	colnames(regOpt.df)[1] <- "bestPenaltyOffset"

#Export if the solution with or without zero clusters give the optimal result
regOpt.df$withOrigoClust <- colnames(meanOptimDf)[which(meanOptimDf[,1:2]==min(meanOptimDf[,1:2]), arr.ind=TRUE, useNames=TRUE)[2]]


#In the exceptional event that one solution is as good for both with and without a origo cluster, this argument is added that always prefers the solution without an origo cluster. It also prefers solutions with a higher penalty in cases where two origo-cluster free solutions give identical results
if(nrow(regOpt.df)>1 && length(grep("nClustWOZero", regOpt.df$withOrigoClust))>0){
  
  regOpt.df <- regOpt.df[which(regOpt.df[,2]=="nClustWOZero"),]
}

if(nrow(regOpt.df)>1){
  
  regOpt.df <- regOpt.df[min(which(regOpt.df[,1]==max(regOpt.df[,1]))),]
}

lowestPenalty <- as.numeric(row.names(meanOptimDf[1,]))
highestPenalty <- as.numeric(row.names(meanOptimDf[nrow(meanOptimDf),]))

if(disableWarnings==FALSE){
  
  if(regOpt.df$bestPenaltyOffset==lowestPenalty){
    print("Warning: the lowest penalty was the most optimal in the range. It is suggested to run with a few lower penalty values to make sure that the most optimal has been found")
  }
  if(regOpt.df$bestPenaltyOffset==highestPenalty){
    print("Warning: the highest penalty was the most optimal in the range. It is suggested to run with a few higher penalty values to make sure that the most optimal has been found")
  }
  
}

#Export the used initCenters, as this needs to be used also when running dClust based on the optimizations.
regOpt.df$initCenters <- initCenters

#Here, the optimization is plotted if wanted.
  if(makeGraph==TRUE){
    pdf(graphName)
    par(mar=c(5, 4, 4, 6) + 0.1)
    #Plot the data
    plot(row.names(meanOptimDf), meanOptimDf[[regOpt.df$withOrigoClust]], pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="",
         type="b",col="black", main="Distance between bootstraps as a function of penaltyOffset values")
    axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("Distance between bootstraps",side=2,line=2.5)
    box()
    
    # Draw the penalty axis
    axis(1,pretty(range(as.numeric(row.names(meanOptimDf))), n=10))
    mtext("penaltyOffset values",side=1,col="black",line=2.5)
    
    # Add Legend
    legend("topleft",legend="Distance (low is good)",
           text.col="black",pch=c(16,15),col="black")
    
    dev.off() 
  }


#Dirty solution to change the resulting names in withOrigoClust to something more meaningful
regOpt.df$withOrigoClust <- ifelse(regOpt.df$withOrigoClust=="distWZero", "yes", "no")

#Return the list.
	regOptList <- list(regOpt.df, meanOptimDf)
	return(regOptList)
}
