#' Function to find the optimal penalty value for penalized K means
#'
#'
#' This function is used before the pKMRun function to identify the optimal regVec value for the specific dataset. This value decides the penalization and consequently also the number of clusters that are identified.
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @importFrom Rcpp evalCpp
#' @importFrom graphics box
#' @param inDataFrameScaled A dataframe with the data that will be used to create the clustering. The data in this dataframe should be scaled in a proper way. Empirically, many datasets seem to be clustered in a meaningful way if they are scaled with the quantileScale function.
#' @param kVec Number of starting points for clusters. This essentially means that it is the highest possible number of clusters that can be defined. The higher the number, the greater the precision, but the computing time is also increased with the number of starting points. Default is 30
#' @param iterations As it sounds, the number of bootstrap reiterations that are performed.
#' @param bootstrapObservations The number of observations that are included in each bootstrap subsampling of the data. NB! The algorithm uses resampling, so the same event can be used twice.
#' @param regVecOffset These values are the ones that are evaluated and the ones that decide the penalization. The number of suggested default values are empirically defined and might not be optimal for a specific dataset, but the algorithm will warn if the most optimal values are on the borders of the range. Note that when this offset is 0, there is no penalization, which means that the algorithm runs normal K-means clustering.
#' @seealso \code{\link{pKMRun}}, \code{\link{pKMPredict}}
#' @return A graph showing the performance of the algorithm under the different regVec values and a list with two components:
#' \describe{
#'     \item{regOpt.df}{A dataframe with one row with all the information about which settings that were used to generate the optimal clustering. The "withOrWithoutZeroClust" information tells the user if the solution with or without a cluster in origo gives the most optimal solution. Generally it is the clustering without an origo cluster ("stabWOZero) that is optimal.}
#'     \item{meanOptimDf}{A dataframe with the information about the results with all tested regVec values}
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
#' x_optim <- pKMOptim(x_scaled)
#' @export pKMOptim
#' @useDynLib DepecheR
pKMOptim <- function(inDataFrameScaled, kVec=30, iterations=50, bootstrapObservations=10000, regVecOffset=c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)){

#The constant k is empirically identified by running a large number of regVec values for a few datasets.
k <- ((bootstrapObservations*sqrt(ncol(inDataFrameScaled)))/1450)
regVec <- regVecOffset*k

dataMat<-data.matrix(inDataFrameScaled, rownames.force = NA)

	no_cores <- detectCores() - 1
	cl <- makeCluster(no_cores, type="FORK")

	optim_list <- parLapply(cl,1:iterations,function(x) grid_search(dataMat,kVec,regVec,1,bootstrapObservations,x))
	stopCluster(cl)

	meanOptimList <- list()

	for(i in 1:length(optim_list[[1]])){

		x <- do.call("rbind", lapply(optim_list, "[[", i))
		meanOptimList[[i]] <- apply(x, 2, mean)

	}
	meanOptimDf <- as.data.frame(meanOptimList)
#Here, the regVec values for the real dataset, not the bootstrap subsamples, is used.
realK <- ((nrow(inDataFrameScaled)*sqrt(ncol(inDataFrameScaled)))/1450)
	rownames(meanOptimDf) <- as.integer(regVecOffset*realK)
	colnames(meanOptimDf) <- c("stabWZero", "stabWOZero", "nClustWZero", "nClustWOZero")

#The simplest solution to the problem below is to first exclude all regVec values that are trivial (result in one cluster) and then go on to the itentification of the variant that renders the distance between the bootstraps lowest.
		meanOptimDfNoTrivial <- meanOptimDf[meanOptimDf[,3]>1 & meanOptimDf[,4]>1,]

	regOpt.df <- as.data.frame(as.numeric(row.names(which(meanOptimDfNoTrivial[,1:2]==min(meanOptimDfNoTrivial[,1:2]), arr.ind=TRUE))))

	colnames(regOpt.df)[1] <- "optimalRegularizationValue"

#Export if the solution with or without zero clusters give the optimal result
regOpt.df$withOrWithoutZeroClust <- colnames(meanOptimDfNoTrivial)[which(meanOptimDfNoTrivial[,1:2]==min(meanOptimDfNoTrivial[,1:2]), arr.ind=TRUE, useNames=TRUE)[2]]

lowestRegVec <- as.numeric(row.names(meanOptimDf[1,]))
highestRegVec <- as.numeric(row.names(meanOptimDf[nrow(meanOptimDf),]))

	if(regOpt.df$optimalRegularizationValue==lowestRegVec){
		print("Warning: the lowest regVec was the most optimal in the range. It is suggested to run with a few lower regVec values to make sure that the most optimal has been found")
	}
	if(regOpt.df$optimalRegularizationValue==highestRegVec){
		print("Warning: the highest regVec was the most optimal in the range. It is suggested to run with a few higher regVec values to make sure that the most optimal has been found")
	}
#Export the used kVec, as this needs to be used also when running pKMRun based on the optimizations.
regOpt.df$kVec <- kVec

#Here, the optimization is plotted. It should be integrated into the function.

pdf("Distance as a function of regVec values.pdf")
par(mar=c(5, 4, 4, 6) + 0.1)
## Plot first set of data and draw its axis
plot(row.names(meanOptimDf), meanOptimDf[[regOpt.df$withOrWithoutZeroClust]], pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="",
   type="b",col="black", main="Distance between bootstraps as a function of regVec values")
axis(2, ylim=c(0,1),col="black",las=1)  ## las=1 makes horizontal labels
mtext("Distance between bootstraps (low is good)",side=2,line=2.5)
box()

## Allow a second plot on the same graph
par(new=TRUE)


## Draw the regVec axis
axis(1,pretty(range(as.numeric(row.names(meanOptimDf))), n=10))
mtext("RegVec values",side=1,col="black",line=2.5)

## Add Legend
legend("topleft",legend="Distance (low is good)",
  text.col="black",pch=c(16,15),col="black")

dev.off()

#Return the list.
	regOptList <- list(regOpt.df, meanOptimDf)
	return(regOptList)
}
