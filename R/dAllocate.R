#' Allocation of observations to pre-established cluster centers.
#'
#'
#' Here, observations of a dataset are allocated to a set of preestablished
#' cluster centers. This is intended to be used for the test set in train-test
#' dataset situations.
#' @importFrom moments kurtosis
#' @param inDataFrame A dataset that should be allocated to a set of cluster
#' centers, for example a richer, but less representative dataset, with all
#' datapoints from all donors, instead of only a set number of values from all.
#' @param depModel This is the result of the original application of the
#' depeche function on the associated, more representative dataset.
#' @seealso \code{\link{depeche}}
#' @return A vector with the same length as number of rows in the inDataFrame,
#' where the cluster identity of each observation is noted.
#'
#' @examples
#' # Retrieve some example data
#' data(testData)
#' \dontrun{
#' # Now arbitrarily (for the sake of the example) divide the data into a
#' # training- and a test set.
#' testDataSample <- sample(1:nrow(testData), size = 10000)
#' testDataTrain <- testData[testDataSample, ]
#' testDataTest <- testData[-testDataSample, ]
#'
#' # Run the depeche function for the train set
#'
#' depeche_train <- depeche(testDataTrain[, 2:15],
#'     maxIter = 20,
#'     sampleSize = 1000
#' )
#'
#' # Allocate the test dataset to the centers of the train dataset
#' depeche_test <- dAllocate(testDataTest[, 2:15], depeche_train
#' )
#'
#' # And finally plot the two groups to see how great the overlap was:
#' clustVecList <- list(list("Ids" =testDataTrain$ids,
#'                           "Clusters" = depeche_train$clusterVector),
#'                      list("Ids" =testDataTest$ids,
#'                           "Clusters" = depeche_test))
#' tablePerId <- do.call("rbind", lapply(seq_along(clustVecList), function(x){
#'                                       locDat <- clustVecList[[x]]
#'                                       locRes <- apply(as.matrix(table(
#'                                       locDat$Ids, locDat$Clusters)),
#'                                       1, function(y) y/sum(y))
#'                                       locResLong <- reshape2::melt(locRes)
#'                                       colnames(locResLong) <-
#'                                       c("Cluster", "Donor", "Fraction")
#'                                       locResLong$Group <- x
#'                                       locResLong
#'                                       }))
#' tablePerId$Cluster <- as.factor(tablePerId$Cluster)
#' tablePerId$Group <- as.factor(tablePerId$Group)
#'
#' library(ggplot2)
#' ggplot(data=tablePerId, aes(x=Cluster, y=Fraction,
#'         fill=Group)) + geom_boxplot() + theme_bw()
#' }
dAllocate <- function(inDataFrame, depModel) {
    if (is.matrix(inDataFrame)) {
        inDataFrame <- as.data.frame(inDataFrame)
    }
    clusterCenters <- depModel$clusterCenters
    logCenterSd <- depModel$logCenterSd
    #The very first step is that the need for log transformation is ruled out,
    #as this procedure should preferably be done with the whole dataset, before
    #even entering depeche.
    if(logCenterSd[[1]]){
        stop("dAllocate does not work well with data that has been internally",
             "log-transformed, due to large effects of small differences in ",
             "the extreme negative tail of the data.",
             "It is instead recommended to transform the ",
             "full dataset before entering depeche, which makes the log-",
             "transformation internally superflous.")
    }
    #Step one here is to scale and normalize the data in the best possible way.
    #This is preferably done with pre-defined values.
    if(is.list(logCenterSd)){
        inDataFrameScaled <- depecheLogCenterSd(inDataFrame, log2Off = TRUE,
                                                logCenterSd[[2]],
                                                logCenterSd[[3]])[[1]]
        #Annoyingly, we have to add back the empty clusters here for this all
        #to be straight forward.
        if(ncol(clusterCenters) < length(logCenterSd[[2]]) &
           length(logCenterSd[[2]]) == ncol(inDataFrame) &
           length(which(colnames(clusterCenters) %in%
                        colnames(inDataFrame))) == ncol(clusterCenters)){
            extraCols <- matrix(0, nrow = nrow(clusterCenters),
                                ncol = length(logCenterSd[[2]])-
                                    ncol(clusterCenters))
            colnames(extraCols) <- colnames(inDataFrame)[
                -which(colnames(inDataFrame) %in% colnames(clusterCenters))]
            clusterCentersFull <- cbind(clusterCenters, extraCols)
            clusterCentersUsed <- clusterCentersFull[
                ,match(colnames(inDataFrame), colnames(clusterCentersFull))]
        } else {
            stop("Mismatch between the data to be allocated and the model")
        }
        clusterCentersScaled <-
            as.matrix(depecheLogCenterSd(as.data.frame(clusterCentersUsed),
                                         log2Off = TRUE,
                                         logCenterSd[[2]],
                                         logCenterSd[[3]])[[1]])
        clusterCentersScaled[which(clusterCentersUsed == 0)] <- 0
    } else {
        inDataFrameScaled <- inDataFrame
        clusterCentersScaled <- clusterCenters
    }

    # Here, all variables that do not
    # contribute to defining a single cluster
    # are removed again.
    clusterCentersReduced <-
        clusterCentersScaled[
            which(rowSums(clusterCentersScaled) != 0),
            which(colSums(clusterCentersScaled) != 0)
        ]
    # If some variables have been excluded as
    # they did not contribute to construction
    # of any cluster, they are removed from
    # the inData here. The special case with only one variable is taken
    # into account.
    # There are two different methods here: one for external and one for
    # internal use. In the first case, there are no colnunm names, but the
    # properties of the cluster centers are also more raw and thus informative.
    if (length(colnames(clusterCentersReduced)) > 0) {
        inDataFrameReduced <-
            inDataFrameScaled[, colnames(clusterCentersReduced)]
    } else {
        inDataFrameReduced <-
            inDataFrameScaled[, which(colSums(clusterCentersScaled) != 0)]
    }

    # Here, a specific case, namely that only one variable contains
    # meaningful information, is taken into account.
    if (is.vector(inDataFrameReduced)) {
        clusterCentersReduced <- as.matrix(clusterCentersReduced)
        inDataFrameReduced <- as.matrix(inDataFrameReduced)
    }

    clusterReallocationResult <- allocate_points(
        as.matrix(inDataFrameReduced),
        clusterCentersReduced, 1
    )[[1]]

    # As allocate_points spontaneously likes to throw out a cluster called 0,
    # this behaviour is controlled here
    clusterReallocationResult <- clusterReallocationResult + 1

    return(clusterReallocationResult)
}
