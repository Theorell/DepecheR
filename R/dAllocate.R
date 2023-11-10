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
#' @export dAllocate
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
        clusterCentersScaled <-
            as.matrix(depecheLogCenterSd(as.data.frame(clusterCenters),
                                         log2Off = TRUE,
                                         logCenterSd[[2]],
                                         logCenterSd[[3]])[[1]])
        clusterCentersScaled[which(clusterCenters == 0)] <- 0
    } else {
        inDataFrameScaled <- inDataFrame
        clusterCentersScaled <- clusterCenters
    }

    # Here, all variables that do not
    # contribute to defining a single cluster
    # are removed.
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
    # internal use. In the first case, there are no colnumn names, but the
    # properties of the cluster centers are also more raw and thus informative.
    if (length(colnames(clusterCentersScaled)) > 0) {
        inDataFrameReduced <-
            inDataFrameScaled[, colnames(clusterCentersScaled)]
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
