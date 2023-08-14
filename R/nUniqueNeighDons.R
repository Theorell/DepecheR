#' How many donors are contained among the nearest neighbors?
#'
#'
#' This function constructs a variable that for each event shows the number of
#' donors in its nearest neighbor surroundings. It builds on the same
#' idea as has been put forward in the Sconify package:
#' -Burns TJ (2019). Sconify: A toolkit for performing KNN-based statistics for
#' flow and mass cytometry data. R package version 1.4.0 and
#' -Hart GT, Tran TM, Theorell J, Schlums H, Arora G, Rajagopalan S, et al.
#' Adaptive NK cells in people exposed to Plasmodium falciparum correlate
#' with protection from malaria. J Exp Med. 2019 Jun 3;216(6):1280–90.
#' First, the k nearest neighbors are defined for cell x. Then, the number of
#' donors in the k nearest neighbor cloud is returned.
#' @param donorData The donor information.
#' @param euclidSpaceData The data cloud in which the nearest neighbors for the
#' events should be identified. Can be a vector, matrix or dataframe.
#' @param neighRows The rows in the dataset that correspond to the neighbors
#' of the donorData points. This can be all the donorData points, or a subset,
#' depending on the setup.
#' @param ctrlRows Optionally, a set of control rows that are used to remove
#' background signal from the neighRows data before sending the data back.
#' @param kNeighK The number of nearest neighbors.
#' @param kMeansK The number of clusters in the initial step of the algorithm.
#' A higher number leads to shorter runtime, but potentially lower accuracy.
#' @return An object of the same dimensions as donorData that has been smoothed.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom gmodels fast.prcomp
#' @importFrom FNN knnx.index
#' @examples
#' data(testData)
#' data(testDataSNE)
#' euclidSpaceData <-
#'     testData[, c(
#'         "SYK", "CD16", "CD57", "EAT.2",
#'         "CD8", "NKG2C", "CD2", "CD56"
#'     )]
#' \dontrun{
#' nDonorsVector <- nUniqueNeighDons(
#'     donorData = as.numeric(testData$label),
#'     euclidSpaceData
#' )
#' }
#' @export nUniqueNeighDons
nUniqueNeighDons <- function(donorData, euclidSpaceData,
                        neighRows = seq_len(nrow(as.matrix(donorData))),
                        ctrlRows, kNeighK = max(100, round(nrow(
                            as.matrix(euclidSpaceData)
                        ) / 10000)),
                        kMeansK = max(1, round(nrow(
                            as.matrix(euclidSpaceData)
                        ) / 1000))) {
    if (is.vector(donorData)) {
        donorData <- data.frame(donorData, 0)
        donorDataIsVector <- TRUE
    } else {
        donorDataIsVector <- FALSE
    }

    if (is.vector(euclidSpaceData)) {
        euclidSpaceData <- data.frame(euclidSpaceData, euclidSpaceData)
    } else if (is.matrix(euclidSpaceData)) {
        euclidSpaceData <- as.data.frame(euclidSpaceData)
    }
    # First, the dimensionality is reduced to 10 dimensions for the
    # euclidSpaceData
    if (ncol(euclidSpaceData) > 10) {
        dataRedDim <- fast.prcomp(euclidSpaceData)$x[, seq(1, 10)]
    } else {
        dataRedDim <- euclidSpaceData
    }
    # Now, the cells are clustered according to this analysis
    kMeansResult <- kmeans(dataRedDim, kMeansK, iter.max = 100)
    kMeansCenters <- kMeansResult$centers
    kMeansClusters <- kMeansResult$cluster
    print("Done with k-means")

    # Here, the rows connected to the neighbors, the control neighbors or
    # neither are defined
    groupVec <- rep("none", nrow(dataRedDim))
    groupVec[neighRows] <- "neigh"
    if (missing(ctrlRows) == FALSE) {
        groupVec[ctrlRows] <- "ctrl"
    }

    rowNumbers <- seq_len(nrow(dataRedDim))

    donorDataClustList <- split(donorData, kMeansClusters)
    dataRedDimClustList <- split(as.data.frame(dataRedDim), kMeansClusters)
    groupVecClustList <- split(groupVec, kMeansClusters)
    rowNumList <- split(rowNumbers, kMeansClusters)
    # clusterIds <- split(kMeansClustersData, kMeansClustersData)

    if (nrow(kMeansCenters) > 11) {
        distCenters <- as.list(as.data.frame(
            t(knnx.index(kMeansCenters, kMeansCenters, 11))
        ))
    } else {
        distCenters <- as.list(as.data.frame(
            t(knnx.index(kMeansCenters, kMeansCenters, nrow(kMeansCenters)))
        ))
    }

    print("Now the first bit is done, and the iterative part takes off")
    nCores <- detectCores() - 1

    cl <- parallel::makeCluster(nCores, type = "SOCK")
    registerDoSNOW(cl)

    allClusters <- names(donorDataClustList)
    firstCluster <- 1
    resultList <- list()
    x <- 1
    while (firstCluster <= length(allClusters)) {
        timeBefore <- Sys.time()
        if (((firstCluster + nCores) - 1) < length(allClusters)) {
            clusterRange <- firstCluster:((firstCluster + nCores) - 1)
        } else {
            clusterRange <- firstCluster:length(allClusters)
        }

        # Now, the datasets are constructed from this cluster range
        locDataRedDimClustList <- dataRedDimClustList[clusterRange]

        # Here, the neighbors and the control neighbors are found in the 11
        # closest clusters for each of the focus clusters.
        locDistCenters <- distCenters[clusterRange]

        neighCtrlNeighReturnList <- lapply(locDistCenters, function(y) {
            locNeighList <- lapply(y, function(z) {
                locNeigh <- dataRedDimClustList[[z]]
                locReturn <- donorDataClustList[[z]]
                locGroup <- groupVecClustList[[z]]
                list(
                    locNeigh[which(locGroup == "neigh"), ],
                    locNeigh[which(locGroup == "ctrl"), ],
                    locReturn[which(locGroup == "neigh"), ],
                    locReturn[which(locGroup == "ctrl"), ]
                )
            })
            locNeighNeigh <- as.matrix(do.call(
                "rbind",
                lapply(locNeighList, "[[", 1)
            ))
            locNeighReturn <- as.matrix(do.call(
                "rbind",
                lapply(locNeighList, "[[", 3)
            ))
            locCtrlNeigh <- as.matrix(do.call(
                "rbind",
                lapply(locNeighList, "[[", 2)
            ))
            locCtrlReturn <- as.matrix(do.call(
                "rbind",
                lapply(locNeighList, "[[", 4)
            ))
            list(locNeighNeigh, locNeighReturn, locCtrlNeigh, locCtrlReturn)
        })

        resultList[[x]] <- unlist(foreach(
            i = seq_along(locDataRedDimClustList),
            .packages = "DepecheR"
        ) %dopar%
            microClust(
                dataCenter = as.matrix(locDataRedDimClustList[[i]]),
                dataNeigh = neighCtrlNeighReturnList[[i]][[1]],
                dataReturn = neighCtrlNeighReturnList[[i]][[2]],
                method = "nUnique", k = kNeighK
            ))

        timeAfter <- as.numeric(Sys.time() - timeBefore)
        print(paste0(
            "Clusters ", clusterRange[1], " to ",
            clusterRange[length(clusterRange)],
            " smoothed in ", timeAfter, " ",
            attributes(timeAfter)$units, ". Now, ",
            length(allClusters) - clusterRange[length(clusterRange)],
            " clusters are left."
        ))
        firstCluster <- (clusterRange[length(clusterRange)] + 1)
        x <- x + 1
    }

    parallel::stopCluster(cl)
    fullResult <- unlist(resultList)
    fullResultOrdered <- fullResult[order(unlist(rowNumList))]

    return(fullResultOrdered)
}
