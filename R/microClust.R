#' This function is the core of the groupProbPlot. See the documentation
#' there for details. 
#' @param dataCenter The original data. 
#' @param dataNeigh The data for the neighbors. Often stronly overlapping with
#' the dataCenter, but for internal reasons, this data cloud is larger
#' than the dataCenter cloud. 
#' @param dataReturn The neighbor data that should be aggregated and sent back.
#' @param method Should median or mean be calculated?
#' @param k Number of neighbors. 
#' @param trim If mean of the neighbors is returned, should it be calculated
#' with trimming?
#' @return A dataset with the same shape as dataCenter, filled with 
#' aggregated information from the k nearest neighbors. 
#' @importFrom FNN knnx.index
#' @importFrom robustbase colMedians
microClust <- function(dataCenter, dataNeigh, dataReturn, 
                       method = "median", k = 11, trim = 0) {
    closest10Pos <- knnx.index(dataNeigh, dataCenter, k = k, 
                               algorithm = "cover_tree")
    if (method == "median") {
        if (is.matrix(dataReturn)) {
            closest10Result <- lapply(
                seq_len(nrow(closest10Pos)),
                function(x)
                    colMedians(dataReturn[closest10Pos[x, ], ],
                        hasNA = FALSE
                    )
            )
            resultDf <- do.call("rbind", closest10Result)
        } else {
            closest10Result <- vapply(
                seq_len(nrow(closest10Pos)),
                function(x)
                    median(dataReturn[closest10Pos[x, ]]), 1
            )
        }
    } else {
        if (is.matrix(dataReturn)) {
            closest10Result <- lapply(
                seq_len(nrow(closest10Pos)),
                function(x)
                    colMeans(dataReturn[closest10Pos[x, ], ])
            )
            resultDf <- do.call("rbind", closest10Result)
        } else {
            closest10Result <- vapply(
                seq_len(nrow(closest10Pos)),
                function(x) mean(dataReturn[closest10Pos[x, ]]),
                1
            )
        }
    }
}