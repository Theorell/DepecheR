#' @importFrom moments kurtosis
#This function is used internally in depeche and dAllocate.
#For all options, see depeche or dAllocate, to whom this function is a
#sub-module. It generates a list containing a scaled and centered dataframe, as
#well as a list that clarifies if and how log, center and scale were applied.
depecheLogCenterSd <- function(inDataFrame, log2Off, center, scale){
    logCenterSd <- list(FALSE, FALSE, 1)
    # Here it is checked if the data has very
    # extreme tails, and if so, the data is
    # log2 transformed
    if (log2Off == FALSE && kurtosis(as.vector(as.matrix(inDataFrame))) >
        100) {
        kurtosisValue1 <- kurtosis(as.vector(as.matrix(inDataFrame)))
        # Here, the log transformation is
        # performed. In cases where the lowest
        # value is 0, everything is simple. In
        # other cases, a slightly more
        # complicated formula is needed
        if (min(inDataFrame) >= 0) {
            inDataFrame <- log2(inDataFrame + 1)
            logCenterSd[[1]] <- TRUE
        } else {
            # First, the data needs to be reasonably
            # log transformed to not too extreme
            # values, but still without loosing
            # resolution.
            inDataMatrixLog <- log2(apply(
                inDataFrame, 2,
                function(x) x - min(x)
            ) + 1)
            # Then, the extreme negative values will
            # be replaced by 0, as they give rise to
            # artefacts.
            inDataMatrixLog[which(is.nan(inDataMatrixLog))] <- 0
            inDataFrame <- as.data.frame(inDataMatrixLog)
            logCenterSd[[1]] <- TRUE
        }

        kurtosisValue2 <- kurtosis(as.vector(as.matrix(inDataFrame)))
        message(
            "The data was found to be heavily tailed (kurtosis ",
            kurtosisValue1, "). Therefore, it was log2-transformed, leading to
            a new kurtosis value of ", kurtosisValue2, "."
        )
    }

    # Centering and overall scaling is
    # performed
    if (length(center) == 1 && (center == "peak" ||
        (center == "default" && ncol(inDataFrame) <= 100))) {
        inDataFrameScaleList <- dScale(inDataFrame,
                                       scale = FALSE,
                                       center = "peak", returnCenter = TRUE
        )
        inDataFramePreScaled <- inDataFrameScaleList[[1]]
        logCenterSd[[2]] <- inDataFrameScaleList[[2]]
        if(ncol(inDataFrame) <= 100){
            message("As the dataset has less than 100 columns, ",
                    "peak centering is applied.")
        } else {
            message("Peak centering is applied although the data has ",
                    "more than 100 columns")
        }
    } else if (length(center) == 1 && (center == "mean" ||
               (center == "default" && ncol(inDataFrame) > 100))) {
        inDataFrameScaleList <- dScale(inDataFrame,
                                       scale = FALSE,
                                       center = "mean", returnCenter = TRUE
        )
        inDataFramePreScaled <- inDataFrameScaleList[[1]]
        logCenterSd[[2]] <- inDataFrameScaleList[[2]]
        if(ncol(inDataFrame) <= 100){
            message("Mean centering is applied although the data has ",
                    "less than 100 columns")
        } else {
            message("As the dataset has more than 100 columns, ",
                    "mean centering is applied.")
        }
    } else if (length(center) > 1){
        if(length(center) == length(inDataFrame)){
            inDataFramePreScaled <-
                as.data.frame(do.call(
                    "cbind",lapply(seq_len(ncol(inDataFrame)),
                                   function(x){
                                       inDataFrame[,x]-center[x]
                                       })))
            dimnames(inDataFramePreScaled) <- dimnames(inDataFrame)
            logCenterSd[[2]] <- center
            message("Centering based on the provded values")
        } else {
            stop("Mismatch between the number of columns and the length of",
                 " the centering vector")
        }

    } else if (center == FALSE) {
        message("No centering performed")
        inDataFramePreScaled <- inDataFrame
    }

    # Here, all the data is divided by the
    # standard deviation of the full dataset
    if(is.logical(scale) && scale){
        sdInDataFramePreScaled <- sd(as.matrix(inDataFramePreScaled))
        inDataFrameScaled <-
            as.data.frame(inDataFramePreScaled / sdInDataFramePreScaled)
        logCenterSd[[3]] <- sdInDataFramePreScaled
    } else if(is.numeric(scale)){
        inDataFrameScaled <-
            as.data.frame(inDataFramePreScaled / scale)
        logCenterSd[[3]] <- scale
        message("Scaling value taken directly from input")
    }

    return(list(inDataFrameScaled, logCenterSd))
}
