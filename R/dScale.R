#' Scaling of a vector or a dataframe.
#'
#'
#' This is a scaling function with a number of alternatives. This method for scaling takes the shape of the data into somewhat more of a consideration than minMaxScale does, but still gives less influence of outliers than more conventional scalin alternatives, such as unit variance scaling.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @param x A numeric/integer vector or dataframe
#' @param control A numeric/integer vector or dataframe of values that could be used to define the range. If no control data is present, the function defaults to using the indata as control data.
#' @param scale If scaling should be performed. Three possible values: a vector with two values indicating the low and high threshold quantiles for the scaling, TRUE, which equals the vector 'c(0.001, 0.999)', and FALSE.
#' @param robustVarScale If the data should be scaled to its standard deviation within the quantiles defined by the scale values above. If TRUE (the default), the data is unit variance scaled based on the standard deviation of the data within the range defined by scale.
#' @param center If centering should be performed. Alternatives are mean', 'peak' and FALSE. 'peak' results in centering around the highest peak in the data, which is useful in most cytometry situations. 'mean' results in mean centering.
#' @param truncate If truncation of the most extreme values should be performed. Three possible values: TRUE, FALSE, and a vector with two values indicating the low and high threshold quantiles for truncation.
#' @param multiCore If the algorithm should be performed on multiple cores. This increases speed in situations when very large datasets (eg >1 000 000 rows) are scaled. With smaller datasets, it works, but is slow. Defaults to FALSE.
#' @param multiplicationFactor A value that all values will be multiplied with. Useful e.g. if the results preferrably should be returned as percent. Defaults to FALSE.
#' @param returnCenter Boolean. If center=TRUE, should the value at the center be returned?
#' @return A vector or dataframe with the same size but where all values in the vector or column of the dataframe have been internally scaled. In addition, if returnCenter=TRUE, a value, or a vector if x is a matrix or a data frame
#' @examples
#' # Load some data
#' data(testData)
#' 
#' # Retrieve the first column
#' x <- testData[, 2]
#' 
#' # The maximum and minimum values are
#' max(x)
#' min(x)
#' 
#' # Run the function without mean centering and with the quantiles set to 0 and 1.
#' y <- dScale(x, scale = c(0, 1), robustVarScale = FALSE, center = FALSE)
#' 
#' # And the data has been scaled to the range between 0 and 1.
#' max(y)
#' min(y)
#' 
#' # Now run the default function for a dataframe
#' summary(testData[, 2:15])
#' 
#' y_df <- dScale(testData[, 2:15])
#' 
#' # Here, the data has first been truncated to the default percentiles, then scaled
#' # to the standard deviation in the remaining interval and finally the center has been
#' # placed where the highest peak in the data is present. NB! Here, no truncation has
#' # been performed in the scaling, only to obtain the scaling values.
#' summary(y_df)
#' @export dScale
dScale <- function(x, control, scale = TRUE, 
    robustVarScale = TRUE, center = "peak", 
    truncate = FALSE, multiplicationFactor = 1, 
    multiCore = FALSE, returnCenter = FALSE) {
    if (any(is(x) == "matrix")) {
        x <- as.data.frame(x)
    }
    
    if (any(is(x) == "numeric") == FALSE && 
        any(is(x) == "integer") == FALSE && 
        any(is(x) == "data.frame") == FALSE) {
        stop("The data is incorrectly formatted, as it is not a vector, a 
            matrix or a dataframe. Change this and try again.")
    }
    
    if (missing("control")) {
        control <- x
    } else {
        if (any(is(control) == "matrix")) {
            control <- as.data.frame(control)
        }
        if (any(is(x) == "data.frame")) {
            control <- rbind(x, control)
        } else {
            control <- c(x, control)
        }
    }
    
    if (is.logical(scale) == TRUE && scale == 
        TRUE) {
        scale <- c(0.001, 0.999)
    }
    
    if (is.logical(truncate) == TRUE && truncate == 
        TRUE) {
        if (is.logical(scale) == TRUE) {
            truncate <- c(0.001, 0.999)
        } else {
            truncate <- scale
        }
    }
    
    if (any(is(x) == "data.frame") == FALSE) {
        result <- dScaleCoFunction(x, control = control, 
            scale = scale, robustVarScale = robustVarScale, 
            truncate = truncate, center = center, 
            multiplicationFactor = multiplicationFactor, 
            returnCenter = returnCenter)
    }
    if (any(is(x) == "data.frame")) {
        if (multiCore == TRUE) {
            no_cores <- floor(detectCores()*0.875)
            cl <- makeCluster(no_cores, type = "SOCK")
            registerDoSNOW(cl)
            if (returnCenter == FALSE) {
                i <- 1
                result <- as.data.frame(foreach(i = seq_len(ncol(x)), 
                  .inorder = TRUE) %dopar% 
                  dScaleCoFunction(x[, i], 
                    control = control[, i], 
                    scale = scale, robustVarScale = robustVarScale, 
                    truncate = truncate, 
                    center = center, multiplicationFactor = multiplicationFactor))
            } else {
                resultComplex <- foreach(i = seq_len(ncol(x)), 
                  .inorder = TRUE) %dopar% 
                  dScaleCoFunction(x[, i], 
                    control = control[, i], 
                    scale = scale, robustVarScale = robustVarScale, 
                    truncate = truncate, 
                    center = center, multiplicationFactor = multiplicationFactor, 
                    returnCenter = TRUE)
                resultDf <- as.data.frame(do.call("cbind", 
                  lapply(resultComplex, "[[", 
                    1)))
                resultCenters <- unlist(lapply(resultComplex, 
                  "[[", 2))
                result <- list(resultDf, 
                  resultCenters)
            }
            stopCluster(cl)
            colnames(result) <- colnames(x)
        } else {
            if (returnCenter == FALSE) {
                result <- as.data.frame(mapply(dScaleCoFunction, 
                  x, control, MoreArgs = list(scale = scale, 
                    robustVarScale = robustVarScale, 
                    truncate = truncate, 
                    center = center, multiplicationFactor = multiplicationFactor), 
                  SIMPLIFY = FALSE))
            } else {
                resultComplex <- mapply(dScaleCoFunction, 
                  x, control, MoreArgs = list(scale = scale, 
                    robustVarScale = robustVarScale, 
                    truncate = truncate, 
                    center = center, multiplicationFactor = multiplicationFactor, 
                    returnCenter = TRUE), 
                  SIMPLIFY = FALSE)
                resultDf <- as.data.frame(do.call("cbind", 
                  lapply(resultComplex, "[[", 
                    1)))
                resultCenters <- unlist(lapply(resultComplex, 
                  "[[", 2))
                result <- list(resultDf, 
                  resultCenters)
            }
        }
    }
    return(result)
}
