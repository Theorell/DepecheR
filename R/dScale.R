#' Scaling of a vector or a dataframe.
#'
#'
#' This is a scaling function with a number of alternatives. This method for 
#' scaling takes the shape of the data into somewhat more of a consideration 
#' than minMaxScale does, but still gives less influence of outliers than more 
#' conventional scalin alternatives, such as unit variance scaling.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar% %do%
#' @param x A numeric/integer vector or dataframe
#' @param control A numeric/integer vector or dataframe of values that could be
#' used to define the range. If no control data is present, the function 
#' defaults to using the indata as control data.
#' @param scale If scaling should be performed. Three possible values: a vector
#' with two values indicating the low and high threshold quantiles for the 
#' scaling, TRUE, which equals the vector 'c(0.001, 0.999)', and FALSE.
#' @param robustVarScale If the data should be scaled to its standard deviation 
#' within the quantiles defined by the scale values above. If TRUE 
#' (the default), the data is unit variance scaled based on the standard 
#' deviation of the data within the range defined by scale.
#' @param center If centering should be performed. Alternatives are mean', 
#' 'peak' and FALSE. 'peak' results in centering around the highest peak in the
#' data, which is useful in most cytometry situations. 'mean' results in mean 
#' centering.
#' @param truncate If truncation of the most extreme values should be performed.
#' Three possible values: TRUE, FALSE, and a vector with two values indicating 
#' the low and high threshold quantiles for truncation.
#' @param nCores If the function is run in multicore mode, which it will if the
#' dataset is large (nrow*ncol>10^6), this decides the number of cores. The 
#' default is currently 87.5 percent with a cap on 10 cores, as no speed 
#' increase is generally seen above 10 cores for normal computers to date. 
#' @param multiplicationFactor A value that all values will be multiplied with.
#' Useful e.g. if the results preferrably should be returned as percent. 
#' Defaults to FALSE.
#' @param returnCenter Boolean. If center=TRUE, should the value at the center
#' be returned?
#' @return A vector or dataframe with the same size but where all values in the
#' vector or column of the dataframe have been internally scaled. In addition,
#' if returnCenter=TRUE, a value, or a vector if x is a matrix or a data frame.
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
#' # Run the function without mean centering and with the quantiles set to 0 
#' # and 1.
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
#' # Here, the data has first been truncated to the default percentiles, then 
#' # scaled to the standard deviation in the remaining interval and finally the
#' # center has been placed where the highest peak in the data is present. 
#' # NB! Here, no truncation has been performed in the scaling, only to obtain
#' # the scaling values.
#' 
#' summary(y_df)
#' 
#' @export dScale
dScale <- function(x, control, scale = TRUE, robustVarScale = TRUE, 
                   center = "peak", truncate = FALSE, multiplicationFactor = 1, 
                   returnCenter = FALSE, nCores="default") {
    if (is.matrix(x)) {
        x <- as.data.frame(x)
    }
    
    if (is.numeric(x) == FALSE && 
        is.integer(x) == FALSE && 
        is.data.frame(x) == FALSE) {
        stop("The data is incorrectly formatted, as it is not a numeric vector, 
             a matrix or a dataframe. Change this and try again.")
    }
    
    if (missing("control")) {
        control <- x
    } else {
        if (is.matrix(control)) {
            control <- as.data.frame(control)
        }
        if (is.data.frame(x)) {
            control <- rbind(x, control)
        } else {
            control <- c(x, control)
        }
    }
    
    if (is.logical(scale)) {
        if(scale){
            scale <- c(0.001, 0.999) 
        }
    }
    
    if (is.logical(truncate)) {
        if(truncate){
            if (is.logical(scale)) {
                truncate <- c(0.001, 0.999)
            } else {
                truncate <- scale
            } 
        }
    }
    
    if (is.data.frame(x) == FALSE) {
        result <- dScaleCoFunction(x, control = control, scale = scale, 
                                   robustVarScale = robustVarScale, 
                                   truncate = truncate, center = center, 
                                   multiplicationFactor = multiplicationFactor,
                                   returnCenter = returnCenter)
    }
    if (is.data.frame(x)) {
        
        if (ncol(x)*nrow(x)>1000000) {
            if( nCores=="default"){
                nCores <- floor(detectCores()*0.875) 
                if(nCores>10){
                    nCores <- 10
                }
            }
            cl <- makeCluster(nCores, type = "SOCK")
            registerDoSNOW(cl)
            i <- 1
            result <- foreach(i = seq_len(ncol(x)), .inorder = TRUE) %dopar%
                dScaleCoFunction(x[, i], control = control[, i], 
                                 scale = scale, 
                                 robustVarScale = robustVarScale, 
                                 truncate = truncate, 
                                 center = center, 
                                 multiplicationFactor 
                                 = multiplicationFactor, 
                                 returnCenter = returnCenter)
            stopCluster(cl)
        } else {
            result <- foreach(i = seq_len(ncol(x)), .inorder = TRUE) %do%
                dScaleCoFunction(x[, i], control = control[, i], 
                                 scale = scale, 
                                 robustVarScale = robustVarScale, 
                                 truncate = truncate, 
                                 center = center, 
                                 multiplicationFactor 
                                 = multiplicationFactor, 
                                 returnCenter = returnCenter)
        }
        if(returnCenter==FALSE){
            result <- as.data.frame(result)
            colnames(result) <- colnames(x)
        } else {
            resultDf <- as.data.frame(do.call("cbind", 
                                              lapply(result, "[[", 1)))
            colnames(resultDf) <- colnames(x)
            resultCenters <- unlist(lapply(result,  "[[", 2))
            result <- list(resultDf, resultCenters)
        }
        
    }
        
    return(result)
}
