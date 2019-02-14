dScaleCoFunction <- function(x, control, scale, robustVarScale, truncate, 
                             center, multiplicationFactor, 
                             returnCenter = FALSE) {
    if (is.logical(scale) == TRUE && scale == FALSE) {
        if (is.logical(truncate) == TRUE) {
            responseVector <- multiplicationFactor * x
        }
        if (length(truncate) == 2) {
            xTruncReal <- truncateData(x, 
                lowQuantile = truncate[1], 
                highQuantile = truncate[2])
            responseVector <- multiplicationFactor * xTruncReal
        }
    }
    
    if (length(scale) == 2) {
        # Define quantile
        bottom <- quantile(control, probs = scale[1], se = FALSE, na.rm = TRUE)
        top <- quantile(control, probs = scale[2], se = FALSE, na.rm = TRUE)
        
        if (robustVarScale == FALSE) {
            if (is.logical(truncate) == TRUE) {
                responseVector <- multiplicationFactor * 
                    ((x - bottom)/(top - bottom))
            }
            if (length(truncate) == 2) {
                xTruncReal <- truncateData(x, 
                                           lowQuantile = truncate[1], 
                                           highQuantile = truncate[2])
                responseVector <- multiplicationFactor * 
                    ((xTruncReal - bottom)/(top - bottom))
            }
        }
        
        if (robustVarScale == TRUE) {
            # First truncate the data to the
            # quantiles defined by the quantiles
            xTruncated <- truncateData(x, lowQuantile = scale[1], 
                                       highQuantile = scale[2])
            
            sdxTruncated <- sd(xTruncated)
            
            # Now the data is scaled
            if (is.logical(truncate) == TRUE) {
                responseVector <- multiplicationFactor * x/sdxTruncated
            }
            if (length(truncate) == 2) {
                xTruncReal <- truncateData(x, lowQuantile = truncate[1], 
                                           highQuantile = truncate[2])
                responseVector <- multiplicationFactor * xTruncReal/sdxTruncated
            }
        }
    }
    
    if (center == "mean") {
        meanValue <- mean(responseVector)
        responseVector <- responseVector - meanValue
        if (returnCenter == TRUE) {
            responseList <- list(responseVector, meanValue)
        }
    }
    
    if (center == "peak") {
        # The peak of the data is defined
        if (length(x) < 500) {
            nBreaks <- 10
        } else {
            nBreaks <- length(x)/50
        }
        histdata <- hist(responseVector, 
            breaks = nBreaks, plot = FALSE)
        zeroPosition <- histdata$mids[match(max(histdata$counts), 
            histdata$counts)]
        
        # And the position for this this peak is
        # subtracted from all points
        responseVector <- responseVector - 
            zeroPosition
        if (returnCenter == TRUE) {
            responseList <- list(responseVector, 
                zeroPosition)
        }
    }
    
    
    if (returnCenter == FALSE) {
        return(responseVector)
    } else {
        return(responseList)
    }
}
