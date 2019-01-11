truncateData <- function(x, control, lowQuantile = 0.001, 
    highQuantile = 0.999) {
    if (any(is(x) == "numeric") == FALSE && 
        any(is(x) == "integer") == FALSE && 
        any(is(x) == "data.frame") == FALSE) {
        stop("Data needs to be either a numeric/integer vector or a dataframe. Change to a suitable object and try again.")
    }
    
    if (missing("control")) {
        control <- x
    }
    
    if (identical(colnames(x), colnames(control)) == 
        FALSE) {
        print("Warning. Column names of the x data and the control data are mismatched or are ordered differently, which may affect the result. Consider correcting this.")
    }
    
    
    if (any(is(x) == "data.frame") == FALSE) {
        result <- truncateDataCoFunction(x, 
            control = control, lowQuantile = lowQuantile, 
            highQuantile = highQuantile)
    }
    if (any(is(x) == "data.frame")) {
        result <- as.data.frame(mapply(truncateDataCoFunction, 
            x, control, MoreArgs = list(lowQuantile = lowQuantile, 
                highQuantile = highQuantile), 
            SIMPLIFY = FALSE))
    }
    
    return(result)
}
