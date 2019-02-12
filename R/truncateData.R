truncateData <- function(x, control, lowQuantile = 0.001, 
    highQuantile = 0.999) {
    if (is.numeric(x) == FALSE && 
        is.integer(x) == FALSE && 
        is.data.frame(x) == FALSE) {
        stop("Data needs to be either a numeric/integer vector or a dataframe. 
             Change to a suitable object and try again.")
    }
    
    if (missing("control")) {
        control <- x
    }
    
    if (identical(colnames(x), colnames(control)) == 
        FALSE) {
        warning("Column names of the x data and the control data are mismatched
                or are ordered differently, which may be wrong. 
                Consider correcting this.")
    }
    
    
    if (is.data.frame(x) == FALSE) {
        result <- truncateDataCoFunction(x, 
            control = control, lowQuantile = lowQuantile, 
            highQuantile = highQuantile)
    }
    if (is.data.frame(x)) {
        result <- as.data.frame(mapply(truncateDataCoFunction, 
            x, control, MoreArgs = list(lowQuantile = lowQuantile, 
                highQuantile = highQuantile), 
            SIMPLIFY = FALSE))
    }
    
    return(result)
}
