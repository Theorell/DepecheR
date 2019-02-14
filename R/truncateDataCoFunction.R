truncateDataCoFunction <- function(x, control, lowQuantile = 1e-04, 
                                   highQuantile = 0.9999) {
    high <- quantile(control, highQuantile)
    low <- quantile(control, lowQuantile)
    
    x[x > high] <- high
    x[x < low] <- low
    return(x)
}
