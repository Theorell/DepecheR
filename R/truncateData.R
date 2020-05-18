# This fucntion is used internally in dScaleCoFunction.
# Ther purpose with the truncation is to decrease the influence of extreme
# outliers on variance calculations and on visualizations.
# For information on the different parameters, see dScale.
truncateData <- function(x, control, lowQuantile = 1e-04,
                         highQuantile = 0.9999) {
    high <- quantile(control, highQuantile)
    low <- quantile(control, lowQuantile)

    x[x > high] <- high
    x[x < low] <- low
    return(x)
}
