truncateDataCoFunction <- function(x, control, lowQuantile = 0.0001, highQuantile = 0.9999) {
  high <- quantile(control, highQuantile)
  low <- quantile(control, lowQuantile)

  x[x > high] <- high
  x[x < low] <- low
  return(x)
}
