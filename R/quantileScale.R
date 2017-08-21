#' Scaling of a vector or a dataframe to a certain quantile.
#'
#'
#' This is a robust alternative to minMax scaling. This method for scaling takes the shape of the data into somewhat more of a consideration than minMaxScale does, but still gives less influence of outliers than more conventional scalin alternatives, such as unit variance scaling.
#' @importFrom Hmisc hdquantile
#' @param x A numeric vector
#' @param control A numeric vector or dataframe of values that could be used to define the range. If no control data is present, the function defaults to using the indata as control data.
#' @param lowQuantile The lower border below which the values are treated as outliers and will be outside of the defined scaling range (0-1*multiplicationFactor).
#' @param highQuantile The higher border above which the values are treated as outliers and will be outside of the defined scaling range (0-1*multiplicationFactor).
#' @param center If mean centering should be performed or not. Defaults to FALSE.
#' @param multiplicationFactor A value that all values will be multiplied with. Useful e.g. if the results preferrably should be returned as percent.
#' @seealso \code{\link{minMaxScale}}
#' @return A vector or dataframe with the same size but where all values in the vector or column of the dataframe have been internally scaled.
#' @examples
#' #Generate a random vector
#' x <- rnorm(1000, 55, 10)
#'
#' #The maximum and minimum values are
#' max(x)
#' min(x)
#'
#' #Run the function
#' y <- quantileScale(x, lowQuantile=0.01, highQuantile=0.99)
#'
#' #And the data has been scaled. Note that the most extreme values are outside of the range 0-1.
#' max(y)
#' min(y)
#'
#' #Do the same but with a dataframe.
#' x_df <- data.frame(cbind(rnorm(1000, 55, 10), rnorm(1000, 2, 90), rnorm(1000, 430, 200)))
#' summary(x_df)
#'
#' #Run the function
#' y_df <- quantileScale(x_df, lowQuantile=0.01, highQuantile=0.99, center=TRUE, multiplicationFactor=100)
#'
#' #And the data has been rescaled to the range between the first and the 99:th percentile of the data in the original dataframe columns. The data has also been centered.
#' summary(y_df)
#'
#' #Here, dataframes are used, and control data is included.
#' x_df <- data.frame(cbind(rnorm(1000, 55, 10), rnorm(1000, 2, 90), rnorm(1000, 430, 200)))
#' summary(x_df)
#' control_df <- data.frame(cbind(rnorm(1000, 55, 15), rnorm(1000, 10, 200), rnorm(1000, 450, 350)))
#'
#' #Run the function
#' y_df <- quantileScale(x_df, control=control_df, lowQuantile=0.01, highQuantile=0.99, center=TRUE, multiplicationFactor=100)
#'
#' #And the data has been rescaled to the range between the first and the 99:th percentile of the data in the control dataframe columns. The data has also been centered. Note that as the range of the values in the control data is larger, the data is compressed.
#' summary(y_df)
#' @export quantileScale
quantileScale <- function(x, control, lowQuantile=0.001, highQuantile=0.999, center=FALSE, multiplicationFactor=1){

  if(class(x)!="numeric" && class(x)!="data.frame"){
    stop("Data needs to be either a numeric vector or a dataframe. Change the class and try again.")
  }

  if(missing("control")){
    control <- x
  }

  if(identical(colnames(x), colnames(control))==FALSE){
    print("Warning. Column names of the x data and the control data are mismatched or are ordered differently, which may affect the result. Consider correcting this.")
  }



    if(class(x)=="numeric"){
    result <- quantileScaleCoFunction(x, control=control, lowQuantile=lowQuantile, highQuantile=highQuantile, center=center, multiplicationFactor=multiplicationFactor)
  }
  if(class(x)=="data.frame"){
    result <- as.data.frame(mapply(quantileScaleCoFunction, x, control, MoreArgs=list(lowQuantile=lowQuantile, highQuantile=highQuantile, center=center, multiplicationFactor=multiplicationFactor), SIMPLIFY = FALSE))
  }
 return(result)
}

quantileScaleCoFunction <- function(x, control, lowQuantile, highQuantile, center, multiplicationFactor){

    #Define quartiles using Harrell-Davis Distribution-Free Quantile Estimator for all values in one column
    top <- hdquantile(control, probs = highQuantile, se=FALSE, na.rm=TRUE)
    bottom <- hdquantile(control, probs = lowQuantile, se=FALSE, na.rm=TRUE)

  responseVector <- multiplicationFactor*((x-bottom)/(top-bottom))

  if(center==TRUE){
    responseVector <- responseVector-mean(responseVector)
  }

  return(responseVector)

}
