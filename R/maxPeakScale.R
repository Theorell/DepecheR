#' Scaling and centering around the peak of the data
#'
#'
#' This function scales the data to the standard deviation within a default or user-defined quantile range. It the centers the data aroun the highest peak within it. 
#' @importFrom Hmisc hdquantile
#' @param x A numeric/integer vector or dataframe
#' @param control A numeric/integer vector or dataframe of values that could be used to define the range. If no control data is present, the function defaults to using the indata as control data.
#' @param lowQuantile The lower border below which the values are treated as outliers and will be outside of the defined scaling range (0-1*multiplicationFactor).
#' @param highQuantile The higher border above which the values are treated as outliers and will be outside of the defined scaling range (0-1*multiplicationFactor).
#' @seealso \code{\link{minMaxScale}} \code{\link{quantileScale}}
#' @return A vector or dataframe with the same size but where all values in the vector or column of the dataframe have been internally scaled.
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData()
#'
#' #The maximum and minimum values are
#' summary(x)
#'
#' #Run the function
#' y <- maxPeakScale(x[2:ncol(x)], lowQuantile=0.01, highQuantile=0.99)
#'
#' #And the data has been rescaled to the standard deviation of range between the first and the 99:th percentile 
#' #of the data in the original dataframe columns. The data has also been centered around the peak of each column.
#' summary(y)
#' 
#' @export maxPeakScale
maxPeakScale <- function(x, control, lowQuantile=0.001, highQuantile=0.999){

  if(class(x)!="numeric" && class(x)!="integer" && class(x)!="data.frame"){
    stop("Data needs to be either a numeric/integer vector or a dataframe. Change the class and try again.")
  }
  
  if(missing("control")){
    control <- x
  }
  
  if(identical(colnames(x), colnames(control))==FALSE){
    print("Warning. Column names of the x data and the control data are mismatched or are ordered differently, which may affect the result. Consider correcting this.")
  }

  if(class(x)!="data.frame"){
    result <- maxPeakScaleCoFunction(x, control=control, lowQuantile=lowQuantile, highQuantile=highQuantile)
  }
  if(class(x)=="data.frame"){
    result <- as.data.frame(mapply(maxPeakScaleCoFunction, x, control, MoreArgs=list(lowQuantile=lowQuantile, highQuantile=highQuantile), SIMPLIFY = FALSE))
  }
  return(result)  	

}

maxPeakScaleCoFunction <- function(x, control, lowQuantile, highQuantile){
 
  #First truncate the data to the quantiles defined by the quantiles
  xTruncated <- truncateData(x, lowQuantile=lowQuantile, highQuantile=highQuantile)
  
  #Now the data is scaled
  xScaled <- scale(xTruncated, center=FALSE)
  
  #And finally the peak of the data is defined
  histdata <- hist(xScaled, breaks=200, plot=FALSE)
  zeroPosition <- histdata$mids[match(max(histdata$counts), histdata$counts)]
  
  #And this peak is subtracted from all points
  xCentered <- xScaled-zeroPosition
  
}

