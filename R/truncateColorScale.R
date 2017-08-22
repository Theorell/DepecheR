#' Truncate the extreme ends of numeric data
#'
#'
#' This function reduces the most extremely low and high observations to their less extreme border. Currently, the standard quantile cutoff is a hundredth of a percent. It is done to decrease these observations influence on color scaling for vizualisations. The function accepts vectors and dataframes.
#' @importFrom Hmisc hdquantile
#' @param x A numeric/integer vector or dataframe.
#' @param control A numeric/integer vector or dataframe of values that could be used to define the range. If no control data is present, the function defaults to using x as control data. NB! The control file needs to be the same class, i.e vector or dataframe, as the x, and in the case of a dataframe it needs to have the same number of columns.
#' @param lowQuantile The lower border below which the values are treated as outliers and will be outside of the defined quantile range.
#' @param highQuantile The higher border above which the values are treated as outliers and will be outside of the defined quantile range.
#' @seealso \code{\link{minMaxScale}}
#' @return A vector or dataframe with the same size but where the extreme values have been substituted with the less extreme values.
#' @examples
#' #Generate a random vector
#' x <- rnorm(1000, 55, 10)
#'
#' #The maximum and minimum values are
#' max(x)
#' min(x)
#'
#' #Here, the function is used to truncate the two most extreme percent of the observations.
#' y <- truncateColorScale(x, lowQuantile=0.01, highQuantile=0.99)
#'
#' #The maximum and minimum values are less extreme:
#' max(y)
#' min(y)
#'
#' #Do the same but with a dataframe
#' x_df <- data.frame(cbind(rnorm(1000, 55, 10), rnorm(1000, 2, 90), rnorm(1000, 430, 200)))
#' summary(x_df)
#'
#' #Run the function
#' y_df <- truncateColorScale(x_df, lowQuantile=0.01, highQuantile=0.99)
#'
#' #And the most extreme values have been reduced
#' summary(y_df)
#'
#' @export truncateColorScale

truncateColorScale <- function(x, control, lowQuantile=0.0001, highQuantile=0.9999){

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
    result <- truncateColorScaleCoFunction(x, control=control, lowQuantile=lowQuantile, highQuantile=highQuantile)
  }
  if(class(x)=="data.frame"){
    result <- as.data.frame(mapply(truncateColorScaleCoFunction, x, control, MoreArgs=list(lowQuantile=lowQuantile, highQuantile=highQuantile), SIMPLIFY = FALSE))
  }

  return(result)
}


truncateColorScaleCoFunction <- function(x, control, lowQuantile=0.0001, highQuantile=0.9999){
  high <- hdquantile(control, highQuantile)
  low <- hdquantile(control, lowQuantile)

  x[x > high] <- high
  x[x < low] <- low
  return(x)
}
