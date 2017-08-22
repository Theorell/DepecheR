#' Scaling to the highest and lowest values in a vector or a dataframe.
#'
#'
#' The simples possible scaling procedure, scaling to the internal highest and lowest values in each column. The function accepts vectors and dataframes.
#' @param x A numeric/integer vector or dataframe.
#' @param control A numeric/integer vector or dataframe of values that could be used to define the range. If no control data is present, the function defaults to using x as control data. NB! The control file needs to be the same class, i.e vector or dataframe, as the x, and in the case of a dataframe it needs to have the same number of columns.
#' @param multiplicationFactor A value that all values will be multiplied with. Useful e.g. if the results preferrably should be returned as percent.
#' @seealso \code{\link{quantileScale}}, @seealso \code{\link{truncateColorScale}}
#' @return A vector or dataframe with the same size but where all values in the vector or in each column of the dataframe have been internally scaled.
#' @examples
#' #Generate a random vector
#' x <- rnorm(1000, 55, 10)
#'
#' #The maximum and minimum values are
#' max(x)
#' min(x)
#'
#' #Run the function
#' y <- minMaxScale(x)
#'
#' #And the data has been scaled to the range of the data. 
#' #In this case, the standard 0-1 range is used, as the 
#' #multiplicationFactor has not been changed from the default.
#' max(y)
#' min(y)
#'
#' #Do the same but with a dataframe
#' x_df <- data.frame(cbind(rnorm(1000, 55, 10), rnorm(1000, 2, 90), rnorm(1000, 430, 200)))
#' summary(x_df)
#'
#' #Run the function
#' y_df <- minMaxScale(x_df, multiplicationFactor=100)
#'
#' #And the data has been rescaled to a percentage of the 
#' #range between the minimal and the maximal values in the 
#' #original dataframe columns.
#' summary(y_df)
#'
#' #Here, dataframes are used, and control data is included.
#' x_df <- data.frame(cbind(rnorm(1000, 55, 10), rnorm(1000, 2, 90), rnorm(1000, 430, 200)))
#' summary(x_df)
#' control_df <- data.frame(cbind(rnorm(1000, 55, 20), rnorm(1000, 10, 200), rnorm(1000, 450, 350)))
#'
#' #Run the function
#' y_df <- quantileScale(x_df, control=control_df, lowQuantile=0.01, 
#' highQuantile=0.99, multiplicationFactor=100)
#'
#' #And the data has been rescaled to a percentage of the range 
#' #between the most extreme values of of the data in the control 
#' #dataframe columns. Note that as the range of the values in the 
#' #control data iscol larger, the data is compressed.
#' summary(y_df)
#' @export minMaxScale
minMaxScale <- function(x, control=FALSE, multiplicationFactor=1){

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
   result <-  multiplicationFactor*(x-min(control))/(max(control)-min(control))
  }
  if(class(x)=="data.frame"){
    result <- multiplicationFactor*((x-apply(control, 2, min, na.rm=TRUE)[col(x)])/(apply(control, 2, max, na.rm=TRUE)[col(x)]-apply(control, 2, min, na.rm=TRUE)[col(x)]))
  }

  return(result)

}
