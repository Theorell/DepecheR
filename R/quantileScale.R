#' Scaling of a vector or a dataframe.
#'
#'
#' This is a scaling function with a number of alternatives. This method for scaling takes the shape of the data into somewhat more of a consideration than minMaxScale does, but still gives less influence of outliers than more conventional scalin alternatives, such as unit variance scaling.
#' @importFrom Hmisc hdquantile
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
#' @param x A numeric/integer vector or dataframe
#' @param control A numeric/integer vector or dataframe of values that could be used to define the range. If no control data is present, the function defaults to using the indata as control data.
#' @param lowQuantile The lower border below which the values are treated as outliers and will be outside of the defined scaling range (0-1*multiplicationFactor).
#' @param highQuantile The higher border above which the values are treated as outliers and will be outside of the defined scaling range (0-1*multiplicationFactor).
#' @param robustVarScale If the data should be scaled to its standard deviation within the quantiles defined by the high and low quantile below. If TRUE (the default), the data is first truncated with truncateData to the quantiles and then the standard deviation scaling is performed.
#' @param center If centering should be performed. Alternatives are "mean", "peak" and FALSE. "peak" results in centering around the highest peak in the data, which is useful in most cytometry situations, the reason it is default. "mean" results in mean centering. 
#' @param truncate If truncation of the most extreme values should be performed. Three possible values, where default is FALSE:
#' @param multiCore If the algorithm should be performed on multiple cores. This increases speed in situations when very large datasets (eg >1 000 000 rows) are scaled. With smaller datasets, it works, but is slow. Defaults to FALSE.
#' #' \describe{
#'     \item{TRUE}{The same quantiles are used as for the low and high quantiles.}
#'     \item{FALSE}{No truncation.}
#'     \item{A vector with two values, eg c(0.01, 0.99)}{Data outside of these quantiles will be truncated.}
#' }   
#' @param multiplicationFactor A value that all values will be multiplied with. Useful e.g. if the results preferrably should be returned as percent.
#' @seealso \code{\link{truncateData}} 
#' @return A vector or dataframe with the same size but where all values in the vector or column of the dataframe have been internally scaled.
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData()
#'
#' #Retrieve the first column
#' x2 <- x[,2]
#' 
#' #The maximum and minimum values are
#' max(x2)
#' min(x2)
#'
#' #Run the function without mean centering and with the quantiles set to 0 and 1.
#' y2 <- quantileScale(x2, robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE)
#'
#' #And the data has been scaled to the range between 0 and 1.
#' max(y2)
#' min(y2)
#'
#' #Now run the default function for a dataframe
#' summary(x[,2:ncol(x)])
#'
#' y_df <- quantileScale(x[,2:ncol(x)])
#'
#' #Here, the data has first been truncated to the default percentiles, then scaled 
#' #to the standard deviation in the remaining interval and finally the center has been
#' #placed where the highest peak in the data is present. NB! Here, no truncation has been performed in the scaling, only to obtain the scaling values.
#' summary(y_df)
#' @export quantileScale 
quantileScale <- function(x, control, lowQuantile=0.001, highQuantile=0.999, robustVarScale=TRUE, center="peak", truncate=FALSE, multiplicationFactor=1, multiCore=FALSE){

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
    result <- quantileScaleCoFunction(x, control=control, robustVarScale=robustVarScale, lowQuantile=lowQuantile, highQuantile=highQuantile, truncate=truncate, center=center, multiplicationFactor=multiplicationFactor)
  }
  if(class(x)=="data.frame"){
    if(multiCore==TRUE){
      no_cores <- detectCores() - 1
      cl = parallel::makeCluster(no_cores, type = "SOCK")
      registerDoSNOW(cl)
      result <- foreach(i=1:ncol(x), .inorder=TRUE) %dopar% quantileScaleCoFunction(x[,i], control=control[,i], robustVarScale=robustVarScale, lowQuantile=lowQuantile, highQuantile=highQuantile, truncate=truncate, center=center, multiplicationFactor=multiplicationFactor)
      parallel::stopCluster(cl)
    } else {
      result <- as.data.frame(mapply(quantileScaleCoFunction, x, control, MoreArgs=list(robustVarScale=robustVarScale, lowQuantile=lowQuantile, highQuantile=highQuantile, truncate=truncate, center=center, multiplicationFactor=multiplicationFactor), SIMPLIFY = FALSE))
    }
      }
 return(result)
}

quantileScaleCoFunction <- function(x, control, lowQuantile, highQuantile, robustVarScale, truncate, center, multiplicationFactor){

    #Define quartiles using Harrell-Davis Distribution-Free Quantile Estimator for all values in one column
    top <- Hmisc::hdquantile(control, probs = highQuantile, se=FALSE, na.rm=TRUE)
    bottom <- Hmisc::hdquantile(control, probs = lowQuantile, se=FALSE, na.rm=TRUE)

  if(robustVarScale==FALSE){
  responseVector <- multiplicationFactor*((x-bottom)/(top-bottom))
  }
  
  if(robustVarScale==TRUE){
    #First truncate the data to the quantiles defined by the quantiles
    xTruncated <- truncateData(x, lowQuantile=lowQuantile, highQuantile=highQuantile)
    
    sdxTruncated <- sd(xTruncated)
    
    #Now the data is scaled
    responseVector <- multiplicationFactor*x/sdxTruncated
    
  }
  
  if(truncate==TRUE){
    responseVector <- truncateData(responseVector, lowQuantile=lowQuantile, highQuantile=highQuantile)
  }
  
  if(length(truncate)==2){
      responseVector <- truncateData(responseVector, lowQuantile=truncate[1], highQuantile=truncate[2])
  }
    
    
  if(center=="mean"){
    responseVector <- responseVector-mean(responseVector)
  }
  
  if(center=="peak"){
    #The peak of the data is defined
    histdata <- hist(responseVector, breaks=200, plot=FALSE)
    zeroPosition <- histdata$mids[match(max(histdata$counts), histdata$counts)]
    
    #And the position for this this peak is subtracted from all points
    responseVector <- responseVector-zeroPosition
  }
  

  return(responseVector)

}
