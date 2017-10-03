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
#' @param scale If scaling should be performed. Two possible values: FALSE or a vector with two values indicating the low and high threshold quantiles for the scaling. "c(0.001, 0.999)" is default
#' @param robustVarScale If the data should be scaled to its standard deviation within the quantiles defined by the scale values above. If TRUE (the default), the data is unit variance scaled based on the standard deviation of the data within the range defined by scale.
#' @param center If centering should be performed. Alternatives are "mean", "peak" and FALSE. "peak" results in centering around the highest peak in the data, which is useful in most cytometry situations, the reason it is default. "mean" results in mean centering. 
#' @param truncate If truncation of the most extreme values should be performed. Three possible values: TRUE, FALSE, and a vector with two values indicating the low and high threshold quantiles for truncation. 
#' @param multiCore If the algorithm should be performed on multiple cores. This increases speed in situations when very large datasets (eg >1 000 000 rows) are scaled. With smaller datasets, it works, but is slow. Defaults to FALSE.
#' @param multiplicationFactor A value that all values will be multiplied with. Useful e.g. if the results preferrably should be returned as percent. Defaults to FALSE.
#' @seealso \code{\link{truncateData}} 
#' @return A vector or dataframe with the same size but where all values in the vector or column of the dataframe have been internally scaled.
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateBimodalData()
#'
#' #Retrieve the first column
#' x2 <- x[,2]
#' 
#' #The maximum and minimum values are
#' max(x2)
#' min(x2)
#'
#' #Run the function without mean centering and with the quantiles set to 0 and 1.
#' y2 <- dScale(x2, scale=c(0,1), robustVarScale=FALSE, center=FALSE)
#'
#' #And the data has been scaled to the range between 0 and 1.
#' max(y2)
#' min(y2)
#'
#' #Now run the default function for a dataframe
#' summary(x[,2:ncol(x)])
#'
#' y_df <- dScale(x[,2:ncol(x)])
#'
#' #Here, the data has first been truncated to the default percentiles, then scaled 
#' #to the standard deviation in the remaining interval and finally the center has been
#' #placed where the highest peak in the data is present. NB! Here, no truncation has been performed in the scaling, only to obtain the scaling values.
#' summary(y_df)
#' @export dScale 
dScale <- function(x, control, scale=c(0.001, 0.999), robustVarScale=TRUE, center="peak", truncate=FALSE, multiplicationFactor=1, multiCore=FALSE){

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
    result <- dScaleCoFunction(x, control=control, scale=scale, robustVarScale=robustVarScale, truncate=truncate, center=center, multiplicationFactor=multiplicationFactor)
  }
  if(class(x)=="data.frame"){
    if(multiCore==TRUE){
      no_cores <- detectCores() - 1
      cl = parallel::makeCluster(no_cores, type = "SOCK")
      registerDoSNOW(cl)
      result <- as.data.frame(foreach(i=1:ncol(x), .inorder=TRUE) %dopar% dScaleCoFunction(x[,i], control=control[,i], scale=scale, robustVarScale=robustVarScale, truncate=truncate, center=center, multiplicationFactor=multiplicationFactor))
      parallel::stopCluster(cl)
      colnames(result) <- colnames(x)
    } else {
      result <- as.data.frame(mapply(dScaleCoFunction, x, control, MoreArgs=list(scale=scale, robustVarScale=robustVarScale, truncate=truncate, center=center, multiplicationFactor=multiplicationFactor), SIMPLIFY = FALSE))
    }
      }
 return(result)
}

dScaleCoFunction <- function(x, control, scale, robustVarScale, truncate, center, multiplicationFactor){



  if(length(scale)==2){
    #Define quantiles using Harrell-Davis Distribution-Free Quantile Estimator for all values in one column
    bottom <- Hmisc::hdquantile(control, probs = scale[1], se=FALSE, na.rm=TRUE)
    top <- Hmisc::hdquantile(control, probs = scale[2], se=FALSE, na.rm=TRUE)

    
    if(robustVarScale==FALSE){
      responseVector <- multiplicationFactor*((x-bottom)/(top-bottom))
    }
    
    if(robustVarScale==TRUE){
      #First truncate the data to the quantiles defined by the quantiles
      xTruncated <- truncateData(x, lowQuantile=scale[1], highQuantile=scale[2])
      
      sdxTruncated <- sd(xTruncated)
      
      #Now the data is scaled
      responseVector <- multiplicationFactor*x/sdxTruncated
      
    }
  }

  
  if(truncate==TRUE){
    responseVector <- truncateData(responseVector, lowQuantile=scale[1], highQuantile=scale[2])
  }

  if(length(truncate)==2){
      responseVector <- truncateData(responseVector, lowQuantile=truncate[1], highQuantile=truncate[2])
  }
    
    
  if(center=="mean"){
    responseVector <- responseVector-mean(responseVector)
  }
  
  if(center=="peak"){
    #The peak of the data is defined
    histdata <- hist(responseVector, breaks=length(x)/50, plot=FALSE)
    zeroPosition <- histdata$mids[match(max(histdata$counts), histdata$counts)]
    
    #And the position for this this peak is subtracted from all points
    responseVector <- responseVector-zeroPosition
  }
  

  return(responseVector)

}
