# Scaling of a vector or a dataframe.
#
#
# This is a scaling function with a number of alternatives. This method for scaling takes the shape of the data into somewhat more of a consideration than minMaxScale does, but still gives less influence of outliers than more conventional scalin alternatives, such as unit variance scaling.
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW 
#' @importFrom foreach foreach %dopar%
# @param x A numeric/integer vector or dataframe
# @param control A numeric/integer vector or dataframe of values that could be used to define the range. If no control data is present, the function defaults to using the indata as control data.
# @param scale If scaling should be performed. Three possible values: a vector with two values indicating the low and high threshold quantiles for the scaling, TRUE, which equals the vector "c(0.001, 0.999)", and FALSE.
# @param robustVarScale If the data should be scaled to its standard deviation within the quantiles defined by the scale values above. If TRUE (the default), the data is unit variance scaled based on the standard deviation of the data within the range defined by scale.
# @param center If centering should be performed. Alternatives are mean", "peak" and FALSE. "peak" results in centering around the highest peak in the data, which is useful in most cytometry situations. "mean" results in mean centering. 
# @param truncate If truncation of the most extreme values should be performed. Three possible values: TRUE, FALSE, and a vector with two values indicating the low and high threshold quantiles for truncation. 
# @param multiCore If the algorithm should be performed on multiple cores. This increases speed in situations when very large datasets (eg >1 000 000 rows) are scaled. With smaller datasets, it works, but is slow. Defaults to FALSE.
# @param multiplicationFactor A value that all values will be multiplied with. Useful e.g. if the results preferrably should be returned as percent. Defaults to FALSE.
# @return A vector or dataframe with the same size but where all values in the vector or column of the dataframe have been internally scaled.
# @examples
# #Generate a default size dataframe with bimodally distributed data
# x <- generateBimodalData()
#
# #Retrieve the first column
# x2 <- x[,2]
# 
# #The maximum and minimum values are
# max(x2)
# min(x2)
#
# #Run the function without mean centering and with the quantiles set to 0 and 1.
# y2 <- dScale(x2, scale=c(0,1), robustVarScale=FALSE, center=FALSE)
#
# #And the data has been scaled to the range between 0 and 1.
# max(y2)
# min(y2)
#
# #Now run the default function for a dataframe
# summary(x[,2:ncol(x)])
#
# y_df <- dScale(x[,2:ncol(x)])
#
# #Here, the data has first been truncated to the default percentiles, then scaled 
# #to the standard deviation in the remaining interval and finally the center has been
# #placed where the highest peak in the data is present. NB! Here, no truncation has 
# #been performed in the scaling, only to obtain the scaling values.
# summary(y_df)
dScale <- function(x, control, scale=TRUE, robustVarScale=TRUE, center="peak", truncate=FALSE, multiplicationFactor=1, multiCore=FALSE){

  if(class(x)!="numeric" && class(x)!="integer" && class(x)!="data.frame"){
    stop("The data is incorrectly formatted, as it is not a vector, a matrix or a dataframe. Change this and try again.")
  }

  if(missing("control")){
    control <- x
  } else {
    control <- rbind(x, control)
  }

  if(is.logical(scale)==TRUE && scale==TRUE){
    scale <- c(0.001, 0.999)
  }
  
  if(is.logical(truncate)==TRUE && truncate==TRUE){
    if(is.logical(scale)==TRUE){
      truncate <- c(0.001, 0.999)
    } else {
      truncate <- scale
    }
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