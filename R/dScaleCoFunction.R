dScaleCoFunction <- function(x, control, scale, robustVarScale, truncate, center, multiplicationFactor){

  if(is.logical(scale)==TRUE && scale==FALSE){
    if(is.logical(truncate)==TRUE){
      responseVector <- multiplicationFactor*x
    }
    if(length(truncate)==2){
      xTruncReal <- truncateData(x, lowQuantile=truncate[1], highQuantile=truncate[2])
      responseVector <- multiplicationFactor*xTruncReal
    }
  } 
  
  if(length(scale)==2){
    #Define quantile
    bottom <- quantile(control, probs = scale[1], se=FALSE, na.rm=TRUE)
    top <- quantile(control, probs = scale[2], se=FALSE, na.rm=TRUE)

    if(robustVarScale==FALSE){
      if(is.logical(truncate)==TRUE){
        responseVector <- multiplicationFactor*((x-bottom)/(top-bottom))
      } 
      if(length(truncate)==2){
        xTruncReal <- truncateData(x, lowQuantile=truncate[1], highQuantile=truncate[2])
        responseVector <- multiplicationFactor*((xTruncReal-bottom)/(top-bottom))
      }
    }
    
    if(robustVarScale==TRUE){
      #First truncate the data to the quantiles defined by the quantiles
      xTruncated <- truncateData(x, lowQuantile=scale[1], highQuantile=scale[2])
      
      sdxTruncated <- sd(xTruncated)
      
      #Now the data is scaled
      if(is.logical(truncate)==TRUE){
          responseVector <- multiplicationFactor*x/sdxTruncated
      } 
      if(length(truncate)==2){
        xTruncReal <- truncateData(x, lowQuantile=truncate[1], highQuantile=truncate[2])
        responseVector <- multiplicationFactor*xTruncReal/sdxTruncated
      }
    }
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
