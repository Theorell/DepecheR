# Function to generate synthetic flow cytometry-like data
#
#
# This function creates a dataframe with any chose size with bimodally distributed data in each column. The probability of the sizes of the respective populations in the bimodal distributions is stochastic.
#' @importFrom dplyr bind_rows
# @param samplings Number of individual probability samplings with identical standard deviation samplings, ie a simple simulated parallel to number of subjects.
# @param dataCols Number of columns in the generated dataframe. Default is 5.
# @param observations Number of rows per donor in the generated dataframe. Default is 10000.
# @return A dataframe with the specified numbers of columns and rows.
# @examples
# #Generate a default size dataframe with bimodally distributed data
# x <- generateBimodalData()
#
# #Plot the first two variables in this dataframe
# plot(x[,2], x[,3])
#
# @export generateBimodalData
generateBimodalData <- function(centers, prop =0.3, dataCols=5, observations=10000){
  if(missing(centers)){
    centers<-rbind(runif(dataCols)+50,runif(dataCols)-50)
  } else {
    dataCols<-ncol(centers)
  }
  stopifnot(nrow(centers)==2)
  
  ids <- c(rep(1,floor(prop*observations)), rep(2,observations-floor(prop*observations)))
  rands <-  matrix( rnorm(observations*dataCols,mean=0,sd=1), observations, dataCols) 
  pop1 <- rands[1:floor(prop*observations),] + rep(centers[1,], each = floor(prop*observations))

  pop2 <- tail(rands,n=observations-floor(prop*observations)) + rep(centers[2,], each = observations-floor(prop*observations))
  
  samples <- rbind(pop1,pop2)
  result <- list(samples,ids)
  names(result)<- c('samples','ids')
  
  return(result)
#   
#   stdevs1 <- sample(c(0.5, 0.6, 0.7, 0.8, 0.9, 1), size=dataCols, replace=TRUE)
#   stdevs2 <- sample(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), size=dataCols, replace=TRUE)
#   stdevs <- cbind(stdevs1, stdevs2)
# 
#     sampleList <- list()
#     for(i in 1:samplings){
#     probabilities1 <- sample(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), size=dataCols, replace=TRUE)
#     probabilities <- cbind(probabilities1, 1-probabilities1)
#     result <- data.frame(mapply(generateBimodalDataCoFunction, data.frame(t(probabilities)), data.frame(t(stdevs)), MoreArgs=list(observations=observations)))
#     Sampling <- rep(i, times=observations)
#     result <- cbind(Sampling, result)
#     sampleList[[i]] <- result
#   }
# resultDf <- bind_rows(sampleList)
#   return(resultDf)
}

