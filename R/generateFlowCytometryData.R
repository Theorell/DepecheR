#' Function to generate synthetic flow cytometry-like data
#'
#'
#' This function creates a dataframe with any chose size with bimodally distributed data in each column. The probability of the sizes of the respective populations in the bimodal distributions is stochastic.
#' @importFrom dplyr bind_rows
#' @param samplings Number of individual probability samplings with identical standard deviation samplings, ie a simple simulated parallel to number of subjects.
#' @param ncols Number of columns in the generated dataframe. Default is 5.
#' @param observations Number of rows per donor in the generated dataframe. Default is 10000.
#' @return A dataframe with the specified numbers of columns and rows.
#' @examples
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateFlowCytometryData()
#'
#' #Plot the first two variables in this dataframe
#' plot(x[,2], x[,3])
#'
#' @export generateFlowCytometryData
generateFlowCytometryData <- function(samplings=1, ncols=5, observations=10000){

  stdevs1 <- sample(c(0.5, 0.6, 0.7, 0.8, 0.9, 1), size=ncols, replace=TRUE)
  stdevs2 <- sample(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), size=ncols, replace=TRUE)
  stdevs <- cbind(stdevs1, stdevs2)

    sampleList <- list()
    for(i in 1:samplings){
    probabilities1 <- sample(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), size=ncols, replace=TRUE)
    probabilities <- cbind(probabilities1, 1-probabilities1)
    result <- data.frame(mapply(generateFlowCytometryDataCoFunction, data.frame(t(probabilities)), data.frame(t(stdevs)), MoreArgs=list(observations=observations)))
    Sampling <- rep(i, times=observations)
    result <- cbind(Sampling, result)
    sampleList[[i]] <- result
  }
resultDf <- bind_rows(sampleList)
  return(resultDf)
}

generateFlowCytometryDataCoFunction <- function(probabilities, stdevs, observations){

  components <- sample(1:2,prob=probabilities,size=observations,replace=TRUE)

  samples <- rnorm(n=observations,mean=c(0,10)[components],sd=stdevs[components])
  return(samples)

}
