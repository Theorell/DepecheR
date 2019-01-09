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
generateBimodalData <- function(centers, prop = 0.3, dataCols = 5, observations = 10000) {
  if (missing(centers)) {
    centers <- rbind(runif(dataCols) + 50, runif(dataCols) - 50)
  } else {
    dataCols <- ncol(centers)
  }
  stopifnot(nrow(centers) == 2)

  ids <- c(rep(1, floor(prop * observations)), rep(2, observations - floor(prop * observations)))
  rands <- matrix(rnorm(observations * dataCols, mean = 0, sd = 1), observations, dataCols)
  pop1 <- rands[1:floor(prop * observations), ] + rep(centers[1, ], each = floor(prop * observations))

  pop2 <- tail(rands, n = observations - floor(prop * observations)) + rep(centers[2, ], each = observations - floor(prop * observations))

  samples <- rbind(pop1, pop2)
  result <- list(samples, ids)
  names(result) <- c("samples", "ids")

  return(result)
}

generateSparseData <- function(modeN = 5, dataCols = 100, observations = 10000) {
  # Check if input ok
  if ((observations / modeN) %% 1 != 0) {
    print("observations has to be divisible by modeN")
    stop()
  }
  obsPerMode <- observations / modeN

  # generate the centers numbers

  randInts <- sample(c(-50, 50), dataCols * modeN, replace = TRUE)
  centers <- matrix(randInts, nrow = modeN, byrow = TRUE)
  # put in sparsity
  for (i in 1:modeN) {
    inds <- sample(1:dataCols, i)
    centers[i, inds] <- 0
  }

  # generate the data
  samples <- matrix(0, nrow = observations, ncol = dataCols)
  ids <- matrix(0, nrow = observations, ncol = 1)
  for (i in 1:modeN) {
    temp <- matrix(rnorm(observations * dataCols, mean = 0, sd = 1), obsPerMode, dataCols)
    samples[seq((1 + obsPerMode * (i - 1)), obsPerMode * (i)), ] <- temp + rep(centers[i, ], each = obsPerMode)
    ids[seq((1 + obsPerMode * (i - 1)), obsPerMode * (i)), ] <- i
  }
  result <- list(samples, ids, centers)
  names(result) <- c("samples", "ids", "centers")
  return(result)
}
