# Function that generates synthetic data from a bimodal distribution.
# centers: A list with two items, each item a vector specifying the center of a
# data cluster.
# prop: A scalar determining relative number of data points in the two clusters
# dataCols: A scalar determining the number of dimensions. Superfluous if
# "centers" is provided.
# observations: A scalar determining the number of sampled data points.
generateBimodalData <- function(centers, prop = 0.3, dataCols = 5,
                                observations = 10000) {
    if (missing(centers)) {
        centers <- rbind(runif(dataCols) + 50, runif(dataCols) - 50)
    } else {
        dataCols <- ncol(centers)
    }
    stopifnot(nrow(centers) == 2)

    ids <- c(
        rep(1, floor(prop * observations)),
        rep(2, observations - floor(prop * observations))
    )
    rands <- matrix(
        rnorm(observations * dataCols, mean = 0, sd = 1),
        observations, dataCols
    )
    pop1 <- rands[seq_len(floor(prop * observations)), ] +
        rep(centers[1, ], each = floor(prop * observations))

    pop2 <- tail(rands, n = observations - floor(prop * observations)) +
        rep(centers[2, ], each = observations - floor(prop * observations))

    samples <- rbind(pop1, pop2)
    result <- list(samples, ids)
    names(result) <- c("samples", "ids")

    return(result)
}

# Function that generates synthetic data, divided over a number of cluster with
# increasing number of zero-dimensions (starting from 0). The sampled points are
# generated equally between the clusters.

# modeN: A scalar determining the number of clusters.
# dataCols: A scalar determing the number of dimensions. DataCols must be
# greater than modeN, to avoid clusters with only zeros.
# obsrevations: A scalar determining the number of sampled data points.

generateSparseData <- function(modeN = 5, dataCols = 100,
                               observations = 10000) {
    # Check if input ok
    if ((observations / modeN) %% 1 != 0) {
        stop("Observations has to be divisible by modeN")
    }
    obsPerMode <- observations / modeN

    # generate the centers numbers

    randInts <- sample(c(-50, 50), dataCols * modeN, replace = TRUE)
    centers <- matrix(randInts, nrow = modeN, byrow = TRUE)
    # put in sparsity
    for (i in seq_len(modeN)) {
        inds <- sample(seq_len(dataCols), i)
        centers[i, inds] <- 0
    }

    # generate the data
    samples <- matrix(0, nrow = observations, ncol = dataCols)
    ids <- matrix(0, nrow = observations, ncol = 1)
    for (i in seq_len(modeN)) {
        temp <- matrix(
            rnorm(obsPerMode * dataCols, mean = 0, sd = 1),
            obsPerMode, dataCols
        )
        samples[seq((1 + obsPerMode * (i - 1)), obsPerMode * (i)), ] <- temp +
            rep(centers[i, ], each = obsPerMode)
        ids[seq((1 + obsPerMode * (i - 1)), obsPerMode * (i)), ] <- i
    }
    result <- list(samples, ids, centers)
    names(result) <- c("samples", "ids", "centers")
    return(result)
}
