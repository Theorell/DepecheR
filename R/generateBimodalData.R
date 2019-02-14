generateBimodalData <- function(centers, 
    prop = 0.3, dataCols = 5, observations = 10000) {
    if (missing(centers)) {
        centers <- rbind(runif(dataCols) + 50, runif(dataCols) - 50)
    } else {
        dataCols <- ncol(centers)
    }
    stopifnot(nrow(centers) == 2)
    
    ids <- c(rep(1, floor(prop * observations)), 
             rep(2, observations - floor(prop * observations)))
    rands <- matrix(rnorm(observations * dataCols, mean = 0, sd = 1), 
                    observations, dataCols)
    pop1 <- rands[seq_len(floor(prop * observations)), ] + 
        rep(centers[1, ], each = floor(prop * observations))
    
    pop2 <- tail(rands, n = observations - floor(prop * observations)) + 
        rep(centers[2, ], each = observations - floor(prop * observations))
    
    samples <- rbind(pop1, pop2)
    result <- list(samples, ids)
    names(result) <- c("samples", "ids")
    
    return(result)
}

generateSparseData <- function(modeN = 5, dataCols = 100, 
                               observations = 10000) {
    # Check if input ok
    if ((observations/modeN)%%1 != 0) {
        stop("Observations has to be divisible by modeN")
    }
    obsPerMode <- observations/modeN
    
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
        temp <- matrix(rnorm(observations * dataCols, mean = 0, sd = 1), 
                       obsPerMode, dataCols)
        samples[seq((1 + obsPerMode * (i - 1)), obsPerMode * (i)), ] <- temp + 
            rep(centers[i, ], each = obsPerMode)
        ids[seq((1 + obsPerMode * (i - 1)), obsPerMode * (i)), ] <- i
    }
    result <- list(samples, ids, centers)
    names(result) <- c("samples", "ids", "centers")
    return(result)
}
