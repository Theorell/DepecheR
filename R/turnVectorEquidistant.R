turnVectorEquidistant <- function(x, startValue = 1, 
    newNumbers) {
    originalNumbers <- sort(unique(x))
    
    if (missing(newNumbers)) {
        newNumbers <- c(startValue:(length(originalNumbers) + 
            (startValue - 1)))
    }
    
    result <- x
    for (i in seq_len(length(originalNumbers))) {
        result[(x == originalNumbers[i])] <- newNumbers[i]
    }
    
    return(as.numeric(result))
}
