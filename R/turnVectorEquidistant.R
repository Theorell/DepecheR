##' Turn a vector into a equally ordered, but equidistant vector.
##'
##'
##' Makes a numeric vector with equal distances between each possible number based on another vector with random distributions. For visualization of clusters on SNE field.
##' @param x Any vector, numeric or otherwise.
##' @param startValue The first value in the range of integers that will be the return vector. Default is 1.
##' @return A numeric vector with equal distances between the unique values in the vector. The lowest unique value in the x vector will get the lowest value in this vector. The second lowest value in the x vector will get the second lowest value and so on.
##' @examples
##' #Make random vector with 10 unique values in the range 1-100
##' x <- sample(sample.int(100, 10), 1000, replace=TRUE)
##'
##' #The unique values are
##' sort(unique(x))
##'
##' #Run the function
##' y <- turnVectorEquidistant(x)
##'
##' #Now, the unique values are distributed equally
##' sort(unique(y))
##'
##' #The first occurence of the lowest value still comes in the same place, 
##' #so the structure of the original vector is retained.
##' identical(grep(min(x),x), grep(min(y), y))
##'
turnVectorEquidistant <- function(x, startValue=1){

	originalNumbers <- sort(unique(x))
	newNumbers <- c(startValue:(length(originalNumbers)+(startValue-1)))
  result <- x
	for(i in 1:length(originalNumbers)){
	  result[(x==originalNumbers[i])] <- newNumbers[i]
	}

	return(as.numeric(result))
}
