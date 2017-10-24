turnVectorEquidistant <- function(x, startValue=1, newNumbers){

	originalNumbers <- sort(unique(x))
	
	if(missingNewNumbers==TRUE){
	  newNumbers <- c(startValue:(length(originalNumbers)+(startValue-1)))
	}
	
  result <- x
	for(i in 1:length(originalNumbers)){
	  result[(x==originalNumbers[i])] <- newNumbers[i]
	}

	return(as.numeric(result))
}
