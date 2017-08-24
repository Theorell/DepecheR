#' Creates a vector of colors of the same length as the data
#'
#'
#' This function takes a vector x and a shorter ordering vector with all the unique values of the x vector in the specific order that the colors should be in and returns a vector of RGB colors the same length as the initial x vector.
#' @importFrom gplots rich.colors
#' @param x A vector, in most cases of identities of individuals or clusters ect.
#' @param order The order, folowing a rainbow distribution, that the colors should be in in the output vector. Defaults to the order that the unique values in x occurs.
#' @return A vector, the same length as x with each unique value substitutet with a color.
#' @examples
#' #Generate a dataframe with bimodally distributed data and a few separate subsamplings
#' x <- generateFlowCytometryData(samplings=5, observations=2000)
#'
#' #Scale the data (not actually necessary in this artificial 
#' #example due to the nature of the generated data)
#' x_scaled <- quantileScale(x=x[2:ncol(x)])
#'
#' #Run Barnes Hut tSNE on this. 
#' library(Rtsne.multicore)
#' xSNE <- Rtsne.multicore(x_scaled, pca=FALSE)
#'
#' #Now use our function
#' xColors <- colorVector(x[,1])
#'
#' #Set a reasonable working directory, e.g.
#' setwd("~/Desktop")
#'
#' #Plot all ids together and use rainbowColors
#' depecheDensity(xYData=as.data.frame(xSNE$Y), commonName="All_samplings", 
#' color=xColors, createDirectory=FALSE)
#'
#' @export colorVector
colorVector <- function(x, order=unique(x)){
	
	orderNumbers <- quantileScale(c(1:length(order)), robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE, multiplicationFactor=100)
	
	#This is done to prevent any events being fully at a 100 or at 0.
	orderColorNumbers <- quantileScale(c(1:length(orderNumbers)), control=c(-3, 103), robustVarScale=FALSE, lowQuantile=0, highQuantile=1, center=FALSE, multiplicationFactor=100)
	
	orderColors <- orderColorNumbers
	for(i in 1:length(orderColorNumbers)){
 
		orderColors[i] <- colorRampPalette(orderColorNumbers[i])(1)

	}

  	#Here, a vector with the same length as the x vector is generated, but where the x info has been substituted with a color.
  	colorVector <- x
  		for(i in 1:length(order)){
    	colorVector[x==order[i]] <- 	orderColors[i]
  	}

  return(colorVector)
  	
}


