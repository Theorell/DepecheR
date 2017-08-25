#' Creates a vector of colors of the same length as the data
#'
#'
#' This function takes a vector x and a shorter ordering vector with all the unique values of the x vector in the specific order that the colors should be in and returns a vector of RGB colors the same length as the initial x vector.
#' @importFrom viridis inferno magma plasma viridis
#' @param x A vector, in most cases of identities of individuals or clusters ect.
#' @param order The order, folowing a rainbow distribution, that the colors should be in in the output vector. Defaults to the order that the unique values in x occurs.
#' @param colorScale The color scale. Inherited from the viridis package. Four possible scales: inferno, magma, plasma, viridis.
#' @return A vector, the same length as x with each unique value substitutet with a color.
#' @seealso \code{\link{dDensityPlot}}, \code{\link{dColorPlot}}, \code{\link{dViolins}}
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
#' dDensityPlot(xYData=as.data.frame(xSNE$Y), idsVector=x[,1], commonName="All_samplings", 
#' color=xColors, createDirectory=FALSE)
#'
#' @export colorVector
colorVector <- function(x, order=unique(x), colorScale="inferno"){
	
  if(colorScale=="inferno"){
    orderColors <- inferno(length(order)) 
  }
  if(colorScale=="viridis"){
    orderColors <- viridis(length(order)) 
  }
  if(colorScale=="plasma"){
    orderColors <- plasma(length(order)) 
  }
  if(colorScale=="magma"){
    orderColors <- magma(length(order)) 
  }
 
  	#Here, a vector with the same length as the x vector is generated, but where the x info has been substituted with a color.
  	colorVector <- x
  		for(i in 1:length(order)){
    	colorVector[x==order[i]] <- 	orderColors[i]
  	}

  return(colorVector)
  	
}


