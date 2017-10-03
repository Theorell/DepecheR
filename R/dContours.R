#' Create density contours for two-dimensional data.
#'
#'
#' Here, contour lines for two-dimensional data are construced. It is primarily thought to be used in the context of SNE plots in this package. This function is used both internally in other functions suchas sneFluoroPlot and sneDensityPlot, but also as a standalone function, as it increases speed greatly to generate the density curves only once per overall analysis.
#' @importFrom MASS kde2d
#' @param sneData A dataframe with two columns containing position information for each observation in the dataset. Typically, this is the raw result from the SNE analysis.
#' @param n The number fo grid points. Default is 100.
#' @seealso \code{\link{dColorPlot}}, \code{\link{dDensityPlot}}, \code{\link{dResidualPlot}}, \code{\link{dWilcoxPlot}}
#' @return A list of three components
#' \describe{
#'     \item{x, y}{The x and y coordinates of the grid points, vectors of length n.}
#'     \item{z}{An n[1] by n[2] matrix of the estimated density: rows correspond to the value of x, columns to the value of y.}
#' }
#'
#' @examples
#' #Generate a dataframe with two numeric, normally distributed vectors.
#' #Generate a default size dataframe with bimodally distributed data
#' x <- generateBimodalData(samplings=2)
#'
#' #Scale this datamframe
#' x_scaled <- dScale(x[,2:ncol(x)], center=FALSE)
#' #Run the function
#' contour_result <- dContours(x_scaled[,1:2])
#'
#' #Plot the result
#' contour(x=contour_result$x, y=contour_result$y, z=contour_result$z, 
#' xlim=c(-0.2, 1.2), 
#' ylim=c(-0.2, 1.2), nlevels=10, lwd=2, drawlabels = FALSE, axes=FALSE, 
#' xaxs="i", yaxs="i")
#'
#' @export dContours

dContours <- function(sneData, n=100){

	sneDataNorm <- dScale(x=sneData, scale=c(0,1), robustVarScale=FALSE, center=FALSE)

	V1 <- sneDataNorm[,1]
	V2 <- sneDataNorm[,2]

#Construct the third dimension with smooth kernel density estimate
den3d <- kde2d(V1, V2, n=n, lims=c(-0.05,1.05, -0.05,1.05))

return(den3d)
}


