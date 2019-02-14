#' Create density contours for two-dimensional data.
#'
#'
#' Here, contour lines for two-dimensional data are construced. It is primarily
#' thought to be used in the context of SNE plots in this package. This function
#' is used both internally in other functions suchas sneFluoroPlot and 
#' sneDensityPlot, but also as a standalone function, as it increases speed 
#' greatly to generate the density curves only once per overall analysis.
#' @importFrom MASS kde2d
#' @param xYData A dataframe with two columns containing position information 
#' for each observation in the dataset. Typically, this is the raw result from 
#' the SNE analysis.
#' @param control A numeric/integer vector or dataframe of values that could be
#' used to define the range in the internal dScale. If no control data is 
#' present, the function defaults to using the indata as control data.
#' @param n The number fo grid points. Default is 100.
#' @seealso \code{\link{dColorPlot}}, \code{\link{dDensityPlot}}, 
#' \code{\link{dResidualPlot}}, \code{\link{dWilcox}}
#' @return A list of three components
#' \describe{
#'     \item{x, y}{The x and y coordinates of the grid points, 
#'     vectors of length n.}
#'     \item{z}{An n[1] by n[2] matrix of the estimated density: rows 
#'     correspond to the value of x, columns to the value of y.}
#' }
#'
#' @examples
#' # Load the test SNE data
#' data(testDataSNE)
#'
#'
#' # Run the function
#' \dontrun{
#' contour_result <- dContours(testDataSNE$Y)
#' }
#' @export dContours
dContours <- function(xYData, control, n = 100) {
    if (missing(control) == FALSE) {
        min1 <- min(control[, 1])
        max1 <- max(control[, 1])
        min2 <- min(control[, 2])
        max2 <- max(control[, 2])
    } else {
        min1 <- min(xYData[, 1])
        max1 <- max(xYData[, 1])
        min2 <- min(xYData[, 2])
        max2 <- max(xYData[, 2])
    }
    lims <- c(min1 - abs(min1 * 0.05), max1 + abs(max1 * 0.05), 
              min2 - abs(min2 * 0.05), max2 + abs(max2 * 0.05))
    
    # Construct the third dimension with
    # smooth kernel density estimate
    den3d <- kde2d(xYData[, 1], xYData[, 2], n = n, lims = lims)
    
    return(den3d)
}
