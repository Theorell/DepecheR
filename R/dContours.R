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
#' contour_result <- dContours(testDataSNE$Y)
#' @export dContours
dContours <- function(xYData, control, n = 100) {
    if (missing(control)) {
        control <- xYData
    }

    range1 <- range(control[, 1])
    sideDist1 <- 0.05 * (range1[2] - range1[1])
    range2 <- range(control[, 2])
    sideDist2 <- 0.05 * (range2[2] - range2[1])

    lims <- c(
        range1[1] - sideDist1,
        range1[2] + sideDist1,
        range2[1] - sideDist2,
        range2[2] + sideDist2
    )

    # Construct the third dimension with
    # smooth kernel density estimate
    den3d <- kde2d(xYData[, 1], xYData[, 2], n = n, lims = lims)

    return(den3d)
}
