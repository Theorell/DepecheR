#' Create a vector of colors of the same length as the data
#'
#'
#' This function takes a vector x and a shorter ordering vector with all the 
#' unique values of the x vector in the specific order that the colors should be
#' in and returns a vector of RGB colors the same length as the initial x 
#' vector.
#' @importFrom viridis inferno magma plasma viridis
#' @importFrom gplots rich.colors
#' @importFrom grDevices rainbow
#' @param x A vector, in most cases of identities of individuals or clusters, 
#' etc.
#' @param colorOrder The order, folowing a rainbow distribution, that the colors 
#' should be in in the output vector. Defaults to the order that the unique 
#' values in x occurs.
#' @param colorScale The color scale. Inherited from the viridis, gplots and 
#' grDevices packages (and the package-specific 'dark_rainbow'). Seven possible
#' scales are pre-made: inferno, magma, plasma, viridis, rich_colors, rainbow 
#' and dark_rainbow. User specified vectors of colors 
#' (e.g. c('#FF0033', '#03AF49')) are also accepted.
#' @return A vector, the same length as x with each unique value substitutet 
#' with a color.
#' @seealso \code{\link{dDensityPlot}}, \code{\link{dColorPlot}}, 
#' \code{\link{dViolins}}
#' @examples
#' # Load some data
#' data(testData)
#'
#' testColor <- dColorVector(testData$ids, colorScale = 'plasma')
#'
#' # In this case, each of the 97 individual donors in the dataset has gotten 
#' their own color code:
#' table(testColor)
#' @export dColorVector
dColorVector <- function(x, colorOrder = unique(x), colorScale = "viridis") {
    if (is.factor(x)) {
        x <- as.character(x)
        colorOrder <- as.character(colorOrder)
    }
    nColors <- length(colorOrder)
    if (length(colorScale) > 1) {
        orderColors <- colorRampPalette(colorScale)(nColors)
    } else {
        orderColors <- switch(colorScale,
                              "inferno"=inferno(nColors),
                              "viridis"=viridis(nColors),
                              "plasma"=plasma(nColors), 
                              "magma"=magma(nColors), 
                              "rich_colors"=rich.colors(nColors),
                              "rainbow"=rainbow(nColors),
                              "dark_rainbow"=
                                  colorRampPalette(c("#990000", "#FFCC00", 
                                                     "#336600", "#000066", 
                                                     "#660033")
                                                   )(nColors))
        
    }
    
    # Here, a vector with the same length as
    # the x vector is generated, but where
    # the x info has been substituted with a
    # color.
    colorVector <- x
    for (i in seq(0,length(colorOrder))) {
        colorVector[x == colorOrder[i]] <- orderColors[i]
    }
    
    return(colorVector)
}
