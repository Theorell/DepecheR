######################
context('scaling')
#dScale
x <- DepecheR:::generateBimodalData()
y_df <- DepecheR:::dScale(x[[1]])

######################
#dScaleCoFunction
x <- DepecheR:::generateBimodalData()
result <- DepecheR:::dScaleCoFunction(x[[1]][,2], control=x[[1]][,2], scale=c(0.001, 0.999), robustVarScale=TRUE, truncate=c(0.001, 0.999), center="mean", multiplicationFactor=1)