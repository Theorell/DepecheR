######################
context('scaling')
#dScale
x <- DepecheR:::generateBimodalData()[,-1]
y_df <- DepecheR:::dScale(x)

######################
#dScaleCoFunction
x <- DepecheR:::generateBimodalData()
result <- DepecheR:::dScaleCoFunction(x[,2], control=x[,2], scale=c(0.001, 0.999), robustVarScale=TRUE, truncate=c(0.001, 0.999), center="mean", multiplicationFactor=1)