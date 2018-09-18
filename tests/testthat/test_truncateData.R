######################
#truncateData
context('truncateData')
x <- DepecheR:::generateBimodalData()[,-1]
xTruncated <- DepecheR:::truncateData(x)

######################
#truncateDataCoFunction
x <- DepecheR:::generateBimodalData()
result <- DepecheR:::truncateDataCoFunction(x[,2], control=x[,2], lowQuantile=0.1, highQuantile=0.9)