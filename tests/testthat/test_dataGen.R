######################
#generateBimodalData
context('dataGen')
x <- DepecheR:::generateBimodalData()

######################
#generateBimodalDataCoFunction
result <- DepecheR:::generateBimodalDataCoFunction(c(0.9, 0.7), c(0.3, 0.7), 1000)