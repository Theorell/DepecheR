#dColorPlot
context('plotting')
data(testData)
data(testDataSNE)
setwd("~/Desktop")
dColorPlot(colorData=testData[,2:15], xYData=testDataSNE$Y, drawColorPalette=TRUE)

testDataSubset <- rbind(testData[1:2000,], testData[95001:97000,])
testDataSNESubset <- rbind(testDataSNE$Y[1:2000,], testDataSNE$Y[95001:97000,])
testColor <- dColorVector(testDataSubset$ids, colorScale="plasma")
dColorPlot(colorData=testColor, xYData=testDataSNESubset, 
           names="separate samplings", addLegend=TRUE, idsVector=testDataSubset$ids)

######################
#dColorPlotCoFunction
data(testData)
data(testDataSNE)
testDataSubset <- rbind(testData[1:2000,], testData[95001:97000,])
testDataSNESubset <- rbind(testDataSNE$Y[1:2000,], testDataSNE$Y[95001:97000,])

densContour <- DepecheR:::dContours(as.data.frame(testDataSNESubset))
testDataSNESubsetFraction <- DepecheR:::dScale(as.data.frame(testDataSNESubset), scale=c(0,1), robustVarScale=FALSE, center=FALSE)
testColor <- dColorVector(testDataSubset$ids, colorScale="plasma")
DepecheR:::dColorPlotCoFunction(colorVariable=testColor, name="separate samplings", xYDataFraction=testDataSNESubsetFraction, title=FALSE, densContour=densContour, bandColor="black", dotSize=400/sqrt(nrow(testDataSNESubsetFraction)))

######################
#dColorVector
data(testData)
testColor <- dColorVector(testData$ids, colorScale="plasma")

#####################
#dContours
x <- DepecheR:::generateBimodalData()
xContours <- DepecheR:::dContours(x[,2:ncol(x)])

######################
#dDensityPlot
data(testData)
data(testDataSNE)
setwd("~/Desktop")

dDensityPlot(xYData=testDataSNE$Y, commonName="All_samplings", 
             color="blue", createDirectory=FALSE)

#Alternative usage
testDataSubset <- rbind(testData[1:2000,], testData[95001:97000,])
testDataSNESubset <- rbind(testDataSNE$Y[1:2000,], testDataSNE$Y[95001:97000,])
testColor <- dColorVector(testDataSubset$ids, colorScale="plasma")

dDensityPlot(xYData=testDataSNESubset, color=testColor, plotEachIdSeparately=TRUE, 
             idsVector=testDataSubset$ids, commonName="sampling")

#Alternative usage
dDensityPlot(xYData=testDataSNESubset, color=testColor, idsVector=testDataSubset$ids,
             commonName="all samplings")

######################
#dDensityPlotCoFunction
xYDataScaled <- DepecheR:::dScale(as.data.frame(testDataSNE$Y), as.data.frame(testDataSNE$Y), scale=c(0,1), robustVarScale=FALSE, center=FALSE)
cols <- colorRampPalette(c("black", "grey", "red"))(256)
densContour <- DepecheR:::dContours(as.data.frame(testDataSNE$Y))
DepecheR:::dDensityPlotCoFunction(xYDataScaled=xYDataScaled, cols=cols, name="test", densContour=densContour, bandColor="black", dotSize=1.5, title=FALSE)


######################
#dResidualPlot
data(testData)
data(testDataSNE)
data(testDataDepeche)
dResidualPlot(xYData=testDataSNE$Y, groupVector=testData[,16], clusterVector=testDataDepeche$clusterVector)