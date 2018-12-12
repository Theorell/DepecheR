#dColorPlot
context('plotting')
data(testData)
data(testDataSNE)
dColorPlot(colorData=testData[,2:3], xYData=testDataSNE$Y, drawColorPalette=TRUE, createDirectory=FALSE, createPlot=FALSE)

testDataSubset <- rbind(testData[1:2000,], testData[95001:97000,])
testDataSNESubset <- rbind(testDataSNE$Y[1:2000,], testDataSNE$Y[95001:97000,])
testColor <- dColorVector(testDataSubset$ids, colorScale="plasma")
dColorPlot(colorData=testColor, xYData=testDataSNESubset, 
           names="separate samplings", addLegend=TRUE, idsVector=testDataSubset$ids, createDirectory=FALSE, createPlot=FALSE)

######################
#dColorPlotCoFunction
data(testData)
data(testDataSNE)
testDataSubset <- rbind(testData[1:2000,], testData[95001:97000,])
testDataSNESubset <- rbind(testDataSNE$Y[1:2000,], testDataSNE$Y[95001:97000,])

densContour <- DepecheR:::dContours(as.data.frame(testDataSNESubset))
testDataSNESubsetFraction <- DepecheR:::dScale(as.data.frame(testDataSNESubset), scale=c(0,1), robustVarScale=FALSE, center=FALSE)
testColor <- dColorVector(testDataSubset$ids, colorScale="plasma")
DepecheR:::dColorPlotCoFunction(colorVariable=testColor, name="separate samplings", xYData=testDataSNESubsetFraction, title=FALSE, densContour=densContour, bandColor="black", dotSize=400/sqrt(nrow(testDataSNESubsetFraction)), createPlot=FALSE)

######################
#dColorVector
data(testData)
testColor <- dColorVector(testData$ids, colorScale="plasma")

#####################
#dContours
x <- DepecheR:::generateBimodalData()
xContours <- DepecheR:::dContours(x[[1]][,2:3])

######################
#dDensityPlot
data(testData)
data(testDataSNE)

dDensityPlot(xYData=testDataSNE$Y, commonName="All_samplings", 
             color="blue", createDirectory=FALSE, createPlot=FALSE)

#Alternative usage
testDataSubset <- rbind(testData[1:1000,], testData[96001:97000,])
testDataSNESubset <- rbind(testDataSNE$Y[1:1000,], testDataSNE$Y[96001:97000,])
testColor <- dColorVector(testDataSubset$ids, colorScale="plasma")

dDensityPlot(xYData=testDataSNESubset, color=testColor, plotEachIdSeparately=TRUE, 
             idsVector=testDataSubset$ids, commonName="sampling", createDirectory=FALSE, createPlot=FALSE)

#Alternative usage
dDensityPlot(xYData=testDataSNESubset, color=testColor, idsVector=testDataSubset$ids,
             commonName="all samplings", createDirectory=FALSE, createPlot=FALSE)

######################
#dDensityPlotCoFunction
xYData <- DepecheR:::dScale(as.data.frame(testDataSNE$Y), as.data.frame(testDataSNE$Y), scale=c(0,1), robustVarScale=FALSE, center=FALSE)
cols <- colorRampPalette(c("black", "grey", "red"))(256)
densContour <- DepecheR:::dContours(as.data.frame(testDataSNE$Y))
DepecheR:::dDensityPlotCoFunction(xYData=xYData, cols=cols, name="test", densContour=densContour, bandColor="black", dotSize=1.5, title=FALSE, createPlot=FALSE)


######################
#dResidualPlot
data(testData)
data(testDataSNE)
data(testDataDepeche)
dResidualPlot(xYData=testDataSNE$Y, groupVector=testData[,16], clusterVector=testDataDepeche$clusterVector, createPlot=FALSE)