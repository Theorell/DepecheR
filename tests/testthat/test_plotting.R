# dColorPlot
context("plotting")
data(testData)
data(testDataSNE)
allDataRows <- c(1:200, 1001:1200, 58001:58200, 
                 59001:59200)
testDataSubset <- testData[allDataRows,]
testDataSNESubset <- testDataSNE$Y[allDataRows,]

dColorPlot(colorData = testData[1:100, 2], xYData = testDataSNE$Y[1:100,], 
           drawColorPalette = TRUE, createDirectory = FALSE, 
           createPlot = FALSE)

testColor <- dColorVector(testDataSubset$ids, colorScale = "plasma")
dColorPlot(colorData = testColor, xYData = testDataSNESubset, names = "separate samplings", addLegend = TRUE, idsVector = testDataSubset$ids, 
    createDirectory = FALSE, createPlot = FALSE)

###################### dColorPlotCoFunction
densContour <- DepecheR:::dContours(as.data.frame(testDataSNESubset))
testDataSNESubsetFraction <- DepecheR:::dScale(as.data.frame(testDataSNESubset), scale = c(0, 1), robustVarScale = FALSE, 
    center = FALSE)
testColor <- dColorVector(testDataSubset$ids, colorScale = "plasma")
DepecheR:::dColorPlotCoFunction(colorVariable = testColor, name = "separate samplings", xYData = testDataSNESubsetFraction, 
    title = FALSE, densContour = densContour, bandColor = "black", dotSize = 400/sqrt(nrow(testDataSNESubsetFraction)), 
    createPlot = FALSE)

###################### dColorVector
testColor <- dColorVector(testDataSubset$ids, colorScale = "plasma")

##################### dContours
xContours <- DepecheR:::dContours(testDataSNESubset)

###################### dDensityPlot
dDensityPlot(xYData = testDataSNESubset, createDirectory=FALSE, createPlot = FALSE)

# Alternative usage
dDensityPlot(xYData = testDataSNESubset, color = testColor, 
             plotEachIdSeparately = TRUE, idsVector = testDataSubset$ids, 
             createDirectory = FALSE, createPlot = FALSE)

# Alternative usage
dDensityPlot(xYData = testDataSNESubset, color = testColor, 
             idsVector = testDataSubset$ids, createDirectory = FALSE, 
             createPlot = FALSE)

###################### dDensityPlotCoFunction
xYData <- dScale(testDataSNESubset, testDataSNESubset, scale = c(0, 1), 
                 robustVarScale = FALSE, center = FALSE)
cols <- colorRampPalette(c("black", "grey", "red"))(256)
densContour <- DepecheR:::dContours(testDataSNESubset)
DepecheR:::dDensityPlotCoFunction(xYData = xYData, cols = cols, name = "test", 
                                  densContour = densContour, 
                                  bandColor = "black", dotSize = 1.5, 
                                  title = FALSE, createPlot = FALSE)


###################### dResidualPlot
data(testDataDepeche)
dResidualPlot(xYData = testDataSNESubset, groupVector = testDataSubset$label, 
              clusterVector = testDataDepeche$clusterVector[allDataRows], 
    createPlot = FALSE)
