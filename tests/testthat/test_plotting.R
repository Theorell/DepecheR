# dColorPlot
context("plotting")
data(testData)
data(testDataSNE)
allDataRows <- c(
    1:200, 1001:1200, 10001:10200,
    11001:11200
)
testDataSubset <- testData[allDataRows, ]
testDataSNESubset <- testDataSNE$Y[allDataRows, ]

dColorPlot(
    colorData = testData[1:100, 2], xYData = testDataSNE$Y[1:100, ],
    createOutput = FALSE
)

dColorPlot(
    colorData = testDataSubset$ids, xYData = testDataSNESubset,
    createOutput = FALSE
)

###################### dColorVector
testColor <- dColorVector(testDataSubset$ids, colorScale = "plasma")

##################### dContours
densContour <- DepecheR:::dContours(testDataSNESubset)

###################### dPlotCoFunction
xYData <- dScale(testDataSNESubset, testDataSNESubset,
    scale = c(0, 1),
    robustVarScale = FALSE, center = FALSE
)

DepecheR:::dPlotCoFunction(
    colorVariable = testColor, plotName = "test",
    xYData = xYData, title = FALSE,
    densContour = densContour, bandColor = "black",
    dotSize = 400 / sqrt(nrow(testDataSNESubsetFraction)),
    plotDir = ".", createOutput = FALSE
)

###################### dDensityPlot
dDensityPlot(xYData = testDataSNESubset, createOutput = FALSE)

# Alternative usage
dDensityPlot(
    xYData = testDataSNESubset, idsVector = testDataSubset$ids,
    plotDir = ".", createOutput = FALSE
)

###################### dDensityPlotCoFunction
DepecheR:::dDensityPlotCoFunction(
    xYData = xYData, color = "blue",
    plotName = "test", densContour = densContour,
    bandColor = "black", dotSize = 1.5,
    title = FALSE, plotDir = ".",
    createOutput = FALSE
)

###################### dResidualPlot
data(testDataDepeche)
dResidualPlot(
    xYData = testDataSNESubset, groupVector = testDataSubset$label,
    clusterVector = testDataDepeche$clusterVector[allDataRows],
    createOutput = FALSE
)

###################### dResidualPlot
dataTrans <-
    testDataSubset[
        , c("SYK", "CD16", "CD57", "EAT.2", "CD8", "NKG2C", "CD2", "CD56")
    ]
groupProbPlot(
    xYData = testDataSNESubset, groupVector = testDataSubset$label,
    dataTrans, createOutput = FALSE
)
