######################
# dViolins
context("violin")
data("testData")
data("testDataDepeche")
dViolins(testDataDepeche$clusterVector, testDataDepeche$sparsityMatrix, inDataFrame = testData[, 2:15], createOutput = FALSE)

######################
# dViolinsCoFunction1
DepecheR:::dViolinsCoFunction1(30, 20)

######################
# dViolinsCoFunction2
data("testData")
data("testDataDepeche")

clusterVector <- testDataDepeche$clusterVector
clusterCenters <- testDataDepeche$sparsityMatrix
percentClusterVector <- DepecheR:::dScale(clusterVector, scale = c(0, 1), robustVarScale = FALSE, center = FALSE, multiplicationFactor = 100)
inDataFocused <- subset(testData[, 2:15], select = colnames(clusterCenters))
clustIndicesSpecific <- sapply(clusterVector, DepecheR:::dViolinsCoFunction1, n = 3)
oneClustAllMu <- clusterCenters[rownames(clusterCenters) == 3, ]
clustColorsSpecific <- c("red", "#d3d3d3")
allVarNames <- colnames(inDataFocused)
oneClustAllVarList <- DepecheR:::dViolinsCoFunction2(inDataFocused[, 1], oneClustAllMu[1], allVarNames[1], clust = clustIndicesSpecific, cols = clustColorsSpecific, clustNum = 3)

######################
# dViolinsCoFunction3
data("testData")
data("testDataDepeche")

clusterVector <- testDataDepeche$clusterVector
clusterCenters <- testDataDepeche$sparsityMatrix
percentClusterVector <- DepecheR:::dScale(clusterVector, scale = c(0, 1), robustVarScale = FALSE, center = FALSE, multiplicationFactor = 100)
inDataFocused <- subset(testData[, 2:15], select = colnames(clusterCenters))
clustIndicesSpecific <- sapply(clusterVector, DepecheR:::dViolinsCoFunction1, n = 3)
oneClustAllMu <- clusterCenters[rownames(clusterCenters) == 3, ]
clustColorsSpecific <- c("red", "#d3d3d3")
allVarNames <- colnames(inDataFocused)
oneClustAllVarList <- DepecheR:::dViolinsCoFunction2(inDataFocused[, 1], oneClustAllMu[1], allVarNames[1], clust = clustIndicesSpecific, cols = clustColorsSpecific, clustNum = 3)

DepecheR:::dViolinsCoFunction3(oneClustAllVarList, plotAll = FALSE, createOutput = FALSE)
