###################### dSplsda
context("sPLSDA")
data("testData")
data("testDataDepeche")
data("testDataSNE")
sPLSDAObject <- dSplsda(xYData = testDataSNE$Y, idsVector = testData$ids, groupVector = testData$label, clusterVector = testDataDepeche$clusterVector, 
    createOutput = FALSE)

# Alternative usage
xYDataPaired <- as.data.frame(testDataSNE$Y)[c(1:30000, 67001:nrow(testDataSNE$Y)), ]
testDataPaired <- testData[c(1:30000, 67001:nrow(testDataSNE$Y)), ]
clusterVectorPaired <- testDataDepeche$clusterVector[c(1:30000, 67001:nrow(testDataSNE$Y))]
sPLSDAObject <- dSplsda(xYData = xYDataPaired, idsVector = testDataPaired$ids, groupVector = testDataPaired$label, 
    clusterVector = clusterVectorPaired, paired = TRUE, name = "d_sPLSDAPlot_paired", groupName1 = "Stimulation 1", 
    groupName2 = "Stimulation 2", createOutput = FALSE)

# Alternative usage
subsetVector <- sample(1:nrow(testData), size = 10000)
testDataSNESubset <- testDataSNE$Y[subsetVector, ]
sPLSDAObject <- dSplsda(xYData = testDataSNESubset, idsVector = testData$ids, groupVector = testData$label, clusterVector = testDataDepeche$clusterVector, 
    displayVector = subsetVector, createOutput = FALSE)

# Alternative usage
testDataRows <- sample(1:nrow(testData), size = 48500)
sPLSDAObject <- dSplsda(xYData = testDataSNE$Y, idsVector = testData$ids, groupVector = testData$label, clusterVector = testDataDepeche$clusterVector, 
    testSampleRows = testDataRows, createOutput = FALSE)

###################### dSplsdaPreCalculations
data("testData")
data("testDataDepeche")
data("testDataSNE")
dSplsdaInData <- DepecheR:::dSplsdaPreCalculations(clusterVector = testDataDepeche$clusterVector, idsVector = testData$ids, 
    groupVector = testData$label, groupName1 = "Stimulation 1", groupName2 = "Stimulation 2")
