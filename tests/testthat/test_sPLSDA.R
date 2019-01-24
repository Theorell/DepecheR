###################### dSplsda
context("sPLSDA")
data("testData")
data("testDataDepeche")
data("testDataSNE")
allDataRows <- c(1:200, 1001:1200, 2001:2200, 3001:3200, 4001:4200, 
                 55001:55200, 56001:56200, 57001:57200, 58001:58200, 
                 59001:59200)
xYData = testDataSNE$Y[allDataRows,]
idsVector = testData$ids[allDataRows]
groupVector = testData$label[allDataRows]
clusterVector = testDataDepeche$clusterVector[allDataRows]

sPLSDAObject <- dSplsda(xYData = xYData, idsVector = idsVector, 
                        groupVector = groupVector, 
                        clusterVector = clusterVector, createOutput = FALSE)

# Alternative usage
sPLSDAObject <- dSplsda(
    xYData = xYData, idsVector = idsVector,
    groupVector = groupVector, clusterVector = clusterVector, 
    paired = TRUE, createOutput = FALSE)

# Alternative usage
subsetVector <- sample(1:2000, size = 100)
sPLSDAObject <- dSplsda(xYData = xYData, idsVector = idsVector, 
                        groupVector = groupVector, 
                        clusterVector = clusterVector, 
                        displayVector = subsetVector, createOutput = FALSE)

# Alternative usage
testDataRows <- sample(1:2000, size = 1000)
sPLSDAObject <- dSplsda(xYData = xYData, idsVector = idsVector, 
                        groupVector = groupVector, 
                        clusterVector = clusterVector, 
                        testSampleRows = testDataRows, createOutput = FALSE)

###################### dSplsdaPreCalculations
data("testData")
data("testDataDepeche")
data("testDataSNE")
dSplsdaInData <- DepecheR:::dSplsdaPreCalculations(
    clusterVector = clusterVector, idsVector = idsVector, 
    groupVector = groupVector, groupName1 = "Stim1", groupName2 = "Stim2")
