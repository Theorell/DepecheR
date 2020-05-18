###################### dSplsda
context("sPLSDA")
data("testData")
data("testDataDepeche")
data("testDataSNE")
allDataRows <- c(
    1:200, 1001:1200, 2001:2200, 3001:3200, 4001:4200,
    10001:10200, 11001:11200, 12001:12200, 13001:13200,
    14001:14200
)
xYData <- testDataSNE$Y[allDataRows, ]
idsVector <- testData$ids[allDataRows]
groupVector <- testData$label[allDataRows]
clusterVector <- testDataDepeche$clusterVector[allDataRows]
subsetVector <- sample(1:2000, size = 100)

sPLSDAObject <- dSplsda(
    xYData = xYData, idsVector = idsVector,
    groupVector = groupVector,
    clusterVector = clusterVector,
    displayVector = subsetVector,
    createOutput = FALSE
)

# Alternative usage
sPLSDAObject <- dSplsda(
    xYData = xYData, idsVector = idsVector,
    groupVector = groupVector,
    clusterVector = clusterVector, paired = TRUE,
    createOutput = FALSE
)

# Alternative usage
sPLSDAObject <- dSplsda(
    xYData = xYData, idsVector = idsVector,
    groupVector = groupVector,
    clusterVector = clusterVector,
    testSampleRows = subsetVector, createOutput = FALSE
)

###################### dSplsdaPreCalculations
data("testData")
data("testDataDepeche")
data("testDataSNE")
dSplsdaInData <- DepecheR:::dSplsdaPreCalculations(
    clusterVector = clusterVector, idsVector = idsVector,
    groupVector = groupVector, groupName1 = "Stim1", groupName2 = "Stim2"
)
