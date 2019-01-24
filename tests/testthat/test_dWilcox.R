###################### dWilcox
context("dWilcox")
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

result <- dWilcox(xYData = xYData, idsVector = idsVector, 
                  groupVector = groupVector, 
                  clusterVector = clusterVector, createOutput = FALSE)

# Alternative usage
subsetVector <- sample(1:nrow(testData), size = 1000)

result <- dWilcox(xYData = xYData, idsVector = idsVector, 
                  groupVector = groupVector, clusterVector = clusterVector, 
    displayVector = subsetVector, createOutput = FALSE)
