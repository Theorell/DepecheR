######################
#dWilcox
context('dWilcox')
data("testData")
data("testDataDepeche")
data("testDataSNE")
setwd("~/Desktop")
result <- dWilcox(xYData=testDataSNE$Y, idsVector=testData$ids, 
                  groupVector=testData$label, clusterVector=testDataDepeche$clusterVector)

#Alternative usage
subsetVector <- sample(1:nrow(testData), size=10000)
testDataSNESubset <- testDataSNE$Y[subsetVector,]

result <- dWilcox(xYData=testDataSNESubset, idsVector=testData$ids, groupVector=testData$label, clusterVector=testDataDepeche$clusterVector, displayVector=subsetVector)