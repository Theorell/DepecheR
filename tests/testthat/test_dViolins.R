###################### dViolins
context("dViolins")
data("testData")
data("testDataDepeche")
dViolins(clusterVector=testDataDepeche$clusterVector, inDataFrame=testData,
         plotClusters=1, plotElements=testDataDepeche$essenceElementList, 
         createOutput = FALSE)

