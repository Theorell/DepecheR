###################### dViolins
context("dVviolins")
data("testData")
data("testDataDepeche")
dViolins(testDataDepeche$clusterVector, testDataDepeche$essenceElementList, 
         inDataFrame=testData, createOutput = FALSE)

