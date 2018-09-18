######################
#dSplsda
context('sPLSDA')
data("testData")
data("testDataDepeche")
data("testDataSNE")
setwd("~/Desktop")
sPLSDAObject <- dSplsda(xYData=testDataSNE$Y, idsVector=testData$ids, 
                        groupVector=testData$label, clusterVector=testDataDepeche$clusterVector)

#Alternative usage
pairingVector <- c(rep(rep(1:29, times=1000), times=2))
xYDataPaired <- as.data.frame(testDataSNE$Y)[39001:nrow(testDataSNE$Y),]
testDataPaired <- testData[39001:nrow(testData),]
clusterVectorPaired <- testDataDepeche$clusterVector[39001:length(testDataDepeche$clusterVector)]
sPLSDAObject <- dSplsda(xYData=xYDataPaired, idsVector=testDataPaired$ids, groupVector=testDataPaired$label, 
                        clusterVector=clusterVectorPaired, pairingVector=pairingVector, name="d_sPLSDAPlot_paired", 
                        groupName1="Stimulation 1", groupName2="Stimulation 2")

#Alternative usage
subsetVector <- sample(1:nrow(testData), size=10000)
testDataSNESubset <- testDataSNE$Y[subsetVector,]
sPLSDAObject <- dSplsda(xYData=testDataSNESubset, idsVector=testData$ids, 
                        groupVector=testData$label, clusterVector=testDataDepeche$clusterVector, displayVector=subsetVector)

#Alternative usage
testDataRows <- sample(1:nrow(testData), size=48500)
sPLSDAObject <- dSplsda(xYData=testDataSNE$Y, idsVector=testData$ids, 
                        groupVector=testData$label, clusterVector=testDataDepeche$clusterVector,
                        testSampleRows=testDataRows)

######################
#dSplsdaPreCalculations
data("testData")
data("testDataDepeche")
data("testDataSNE")
setwd("~/Desktop")
dSplsdaInData <- DepecheR:::dSplsdaPreCalculations(clusterVector=testDataDepeche$clusterVector, idsVector=testData$ids, groupVector=testData$label, groupName1="Stimulation 1", groupName2="Stimulation 2")
