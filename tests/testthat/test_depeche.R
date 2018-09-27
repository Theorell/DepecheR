######################
#depeche
context('depeche')
setwd("~/Desktop")
x<- DepecheR:::generateBimodalData(observations = 20)
out <- depeche(x$samples, maxIter=100)
#testDataDepecheResult <- depeche(testData[,2:15], sampleSizes=500, maxIter=20)

#Alternative usage
#testDataDepecheResultDual <- depeche(testData[,2:15], dualDepecheSetup=data.frame(rep(1:2, each=7), 
 #                                                                                 colnames(testData[,2:15])), penalties=c(64, 128), sampleSizes=500, selectionSampleSize=500, 
#                                     maxIter=20, ids=testData$ids)

######################
#depecheAllData
#x <- DepecheR:::generateBimodalData()[,-1]
#x_depeche_all <- DepecheR:::depecheAllData(x, 1)

######################
#depecheCoFunction
#x <- DepecheR:::generateBimodalData()[,-1]
#x_scaled <- DepecheR:::dScale(x)
#depecheCoResult <- DepecheR:::depecheCoFunction(inDataFrameScaled=x_scaled, firstClusterNumber=1, penalties=c(2,4), sampleSizes=500, selectionSampleSize=500, k=30, minARIImprovement=0.01, minARI=0.95, maxIter=20)
