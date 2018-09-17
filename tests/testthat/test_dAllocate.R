#dAllocate
context("dAllocate")

  data(testData)
  testDataSample <- sample(1:nrow(testData), size=48500)
  testDataTrain <- testData[testDataSample,]
  testDataTest <- testData[-testDataSample,]
  
  setwd("~/Desktop")
  x_depeche_train <- depeche(testDataTrain[,2:15], maxIter=20, sampleSizes=1000, ids=testDataTrain$ids)
  
  inDataFramePreScaled <- DepecheR:::dScale(testDataTest[,2:15], scale=FALSE, center="peak")
  #Here, all the data is divided by the standard deviation of the full dataset
  sdInDataFramePreScaled <- sd(as.matrix(inDataFramePreScaled))
  inDataFrameScaled <- inDataFramePreScaled/sdInDataFramePreScaled
  
  x_depeche_test <- DepecheR:::dAllocate(inDataFrameScaled, clusterCenters=x_depeche_train$clusterCenters, ids=testDataTest$ids)
  
#test_that("dAllocate expected output"{})