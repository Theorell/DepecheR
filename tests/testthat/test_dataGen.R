######################
#generateBimodalData
context('dataGen')

centers <- rbind(c(0,0),c(10,10))
n <- 1000
prop<-0.1
x <- DepecheR:::generateBimodalData(centers,observations = n, prop = prop)
test_that("generateBimodalData expected output simple centers", {
  expect_equal(nrow(x$samples),n)
  expect_equal(ncol(x$samples),2)
  expect_true(all(x$ids[1:n*prop]==1))
  expect_true(all(tail(x$ids,n-n*prop)==2))
  expect_true(all(colMeans(x$samples)<centers[2,]*(1-prop)+0.1))
  expect_true(all(colMeans(x$samples)>centers[2,]*(1-prop)-0.1))
})

x <- DepecheR:::generateBimodalData(observations = n)
test_that("generateBimodalData expected output default centers", {
  expect_equal(nrow(x$samples[x$ids==1,]),n*0.3)
  expect_equal(nrow(x$samples[x$ids==2,]),n*0.7)
  expect_true(mean(x$samples[x$ids==1,])<51)
  expect_true(mean(x$samples[x$ids==1,])>49)
})

x <- DepecheR:::generateSparseData(modeN=3, observations=30)
test_that("generateSparseData expected output default centers", {
  binSamples <- abs(x$samples) < 25
  groundTruth <- t(cbind(replicate(10,x$centers[1,]==0),replicate(10,x$centers[2,]==0),replicate(10,x$centers[3,]==0)))
  expect_true(all(binSamples == groundTruth))
})