######################
#depeche
context('depeche')
x<- DepecheR:::generateBimodalData(observations = 20,dataCols = 150)
out <- depeche(x$samples, maxIter=8, createOutput=FALSE)
test_that("depeche expected output", {
  expect_true(max(out$clusterVector)==2)
  expect_true(all(x$ids ==out$clusterVector) ||all(x$ids %% 2 +1 ==out$clusterVector))
})
