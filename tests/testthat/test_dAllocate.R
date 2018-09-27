#dAllocate

#IDEA: TAKE A DATASET THAT WE ALLOCATE CORRECTLY TO 99%, THEN CHECK THAT THIS HAPPENS WITH THE RAND INDEX!!!
context("dAllocate")
cols <- 150
x<- DepecheR:::generateBimodalData(observations = 20,dataCols = cols)
centers<-as.data.frame(rbind(runif(cols)+50,runif(cols)-50))
colnames(centers) <- paste('V',1:cols, sep="")
colnames(centers)
out <-DepecheR:::dAllocate(as.data.frame(x$samples), out$clusterCenters, ids=x$ids)
test_that("dAllocate expected output", {
  expect_true(all(out$realloClusterVector==x$ids) )
})
