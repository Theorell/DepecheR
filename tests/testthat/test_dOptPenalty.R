######################
#dOptPenalty
context('dOptPenatly')
x <- DepecheR:::generateBimodalData(observations = 100,dataCols = 10)
x_scaled <- DepecheR:::dScale(as.data.frame(x$samples), center = FALSE)
dOPR <- DepecheR:::dOptPenalty(x_scaled, k=30, maxIter=10, disableWarnings=TRUE, minARI=0.95, makeGraph=FALSE)

ARI <- dOPR[[2]][,1]
nClust <- dOPR[[2]][,2]
bestPenalty <- dOPR[[1]]$bestPenalty
rowNamesdOpr <- row.names(dOPR[[2]])
nClust[which(rowNamesdOpr==bestPenalty)]
ARI[which(rowNamesdOpr==bestPenalty)]

test_that("dOptPenalty expected output", {
  expect_equal(nClust[which(rowNamesdOpr==bestPenalty)],2)
  expect_equal(ARI[which(rowNamesdOpr==bestPenalty)],1)
})