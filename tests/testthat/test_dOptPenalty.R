######################
# dOptPenalty
context("dOptPenatly")

x <- DepecheR:::generateBimodalData(observations = 100, dataCols = 10)
x_scaled <- DepecheR:::dScale(as.data.frame(x$samples), center = FALSE)
dOPR <- DepecheR:::dOptPenalty(x_scaled, sampleSize = 100, k = 30, penalties = c(0, 2, 4, 8, 16, 32, 128), maxIter = 10, minARIImprovement = 0.01, optimARI = 0.95, makeGraph = FALSE)

ARI <- dOPR[[2]][, 1]
nClust <- dOPR[[2]][, 2]
bestPenalty <- dOPR[[1]]$bestPenalty
rowNamesdOpr <- row.names(dOPR[[2]])
nClust[which(rowNamesdOpr == bestPenalty)]
ARI[which(rowNamesdOpr == bestPenalty)]

test_that("dOptPenalty expected output", {
  expect_equal(nClust[which(rowNamesdOpr == bestPenalty)], 2)
  expect_equal(ARI[which(rowNamesdOpr == bestPenalty)], 1)
})
