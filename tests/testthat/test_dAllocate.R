# dAllocate

context("dAllocate")
cols <- 150
set.seed(191)
x <- DepecheR:::generateBimodalData(observations = 20, dataCols = cols)
centers <- rbind(runif(cols) + 50, runif(cols) - 50)
out <- DepecheR:::dAllocate(x$samples, list("clusterCenters" =centers,
                                 "logCenterSd" = FALSE))
test_that("dAllocate expected output", {
    expect_true(all(out == x$ids))
})

