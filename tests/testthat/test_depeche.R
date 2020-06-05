###################### depeche

context("depeche")
x <- DepecheR:::generateBimodalData(observations = 100, dataCols = 150)
out <- depeche(x$samples, maxIter = 8, nCores = 2, createOutput = FALSE)
test_that("depeche expected output", {
    expect_true(max(out$clusterVector) == 2)
    expect_true(all(x$ids == out$clusterVector) ||
        all(x$ids %% 2 + 1 == out$clusterVector))
})

# sparsity test

x <- DepecheR:::generateSparseData(observations = 500)
out <- depeche(x$samples,
    maxIter = 8, createOutput = FALSE, nCores = 2,
    center = FALSE
)
binCenters <- x$centers == 0
binClCenters <- out$clusterCenters == 0
test_that("depeche expected sparsity", {
    for (i in 1:5) {
        expect_true(all(binCenters[i, ] ==
            binClCenters[which(rowSums(binClCenters) == i), ]))
    }
})
