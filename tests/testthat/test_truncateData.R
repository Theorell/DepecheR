######################
# truncateData
context("truncateData")


# First normal version
low <- 0.2
x <- c(1, 2, -5, 10, 20)
xTruncated <- DepecheR:::truncateData(x, lowQuantile = low, highQuantile = 1)
quant <- quantile(x, low)
x[3] <- quant
test_that("truncateData expected output", {
  expect_equal(x, xTruncated)
})

######################
# truncateDataCoFunction
x <- c(1000, 4000, -9, 100, 101)
xTruncated <- DepecheR:::truncateDataCoFunction(x, control = c(0, 100), lowQuantile = 0, highQuantile = 1)
x[3] <- 1
test_that("truncateData expected output", {
  expect_equal(c(100, 100, 0, 100, 100), xTruncated)
})
