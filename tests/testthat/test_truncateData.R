###################### truncateData
context("truncateData")
x <- c(1000, 4000, -9, 100, 101)
xTruncated <- DepecheR:::truncateData(x, control = c(0, 100), lowQuantile = 0, 
                                      highQuantile = 1)
x[3] <- 1
test_that("truncateData expected output", {
    expect_equal(c(100, 100, 0, 100, 100), xTruncated)
})
