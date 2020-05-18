context("neighSmooth")
focDat <- c(rep(0, 75), rep(1, 25), rep(0, 50), rep(1, 50))
euSpaceDat <- c(rep(0, 100), rep(1, 100))
resDat <- neighSmooth(focusData = focDat, euclidSpaceData = euSpaceDat)
test_that("neighSmooth expected output", {
    expect_true(all(resDat[seq(1, 100)] == 0.25))
})
