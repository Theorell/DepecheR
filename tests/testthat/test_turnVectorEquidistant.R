###################### turnVectorEquidistant
context("turnVectorEquidistant")

x <- c(1, 3, 3, -3, 1000, 0.2)
solution <- c(3, 4, 4, 1, 5, 2)
translation <- 5
test_that("turnVectorEquidistant expected output", {
    expect_equal(solution, DepecheR:::turnVectorEquidistant(x))
    expect_equal(solution + translation - 1, DepecheR:::turnVectorEquidistant(x, startValue = translation))
})
