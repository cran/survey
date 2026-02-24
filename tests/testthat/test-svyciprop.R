library(testthat)

test_that("svyciprop's fancy confidence intervals are indeed asymmetric", {
  input <- mtcars
  input$carb <- factor(input$carb)
  svy_mtcars <- svydesign(ids = ~0, weights = NULL, data = input)
  
  for (fancy_method in c("logit", "likelihood", "asin", "beta", "xlogit", "wilson")) {
    this_ci <- svyciprop(~am, design=svy_mtcars, method=fancy_method)
    if (this_ci[1]<0.5) {
      # expect the left tail to be shorter than the right
      expect_lt(
        this_ci[1] - confint(this_ci)[1], confint(this_ci)[2] - this_ci[1]
      )
    }
  }
  
  plain_method <- "mean"
  this_ci <- svyciprop(~am, design=svy_mtcars, method=plain_method)
  expect_equal(
    this_ci[1] - confint(this_ci)[1], confint(this_ci)[2] - this_ci[1]
  )

})


