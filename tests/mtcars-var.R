library(survey)

# don't throw an error on domains of size 1, just return NA
input <- mtcars
input$carb <- factor(input$carb)
design <- svydesign(ids = ~0, weights = NULL, data = input)
svyby(
  ~mpg,
  ~carb,
  design,
  svyvar
)


## same n with na.rm=TRUE as subset(, !is.na)
input$mpg[1]<-NA
design <- svydesign(ids = ~0, weights = NULL, data = input)
stopifnot(all.equal(svyvar(~mpg, design, na.rm=TRUE),
          svyvar(~mpg, subset(design, !is.na(mpg)))))
