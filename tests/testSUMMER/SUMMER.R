## CRAN won't allow this test unless INLA is a listed dependency
## which is impossible, so the test is hidden away here.

library(survey)

data("DemoData2",package="SUMMER")
data("DemoMap2", package="SUMMER")

if(require("INLA",quietly=TRUE)){
INLA::inla.setOption(num.threads="1:1")
 library(survey)
 des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
                  weights = ~weights, data = DemoData2, nest = TRUE)
 Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)

 cts.cov.res <- svysmoothArea(tobacco.use ~ age, 
                          domain = ~region,
                          design = des0,
                          adj.mat = DemoMap2$Amat, 
                          X.domain = Xmat,
                          pc.u = 1,
                          pc.alpha = 0.01,
                          pc.u.phi = 0.5,
                          pc.alpha.phi = 2/3)
print(cts.cov.res)
plot(cts.cov.res)
summary(cts.cov.res)
}
