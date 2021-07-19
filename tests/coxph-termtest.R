library(survey)
library(survival)
set.seed(2021-6-25)
test1 <- list(time=c(4,3,1,1,2,2,3),
              status=c(1,1,1,0,1,1,0),
              x1=as.factor(rbinom(7, 2, 0.5)),
              x=c(0,2,1,1,1,0,0))
# Fit a stratified model
mod_c <- coxph(Surv(time, status) ~ x1 + x, test1)
mod_d <- coxph(Surv(time, status) ~ x + x1, test1)
stopifnot(all.equal(regTermTest(mod_c, ~x1, df = Inf)[c("chisq","df","test.terms","p")],
regTermTest(mod_d, ~x1, df = Inf)[c("chisq","df","test.terms","p")]))
