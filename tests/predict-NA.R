## predict.svyglm with rank-deficient model
## from Terry Therneau, 6/24

load("sdata2.rda") 
library(survival)
library(survey)

pfit2 <- svyglm(pseudo ~ trt + flt3 + factor(time), design= sdata2)
pfit3 <- svyglm(pseudo ~ trt + flt3 + factor(time), design=sdata2,
                family= gaussian(link= blog()))
pfit2b <- update(pfit2, . ~ . + flt3:time)
pfit3b <- update(pfit3, . ~ . + flt3:time)
dummy <- expand.grid(trt=c("A", "B"), flt3= LETTERS[1:3], time=24* 1:3)

phat2b <- predict(pfit2b, newdata=dummy, type='response',se.fit=FALSE)
phat3b <- predict(pfit3b, newdata=dummy, type='response',se.fit=TRUE)

stopifnot(class(phat2b)=="numeric")
stopifnot(class(phat3b)=="svystat")

stopifnot(all.equal(as.vector(coef(phat3b))[1],12.3296615700191))
stopifnot(all.equal(as.vector(phat2b)[1], 11.4568358310328))
