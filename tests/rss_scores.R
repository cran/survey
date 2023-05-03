## Example from Rao, Scott, and Skinner 1998 Statistica Sinica
library(survey)
data(myco)
dmyco<-svydesign(id=~1, strata=~interaction(Age,leprosy),weights=~wt,data=myco)
m_full<-svyglm(leprosy~I((Age+7.5)^-2)+Scar, family=quasibinomial, design=dmyco)
m_null<-svyglm(leprosy~I((Age+7.5)^-2), family=quasibinomial, design=dmyco)

stopifnot(isTRUE(all.equal(coef(m_null), c(`(Intercept)`=-4.6, `I((Age + 7.5)^-2)`=-427),tol=1e-2)))

s<-svyscoretest(m_full, ~Scar)
stopifnot(abs(s[1]-10.73)<0.05)

t<-svyscoretest(m_full,~Scar,method="individual")
stopifnot(abs(coef(t)- -32.61)<0.1)
stopifnot(abs(vcov(t)-99.1)<0.1)
 
