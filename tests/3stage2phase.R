## simulated data with three-stage sample at phase 1, SRS at phase 2
## motivated by dietary biomarker substudy in HCHS

library(survey)

load("simdata1.RData")

twophase.full = twophase(id=list(~block+house+ind,~1),
                         strata=list(~strat,NULL),
                         probs=list(~P.block+P.house+P.ind,NULL),
                         subset=~phase2,
                         data=simdata1,method='full')

twophase.approx = twophase(id=list(~block+house+ind,~1),
                         strata=list(~strat,NULL),
                         probs=list(~P.block+P.house+P.ind,NULL),
                         subset=~phase2,
                         data=simdata1,method='approx')

twophase.rep = twophase(id=list(~block,~1),
                         strata=list(~strat,NULL),
                         probs=list(~I(P.block*P.house*P.ind),NULL),
                         subset=~phase2,
                         data=simdata1,method='full')


twophase.repapprox = twophase(id=list(~block,~1),
                         strata=list(~strat,NULL),
                         probs=list(~I(P.block*P.house*P.ind),NULL),
                         subset=~phase2,
                         data=simdata1,method='approx')


svymean(~age, twophase.full)
svymean(~age, twophase.approx)
svymean(~age, twophase.rep)
svymean(~age, twophase.repapprox)
