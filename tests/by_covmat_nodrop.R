library(survey)
load("db.rda")
d<-svydesign(ids=~psu,data=db,weights=~pweight,fpc=~psu_size)
pop.sex<-data.frame(sex=c("男","女"),Freq=c(4592,4190))
pop.hos<-data.frame(hospital=c("hos1","hos2",
                               "hos3","hos4",
                               "hos5","hos6",
                               "hos7"),
                    Freq=c(1796+383,2917,1805,316,420,389,756))
d<-rake(d,sample.margins=list(~sex,~hospital),population.margins=list(pop.sex,pop.hos))
## with option covmat = TRUE, 错误于inflmats[[i]][idxs[[i]], ] <- infs[[i]]: 被替换的项目不是替换值长度的倍数
svyby(formula = ~FT_num+DMFT+ft_num+dmft, by=~sex,design = d,
      FUN=svytotal,na.rm=TRUE,covmat = TRUE,deff=TRUE)


## now check correctness
data(api) 
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
calstrat<-calibrate(dstrat,~stype+0, pop=c(4421,755,1018))
svyby(~enroll, ~comp.imp, FUN=svytotal, covmat=TRUE,deff=TRUE,dstrat)->a
svyby(~enroll, ~comp.imp, FUN=svytotal, covmat=TRUE,deff=TRUE,calstrat)->b

svycontrast(a,c(1,-1))->delta_a
svycontrast(b,c(1,-1))->delta_b

stopifnot(isTRUE(all.equal(coef(delta_a),coef(delta_b),tol=1e-7)))
stopifnot(isTRUE(all.equal(SE(delta_a),SE(delta_b),tol=1e-7)))
