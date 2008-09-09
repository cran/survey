
library(survey)
library(RSQLite)

data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, fpc=~fpc,data=apiclus1)
dbclus1<-svydesign(id=~dnum, weights=~pw, fpc=~fpc,
data="apiclus1",dbtype="SQLite", dbname=system.file("api.db",package="survey"))

m<-svymean(~api00+stype,dclus1)
m.db<-svymean(~api00+stype, dbclus1)
all.equal(coef(m),coef(m.db))
all.equal(vcov(m), vcov(m.db))

r<-svyratio(~api.stu, ~enroll, design=dclus1)
r.db<-svyratio(~api_stu, ~enroll, design=dbclus1)
all.equal(coef(r), coef(r.db))
all.equal(SE(r), SE(r.db))

b<-svyby(~api99+api00,~stype, design=dclus1, svymean, deff=TRUE)
b.db<-svyby(~api99+api00,~stype, design=dbclus1,svymean, deff=TRUE)
all.equal(coef(b), coef(b.db))
all.equal(SE(b), SE(b.db))
all.equal(deff(b), deff(b.db))

l<-svyglm(api00~api99+mobility, design=dclus1)
l.db<-svyglm(api00~api99+mobility, design=dbclus1)
all.equal(coef(l),coef(l.db))
all.equal(vcov(l), vcov(l.db))

