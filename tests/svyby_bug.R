library(survey)
options(warn=2)

## Caused warnings and unhelpful results in 4.1_1 (Guilherme Jacob)
data(api)
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
svyby(~api99, ~stype, dclus1, svymean )

set.seed(123)
apiclus1$api99[ sample.int( nrow(apiclus1) , 5 ) ] <- NA
dclus1.na <-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

# subsetting w/ na.rm = FALSE...
svymean( ~api99 , subset( dclus1.na , stype == "E" ) , na.rm = FALSE )
svymean( ~api99 , subset( dclus1.na , stype == "H" ) , na.rm = FALSE )
svymean( ~api99 , subset( dclus1.na , stype == "M" ) , na.rm = FALSE )

# ... looks like this:
svyby(~api99, ~stype, dclus1.na , svymean )

# subsetting w/ na.rm = TRUE...
svymean( ~api99 , subset( dclus1.na , stype == "E" ) , na.rm = TRUE )
svymean( ~api99 , subset( dclus1.na , stype == "H" ) , na.rm = TRUE )
svymean( ~api99 , subset( dclus1.na , stype == "M" ) , na.rm = TRUE )

# ... looks like this
svyby(~api99, ~stype, dclus1.na , svymean , na.rm = TRUE )

# Without missing values, this works:
svyby(~api99, ~stype, dclus1 , svymean , na.rm = TRUE , covmat = TRUE )

# ... but this breaks!
svyby(~api99, ~stype, dclus1.na , svymean , na.rm = TRUE , covmat = TRUE )

# ... and i don't think this is the expected behavior
svyby( ~api99, ~stype, dclus1.na , svymean , na.rm.all = TRUE , covmat = TRUE )
svyby( ~api99, ~stype, dclus1.na , svymean , na.rm.all = TRUE , na.rm = TRUE , covmat = TRUE )


## <TL> Now some more as tests
svyby(~api99, ~stype, dclus1.na , svytotal , na.rm = TRUE , covmat = TRUE )
svyby(~api99, ~stype, dclus1.na , svyratio , na.rm = TRUE , denominator=~api00, covmat = TRUE )

ff<-function(f,d,...,na.rm=TRUE) svyglm(f,d,...)
svyby(api99~1, ~stype, dclus1.na , ff , na.rm = TRUE , covmat = TRUE )
