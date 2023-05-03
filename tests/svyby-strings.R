
## check that stringsAsFactors fixes the string levels problem
data(api, package = "survey")
des <- survey::svydesign(id = ~dnum, weights = ~pw, data = apiclus1, fpc = ~fpc)
est0 <-
  survey::svyby(design=des, formula=~cname, by=~both, FUN=survey::svymean, keep.var=TRUE, stringsAsFactors=TRUE) 

stopifnot(isTRUE(all(dim(est0)==c(2,23))))
