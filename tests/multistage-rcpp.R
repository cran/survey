## temporarily disable while it's being fixed
## q("no")

## 
## Check that multistage() and multistage_rcpp() give same results
## for different options.
##
library(survey)

# Check for a simple random sample ----

  data('api', package = 'survey')
  api_srs_design <- svydesign(
    data = apisrs,
    ids = ~ 1,
    weights = ~ 1
  )
  
  x <- as.matrix(api_srs_design$variables[,c('api00','api99')] / api_srs_design$prob)
  
  base_r_result <- survey:::multistage(x = x,
                                       clusters = api_srs_design$cluster,
                                       stratas = api_srs_design$strata,
                                       nPSUs = api_srs_design$fpc$sampsize, fpcs = api_srs_design$fpc$popsize,
                                       lonely.psu=getOption("survey.lonely.psu"),
                                       one.stage=TRUE, stage = 1, cal = NULL)
  
  rcpp_result <- survey:::multistage_rcpp(x = x,
                                          clusters = api_srs_design$cluster,
                                          stratas = api_srs_design$strata,
                                          nPSUs = api_srs_design$fpc$sampsize, fpcs = api_srs_design$fpc$popsize,
                                          lonely.psu=getOption("survey.lonely.psu"),
                                          one.stage=TRUE, stage = 1, cal = NULL)
  
  if (!isTRUE(all.equal(base_r_result, rcpp_result))) {
    stop("Differences between `multistage()` and `multistage_rcpp()` for SRS")
  }
  
# Check for a stratified simple random sample ----
  
  apistrat_design <- svydesign(
    data = apistrat,
    id =~ 1,
    strata =~ stype,
    weights =~ pw,
    fpc =~ fpc
  )
  
  x <- as.matrix(apistrat_design$variables[,c('api00','api99')] / apistrat_design$prob)
  
  base_r_result <- survey:::multistage(x = x,
                                       clusters = apistrat_design$cluster,
                                       stratas = apistrat_design$strata,
                                       nPSUs = apistrat_design$fpc$sampsize, fpcs = apistrat_design$fpc$popsize,
                                       lonely.psu=getOption("survey.lonely.psu"),
                                       one.stage=TRUE, stage = 1, cal = NULL)
  
  rcpp_result <- survey:::multistage_rcpp(x = x,
                                          clusters = apistrat_design$cluster,
                                          stratas = apistrat_design$strata,
                                          nPSUs = apistrat_design$fpc$sampsize, fpcs = apistrat_design$fpc$popsize,
                                          lonely.psu=getOption("survey.lonely.psu"),
                                          one.stage=TRUE, stage = 1, cal = NULL)
  
  if (!isTRUE(all.equal(base_r_result, rcpp_result))) {
    stop("Differences between `multistage()` and `multistage_rcpp()` for stratified sample")
  }
  
  ##_ Check whether expected true zeroes are actually computed as zeroes ----
  
  x <- as.matrix(model.matrix(~ -1 + stype, data = apistrat_design$variables) / (sum(apistrat_design$prob)/apistrat_design$prob))
  
  base_r_result <- survey:::multistage(x = x,
                                       clusters = apistrat_design$cluster,
                                       stratas = apistrat_design$strata,
                                       nPSUs = apistrat_design$fpc$sampsize, fpcs = apistrat_design$fpc$popsize,
                                       lonely.psu=getOption("survey.lonely.psu"),
                                       one.stage=TRUE, stage = 1, cal = NULL)
  
  rcpp_result <- survey:::multistage_rcpp(x = x,
                                          clusters = apistrat_design$cluster,
                                          stratas = apistrat_design$strata,
                                          nPSUs = apistrat_design$fpc$sampsize, fpcs = apistrat_design$fpc$popsize,
                                          lonely.psu=getOption("survey.lonely.psu"),
                                          one.stage=TRUE, stage = 1, cal = NULL)
  
  if (!all(rcpp_result == 0) | !all(diag(rcpp_result) >= 0)) {
    stop("Computations which should equal zero do not.")
  }

# Create a stratified, multistage design ----
  
  data(mu284)
  
  ## Create three strata, the third of which has only one PSU
  mu284_stratified <- rbind(
    transform(mu284, stratum = 1,
              y2 = y1 + c(0, 2, 1, 0, -1, 2, 0, 3, 0, 0, -1, -1, 1, -1, 0)),
    transform(mu284, stratum = 2,
              y2 = y1 + c(-1, 0, -1, 0, -1, 0, 0, 1, 2, -1, 1, 1, 0, 1, 0)),
    transform(mu284[1,], stratum = 3,
              y2 = y1 + 2)
  )
  
  ## Create domain variables which yield a lonely PSU or a lonely SSU
  mu284_stratified[['DOMAIN_W_LONELY_PSU']] <- mu284_stratified$stratum == 1 | (mu284_stratified$stratum == 2 & mu284_stratified$id1 == 19)
  mu284_stratified[['DOMAIN_W_LONELY_SSU']] <- mu284_stratified$stratum == 1 | (mu284_stratified$stratum == 2 & mu284_stratified$id1 == 19 & mu284_stratified$id2 == 1)
  
  ## Create a survey design object with no lonely PSUs
  dmu284_strat <- svydesign(data = subset(mu284_stratified, stratum != 3),
                            strata = ~ stratum, nest = TRUE,
                            id = ~ id1 + id2, fpc = ~ n1 + n2)
  ## Create a survey design with one lonely PSU
  dmu284_strat_w_lonely <- svydesign(data = mu284_stratified,
                                     strata = ~ stratum, nest = TRUE,
                                     id = ~ id1 + id2, fpc = ~ n1 + n2)
  ## Create subsetted design objects, subsetted to domain with a lonely PSU or SSU
  dmu284_strat_domain_lonely_psu <- subset(dmu284_strat, DOMAIN_W_LONELY_PSU)
  dmu284_strat_domain_lonely_ssu <- subset(dmu284_strat, DOMAIN_W_LONELY_SSU)
  
  
# Check same results for different values of 'lonely PSU' options ----
  
  lonely_psu_options <- c(certainty = 'certainty',
                          remove    = 'remove',
                          average   = 'average',
                          adjust    = 'adjust')
  one_stage_options <- c(TRUE, FALSE)
  
  design_lonely_psu_comparisons <- expand.grid('survey.lonely.psu' = lonely_psu_options,
                                               'one.stage' = one_stage_options,
                                               stringsAsFactors = FALSE)
  design_lonely_psu_comparisons[['results_match']] <- NA
  
  ##_ Check for lonely PSUs caused by design ----
  for (i in seq_len(nrow(design_lonely_psu_comparisons))) {
    survey.lonely.psu <- design_lonely_psu_comparisons[['survey.lonely.psu']][i]
    one.stage <- design_lonely_psu_comparisons[['one.stage']][i]
    options('survey.lonely.psu' = survey.lonely.psu)
    
    x <- as.matrix(dmu284_strat_w_lonely$variables[,c('y1','y2')] / dmu284_strat_w_lonely$prob)
    
    base_r_result <- survey:::multistage(x = x,
                                         clusters = dmu284_strat_w_lonely$cluster,
                                         stratas = dmu284_strat_w_lonely$strata,
                                         nPSUs = dmu284_strat_w_lonely$fpc$sampsize, fpcs = dmu284_strat_w_lonely$fpc$popsize,
                                         lonely.psu=getOption("survey.lonely.psu"),
                                         one.stage=one.stage, stage = 1, cal = NULL)
    
    rcpp_result <- survey:::multistage_rcpp(x = x,
                                            clusters = dmu284_strat_w_lonely$cluster,
                                            stratas = dmu284_strat_w_lonely$strata,
                                            nPSUs = dmu284_strat_w_lonely$fpc$sampsize, fpcs = dmu284_strat_w_lonely$fpc$popsize,
                                            lonely.psu=getOption("survey.lonely.psu"),
                                            one.stage=one.stage, stage = 1, cal = NULL)
    
      design_lonely_psu_comparisons[['results_match']][i] <- isTRUE(all.equal(base_r_result, rcpp_result))
  }
  
  print(design_lonely_psu_comparisons)
  if (!all(design_lonely_psu_comparisons$results_match)) {
    stop("Results for design lonely PSUs differ between base R and Rcpp implementations.")
  }

  ##_ Check for domain lonely PSUs caused by subsetting ----
  
  domain_lonely_psu_comparisons <- expand.grid('survey.lonely.psu' = lonely_psu_options,
                                               'survey.adjust.domain.lonely' = c(FALSE, TRUE),
                                               'one.stage' = one_stage_options,
                                               stringsAsFactors = FALSE)
  domain_lonely_psu_comparisons[['results_match']] <- NA
  
  for (i in seq_len(nrow(domain_lonely_psu_comparisons))) {
      
    options('survey.lonely.psu' = domain_lonely_psu_comparisons[['survey.lonely.psu']][i])
    options('survey.adjust.domain.lonely' = domain_lonely_psu_comparisons[['survey.adjust.domain.lonely']][i])
    one.stage <- domain_lonely_psu_comparisons[['one.stage']][i]
    
    
    x <- as.matrix(dmu284_strat_domain_lonely_psu$variables[,c('y1','y2')] / dmu284_strat_domain_lonely_psu$prob)
    suppressWarnings({
      base_r_result <- survey:::multistage(x = x,
                                           clusters = dmu284_strat_domain_lonely_psu$cluster,
                                           stratas = dmu284_strat_domain_lonely_psu$strata,
                                           nPSUs = dmu284_strat_domain_lonely_psu$fpc$sampsize,
                                           fpcs = dmu284_strat_domain_lonely_psu$fpc$popsize,
                                           lonely.psu=getOption("survey.lonely.psu"),
                                           one.stage=one.stage, stage = 1, cal = NULL)
    })
    
    rcpp_result <- survey:::multistage_rcpp(x = x,
                                            clusters = dmu284_strat_domain_lonely_psu$cluster,
                                            stratas = dmu284_strat_domain_lonely_psu$strata,
                                            nPSUs = dmu284_strat_domain_lonely_psu$fpc$sampsize,
                                            fpcs = dmu284_strat_domain_lonely_psu$fpc$popsize,
                                            lonely.psu=getOption("survey.lonely.psu"),
                                            one.stage=one.stage, stage = 1, cal = NULL)
    
    domain_lonely_psu_comparisons[['results_match']][i] <- isTRUE(all.equal(base_r_result, rcpp_result))
  }
  
  print(domain_lonely_psu_comparisons)
  if (!all(domain_lonely_psu_comparisons$results_match)) {
    stop("Results for domain lonely PSUs differ between base R and Rcpp implementations.")
  }
  
# Checks for edge cases ----
  
  ## Zero records ----
  
  zero_row_design <- dmu284_strat |> subset(id1 == "Nonexistent")
  
  zero_row_result <- survey:::multistage_rcpp(
    x = as.matrix(zero_row_design$variables[,c('y1', 'y2')]),
    clusters = zero_row_design$cluster,
    stratas = zero_row_design$strata,
    nPSUs = zero_row_design$fpc$sampsize, 
    fpcs = zero_row_design$fpc$popsize,
    lonely.psu=getOption("survey.lonely.psu"),
    one.stage=one.stage, stage = 1, cal = NULL
  )
  
  if (any(is.na(zero_row_result)) || any(zero_row_result != 0)) {
    stop("`multistage_rcpp()` should return zeros for zero-row inputs.")
  }
  
