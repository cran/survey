## From Ben Schneider, https://github.com/bschneidr/r-forge-survey-mirror/pull/7
library(survey)
data('api', package = 'survey')
  
  boot_design <- svydesign(
    ids = ~ 1, strata = ~ stype,
    weights = ~ pw,
    data = apistrat,
  ) |> as.svrepdesign(type = "boot")
  
# Attempt to estimate variance of quantile using direct replication ----
  new <- svyquantile(
    x = ~ api00 + api99,
    quantiles = c(0.25, 0.75),
    design = boot_design,
    interval.type = "quantile",
    return.replicates = TRUE
  )
  
print(new)

  old <- oldsvyquantile(
    x = ~ api00 + api99,
    quantiles = c(0.25, 0.75),
    design = boot_design,
    interval.type = "quantile",
    return.replicates = TRUE
  )
  
print(old)
  
confint(new)
confint(old)
