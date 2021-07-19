library(survey)

# Dummy data for testing nonresponse adjustments
toy <- data.frame(id = 1:10,
                  dummy = 1,
                  class = c(rep("A", 5), rep("B", 5)),
                  responded = c(rep(TRUE, 8), FALSE, FALSE),
                  weight = rep(100, 10))

# With jackknife replicate weights
toy_repweights <- matrix(rep(1000/9, 100), nrow = 10)
diag(toy_repweights) <- rep(0, 10)

# Scramble up which person is in which jackknife group
toy_repweights <- toy_repweights[sample(1:10, size = 10), ]

toy_design <- svrepdesign(variables = toy[, 1:4],
                          weights = toy$weight,
                          repweights = toy_repweights,
                          type = "JK1",
                          scale = 0.9)

# Get the sum of the weights for the full sample
poptotals <- as.data.frame(svyby(formula = ~dummy,
                                 by = ~class,
                                 design = toy_design,
                                 FUN = survey::svytotal))

poptotals$Freq <- poptotals$dummy
poptotals$dummy <- NULL
poptotals$se <- NULL

# Adjust the weights of the responding sample to match to the full sample
adjusted <-  postStratify(design = subset(toy_design, responded),
                          strata = ~class,
                          population = poptotals)

# This works for the weights...
svyby(formula = ~dummy,
      by = ~class,
      design = adjusted,
      FUN = survey::svytotal)

# ...and some, but not all, the replicate weights
stopifnot(all.equal(colSums(toy_design$repweights), colSums(adjusted$repweights)))


