library(survey)

data_no_strat <- data.frame(x = 1:10, id = 1:10, fpc = 10, probs = 1)

design_no_strat <- svydesign(id = ~id, 
                             probs = ~probs, 
                             data = data_no_strat, 
                             fpc = ~fpc)

quantiles <- c(0.01, 0.05, 0.1, 0.15, seq(21, 81, 10)*0.01, 0.85, 0.9, 0.95, 0.99)

res <- svyquantile(~x, design_no_strat, quantiles, ci = TRUE,
            interval.type = "mean", qrule = "hf1")

stopifnot(all(diff(res$x[,"quantile"])>=0))


