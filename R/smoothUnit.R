svysmoothUnit <- function(formula,
                       domain,
                       design,
                       family = c("gaussian", "binomial"),
                       X.pop = NULL,
                       adj.mat = NULL,
                       domain.size = NULL,
                       pc.u = 1,
                       pc.alpha = 0.01,
                       pc.u.phi = 0.5,
                       pc.alpha.phi = 2/3,
                       level = .95, n.sample = 250,
                       return.samples = FALSE,
                       X.pop.weights = NULL,...) {
    UseMethod("svysmoothUnit",design)
}

svysmoothUnit.survey.design <- function(formula,
                       domain,
                       design,
                       family = c("gaussian", "binomial"),
                       X.pop = NULL,
                       adj.mat = NULL,
                       domain.size = NULL,
                       pc.u = 1,
                       pc.alpha = 0.01,
                       pc.u.phi = 0.5,
                       pc.alpha.phi = 2/3,
                       level = .95, n.sample = 250,
                       return.samples = FALSE,
                       X.pop.weights = NULL,...) {
  if (!requireNamespace("SUMMER", quietly = TRUE)) {
    stop(
      "Package \"SUMMER\" must be installed to use this function.",
      call. = FALSE
    )
  }
  family<-match.arg(family)
  rval<-SUMMER::smoothUnit(formula = formula,
                            domain = domain,
                            design = design,
                            family = family,
                            X.pop = X.pop,
                            adj.mat = adj.mat,
                            domain.size = domain.size,
                            pc.u = pc.u,
                            pc.alpha = pc.alpha,
                            pc.u.phi = pc.u.phi,
                            pc.alpha.phi = pc.alpha.phi,
                            level = level,
                            n.sample = n.sample,
                            return.samples = return.samples,
                           X.pop.weights = X.pop.weights)
  rval$call<-sys.call(-1)
  rval
}


svysmoothUnit.default <- function(formula,
                       domain,
                       design,
                       family = c("gaussian", "binomial"),
                       X.pop = NULL,
                       adj.mat = NULL,
                       domain.size = NULL,
                       pc.u = 1,
                       pc.alpha = 0.01,
                       pc.u.phi = 0.5,
                       pc.alpha.phi = 2/3,
                       level = .95, n.sample = 250,
                       return.samples = FALSE,
                       X.pop.weights = NULL,...) {

stop("svysmoothUnit is only available for survey.design objects")


}
