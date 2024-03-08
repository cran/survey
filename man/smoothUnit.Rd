\name{smoothUnit}
\alias{svysmoothUnit}
\title{Smooth via basic unit level model}
\usage{
svysmoothUnit(
  formula,
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
  level = 0.95,
  n.sample = 250,
  return.samples = FALSE,
  X.pop.weights = NULL,...
)
}
\arguments{
\item{formula}{An object of class 'formula' describing the model to be fitted.}

\item{domain}{One-sided formula specifying factors containing domain labels}

\item{design}{An object of class "survey.design" containing the data for the model}

\item{family}{of the response variable, currently supports 'binomial' (default with logit link function) or 'gaussian'.}

\item{X.pop}{Data frame of population unit-level covariates. One of the column name needs to match the domain specified, in order to be linked to the data input. Currently only supporting time-invariant covariates.}

\item{adj.mat}{Adjacency matrix with rownames matching the domain labels. If set to NULL, the IID spatial effect will be used.}

\item{domain.size}{Data frame of domain sizes. One of the column names needs to match the name of the domain variable, in order to be linked to the data input and there must be a column names 'size' containing domain sizes. The default option is no transformation, but logit and log are implemented.}

\item{pc.u}{Hyperparameter U for the PC prior on precisions. See the INLA documentation for more details on the parameterization.}

\item{pc.alpha}{Hyperparameter alpha for the PC prior on precisions.}

\item{pc.u.phi}{Hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.}

\item{pc.alpha.phi}{Hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.}

\item{level}{The specified level for the posterior credible intervals}

\item{n.sample}{Number of draws from posterior used to compute
  summaries}

\item{return.samples}{If TRUE, return matrix of posterior samples of area level quantities}

\item{X.pop.weights}{Optional vector of weights to use when aggregating
  unit level predictions}
\item{\dots}{for future expansion}
}
\value{
A \code{svysae} object
}
\description{
Generates small area estimates by smoothing direct estimates using a basic
unit level model. This model assumes sampling is ignorable (no selection
  bias). It's a Bayesian linear (\code{family="gaussian"}) or generalised linear
  (\code{family="binomial"}) mixed model for the unit-level data with
  individual-level covariates and area-level random effects.
}
\seealso{The \code{survey-sae} vignette}
\references{Battese, G. E., Harter, R. M., & Fuller, W. A. (1988). An Error-Components Model for Prediction of County Crop Areas Using Survey and Satellite Data. Journal of the American Statistical Association, 83(401), 28-36.
}