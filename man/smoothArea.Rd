\name{smoothArea}
\alias{svysmoothArea}
\title{Small area estimation via basic area level model}
\usage{
svysmoothArea(
  formula,
  domain,
  design = NULL,
  adj.mat = NULL,
  X.domain = NULL,
  direct.est = NULL,
  domain.size = NULL,
  transform = c("identity", "logit", "log"),
  pc.u = 1,
  pc.alpha = 0.01,
  pc.u.phi = 0.5,
  pc.alpha.phi = 2/3,
  level = 0.95,
  n.sample = 250,
  var.tol = 1e-10,
  return.samples = FALSE,...
)
}
\arguments{
\item{formula}{An object of class 'formula' describing the model to be fitted.
If direct.est is specified, the right hand side of the formula is not necessary.}

\item{domain}{One-sided formula specifying factors containing domain labels}

\item{design}{An object of class "svydesign" containing the data for the model}

\item{adj.mat}{Adjacency matrix with rownames matching the domain labels. If set to \code{NULL}, the IID spatial effect will be used.}

\item{X.domain}{Data frame of areal covariates. One of the column names needs to match the name of the domain variable, in order to be linked to the data input. Currently only supporting time-invariant covariates.}

\item{direct.est}{Data frame of direct estimates, with first column containing the domain variable, second column containing direct estimate, and third column containing the variance of direct estimate.}

\item{domain.size}{Data frame of domain sizes. One of the column names needs to match the name of the \code{domain} variable, in order to be linked to the data input and there must be a column names 'size' containing domain sizes.}

\item{transform}{Optional transformation applied to the direct estimates before fitting area level model. The default option is no transformation, but logit and log are implemented.}

\item{pc.u}{Hyperparameter U for the PC prior on precisions. See the INLA documentation for more details on the parameterization.}

\item{pc.alpha}{Hyperparameter alpha for the PC prior on precisions.}

\item{pc.u.phi}{Hyperparameter U for the PC prior on the mixture probability phi in BYM2 model.}

\item{pc.alpha.phi}{Hyperparameter alpha for the PC prior on the mixture probability phi in BYM2 model.}

\item{level}{The specified level for the posterior credible intervals}

\item{n.sample}{Number of draws from posterior used to compute summaries}

\item{var.tol}{Tolerance parameter; if variance of an area's direct estimator is below this value, that direct estimator is dropped from model}

\item{return.samples}{If TRUE, return matrix of posterior samples of
  area level quantities}
\item{\dots}{for future methods}
}
\value{
A \code{svysae} object
}
\description{
Generates small area estimates  by smoothing direct estimates using an
area level model
}
\details{
The basic area level model is a Bayesian version of the Fay-Herriot
model (Fay & Herriot,1979). It treats direct estimates of small area quantities as response
data and explicitly models differences between areas using covariate
information and random effects. The Fay-Herriot model can be viewed as a
two-stage model: in the first stage, a sampling model represents the
sampling variability of a direct estimator and in the second stage, a
linking model describes the between area differences in small area
quantities. More detail is given in section 4 of Mercer et al (2015).
  }
  \seealso{The \code{survey-sae} vignette}
  \references{
Fay, Robert E., and Roger A. Herriot. (1979). Estimates of Income for
  Small Places: An Application of James-Stein Procedures to Census
  Data. Journal of the American Statistical Association 74 (366a):
  269-77. 

 Mercer LD, Wakefield J, Pantazis A, Lutambi AM, Masanja H, Clark
    S. Space-Time Smoothing of Complex Survey Data: Small Area
    Estimation for Child Mortality. Ann Appl Stat. 2015
    Dec;9(4):1889-1905. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4959836/}
  }

  \examples{
## artificial data from SUMMER package
## Uses too many cores for a CRAN example

\dontrun{
 hasSUMMER<-tryCatch({
   data("DemoData2",package="SUMMER")
   data("DemoMap2", package="SUMMER")
  }, error=function(e) FALSE)

if (!isFALSE(hasSUMMER)){
 library(survey)
 des0 <- svydesign(ids = ~clustid+id, strata = ~strata,
                  weights = ~weights, data = DemoData2, nest = TRUE)
 Xmat <- aggregate(age~region, data = DemoData2, FUN = mean)

 cts.cov.res <- svysmoothArea(tobacco.use ~ age, 
                          domain = ~region,
                          design = des0,
                          adj.mat = DemoMap2$Amat, 
                          X.domain = Xmat,
                          pc.u = 1,
                          pc.alpha = 0.01,
                          pc.u.phi = 0.5,
                          pc.alpha.phi = 2/3)
 print(cts.cov.res)
 plot(cts.cov.res)
}
}
}  