\name{twophase}
\alias{twophase}
\alias{twophasevar}
\alias{twophase2var}
\alias{[.twophase}
\alias{subset.twophase}
\alias{print.twophase}
\alias{summary.twophase}
\alias{print.summary.twophase}
\alias{model.frame.twophase}
\alias{na.fail.twophase}
\alias{na.omit.twophase}
\alias{na.exclude.twophase}
\alias{svyrecvar.phase1}
\alias{multistage.phase1}
\alias{onestage.phase1}
\alias{onestrat.phase1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Two-phase designs}
\description{
In a two-phase design a sample is taken from a population and a
subsample taken from the sample, typically stratified by variables not
known for the whole population.  The second phase can use any design 
supported for single-phase sampling. The first phase must currently
be one-stage element  or cluster sampling
}
\usage{
twophase(id, strata = NULL, probs = NULL, weights = NULL, fpc = NULL,
subset, data, method=c("full","approx","simple"))
twophasevar(x,design)
twophase2var(x,design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{list of two formulas for sampling unit identifiers}
  \item{strata}{list of two formulas (or \code{NULL}s) for stratum identifies}
  \item{probs}{ list of two formulas (or \code{NULL}s) for sampling probabilities}
  \item{weights}{Only for \code{method="approx"}, list of two formulas (or \code{NULL}s) for sampling weights}
  \item{fpc}{list of two formulas (or \code{NULL}s) for finite
    population corrections}
  \item{subset}{formula specifying which observations are selected in
    phase 2}
  \item{data}{Data frame will all data for phase 1 and 2}
  \item{method}{\code{"full"} requires (much) more memory, but gives unbiased
    variance estimates for general multistage designs at both phases.
    \code{"simple"} or \code{"approx"} uses the standard error calculation from
    version 3.14 and earlier,  which uses much less memory and is correct for designs with simple
    random sampling at phase one and stratified random sampling at phase two.
  }
  \item{x}{probability-weighted estimating functions}
  \item{design}{two-phase design}
}
\details{
  The population for the second phase is the first-phase sample. If the
  second phase sample uses stratified (multistage cluster) sampling
  without replacement and all the stratum and sampling unit identifier
  variables are available for the whole first-phase sample it is
  possible to estimate the sampling probabilities/weights and the
  finite population correction. These would then  be specified as
  \code{NULL}.

  Two-phase case-control and case-cohort studies in biostatistics will
  typically have simple random sampling with replacement as the first
  stage. Variances given here may differ slightly from those in the
  biostatistics literature where a model-based estimator of the
  first-stage variance would typically be used.
  
  Variance computations are based on the conditioning argument in
  Section 9.3 of Sarndal et al. Method \code{"full"} corresponds exactly
  to the formulas in that reference. Method \code{"simple"} or
  \code{"approx"} (the two are the same) uses less time and memory but
  is exact only for some special cases. The most important special case
  is the two-phase epidemiologic designs where phase 1 is simple random
  sampling from an infinite population and phase 2 is stratified random
  sampling.  See the \code{tests} directory for a worked example. The
  only disadvantage of method="simple" in these cases is that
  standardization of margins (\code{\link{marginpred}}) is not available.


  For \code{method="full"}, sampling probabilities must be available for
  each stage of sampling, within each phase.  For multistage sampling
  this requires specifying either \code{fpc} or \code{probs} as a
  formula with a term for each stage of sampling.  If no \code{fpc} or
  \code{probs} are specified at phase 1 it is treated as simple random
  sampling from an infinite population, and population totals will not
  be correctly estimated, but means, quantiles, and regression models
  will be correct.
  }

\value{

  \code{twophase} returns an object of class \code{twophase2} (for
  \code{method="full"}) or \code{twophase}.  The structure of
  \code{twophase2} objects may change as unnecessary components are removed.

  \code{twophase2var} and \code{twophasevar} return a variance matrix with an attribute
  containing the separate phase 1 and phase 2 contributions to the variance.
}
\references{
Sarndal CE, Swensson B, Wretman J (1992) "Model Assisted Survey Sampling"
Springer.


Breslow NE and Chatterjee N, Design and analysis of two-phase
studies with binary outcome applied to Wilms tumour prognosis.  "Applied
Statistics"  48:457-68, 1999

Breslow N, Lumley T, Ballantyne CM, Chambless LE, Kulick M. (2009)
Improved Horvitz-Thompson estimation of model parameters from two-phase
stratified samples: applications in epidemiology. Statistics in
Biosciences. doi 10.1007/s12561-009-9001-6


Lin, DY and Ying, Z (1993). Cox regression with incomplete covariate measurements.
"Journal of the American Statistical Association" 88: 1341-1349.
}

\seealso{\code{\link{svydesign}}, \code{\link{svyrecvar}} for multi*stage*
  sampling

  \code{\link{calibrate}} for calibration (GREG) estimators.

  \code{\link{estWeights}} for two-phase designs for missing data.

  The "epi" and "phase1" vignettes for examples and technical details.
}
\examples{
 ## two-phase simple random sampling.
 data(pbc, package="survival")
 pbc$randomized<-with(pbc, !is.na(trt) & trt>0)
 pbc$id<-1:nrow(pbc)
 d2pbc<-twophase(id=list(~id,~id), data=pbc, subset=~randomized)
 svymean(~bili, d2pbc)

 ## two-stage sampling as two-phase
 data(mu284)
 ii<-with(mu284, c(1:15, rep(1:5,n2[1:5]-3)))
 mu284.1<-mu284[ii,]
 mu284.1$id<-1:nrow(mu284.1)
 mu284.1$sub<-rep(c(TRUE,FALSE),c(15,34-15))
 dmu284<-svydesign(id=~id1+id2,fpc=~n1+n2, data=mu284)
 ## first phase cluster sample, second phase stratified within cluster
 d2mu284<-twophase(id=list(~id1,~id),strata=list(NULL,~id1),
                     fpc=list(~n1,NULL),data=mu284.1,subset=~sub)
 svytotal(~y1, dmu284)
 svytotal(~y1, d2mu284)
 svymean(~y1, dmu284)
 svymean(~y1, d2mu284)

 ## case-cohort design: this example requires R 2.2.0 or later
 library("survival")
 data(nwtco)

 ## stratified on case status
 dcchs<-twophase(id=list(~seqno,~seqno), strata=list(NULL,~rel),
         subset=~I(in.subcohort | rel), data=nwtco)
 svycoxph(Surv(edrel,rel)~factor(stage)+factor(histol)+I(age/12), design=dcchs)

 ## Using survival::cch 
 subcoh <- nwtco$in.subcohort
 selccoh <- with(nwtco, rel==1|subcoh==1)
 ccoh.data <- nwtco[selccoh,]
 ccoh.data$subcohort <- subcoh[selccoh]
 cch(Surv(edrel, rel) ~ factor(stage) + factor(histol) + I(age/12), data =ccoh.data,
        subcoh = ~subcohort, id=~seqno, cohort.size=4028, method="LinYing")


 ## two-phase case-control
 ## Similar to Breslow & Chatterjee, Applied Statistics (1999) but with
 ## a slightly different version of the data set
 
 nwtco$incc2<-as.logical(with(nwtco, ifelse(rel | instit==2,1,rbinom(nrow(nwtco),1,.1))))
 dccs2<-twophase(id=list(~seqno,~seqno),strata=list(NULL,~interaction(rel,instit)),
    data=nwtco, subset=~incc2)
 dccs8<-twophase(id=list(~seqno,~seqno),strata=list(NULL,~interaction(rel,stage,instit)),
    data=nwtco, subset=~incc2)
 summary(glm(rel~factor(stage)*factor(histol),data=nwtco,family=binomial()))
 summary(svyglm(rel~factor(stage)*factor(histol),design=dccs2,family=quasibinomial()))
 summary(svyglm(rel~factor(stage)*factor(histol),design=dccs8,family=quasibinomial()))

 ## Stratification on stage is really post-stratification, so we should use calibrate()
 gccs8<-calibrate(dccs2, phase=2, formula=~interaction(rel,stage,instit))
 summary(svyglm(rel~factor(stage)*factor(histol),design=gccs8,family=quasibinomial()))

 ## For this saturated model calibration is equivalent to estimating weights.
 pccs8<-calibrate(dccs2, phase=2,formula=~interaction(rel,stage,instit), calfun="rrz")
 summary(svyglm(rel~factor(stage)*factor(histol),design=pccs8,family=quasibinomial()))

 ## Since sampling is SRS at phase 1 and stratified RS at phase 2, we
 ## can use method="simple" to save memory.
 dccs8_simple<-twophase(id=list(~seqno,~seqno),strata=list(NULL,~interaction(rel,stage,instit)),
    data=nwtco, subset=~incc2,method="simple")
 summary(svyglm(rel~factor(stage)*factor(histol),design=dccs8_simple,family=quasibinomial()))



}
\keyword{survey}% at least one, from doc/KEYWORDS

