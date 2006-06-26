\name{twophase}
\alias{twophase}
\alias{twophasevar}
\alias{[.twophase}
\alias{dim.twophase}
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
The internal structure of \code{twophase} objects may change in the future.
}
\usage{
twophase(id, strata = NULL, probs = NULL, weights = NULL, fpc = NULL,
subset, data)
twophasevar(x,design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{list of two formulas for sampling unit identifiers}
  \item{strata}{list of two formulas (or \code{NULL}s) for stratum identifies}
  \item{probs}{ list of two formulas (or \code{NULL}s) for sampling probabilities}
  \item{weights}{list of two formulas (or \code{NULL}s) for sampling weights}
  \item{fpc}{list of two formulas (or \code{NULL}s) for finite
    population corrections}
  \item{subset}{formula specifying which observations are selected in
    phase 2}
  \item{data}{Data frame will all data for phase 1 and 2}
  \item{x}{probability-weighted estimating functions}
  \item{design}{two-phase design}
}
\details{
  The population for the second phase is the first-phase sample. If the
  second phase sample uses stratified (multistage cluster) sampling
  without replacement and all the stratum and sampling unit identifier
  variables are available for the whole first-phase sample it is
  possible to estimate the sampling probabilities/weights and the
  finite population correction. These would then  be specified as \code{NULL}.
  
  Two-phase case-control and case-cohort studies in biostatistics will
  typically have simple random sampling with replacement as the first
  stage. Variances given here may differ slightly from those in the
  biostatistics literature where a model-based estimator of the
  first-stage variance would be used (cf Therneau & Li, 1999)
  
  Variance computations are based on the conditioning argument in
  Section 9.3 of Sarndal et al, but do not use their formulas. In particular,
  the phase-1 variance uses a formula that has bias O(1/phase 2 sample size) but
  is easier to extend to multistage sampling designs. See the \code{tests} 
  directory for a worked example.
  }
\note{
  The variance computations were derived under the assumption that the
  PSUs in the second stage are nested in the ultimate sampling units at
  the first stage. 
}
\value{

  \code{twophase} returns an object of class \code{twophase}, whose
  structure is liable to change without notice. Currently it contains
  three survey design objects for the full phase 1 sample, the subset of
  the phase 1 sample that is in phase 2, and the phase 2 sample.

  \code{twophasevar} returns a variance matrix with an attribute
  containing the separate phase 1 and phase 2 contributions to the variance.
}
\references{
Sarndal CE, Swensson B, Wretman J (1992) "Model Assisted Survey Sampling"
Springer.

Barlow, WE (1994). Robust variance estimation for the case-cohort design. "Biometrics"
50: 1064-1072

Breslow NW and Chatterjee N, Design and analysis of two-phase
studies with binary outcome applied to Wilms tumour prognosis.  "Applied
Statistics"  48:457-68, 1999

Therneau TM and Li H., Computing the Cox model for case-cohort designs. "Lifetime Data Analysis"  5:99-112, 1999

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
 pbc$id<-1:nrow(pbc)
 d2pbc<-twophase(id=list(~id,~id), data=pbc, subset=~I(trt==-9))
 svymean(~bili, d2pbc)
 with(pbc, mean(bili[trt==-9]))
 with(pbc, sd(bili[trt==-9])/sqrt(sum(trt==-9)))

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


}
\keyword{survey}% at least one, from doc/KEYWORDS

