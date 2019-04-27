\name{withPV.survey.design}
\alias{withPV.survey.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Analyse plausible values in surveys
}
\description{
Repeats an analysis for each of a set of 'plausible values' in a survey data set, returning a list suitable for \code{mitools::MIcombine}. The default method works for both standard and replicate-weight designs but not for two-phase designs. 
}
\usage{
\S3method{withPV}{survey.design}(mapping, data, action, rewrite=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mapping}{
A formula or list of formulas describing each variable in the analysis that has plausible values. The left-hand side of the formula is the name to use in the analysis; the right-hand side gives the names in the dataset.
}
  \item{data}{
A survey design object, as created by \code{svydesign} or \code{svrepdesign}
}
  \item{action}{
With \code{rewrite=TRUE}, a function taking a survey design object as
its only argument, or a quoted expression.  With \code{rewrite=TRUE}
a function taking a survey design object as its only argument, or a
quoted expression with \code{.DESIGN} referring to the survey design object to be used.
}
\item{rewrite}{
Rewrite \code{action} before evaluating it (versus constructing new data
sets)
  }
  \item{\dots}{
For methods
}
}
\value{
A list of the results returned by each evaluation of \code{action}, with the call as an attribute.
}


\seealso{
\code{\link{with.svyimputationList}}
}
\examples{
if(require(mitools)){
data(pisamaths, package="mitools")
des<-svydesign(id=~SCHOOLID+STIDSTD, strata=~STRATUM, nest=TRUE,
	weights=~W_FSCHWT+condwt, data=pisamaths)

oo<-options(survey.lonely.psu="remove")

results<-withPV(list(maths~PV1MATH+PV2MATH+PV3MATH+PV4MATH+PV5MATH),
   data=des,
   action=quote(svyglm(maths~ST04Q01*(PCGIRLS+SMRATIO)+MATHEFF+OPENPS, design=des)),
   rewrite=TRUE)

summary(MIcombine(results))
options(oo)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}% use one of  RShowDoc("KEYWORDS")