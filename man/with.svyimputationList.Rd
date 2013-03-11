\name{with.svyimputationList}
\alias{with.svyimputationList}
\alias{subset.svyimputationList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analyse multiple imputations}
\description{
Performs a survey analysis on each of the designs in a
\code{svyimputationList} objects and returns a list of results suitable
for \code{MIcombine}. The analysis may be specified as an expression or
as a function.
}
\usage{
\method{with}{svyimputationList}(data, expr, fun, ...,multicore=getOption("survey.multicore"))
\method{subset}{svyimputationList}(x, subset,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data,x}{A \code{svyimputationList} object }
  \item{expr}{An expression giving a survey analysis}
  \item{fun}{A function taking a survey design object as its argument }
  \item{\dots}{for future expansion }
  \item{multicore}{Use \code{multicore} package to distribute imputed data sets over multiple processors?}
  \item{subset}{An logical expression specifying the subset}
 
}

\value{
A list of the results from applying the analysis to each design object.
}

\seealso{\code{MIcombine}, in the \code{mitools} package }
\examples{
library(mitools)
data.dir<-system.file("dta",package="mitools")
files.men<-list.files(data.dir,pattern="m.\\\\.dta$",full=TRUE)
men<-imputationList(lapply(files.men, foreign::read.dta))
files.women<-list.files(data.dir,pattern="f.\\\\.dta$",full=TRUE)
women<-imputationList(lapply(files.women, foreign::read.dta))
men<-update(men, sex=1)
women<-update(women,sex=0)
all<-rbind(men,women)

designs<-svydesign(id=~id, strata=~sex, data=all)
designs

results<-with(designs, svymean(~drkfre))

MIcombine(results)

summary(MIcombine(results))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey }% __ONLY ONE__ keyword per line
