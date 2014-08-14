\name{svydesign}
\alias{svydesign}
\alias{svydesign.default}
\alias{svydesign.imputationList}
\alias{svydesign.character}
\alias{na.omit.survey.design}
\alias{na.exclude.survey.design}
\alias{na.fail.survey.design}

\title{Survey sample analysis.}
\description{
  Specify a complex survey design.
}
\usage{
svydesign(ids, probs=NULL, strata = NULL, variables = NULL, fpc=NULL,
data = NULL, nest = FALSE, check.strata = !nest, weights=NULL,pps=FALSE,...)
\method{svydesign}{default}(ids, probs=NULL, strata = NULL, variables = NULL,
  fpc=NULL,data = NULL, nest = FALSE, check.strata = !nest, weights=NULL,
  pps=FALSE,variance=c("HT","YG"),...)
\method{svydesign}{imputationList}(ids, probs = NULL, strata = NULL, variables = NULL, 
    fpc = NULL, data, nest = FALSE, check.strata = !nest, weights = NULL, pps=FALSE,
     ...)
\method{svydesign}{character}(ids, probs = NULL, strata = NULL, variables = NULL, 
    fpc = NULL, data, nest = FALSE, check.strata = !nest, weights = NULL, pps=FALSE,
    dbtype = "SQLite", dbname, ...)
}
%- maybe also `usage' for other objects documented here.'`
\arguments{
  \item{ids}{Formula or data frame specifying cluster ids from largest
    level to smallest level, \code{~0} or \code{~1} is a formula for no clusters.}
  \item{probs}{Formula or data frame specifying cluster sampling probabilities}
  \item{strata}{Formula or vector specifying strata, use \code{NULL} for no strata}
  \item{variables}{Formula or data frame specifying the variables
    measured in the survey. If \code{NULL}, the \code{data} argument is
    used.}
  \item{fpc}{Finite population correction: see Details below}
  \item{weights}{Formula or vector specifying sampling weights as an
    alternative to \code{prob}}
  \item{data}{Data frame to look up variables in the formula
    arguments, or database table name, or \code{imputationList} object, see below}
  \item{nest}{If \code{TRUE}, relabel cluster ids to enforce nesting
    within strata}
  \item{check.strata}{If \code{TRUE}, check that clusters are nested in
    strata}.
  \item{pps}{\code{"brewer"} to use Brewer's approximation for PPS 
    sampling without replacement. \code{"overton"} to use
    Overton's approximation. An object of class \code{\link{HR}} to use the Hartley-Rao approximation. An
    object of class \code{\link{ppsmat}} to use the Horvitz-Thompson estimator.}
  \item{dbtype}{name of database driver to pass to \code{dbDriver}}
  \item{dbname}{name of database (eg file name for SQLite)}
   \item{variance}{For \code{pps} without replacement, use \code{variance="YG"} for the Yates-Grundy estimator instead of the Horvitz-Thompson estimator}
  \item{\dots}{for future expansion}
}
\details{
   The \code{svydesign} object combines a data frame and all the survey
   design information needed to analyse it.  These objects are used by
   the survey modelling and summary functions.  The
   \code{id} argument is always required, the \code{strata},
   \code{fpc}, \code{weights} and \code{probs} arguments are
   optional. If these variables are specified they must not have any
   missing values.
   
   By default, \code{svydesign} assumes that all PSUs, even those in
   different strata, have a unique value of the \code{id}
   variable. This allows some data errors to be detected. If your PSUs
   reuse the same identifiers across strata then set \code{nest=TRUE}.
   
   
   The finite population correction (fpc) is used to reduce the variance when
   a substantial fraction of the total population of interest has been
   sampled. It may not be appropriate if the target of inference is the
   process generating the data rather than the statistics of a
   particular finite population.
   
   The finite population correction can be specified either as the total
   population size in each stratum or as the fraction of the total
   population that has been sampled. In either case the relevant
   population size is the sampling units.  That is, sampling 100 units
   from a population stratum of size 500 can be specified as 500 or as
   100/500=0.2.  The exception is for PPS sampling without replacement, where the
   sampling probability (which will be different for each PSU) must be used.
   
   If population sizes are specified but not sampling probabilities or
   weights, the sampling probabilities will be computed from the
   population sizes assuming simple random sampling within strata. 
   
   For multistage sampling the \code{id} argument should specify a
   formula with the cluster identifiers at each stage.  If subsequent
   stages are stratified \code{strata} should also be specified as a
   formula with stratum identifiers at each stage.  The population size
   for each level of sampling should also be specified in \code{fpc}.
   If \code{fpc} is not specified then sampling is assumed to be with
   replacement at the top level and only the first stage of cluster is
   used in computing variances. If \code{fpc} is specified but for fewer
   stages than \code{id}, sampling is assumed to be complete for
   subsequent stages.   The variance calculations for
   multistage sampling assume simple or stratified random sampling
   within clusters at each stage except possibly the last.

   For PPS sampling without replacement it is necessary to specify the
   probabilities for each stage of sampling using the \code{fpc}
   arguments, and an overall \code{weight} argument should not be
   given. At the moment, multistage or stratified PPS sampling without
   replacement is supported only with \code{pps="brewer"}, or by
   giving the full joint probability matrix using
   \code{\link{ppsmat}}. [Cluster sampling is supported by all
   methods, but not subsampling within clusters].  
   
   The \code{dim}, \code{"["}, \code{"[<-"} and na.action methods for
   \code{survey.design} objects operate on the dataframe specified by
   \code{variables} and ensure that the design information is properly
   updated to correspond to the new data frame.  With the \code{"[<-"}
   method the new value can be a \code{survey.design} object instead of a
   data frame, but only the data frame is used. See also
   \code{\link{subset.survey.design}} for a simple way to select
   subpopulations.

The \code{model.frame} method extracts the observed data.


 If the strata with only one PSU are not self-representing (or they are,
but \code{svydesign} cannot tell based on \code{fpc}) then the handling
of these strata for variance computation is determined by
\code{options("survey.lonely.psu")}.  See \code{\link{svyCprod}} for
details.

\code{data} may be a character string giving the name of a table or view
in a relational database that can be accessed through the \code{DBI} or \code{ODBC}
interfaces. For DBI interfaces \code{dbtype} should be the name of the database
driver and \code{dbname} should be the name by which the driver identifies
the specific database (eg file name for SQLite). For ODBC databases
\code{dbtype} should be \code{"ODBC"} and \code{dbname} should be the
registed DSN for the database. On the Windows GUI, \code{dbname=""} will
produce a dialog box for interactive selection. 

The appropriate database interface package must already be loaded (eg
\code{RSQLite} for SQLite, \code{RODBC} for ODBC).  The survey design
object will contain only the design meta-data, and actual variables will
be loaded from the database as needed.  Use
\code{\link[=close.DBIsvydesign]{close}} to close the database connection and
\code{\link[=open.DBIsvydesign]{open}} to reopen the connection, eg, after
loading a saved object.

The database interface does not attempt to modify the underlying
database and so can be used with read-only permissions on the database.
 

If \code{data} is an \code{imputationList} object (from the "mitools"
package), \code{svydesign} will return a \code{svyimputationList} object
containing a set of designs. Use \code{\link{with.svyimputationList}} to
do analyses on these designs and \code{MIcombine} to combine the results.

}

\value{
An object of class \code{survey.design}.
}
\author{Thomas Lumley}


\seealso{
  \code{\link{as.svrepdesign}} for converting to replicate weight designs,
  \code{\link{subset.survey.design}} for domain estimates,
  \code{\link{update.survey.design}} to add variables.

  \code{mitools} package for using multiple imputations

  \code{\link{svyrecvar}} and \code{\link{svyCprod}} for details of
  variance estimation

  \code{\link{election}} for examples of PPS sampling without replacement.

  \url{http://faculty.washington.edu/tlumley/survey/} for examples of
  database-backed objects.
}



\examples{
  data(api)
# stratified sample
dstrat<-svydesign(id=~1,strata=~stype, weights=~pw, data=apistrat, fpc=~fpc)
# one-stage cluster sample
dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
# two-stage cluster sample: weights computed from population sizes.
dclus2<-svydesign(id=~dnum+snum, fpc=~fpc1+fpc2, data=apiclus2)

## multistage sampling has no effect when fpc is not given, so
## these are equivalent.
dclus2wr<-svydesign(id=~dnum+snum, weights=weights(dclus2), data=apiclus2)
dclus2wr2<-svydesign(id=~dnum, weights=weights(dclus2), data=apiclus2)

## syntax for stratified cluster sample
##(though the data weren't really sampled this way)
svydesign(id=~dnum, strata=~stype, weights=~pw, data=apistrat,
nest=TRUE)

## PPS sampling without replacement
data(election)
dpps<- svydesign(id=~1, fpc=~p, data=election_pps, pps="brewer")

##database example: requires RSQLite
\dontrun{
library(RSQLite)
dbclus1<-svydesign(id=~dnum, weights=~pw, fpc=~fpc,
data="apiclus1",dbtype="SQLite", dbname=system.file("api.db",package="survey"))

}
}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{univar}% __ONLY ONE__ keyword per line
\keyword{manip}% __ONLY ONE__ keyword per line
