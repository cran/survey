\name{svytable}
\alias{svreptable}
\alias{svytable}
\alias{svytable.svyrep.design}
\alias{svytable.survey.design}
\alias{svychisq}
\alias{svychisq.survey.design}
\alias{svychisq.svyrep.design}
\alias{summary.svytable}
\alias{print.summary.svytable}
\alias{summary.svreptable}
\alias{degf}
\alias{degf.svyrep.design}
\alias{degf.survey.design2}
\alias{degf.twophase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Contingency tables for survey data}
\description{
  Contingency tables and chisquared tests of association for survey data.
}
\usage{
\method{svytable}{survey.design}(formula, design, Ntotal = NULL, round = FALSE,...)
\method{svytable}{svyrep.design}(formula, design, Ntotal = sum(weights(design, "sampling")), round = FALSE,...)
\method{svychisq}{survey.design}(formula, design, 
   statistic = c("F",  "Chisq","Wald","adjWald","lincom","saddlepoint"),na.rm=TRUE,...)
\method{svychisq}{svyrep.design}(formula, design, 
   statistic = c("F",  "Chisq","Wald","adjWald","lincom","saddlepoint"),na.rm=TRUE,...)
\method{summary}{svytable}(object, 
   statistic = c("F","Chisq","Wald","adjWald","lincom","saddlepoint"),...)
degf(design, ...)
\method{degf}{survey.design2}(design, ...)
\method{degf}{svyrep.design}(design, tol=1e-5,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula specifying margins for the table (using \code{+} only)}
  \item{design}{survey object}
  \item{statistic}{See Details below}
  \item{Ntotal}{A population total or set of population stratum totals
    to normalise to.}
  \item{round}{Should the table entries be rounded to the nearest
    integer?}
  \item{na.rm}{Remove missing values}
  \item{object}{Output from \code{svytable}}
  \item{...}{For \code{svytable} these are passed to \code{xtabs}. Use
    \code{exclude=NULL}, \code{na.action=na.pass} to include \code{NA}s
    in the table}
  \item{tol}{Tolerance for \code{\link{qr}} in computing the matrix rank}
 }
\details{
 
The \code{svytable} function computes a weighted crosstabulation. This
is especially useful for producing graphics. It is sometimes easier
to use \code{\link{svytotal}} or \code{\link{svymean}}, which also
produce standard errors, design effects, etc.

The frequencies in the table can be normalised to some convenient total
such as 100 or 1.0 by specifying the \code{Ntotal} argument.  If the
formula has a left-hand side the mean or sum of this variable rather
than the frequency is tabulated.

The \code{Ntotal} argument can be either a single number or a data
frame whose first column gives the (first-stage) sampling strata and
second column the population size in each stratum.  In this second case
the \code{svytable} command performs `post-stratification': tabulating
and scaling to the population within strata and then adding up the
strata.

As with other \code{xtabs} objects, the output of \code{svytable} can be
processed by \code{ftable} for more attractive display. The
\code{summary} method for \code{svytable} objects calls \code{svychisq}
for a test of independence.

\code{svychisq} computes first and second-order Rao-Scott corrections to
the Pearson chisquared test, and two Wald-type tests.

The default (\code{statistic="F"}) is the Rao-Scott second-order
correction.  The p-values are computed with a Satterthwaite
approximation to the distribution and with denominator degrees of
freedom as recommended by Thomas and Rao (1990). The alternative
\code{statistic="Chisq"} adjusts the Pearson chisquared statistic by a
design effect estimate and then compares it to the chisquared
distribution it would have under simple random sampling.

The \code{statistic="Wald"} test is that proposed by Koch et al (1975)
and used by the SUDAAN software package. It is a Wald test based on the
differences between the observed cells counts and those expected under
independence. The adjustment given by \code{statistic="adjWald"} reduces
the statistic when the number of PSUs is small compared to the number of
degrees of freedom of the test. Thomas and Rao (1990) compare these
tests and find the adjustment benefical.

\code{statistic="lincom"} replaces the numerator of the Rao-Scott F with
the exact asymptotic distribution, which is a linear combination of
chi-squared variables (see \code{\link{pchisqsum}}, and
\code{statistic="saddlepoint"} uses a saddlepoint approximation to this
distribution.  The \code{CompQuadForm} package is needed for
\code{statistic="lincom"} but not for
\code{statistic="saddlepoint"}. The saddlepoint approximation is
especially useful when the p-value is very small (as in large-scale
multiple testing problems).

For designs using replicate weights the code is essentially the same as
for designs with sampling structure, since the necessary variance
computations are done by the appropriate methods of
\code{\link{svytotal}} and \code{\link{svymean}}.  The exception is that
the degrees of freedom is computed as one less than the rank of the
matrix of replicate weights (by \code{degf}).


At the moment, \code{svychisq} works only for 2-dimensional tables.
}
\value{
  The table commands return an \code{xtabs} object, \code{svychisq}
  returns a \code{htest} object.
}
\references{
Davies RB (1973). "Numerical inversion of a characteristic function"
Biometrika 60:415-7

P. Duchesne, P. Lafaye de Micheaux (2010) "Computing the distribution of
quadratic forms: Further comparisons between the Liu-Tang-Zhang
approximation and exact methods", Computational Statistics and Data
Analysis, Volume 54,  858-862

Koch, GG, Freeman, DH, Freeman, JL (1975) "Strategies in the
multivariate analysis of data from complex surveys" International
Statistical Review 43: 59-78
  
Rao, JNK, Scott, AJ (1984) "On Chi-squared Tests For Multiway
Contigency Tables with Proportions Estimated From Survey Data"  Annals
of Statistics 12:46-60.

Sribney WM (1998) "Two-way contingency tables for survey or clustered
data" Stata Technical Bulletin 45:33-49.

Thomas, DR, Rao, JNK (1990) "Small-sample comparison of level and power
for simple goodness-of-fit statistics under cluster sampling" JASA 82:630-636
}

\note{Rao and Scott (1984) leave open one computational issue. In
  computing `generalised design effects' for these tests, should the
  variance under simple random sampling be estimated using the observed
  proportions or the the predicted proportions under the null
  hypothesis? \code{svychisq} uses the observed proportions, following
  simulations by Sribney (1998), and the choices made in Stata}


\seealso{\code{\link{svytotal}} and \code{\link{svymean}} report totals
  and proportions by category for factor variables.
  
  See \code{\link{svyby}} and \code{\link{ftable.svystat}} to construct
  more complex tables of summary statistics.

  See \code{\link{svyloglin}} for loglinear models.

  See \code{\link{regTermTest}} for Rao-Scott tests in regression models.
}
\examples{
  data(api)
  xtabs(~sch.wide+stype, data=apipop)

  dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
  summary(dclus1)

  (tbl <- svytable(~sch.wide+stype, dclus1))
  plot(tbl)
  fourfoldplot(svytable(~sch.wide+comp.imp+stype,design=dclus1,round=TRUE), conf.level=0)

  svychisq(~sch.wide+stype, dclus1)
  summary(tbl, statistic="Chisq")
  svychisq(~sch.wide+stype, dclus1, statistic="adjWald")

  rclus1 <- as.svrepdesign(dclus1)
  summary(svytable(~sch.wide+stype, rclus1))
  svychisq(~sch.wide+stype, rclus1, statistic="adjWald")

}
\keyword{survey}% at least one, from doc/KEYWORDS
\keyword{category}% __ONLY ONE__ keyword per line
\keyword{htest}% __ONLY ONE__ keyword per line
