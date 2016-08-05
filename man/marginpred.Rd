\name{marginpred}
\Rdversion{1.1}
\alias{marginpred}
\alias{marginpred.svycoxph}
\alias{marginpred.svykmlist}
\alias{marginpred.svyglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Standardised predictions (predictive margins) for regression models.
}
\description{
Reweights the design (using \code{\link{calibrate}}) so that the adjustment variables are uncorrelated
with the variables in the model, and then performs predictions by
calling \code{predict}.  When the adjustment model is saturated this is
equivalent to direct standardization on the adjustment variables. 

The \code{svycoxph} and \code{svykmlist} methods return survival curves.
}
\usage{
marginpred(model, adjustfor, predictat, ...)
\method{marginpred}{svycoxph}(model, adjustfor, predictat, se=FALSE, ...)
\method{marginpred}{svykmlist}(model, adjustfor, predictat, se=FALSE, ...)
\method{marginpred}{svyglm}(model, adjustfor, predictat,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    A regression model object of a class that has a \code{marginpred} method
  }
  \item{adjustfor}{
    Model formula specifying adjustment variables, which must be in the
    design object of the model
  }
  \item{predictat}{
    A data frame giving values of the variables in \code{model} to
    predict at}
  \item{se}{Estimate standard errors for the survival curve (uses a lot
    of memory if the sample size is large)}
  \item{\dots}{Extra arguments, passed to the \code{predict} method for \code{model}}
}

\seealso{
  \code{\link{svypredmeans}} for the method of Graubard and Korn implemented in SUDAAN.

  \code{\link{calibrate}}

  \code{\link{predict.svycoxph}}
}
\examples{
## generate data with apparent group effect from confounding
set.seed(42)
df<-data.frame(x=rnorm(100))
df$time<-rexp(100)*exp(df$x-1)
df$status<-1
df$group<-(df$x+rnorm(100))>0
des<-svydesign(id=~1,data=df)
newdf<-data.frame(group=c(FALSE,TRUE), x=c(0,0))

## Cox model
m0<-svycoxph(Surv(time,status)~group,design=des)
m1<-svycoxph(Surv(time,status)~group+x,design=des)
## conditional predictions, unadjusted and adjusted
cpred0<-predict(m0, type="curve", newdata=newdf, se=TRUE)
cpred1<-predict(m1, type="curve", newdata=newdf, se=TRUE)
## adjusted marginal prediction
mpred<-marginpred(m0, adjustfor=~x, predictat=newdf, se=TRUE)

plot(cpred0)
lines(cpred1[[1]],col="red")
lines(cpred1[[2]],col="red")
lines(mpred[[1]],col="blue")
lines(mpred[[2]],col="blue")

## Kaplan--Meier
s2<-svykm(Surv(time,status>0)~group, design=des)
p2<-marginpred(s2, adjustfor=~x, predictat=newdf,se=TRUE)
plot(s2)
lines(p2[[1]],col="green")
lines(p2[[2]],col="green")

## logistic regression
logisticm <- svyglm(group~time, family=quasibinomial, design=des)
newdf$time<-c(0.1,0.8)
logisticpred <- marginpred(logisticm, adjustfor=~x, predictat=newdf)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}

