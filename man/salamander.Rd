\name{salamander}
\docType{data}
\alias{salamander}
\concept{generalized linear mixed model}
\title{Salamander mating data set from McCullagh and Nelder (1989)}
\description{This data set presents the outcome of three experiments
  conducted at the University of Chicago in 1986 to study interbreeding
  between populations of mountain dusky salamanders (McCullagh and
  Nelder, 1989, Section 14.5). The analysis here is from Lumley (1998,
  section 5.3)}

\usage{data(salamander)}
\format{
  A data frame with the following columns:

   \describe{
\item{Mate}{Whether the salamanders mated (1) or did not mate (0).}
\item{Cross}{Cross between female and male type. A factor with four levels: \code{R/R},\code{R/W},\code{W/R}, and \code{W/W}. The type of the female salamander is listed first and the male is listed second. Rough Butt is represented by R and White Side is represented by W. For example, \code{Cross=W/R} indicates a White Side female was crossed with a Rough Butt male.}
\item{Male}{Identification number of the male salamander. A factor.}
\item{Female}{Identification number of the female salamander. A factor.}

}
}


\references{McCullagh P. and Nelder, J. A. (1989)  \emph{Generalized
    Linear Models}.  Chapman and Hall/CRC.
Lumley T (1998) PhD thesis, University of Washington

}


\examples{
data(salamander)
salamander$mixed<-with(salamander, Cross=="W/R" | Cross=="R/W")
salamander$RWvsWR<-with(salamander,  ifelse(mixed,
          ((Cross=="R/W")-(Cross=="W/R"))/2,
          0))
xsalamander<-xdesign(id=list(~Male, ~Female), data=salamander,
    overlap="unbiased")

## Adjacency matrix
## Blocks 1 and 2 are actually the same salamanders, but
## it's traditional to pretend they are independent.
image(xsalamander$adjacency)

## R doesn't allow family=binomial(identity)
success <- svyglm(Mate~mixed+RWvsWR, design=xsalamander,
    family=quasi(link="identity", variance="mu(1-mu)"))
summary(success)
}
\keyword{datasets}

