\name{wald.strata}

\alias{wald.strata}

\title{ Wald Test for Stratification Factors }

\description{
  It performs a Wald test to test the hypothesis of equal spline coefficients among strata 
  in the approximation of baseline risk function.
}

\usage{
wald.strata(fit)
} 

\arguments{
  \item{fit}{an object of class \code{jointModel} with \code{method = "spline-PH-GH"}
  and with a strata specification in the survival part.}
}

\value{
  an object of class \code{wald.strata} with components:
  \item{alternative}{a character string naming the alternative.}
  \item{Result}{a numeric matrix with the results of the Wald test.}
}

\note{
This test is valid when the same knots have been used across strata.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
\dontrun{
fitLME <- lme(log(serBilir) ~ drug * year - drug, random = ~ year | id, 
    data = pbc2)
fitSURV <- coxph(Surv(years, status2) ~ drug + strata(hepatomegaly), 
    data = pbc2.id, x = TRUE)
fit.pbc <- jointModel(fitLME, fitSURV, timeVar = "year", method = "spline-PH-GH")
wald.strata(fit.pbc)
}
}

\keyword{multivariate}
\keyword{regression}

