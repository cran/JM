\name{ranef}

\alias{ranef.jointModel}

\title{Random Effects Estimates for Joint Models}

\description{
Extracts the random effects estimates from a fitted joint model.
}

\usage{
\method{ranef}{jointModel}(object, postVar = FALSE, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{postVar}{logical; if \code{TRUE} the variance-covariance matrix of these estimates is also returned. 
  See \bold{Details} for more info.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  For methods \code{"ch-GH"}, \code{"weibull-GH"}, \code{"ph-GH"} the \code{ranef()} method returns the 
  posterior means (and posterior variances if \code{postVar = TRUE}), whereas for method \code{"ch-Laplace"} it returns the
  posterior modes (and the inverse of the negative of the second-order derivative of the log-posterior with respect to the 
  random effects, if \code{postVar = TRUE}). 
}

\value{
  a numeric matrix with rows denoting the individuals and columns the random effects (e.g., intercepts, slopes, etc.). 
  If \code{postVar = TRUE}, the numeric matrix has an extra attribute ``postVar".
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{coef.jointModel}}, \code{\link{fixef.jointModel}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ 1 | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)

# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")
ranef(fitJOINT)
}
}

\keyword{methods}
