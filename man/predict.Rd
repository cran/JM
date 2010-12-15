\name{predict}

\alias{predict.jointModel}

\title{Predictions for Joint Models}

\description{
  Calculates predicted values for the longitudinal part of a joint model.
}

\usage{
\method{predict}{jointModel}(object, newdata, se.fit = FALSE, level = 0.95,
    returnData = FALSE, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{newdata}{a data frame in which to look for variables with which to predict.}
  \item{se.fit}{logical if \code{TRUE} the standard error of the predicted
    values along with a \code{level} confidence interval are
    returned.}
  \item{level}{a numeric value denoting the tolerance/confidence level.}
  \item{returnData}{logical; if \code{TRUE} the data frame supplied in
    \code{newdata} is returned augmented with the outputs of the function.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  This function computed predicted values only for the longitudinal outcome and
  only for the fixed-effects part of the longitudinal submode. In particular,
  let \eqn{X} denote the fixed-effects design matrix calculated using
  \code{newdata}. The \code{predict()} calculates \eqn{\hat{y} = X \hat{\beta}},
  and if \code{se.fit = TRUE}, \eqn{var(\hat{y}) = X V X^t}, with \eqn{V}
  denoting the covariance matrix of \eqn{\hat{\beta}}. The confidence interval
  is based on a normal approximation.
}

\value{
  If \code{se.fit = FALSE} a numeric vector of predicted values, otherwise a
  list with components \code{pred} the predicted values, \code{se.fit} the
  standard error for the fitted values, and \code{low} and \code{upp} the lower
  and upper limits of the confidence interval. If \code{returnData = TRUE}, it
  returns the data frame \code{newdata} with the previously mentioned components
  added.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{jointModel}}}

\examples{
\dontrun{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, 
    random = ~ year | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, 
    data = pbc2.id, x = TRUE)
# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")

DF <- with(pbc2, expand.grid(drug = levels(drug),
    year = seq(min(year), max(year), len = 100)))
Ps <- predict(fitJOINT, DF, se.fit = TRUE, return = TRUE)
require(lattice)
xyplot(pred + low + upp ~ year | drug, data = Ps,
    type = "l", col = c(2,1,1), lty = c(1,2,2), lwd = 2,
    ylab = "Average log serum Bilirubin")
}
}

\keyword{methods}
