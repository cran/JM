\name{residuals}

\alias{residuals.jointModel}

\title{Residuals for Joint Models}

\description{
  Calculates residuals for joint models.
}

\usage{
\method{residuals}{jointModel}(object, process = c("Longitudinal", "Event"), 
  type = c("Marginal", "Subject", "stand-Marginal", 
    "stand-Subject"), \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{process}{for which model (i.e., linear mixed model or survival model) to calculate residuals.}
  \item{type}{what type of residuals to calculate. See \bold{Details}.}
  \item{\dots}{additional arguments; currently none is used.}
}

\details{
  Currently residuals are calculated for the longitudinal outcome. If \code{type = "Marginal"} these are 
  \eqn{e_{ij} = y_{ij} - x_{ij}^T \hat{\beta}}, whereas for \code{type = "Subject"}, \eqn{e_{ij} = y_{ij}
   - x_{ij}^T \hat{\beta} - z_{ij}^T b_i}, where \eqn{i} denotes the subject and \eqn{j} the measurement,
   \eqn{y_{ij}} the longitudinal responses, \eqn{x_{ij}^T} and \eqn{z_{ij}^T} the corresponding rows of the
   fixed and random effects design matrices, respectively, and \eqn{\beta} and \eqn{b_i} denote the fixed effects 
   and random effects components. 
  
  If \code{type = "stand-Marginal"} or \code{type = "stand-Subject"}, the above defined residuals are divided by
  the estimated standard deviation of the error term.
}

\value{
  a numeric vector of residual values.
}

\author{Dimitris Rizopoulos \email{dimitris.rizopoulos@med.kuleuven.be}}

\note{
 Residuals are not currently computed for the survival outcome.
}

\seealso{\code{\link{fitted.jointModel}}}

\examples{
# linear mixed model fit
fitLME <- lme(sqrt(CD4) ~ obstime * drug - drug, 
    random = ~ 1 | patient, data = aids)
# cox model fit
fitCOX <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)

# joint model fit, under the additive log cumulative hazard model
fitJOINT <- jointModel(fitLME, fitCOX, 
    timeVar = "obstime", method = "ch-GH")

head(cbind(
    "Marginal" = residuals(fitJOINT),
    "std-Marginal" = residuals(fitJOINT, type = "stand-Marginal"),
    "Subject" = residuals(fitJOINT, type = "Subject"),
    "std-Subject" = residuals(fitJOINT, type = "stand-Subject")
))
}

\keyword{methods}
