\name{plot}

\alias{plot.jointModel}

\title{Plot Diagnostics for Joint Models}

\description{
  Produces a variety of plots for fitted joint models.
}

\usage{
\method{plot}{jointModel}(x, which = 1:4, caption = c("Residuals vs Fitted", 
  "Normal Q-Q", "Marginal Survival", "Marginal Cumulative Hazard", 
  "Marginal log Cumulative Hazard", "Baseline Hazard", 
  "Cumulative Baseline Hazard", "Subject-specific Survival", 
  "Subject-specific Cumulative Hazard", 
  "Subject-specific log Cumulative Hazard"), survTimes = NULL, 
  main = "", 
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
  \dots, ids = NULL, add.smooth = getOption("add.smooth"), 
  add.qqline = TRUE, add.KM = FALSE, cex.caption = 1)
}

\arguments{
  \item{x}{an object inheriting from class \code{jointModel}.}
  \item{which}{which types of plots to produce, specify a subset of the numbers 1:10.}
  \item{caption}{captions to appear above the plots defined by argument \code{which}.}
  \item{survTimes}{a vector of survival times for which the survival, cumulative hazard or 
    log cumulative hazard will be computed. Default is \code{seq(minT, maxT, length = 15)}, where \code{minT} and 
    \code{maxT} are the minimum and maximum observed survival times, respectively.}
  \item{main}{a character string specifying the title in the plot.}
  \item{ask}{logical; if \code{TRUE}, the user is asked before each plot, see \code{par(ask=.)}.}
  \item{\dots}{other parameters to be passed through to plotting functions.}
  \item{ids}{a numeric vector specifying which subjects, the subject-specific plots will include; 
    default is all subjects.}
  \item{add.smooth}{logical; if \code{TRUE} a smooth line is superimposed in the "Residuals vs Fitted" plot.}
  \item{add.qqline}{logical; if \code{TRUE} a qq-line is superimposed in the "Normal Q-Q" plot.}
  \item{add.KM}{logical; if \code{TRUE} the Kaplan-Meier estimate of the survival function is superimposed in the 
    "Marginal Survival" plot.}
  \item{cex.caption}{magnification of captions.}
}

\note{
  The plots of the baseline hazard and the cumulative baseline hazard are only produced when the joint model has
  been fitted using \code{method = "ph-GH"}.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\seealso{\code{\link{jointModel}}}

\examples{
# linear mixed model fit
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ 1 | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")

plot(fitJOINT, 3, add.KM = TRUE, col = "red", lwd = 2)

par(mfrow = c(2, 2))
plot(fitJOINT)

}

\keyword{methods}
