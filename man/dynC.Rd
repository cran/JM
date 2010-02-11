\name{dynC}
\alias{dynC}
\alias{dynC.jointModel}

\title{ Dynamic Discrimination Index for Joint Models }

\description{
  It computes the dynamic discrimination index presented in Rizopoulos (2010).
}

\usage{
dynC(object, \dots)

\method{dynC}{jointModel}(object, dt, data, idVar = "id", 
    times = NULL, nt = 10, simulate = FALSE, estimator = c("median", "mean"),
    M = 10, validate = FALSE, B = 10, verbose = FALSE, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jointModel}.}
  \item{dt}{a numeric value indicating the length of the interval of primary interest within which we want to 
        distinguish between subjects who died within the interval from subjects who survived longer than that.}
  \item{data}{a data frame that contains all variables use in the fit of the joint model. 
        This should include both the variables for the longitudinal and survival parts.}
  \item{idVar}{a character string identifying the variable in \code{data} that distinguishes between subjects.}
  \item{times}{at which time points to compute the AUCs in order to obtain the dynamic C index.}
  \item{nt}{the number of time points to compute the AUCs, if the number of unique visit times in the sample is greater than \code{nt}.}
  \item{simulate}{the \code{simulate} argument of \code{\link{survfitJM}}.}
  \item{estimator}{the \code{estimator} argument of \code{\link{survfitJM}}.}
  \item{M}{the \code{M} argument of \code{\link{survfitJM}}.}
  \item{validate}{logical; if \code{TRUE} the Bootstrap validated index is produced -- 
    \bold{this feature is not yet implemented}.}
  \item{B}{numeric value indicating how many Bootstrap samples to use -- \bold{this feature is not yet implemented}.}
  \item{verbose}{logical; if \code{TRUE} information is printed after each Bootstrap sample -- \bold{this feature is not yet implemented}.}
  \item{\dots}{extra arguments; currently none is used.}
}

\details{
More details can be found in Rizopoulos (2010).
}

\value{
  An object of class \code{dynC} is a list with components,
  \item{AUCt}{A numeric matrix with the values of the Area Under the ROC Curve for the different time points.}
  \item{dynC}{the value of the dynamic C index.}
}

\references{
Antolini, L., Boracchi, P., and Biganzoli, E. (2005). A time-dependent discrimination index
for survival data. \emph{Statistics in Medicine} \bold{24}, 3927--3944.

Heagerty, P. and Zheng, Y. (2005). Survival model predictive accuracy and ROC curves.
\emph{Biometrics} \bold{61}, 92--105.

Rizopoulos, D. (2010). Dynamic prediction in joint models for longitudinal and time-to-event data. 
\emph{Submitted}.
}

\author{
    Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
    \code{\link{jointModel}},
    \code{\link{survfitJM}}
}

\examples{
\dontrun{
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:(drug + AZT + prevOI + gender), 
    random = ~ obstime | patient, data = aids)
fitSURV <- coxph(Surv(Time, death) ~ drug + AZT + prevOI + gender, 
    data = aids.id, x = TRUE)
fit.aids <- jointModel(fitLME, fitSURV, timeVar = "obstime", method = "piecewise-PH-GH")

# the following will take some time to execute...
C.JMaids1 <- dynC(fit.aids, dt = 2, data = aids, idVar = "patient")
C.JMaids1

C.JMaids2 <- dynC(fit.aids, dt = 6, data = aids, idVar = "patient")
C.JMaids2

C.JMaids3 <- dynC(fit.aids, dt = 8, data = aids, idVar = "patient")
C.JMaids3
}
}

\keyword{methods}
