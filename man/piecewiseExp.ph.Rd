\name{piecewiseExp.ph}

\alias{piecewiseExp.ph}

\title{ Proportional Hazards Models with Piecewise Constant Baseline Hazard Function }

\description{
  Based on a fitted Cox model this function fits the corresponding relative risk model with a 
	piecewise constant baseline hazard using the Poisson regression equivalence
}

\usage{
piecewiseExp.ph(coxObject, knots = NULL, length.knots = 6)
} 

\arguments{
  \item{coxObject}{ an object of class \code{coxph}.}
  \item{knots}{A numeric vector denoting the internal knots (cut points) defining the intervals in which the baseline hazard is assumed constant.}
  \item{length.knots}{a numeric value denoting the number of internal knots to use in the fit. 
		Used when \code{knots = NULL}.}
}

\value{
  an object of class \code{glm}.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\examples{
coxFit <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)
piecewiseExp.ph(coxFit)
}

\keyword{multivariate}
\keyword{regression}

