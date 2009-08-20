\name{JM}

\alias{JM-package}
\alias{JM}

\docType{package}

\title{
Joint Modelling of Longitudinal and Time-to-Event Data in R
}

\description{
This package fits shared parameter models for the joint modelling of normal longitudinal responses and event times 
under a maximum likelihood approach. Various options for the survival model and optimization/integration algorithms are 
provided.
}

\details{
\tabular{ll}{
Package: \tab JM\cr
Type: \tab Package\cr
Version: \tab 0.4-0\cr
Date: \tab 2009-08-20\cr
License: \tab \acronym{GPL} \cr
}

The package has a single model-fitting function called \code{\link{jointModel}}, which accepts as main arguments a linear
mixed effects object fit returned by function \code{lme()} of package \pkg{nlme}, and a survival object fit returned 
by either function \code{coxph()} or function \code{survreg()} of package \pkg{survival}. In addition, the \code{method}
argument of \code{jointModel()} specifies the type of the survival submodel to be fitted and the type of the numerical 
integration technique; available options are:
\describe{
    \item{\code{"ph-GH"}}{the time-dependent version of a proportional hazards model with unspecified baseline hazard 
        function. The Gauss-Hermite integration rule is used to approximate the required integrals. (This option 
        corresponds to the joint model proposed by Wulfsohn and Tsiatis, 1997)}
    \item{\code{"weibull-PH-GH"}}{the Weibull model under the proportional hazards formulation. The Gauss-Hermite 
        integration rule is used to approximate the required integrals.}
    \item{\code{"weibull-AFT-GH"}}{the Weibull model under the accelerated failure time formulation. The Gauss-Hermite 
        integration rule is used to approximate the required integrals.}
    \item{\code{"piecewise-PH-GH"}}{a proportional hazards model with a piecewise constant baseline risk function. The Gauss-Hermite 
        integration rule is used to approximate the required integrals.}    
    \item{\code{"ch-GH"}}{an additive log cumulative hazard model, in which the log cumulative baseline hazard is 
        approximated using B-splines. The Gauss-Hermite integration rule is used to approximate the required integrals.}
    \item{\code{"ch-Laplace"}}{an additive log cumulative hazard model, in which the log cumulative baseline hazard is 
        approximated using B-splines. A fully exponential Laplace approximation method is used to approximate the 
        required integrals (Rizopoulos et al., 2009).}
}
}

\author{
Dimitris Rizopoulos

Maintainer: Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>
}

\references{
Henderson, R., Diggle, P. and Dobson, A. (2000) Joint modelling of longitudinal measurements and event time data.
\emph{Biostatistics} \bold{1}, 465--480.

Rizopoulos, D., Verbeke, G. and Lesaffre, E. (2009) Fully exponential Laplace approximation for the joint
modelling of survival and longitudinal data. \emph{Journal of the Royal Statistical Society, Series B} \bold{71}, 637--654.

Rizopoulos, D., Verbeke, G. and Molenberghs, G. (2009) Multiple-imputation-based residuals and diagnostic 
plots for joint models of longitudinal and survival outcomes. \emph{Biometrics}, to appear (doi: 10.1111/j.1541-0420.2009.01273.x).

Tsiatis, A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. 
\emph{Statistica Sinica} \bold{14}, 809--834.

Wulfsohn, M. and Tsiatis, A. (1997) A joint model for survival and longitudinal data measured with error. 
\emph{Biometrics} \bold{53}, 330--339.
}

\seealso{\code{\link{jointModel}}}

\keyword{package}
\keyword{multivariate}
