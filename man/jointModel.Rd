\name{jointModel}

\alias{jointModel}

\title{Joint Models for Longitudinal and Survival Data}

\description{
This function fits shared parameter models for the joint modelling of normal longitudinal responses and time-to-event
data under a maximum likelihood approach. Various options for the survival model are available.
}

\usage{
jointModel(lmeObject, survObject, timeVar, 
    method = c("weibull-AFT-GH", "weibull-PH-GH", 
    "piecewise-PH-GH", "Cox-PH-GH", "spline-PH-GH", "ch-Laplace"), 
    lag = 0, scaleWB = NULL, init = NULL, control = list(), \dots)
}

\arguments{
  \item{lmeObject}{an object inheriting from class \code{lme} (see also \bold{Note}).}
  \item{survObject}{an object inheriting from class \code{coxph} or class \code{survreg}. In the call to \code{coxph()}
    or \code{survreg()}, you need to specify the argument \code{x = TRUE} such that the design matrix is contained in 
    the object fit. See \bold{Examples}.}
  \item{timeVar}{a character string indicating the time variable in the linear mixed effects model.}
  \item{method}{a character string specifying the type of joint model to fit. See \bold{Details}.}
  \item{lag}{a numeric value denoting a lag effect in the time-dependent covariate represented by the mixed model; default is 0.}
  \item{scaleWB}{a numeric scalar denoting a fixed value for the scale parameter of the Weibull hazard; used only when 
    \code{method = "weibull-AFT-GH"} or \code{method = "weibull-PH-GH"}. The default \code{NULL} means that the scale parameter is estimated.}
  \item{init}{a list of user-specified initial values. The initial values list must have the following components:
    \describe{
        \item{betas}{the vector of fixed effects for the linear mixed effects model.}
        \item{sigma}{the measurement error standard deviation for the linear mixed effects model.}
        \item{D}{the variance-covariance matrix of the random effects.}
        \item{gammas}{the vector of baseline covariates for the survival model. For \code{method = "ch-Laplace"} this vector should 
            first contain initial values for the sorted B-spline coefficients used to model the log cumulative baseline hazard.}
        \item{alpha}{the association parameter.}
        \item{xi}{the vector of baseline risk function values within the intervals specified by the knots; specified only 
            when \code{method = "piecewise-PH-GH"}.}
        \item{gammas.bs}{the vector of spline coefficients; specified only when \code{method = "spline-PH-GH"}.}
        \item{sigma.t}{the scale parameter for the Weibull baseline risk function; specified only when 
            \code{method = "weibull-AFT-GH"} or \code{method = "weibull-PH-GH"}.}
        \item{lambda0}{a vector of the baseline hazard values at the sorted unique event times; specified only when 
            \code{method = "Cox-PH-GH"}.}
    }
    If the user-specified list of initial values does not contain some of these components or contains components 
    not of the appropriate length, then the default initial values are used instead.
  }
  \item{control}{a list of control values with components:
    \describe{
        \item{only.EM}{logical; if \code{TRUE} only the EM algorithm is used in the optimization, otherwise if 
            convergence has not been achieved a quasi-Newton algorithm is initiated. Default is \code{FALSE} except for 
            \code{method = "Cox-PH-GH"} for which only the EM algorithm is available.}
        \item{iter.EM}{the number of EM iterations. Default is 50 except for \code{method = "Cox-PH-GH"} for which
            the default is 200.}
        \item{iter.qN}{the number of quasi-Newton iterations. Default is 150.}
        \item{optimizer}{a character string indicating which optimizer to use; options are "optim" (default) and
            "nlminb".}
        \item{tol1}{tolerance value for convergence in the parameters; see \bold{Details}. Default is 1e-03.}
        \item{tol2}{tolerance value for convergence in the parameters; see \bold{Details}. Default is 1e-04.}
        \item{tol3}{tolerance value for convergence in the log-likelihood; see \bold{Details}. Default is 
            \code{sqrt(.Machine$double.eps)}.}
        \item{numeriDeriv}{a character string indicating which type of numerical derivative to use to compute the
            Hessian matrix; options are "fd" (default) denoting the forward difference approximation, and "cd" 
            denoting the central difference approximation.}
        \item{eps.Hes}{tolerance value used in the numerical derivative method. Default is 1e-06; if you choose 
            \code{numeriDeriv = "cd"} a larger value (e.g., 1e-04) is suggested.}
        \item{parscale}{the \code{parscale} control argument for \code{optim()}, or the \code{scale} argument for 
            \code{nlminb()}. It should be a numeric vector of length equal to the number of parameters. Default is 0.01
            for all parameters.}
        \item{step.max}{tolerance value for the maximum step size in the Newton-Raphson algorithm used to update the 
            parameters of the survival submodel for \code{method = "ch-Laplace"}. Default is 0.1.}
        \item{backtrackSteps}{the number of backtrack steps to use when updating the parameters of the survival submodel
            under \code{method = "ch-Laplace"}.}
        \item{knots}{a numeric vector of the knots positions for the piecewise constant baseline risk function of for 
            the log times used in the B-splines approximation of the log cumulative baseline hazard; therefore, this argument 
            is relevant only when \code{method = "piecewise-PH-GH"}, \code{method = "spline-PH-GH"} or \code{method = "ch-Laplace"}. 
            The default is to place equally-spaced \code{lng.in.kn} knots in the quantiles of the observed event times. For stratified models
            fitted with \code{method = "spline-PH-GH"} this should be a list with elements numeric vectors of knots positions for each strata.}
        \item{lng.in.kn}{the number of internal knots; relevant only when when \code{method = "piecewise-PH-GH"} where it 
            denotes the number of internal knots for the piecewise constant baseline risk function or when \code{method = "spline-PH-GH"} 
            or \code{method = "ch-Laplace"} where it denotes the number of internal knots for B-splines approximation of the log 
            baseline hazard. Default is 6 when \code{method = "piecewise-PH-GH"} and 5 otherwise.}
        \item{equal.strata.knots}{logical; if \code{TRUE} (the default), then the same knots are used in the approximation of the 
            baseline risk function in different strata when \code{method = "spline-PH-GH"}.}
        \item{ord}{a positive integer denoting the order of the B-splines used to approximate the log cumulative hazard
            (default is 4); relevant only when \code{method = "spline-PH-GH"} or \code{method = "ch-Laplace"}.}
        \item{GHk}{the number of Gauss-Hermite quadrature points used to approximate the integrals over the random 
            effects. The default is 15 for one- or two-dimensional integration and 9 otherwise.}
        \item{GKk}{the number of Gauss-Kronrod points used to approximate the integral involved in the calculation of the survival function. Two 
            options are available, namely 7 or 15. For \code{method = "weibull-PH-GH"} and \code{method = "weibull-AFT-GH"} 15 are used, 
            whereas for \code{method = "piecewise-PH-GH"} 7.}
        \item{verbose}{logical; if \code{TRUE}, the parameter estimates and the log-likelihood value are printed during
            the optimization procedure. Default is \code{FALSE}.}
    }
  }
  \item{\dots}{options passed to the \code{control} argument.}
}

\details{
Function \code{jointModel} fits joint models for longitudinal and survival data (more detailed information about the formulation of these
models can be found in Rizopoulos (2010)). For the longitudinal responses the linear mixed effects model represented by the \code{lmeObject} is 
assumed. For the survival times three options are available. In particular, let \eqn{T_i} denote the time-to-event for the \eqn{i}th sample unit, 
\eqn{w_i} denotes the vector of baseline covariates in \code{survObject}, with associated parameter vector \eqn{\gamma}, \eqn{m_i(t)} the value 
of the longitudinal outcome at time point \eqn{t} as approximated by the linear mixed model (i.e., \eqn{m_i(t)} equals the fixed effects part 
\code{+} random effects part of the linear mixed effects model for sample unit \eqn{i}), and \eqn{\alpha} the association parameter. Then, 
for \code{method = "weibull-AFT-GH"} a time-dependent Weibull model under the accelerated failure time formulation is assumed. For 
\code{method = "weibull-PH-GH"} a time-dependent relative risk model is postulated with a Weibull baseline risk function. For 
\code{method = "piecewise-PH-GH"} a time-dependent relative risk model is postulated with a piecewise constant baseline risk function.
For \code{method = "spline-PH-GH"} a time-dependent relative risk model is assumed in which the log baseline risk function is approximated
using B-splines. For \code{method = "ch-Laplace"} an additive model on the log cumulative hazard scale is assumed 
(see Rizopoulos et al., 2009 for more info). Finally, for \code{method = "Cox-PH-GH"} a time-dependent relative risk model 
is assumed where the baseline risk function is left unspecified (Wulfsohn and Tsiatis, 1997).

For \code{method = "spline-PH-GH"} it is also allowed to include stratification factors. These should be included in the specification of
the \code{survObject} using function \code{strata()}. Note that in this case \code{survObject} must only be a 'coxph' object.

For all survival models except for the time-dependent proportional hazards model, the optimization algorithm starts 
with EM iterations, and if convergence is not achieved, it switches to quasi-Newton iterations (i.e., BFGS in 
\code{optim()} or \code{nlminb()}, depending on the value of the \code{optimizer} control argument). For \code{method = "Cox-PH-GH"} only the 
EM algorithm is used. During the EM iterations, convergence is declared if either of the following two conditions is satisfied: (i) \eqn{ L(\theta^{it}) - L(\theta^{it - 1}) < tol_3 
\{ | L(\theta^{it - 1}) | + tol_3 \} }, or (ii) \eqn{ \max \{ | \theta^{it} - \theta^{it - 1} | / ( | \theta^{it - 1} | 
+ tol_1) \} < tol_2}, where \eqn{\theta^{it}} and \eqn{\theta^{it - 1}} is the vector of parameter values at the current 
and previous iterations, respectively, and \eqn{L(.)} is the log-likelihood function. The values for \eqn{tol_1}, 
\eqn{tol_2} and \eqn{tol_3} are specified via the \code{control} argument. During the quasi-Newton iterations, the 
default convergence criteria of either \code{optim()} or \code{nlminb()} are used.

The required integrals are approximated using the Gauss-Hermite quadrature rule for \code{method = "weibull-AFT-GH"}, 
\code{method = "weibull-PH-GH"}, \code{method = "piecewise-PH-GH"}, \code{method = "spline-PH-GH"} and \code{method = "Cox-PH-GH"}, 
whereas for \code{method = "ch-Laplace"} the fully exponential Laplace approximation described in Rizopoulos et al. (2009) is used. 
This last option is more suitable when high-dimensional random effects vectors are considered (e.g., when modelling nonlinear 
subject-specific trajectories with splines or high-order polynomials).
}

\value{
See \code{\link{jointModelObject}} for the components of the fit.
}

\references{
Henderson, R., Diggle, P. and Dobson, A. (2000) Joint modelling of longitudinal measurements and event time data.
\emph{Biostatistics} \bold{1}, 465--480.

Hsieh, F., Tseng, Y.-K. and Wang, J.-L. (2006) Joint modeling of survival and longitudinal data: Likelihood 
approach revisited. \emph{Biometrics} \bold{62}, 1037--1043.

Rizopoulos, D. (2010) JM: an R package for the joint modelling of longitudinal and time-to-event data. Submitted.
(\url{http://rwiki.sciviews.org/doku.php?id=packages:cran:jm#helpful_material})

Rizopoulos, D., Verbeke, G. and Lesaffre, E. (2009) Fully exponential Laplace approximation for the joint
modelling of survival and longitudinal data. \emph{Journal of the Royal Statistical Society, Series B} \bold{71}, 637--654.

Rizopoulos, D., Verbeke, G. and Molenberghs, G. (2009) Multiple-imputation-based residuals and diagnostic plots
for joint models of longitudinal and survival outcomes. \emph{Biometrics} \bold{66}, 20--29.

Tsiatis, A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. 
\emph{Statistica Sinica} \bold{14}, 809--834.

Wulfsohn, M. and Tsiatis, A. (1997) A joint model for survival and longitudinal data measured with error. 
\emph{Biometrics} \bold{53}, 330--339.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}

\note{
1. The \code{lmeObject} should represent a linear mixed model with a simple random-effects structure, i.e., only the 
\code{pdDiag()} class is currently allowed.

2. The \code{lmeObject} should not contain any within-group correlation structure (i.e., \code{correlation} argument of
\code{lme()}) or within-group heteroscedasticity structure (i.e., \code{weights} argument of \code{lme()}).

3. It is assumed that the linear mixed effects model \code{lmeObject} and the survival model \code{survObject} have been 
fitted to the same subjects. Moreover, it is assumed that the ordering of the subjects is the same for both 
\code{lmeObject} and \code{survObject}, i.e., that the first line in the data frame containing the event times 
corresponds to the first set of lines identified by the grouping variable in the data frame containing the repeated 
measurements, and so on.

4. In the \code{print} and \code{summary} generic functions for class \code{jointModel}, the estimated coefficients (and
standard errors for the \code{summary} generic) for the event process are augmented with the element "Assoct" that 
corresponds to the association parameter \eqn{\alpha} (see \bold{Details}).

5. The standard errors returned by the \code{summary} generic function for class \code{jointModel} when 
\code{method = "Cox-PH-GH"} are based on the profile score vector (i.e., given the NPMLE for the unspecified baseline 
hazard). Hsieh et al. (2006) have noted that these standard errors are underestimated.

6. For more examples of \code{jointModel()} check \url{http://rwiki.sciviews.org/doku.php?id=packages:cran:jm}.
}

\seealso{
    \code{\link{jointModelObject}},
    \code{\link{anova.jointModel}},
    \code{\link{coef.jointModel}},
    \code{\link{fixef.jointModel}},
    \code{\link{ranef.jointModel}},
    \code{\link{fitted.jointModel}},
    \code{\link{residuals.jointModel}},
    \code{\link{plot.jointModel}},
    \code{\link{survfitJM}},
    \code{\link{dynC}}
}

\examples{
# linear mixed model fit (random intercepts)
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ 1 | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")
fitJOINT
summary(fitJOINT)

\dontrun{
# linear mixed model fit (random intercepts + random slopes)
fitLME <- lme(log(serBilir) ~ drug * year, random = ~ year | id, data = pbc2)
# survival regression fit
fitSURV <- survreg(Surv(years, status2) ~ drug, data = pbc2.id, x = TRUE)
# joint model fit, under the (default) Weibull model
fitJOINT <- jointModel(fitLME, fitSURV, timeVar = "year")
fitJOINT
summary(fitJOINT)
}

# linear mixed model fit
fitLME <- lme(sqrt(CD4) ~ obstime * drug - drug, 
    random = ~ 1 | patient, data = aids)
# cox model fit
fitCOX <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)
# joint model fit
fitJOINT <- jointModel(fitLME, fitCOX, 
    timeVar = "obstime", method = "spline-PH-GH")
fitJOINT
summary(fitJOINT)
}

\keyword{multivariate}
\keyword{regression}
