`fixef.jointModel` <-
function (object, process = c("Longitudinal", "Event"), include.splineCoefs = FALSE, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    process <- match.arg(process)
    if (process == "Longitudinal") {
        object$coefficients$betas
    } else {
        gammas <- object$coefficients$gammas
        if (object$method %in% c("ch-GH", "ch-Laplace") && !include.splineCoefs) {
            ng <- length(gammas)
            nw <- ncol(object$x$W)
            gammas <- if (is.null(nw)) NULL else gammas[seq(ng - nw + 1, ng)]
        }
        c(gammas, "Assoct" = as.vector(object$coefficients$alpha))
    }
}

