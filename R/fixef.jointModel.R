fixef.jointModel <-
function (object, process = c("Longitudinal", "Event"), include.splineCoefs = FALSE, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    process <- match.arg(process)
    if (process == "Longitudinal") {
        object$coefficients$betas
    } else {
        gammas <- object$coefficients$gammas
        if (object$method == "ch-Laplace" && !include.splineCoefs) {
            ng <- length(gammas)
            nw <- ncol(object$x$W)
            gammas <- if (is.null(nw)) NULL else gammas[seq(ng - nw + 1, ng)]
        }
        out <- c(gammas, "Assoct" = as.vector(object$coefficients$alpha))
        if (object$method == "weibull-AFT-GH")
            out <- - out
        if ((lag <- object$y$lag) > 0) {
            ii <- names(out) == "Assoct"
            names(out)[ii] <- paste("Assoct(lag=", lag, ")", sep = "")
        }
        out
    }
}

