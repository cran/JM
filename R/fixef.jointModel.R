fixef.jointModel <-
function (object, process = c("Longitudinal", "Event"), 
        include.splineCoefs = FALSE, ...) {
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
        out <- c(gammas, "Assoct" = as.vector(object$coefficients$alpha), 
            "Assoct.s" = as.vector(object$coefficients$Dalpha))
        if (object$method == "weibull-AFT-GH")
            out <- - out
        if ((lag <- object$y$lag) > 0) {        
            if (object$parameterization %in% c("value", "both")) {
                ii <- names(gammas) == "Assoct"
                names(gammas)[ii] <- paste("Assoct(lag=", lag, ")", sep = "")
            }
            if (object$parameterization %in% c("slope", "both")) { 
                jj <- names(gammas) == "Assoct.s"
                names(gammas)[jj] <- paste("Assoct.s(lag=", lag, ")", sep = "")
            }
        }
        out
    }
}

