summary.jointModel <-
function (object, ...) {
    if (!inherits(object, "jointModel"))
        stop("Use only with 'jointModel' objects.\n")
    VarCov <- vcov(object)
    betas <- object$coefficients$betas
    indY <- seq(1, length(betas))
    sds <- sqrt(diag(VarCov[indY, indY]))
    coefsY <- cbind("Value" = betas, "Std.Err" = sds, "z-value" = betas / sds, 
        "p-value" = 2 * pnorm(abs(betas / sds), lower.tail = FALSE))
    if (object$method == "ph-GH") {
        gammas <- c(object$coefficients$gammas, "Assoct" = as.vector(object$coefficients$alpha))
        indT <- seq(length(betas) + 2, length(betas) + length(gammas) + 1)
    } else if (object$method == "weibull-PH-GH") {
        gammas <- c(object$coefficients$gammas, "Assoct" = as.vector(object$coefficients$alpha),
            "log(scale)" = log(as.vector(object$coefficients$sigma.t)))
        indT <- seq(length(betas) + 2, length(betas) + length(gammas) + 1)
    } else if (object$method == "weibull-AFT-GH") {
        gammas <- c(-object$coefficients$gammas, "Assoct" = -as.vector(object$coefficients$alpha),
            "log(scale)" = log(as.vector(object$coefficients$sigma.t)))
        indT <- seq(length(betas) + 2, length(betas) + length(gammas) + 1)
    } else if (object$method == "piecewise-PH-GH") {
        gammas <- c(object$coefficients$gammas, "Assoct" = as.vector(object$coefficients$alpha),
            log(as.vector(object$coefficients$xi)))
        names(gammas)[seq(length(gammas) - object$x$Q + 1, length(gammas))] <- paste("log(xi.", seq_len(object$x$Q), ")", sep = "")
        indT <- seq(length(betas) + 2, length(betas) + length(gammas) + 1)
    } else {
        gms <- object$coefficients$gammas
        ng <- length(gms)
        nw <- ncol(object$x$W)
        if (is.null(nw))
            nw <- 0
        gms <- gms[- seq(1, ng - nw)]
        gammas <- c(gms, "Assoct" = as.vector(object$coefficients$alpha))
        indT <- seq(length(betas) + 2 + ng - nw, length(betas) + ng + 2)
    }
    sds <- if (length(indT) > 1) sqrt(diag(VarCov[indT, indT])) else sqrt(VarCov[indT, indT])
    coefsT <- cbind("Value" = gammas, "Std.Err" = sds, "z-value" = gammas / sds,
        "p-value" = 2 * pnorm(abs(gammas / sds), lower.tail = FALSE))
    out <- list("CoefTable-Long" = coefsY, "CoefTable-Event" = coefsT, D = object$coefficients$D, 
        sigma = object$coefficients$sigma, logLik = as.vector(logLik(object)), AIC = AIC(object), 
        BIC = AIC(object, k = log(object$N)))
    out$N <- object$N
    out$n <- object$n
    out$d <- object$d
    out$id <- object$id
    out$method <- object$method
    out$control <- object$control
    out$knots <- unique(object$knots)
    out$conv <- object$conv
    out$call <- object$call
    class(out) <- "summary.jointModel"
    out
}

