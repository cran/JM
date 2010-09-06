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
    if (object$method == "Cox-PH-GH") {
        gammas <- c(object$coefficients$gammas, "Assoct" = as.vector(object$coefficients$alpha),
            "Assoct.s" = as.vector(object$coefficients$Dalpha))
        indT <- grep("T.", colnames(VarCov), fixed = TRUE)
    } else if (object$method == "weibull-PH-GH") {
        gammas <- c(object$coefficients$gammas, "Assoct" = as.vector(object$coefficients$alpha),
            "Assoct.s" = as.vector(object$coefficients$Dalpha), "log(scale)" = log(as.vector(object$coefficients$sigma.t)))
        indT <- grep("T.", colnames(VarCov), fixed = TRUE)
    } else if (object$method == "weibull-AFT-GH") {
        gammas <- c(object$coefficients$gammas, "Assoct" = as.vector(object$coefficients$alpha),
            "Assoct.s" = as.vector(object$coefficients$Dalpha), "log(scale)" = log(as.vector(object$coefficients$sigma.t)))
        gammas[seq(1, length(gammas) - 1)] <- - gammas[seq(1, length(gammas) - 1)]
        indT <- grep("T.", colnames(VarCov), fixed = TRUE)
    } else if (object$method == "piecewise-PH-GH") {
        gammas <- c(object$coefficients$gammas, "Assoct" = as.vector(object$coefficients$alpha),
            "Assoct.s" = as.vector(object$coefficients$Dalpha), log(as.vector(object$coefficients$xi)))
        names(gammas)[seq(length(gammas) - object$x$Q + 1, length(gammas))] <- paste("log(xi.", seq_len(object$x$Q), ")", sep = "")
        indT <- grep("T.", colnames(VarCov), fixed = TRUE)
    } else if (object$method == "spline-PH-GH") {
        gammas <- c(object$coefficients$gammas, "Assoct" = as.vector(object$coefficients$alpha),
            "Assoct.s" = as.vector(object$coefficients$Dalpha), object$coefficients$gammas.bs)
        indT <- grep("T.", colnames(VarCov), fixed = TRUE)
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
    sds <- if (length(indT) > 1) sqrt(diag(VarCov[indT, indT])) else sqrt(VarCov[indT, indT])
    if (!is.null(object$scaleWB))
        sds <- c(sds, NA)
    coefsT <- cbind("Value" = gammas, "Std.Err" = sds, "z-value" = gammas / sds,
        "p-value" = 2 * pnorm(abs(gammas / sds), lower.tail = FALSE))
    out <- list("CoefTable-Long" = coefsY, "CoefTable-Event" = coefsT, D = object$coefficients$D, 
        sigma = object$coefficients$sigma, logLik = as.vector(logLik(object)), AIC = AIC(object), 
        BIC = AIC(object, k = log(object$n)))
    out$N <- object$N
    out$n <- object$n
    out$d <- object$d
    out$id <- object$id
    out$method <- object$method
    out$control <- object$control
    out$knots <- unique(object$knots)
    out$conv <- object$conv
    out$parameterization <- object$parameterization
    out$call <- object$call
    class(out) <- "summary.jointModel"
    out
}

